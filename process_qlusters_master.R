###########
#  Process qluster master code. 
# Based in the code at: ./codes/process_clusters/process_qlusters_chim_jan19.R
# Last modification - June 2019 :  Generalized version compacted to a single function. 
#  Author: Diana Vera Cruz
###########

##################
######  STEPS
##################

# 1. Read replicas datasets.
# 2. Shared haplotypes.
# 3. Merge non-shared haplotypes to their closest shared haplotype. 
# 4. Relative frequency cutoff -> 0.5% from Synthetic plasmids -> parameter to adjust it. 
# 5. Take out chimeras :  Partial alignment  + Parents existence in sample. 
# 6. Store results:  Tab & seqs. 

##################
###### MAIN FUNCTION
##################

## PARAMETERS: 
#### w_dir = working directory that contains all the study/dataset information.
#### subset_dir = Subfirectory name / adresss from w_dir containes qluster results. 
#### gene = Gene name that also identifies reference sequences and directories for the haplotypes per sample. 
#### cut_off = Minimum threshold for haplotypes shared -> default is 0.5% 
#### chim = T/F logical, chimeras need to be detected and labeled or not.
#### collapse_non_shared = T/F logical. Collapse sequences with haplotype frequencies below cut-off to their closest valid haplotype.
#### s_names = NULL -> replicates names/IDs to be excluded from the analysis.
#### summ_dir -> Summary directory to output results. 
         ## w_dir/summ_dir
         ## w_dir/summ_dir/gene ->  summary for all the samples of such gene. 
          ## w_dir/summ_dir/gene/(pq/graphs) ->   pq has csv for each sample and graphs pdfs of correlation/slope graphs. 

## Process clusters main function. 
process_qlusters = function(w_dir, subset_dir, gene, summ_dir = 'summary', cut_off=0.005, chim = T, collapse_non_shared = T, s_names = NULL){
  require(ShortRead); require(ape); require(magrittr); require(reshape2); require(ggplot2); require(gridExtra)
  
  ################
  ##  0. FUNCTIONS 
  ################
  ## Function to retrieve all the files from Seekdeep qluster ->  2 levels:  Samples - Replicas & check that qluster is done. 
  #### w_dir -> Directory containing the samples directories for SeekDeep qluster results. (Single dir)
  ## Returns a dataframe with samples IDs, replicas and directories. (3 cols, rows = num total replicas)
  process_qlusters.sam_table = function(w_dir){  # w_dir = './haplotypes2'
    ##  Retrieve  samples/replicas dirs. 
    d = list.files( w_dir, full.names=T)
    d = sapply(d, list.files, full.names=T) %>% unlist %>% as.vector##  List all the samples directories  within each id sample.(Replicas directories) 
    sam_tab = data.frame( (strsplit(gsub(w_dir, '', d), '/') %>% unlist %>% matrix(ncol=3, nrow=length(d), byrow=T)) [,-1] , d , stringsAsFactors = F)
    names(sam_tab) = c('ix', 's_name', 'dir')
    ## Restrict dataset to samples with qluster data. 
    sam_tab = sam_tab[sapply(sam_tab$dir, function(i) file.exists(paste(i, 'outputInfo.tab.txt', sep='/'))),]
    return(sam_tab)
  }
  
  ## Function to get data froma  given replica. 
  #### replica ->  full path to replica directory.
  #### check_primers -> useful to compare seq errors using as benchmark the primer regions. default is F
  #### off_indels -> remove clusters/haplotypes with seq length diff to the median length. 
  ## Returns a list with two dataframes, freqs and seqs. 
  process_qlusters.replica = function(replica, check_primers = F, off_indels = T){ 
    ## Load data: Sequences + tab 
    tmp_freqs = data.frame(s_name= gsub('.+\\/(\\w+)$', '\\1', replica), read.table(paste(replica, 'outputInfo.tab.txt', sep='/'), header=T) )[,1:4]
    tmp_seqs = readFastq(paste(replica, 'output.fastq', sep='/')) %>% sread
    
    ##  Correct s_name and shorter Cluster id. 
    tmp_freqs$s_name = gsub(pattern = '\\S+\\/', rep = '', x = as.character(tmp_freqs$s_name[1]), perl=TRUE)
    #tmp_freqs$run = gsub('\\d+', '', tmp_freqs$s_name[1], perl=T)
    tmp_freqs$ClusterId = gsub(pattern = '_t\\d+', rep = '', x = tmp_freqs$ClusterId, perl=TRUE)
    ## Check it contain the name of the run id, rather than other weird name. 
    tmp_freqs$ClusterId = gsub(pattern = '\\w+\\.', rep = paste(tmp_freqs$s_name[1], '.', sep=''), x = tmp_freqs$ClusterId, perl=TRUE)
    names(tmp_seqs) = tmp_freqs$ClusterId ## Correspondence btw freqs and seqs. 
    
    ln = width(tmp_seqs)[1:10] %>% median  ##  Expected length. 
    
    ## Filter by sequence length.  NO INDELS # Assuming both gB & gL are esential, no indels should be there. 
    if(sum(width(tmp_seqs) != ln) > 0 & off_indels == T){
      tmp_ix = which(width(tmp_seqs) != ln ) #& tmp_freqs$ClusterSize < 100
      tmp_freqs = tmp_freqs[- tmp_ix, ] 
      tmp_seqs = tmp_seqs[- tmp_ix ]
    }
    
    if(check_primers == T & off_indels == T){
      ##  Primers regions have to be equal by construction , otherwise are errors during PCR / Sequencing. 
      ##  Retrieve primer regions and replace by canonical primer seq.
      cm_primers = (paste(subseq(tmp_seqs,1,15), subseq(tmp_seqs, ln-14, ln), sep='') %>% DNAStringSet %>% consensusMatrix)[1:4,]
      primers = sapply(1:ncol(cm_primers), function(i) which(cm_primers[,i] == max( cm_primers[,i] ) )) %>% names %>% paste(collapse='') %>% DNAString
      tmp_seqs = paste(subseq(primers,1,15), subseq(tmp_seqs, 16, ln-15), subseq(primers,16,30),sep='') %>% DNAStringSet
    }
    
    ## Eliminate redundant sequences : (Sometimes from SK qluster results, also because of primer errors.) 
    tmp2 = factor(as.character(tmp_seqs))
    if(length(levels(tmp2)) < length(tmp_seqs)){
      tmp_freqs$seq = tmp2
      q_size = tapply(tmp_freqs$ClusterSize, tmp_freqs$seq, sum)  ##  New cluster size
      tmp_freqs = tmp_freqs[tapply(1:nrow(tmp_freqs), tmp_freqs$seq, function(i) i[1]),]  ## Reduce set
      tmp_freqs$ClusterSize = q_size
      tmp_seqs = as.character(tmp_freqs$seq) %>% DNAStringSet  ## Modify seqs set. 
      names(tmp_seqs) = tmp_freqs$ClusterId
      tmp_freqs = tmp_freqs[,-which(names(tmp_freqs) == 'seq')]
    }
    return(list(freqs=tmp_freqs, seqs=tmp_seqs))
  }
  
  ## Function to plot correlation by pairs of replicas of a sample.  
  #### tmp -> sample matrix of (m,2) with haplotype frequency, m = each unique haplotype; n = replica, only 2. 
  #### file_graph -> file to store the graph.
  ## Returns stats info about correlation / slope between pairs of replicas.
  process_qlusters.graphs = function(tmp, file_graph = NULL){
    if(ncol(tmp) == 2 & nrow(tmp) > 0){
      reads_rep = colSums(tmp)
      tmp2 = apply(tmp, 2, function(i) i/sum(i))
      tmp_log = log10(tmp2)
      tmp_stats = list()
      ##  Calculate cor & lm 
      tmp_stats$cor = cor(tmp_log[,1], tmp_log[,2])
      tmp_stats$lm = lm(tmp_log[,1] ~ 0 + tmp_log[,2])
      slope = tmp_stats$lm$coefficients
      if(abs(slope - 1) > abs((1/slope) - 1)){  ## Check for best slope. 
        tmp_stats$lm = lm(tmp_log[,2] ~ 0 + tmp_log[,1])
        tmp2 = tmp2[,2:1];  tmp_log = tmp_log[,2:1]
      }
      tmp_stats$order = colnames(tmp2)
      ## Create graph.
      g_lab = c(tmp_stats$lm$coefficients, tmp_stats$cor) %>% round(digits = 3)
      g_lab = paste(c('m', 'cor'), g_lab, sep = ' = ')
      coff_ix = which( apply(tmp, 1, sum) >= (cut_off * sum(reads_rep))  )
      hap_coff = sum(tmp[coff_ix,]) / sum(reads_rep)
      if(length(coff_ix) > 1) hap_coff = rowSums(tmp[coff_ix,]) / sum(reads_rep)
      
      if(!is.null(file_graph)) pdf(file_graph, width = 4, height = 8)
      par(mfrow=c(2,1), mai = c(1, 1, 0.5, 0.5))
      plot(c(-5, 0), c(-5, 0), col='red', type='l', xlab = colnames(tmp)[1], ylab = colnames(tmp)[2], main = paste(names(ix)[s], gene), xlim=c(-5, 0), ylim=c(-5, 0)
      );  abline(h = log10(cut_off), col = 'orange'
      ); abline(v = log10(cut_off), col = 'orange'
      ); text(x = c(-4, -4), y = c(-0.5, -1.5), labels = g_lab, col = 'purple'
      ); points(tmp_log[,1], tmp_log[,2], pch = 16, cex = 0.7
      ); points(log10(hap_coff), log10(hap_coff), pch = 16, cex = 0.55, col = 'red')     
      
      plot(c(0,1), c(0,1), col='red', type='l', xlab = paste(colnames(tmp)[1], reads_rep[1], sep='  R= '), ylab = paste(colnames(tmp)[2], reads_rep[2], sep='  R= '), xlim=c(0,1), ylim=c(0,1)
      );  points((tmp2[,1]), (tmp2[,2]), pch = 16, cex = 0.6
      ); text(x = 0.15, y = 0.9, labels = paste('# Hap =', nrow(tmp)))
      if(!is.null(file_graph)) dev.off()
      
      return(tmp_stats)  ## Return stats  data. 
    }else return(NULL)
  }
  
  ## Function to find potential quimeric haplotypes sequences. 
  #### h_seq -> DNAStringSet with sequences to investigate -> SAME LENGTH.
  #### refs -> Names of the sequences in h_seq to use as references. 
  ## Returns potential references pairs that can explain haplotype as chimeric instead of de novo mutations. 
  find_chimeras = function(h_seq, refs){
    d = as.matrix(ape::dist.dna(as.DNAbin(h_seq), 'raw'))*width(h_seq)[1]
    #print(d[refs, refs])
    pairs = data.frame()
    novo = dif_ref = dif_alt = list()
    ##  Found possible recombinants per each of the sequences used for the query. 
    for(h in setdiff(rownames(d),refs)){
      ix = order(d[h, refs])
      ## Distance between haplotype h to the other parental strains = distance of parentakl strain of h  to other parental seqs + distance of h to its parental.
      #d[h, refs[ix[-1]]] == (d[refs[ix[1]], refs[ix[-1]]] + d[h, refs[ix[1]]])   c = a + b  
      ## If equal: No recombination.  # If different + Recombination or new mutation.
      a = round( d[h, refs[ix[1]]] )  # distance from hap h to parental reference.
      b = round( d[refs[ix[1]], refs[ix[-1]]] )  # distance from parental ref to other non parental references (alt).
      c = round( d[h, refs[ix[-1]]] ) # distance from hap h to other non parental references.
      ix2 = ix[-1][c != (a+b)]  ## Number of other references non parental where the distance doesn't follow c = a + b
      if (length(ix2)>0){
        ##Either solving the system of equations or obtaining the variable sites, the results will be the same, so let's use sets XD
        # z = (b+c-a)/2; x = b-z; y = c-z
        ## Applying sets.
        for (alt in refs[ix2]){
          temp = lapply(h_seq[c(h, refs[ix[1]], alt)], function(i) consensusMatrix(DNAMultipleAlignment(i))[1:4,])
          #which(apply(temp[[1]]+temp[[2]]+temp[[3]], 2, function(i) sum(i>0)) > 1) # total var sites
          a = which(apply(temp[[1]]+temp[[2]], 2, function(i) sum(i>0)) > 1) # var sites between h and ref
          c = which(apply(temp[[1]]+temp[[3]], 2, function(i) sum(i>0)) > 1) # var sites between h and alt
          a2 = sapply(a, function(i) paste(rownames(temp[[1]])[temp[[1]][,i]==1], i, rownames(temp[[2]])[temp[[2]][,i]==1], sep=':') )  ## Adding the nt
          c2 = sapply(c, function(i) paste(rownames(temp[[1]])[temp[[1]][,i]==1], i, rownames(temp[[3]])[temp[[3]][,i]==1], sep=':') )  ## Adding the nt
          y2 = intersect(a2,c2) # de novo mutations
          y = as.numeric(gsub('A|T|C|G|:','',y2))
          x = sort(setdiff(a, y))
          z = sort(setdiff(c, y))
          or = order(c(min(x), max(x),min(z), max(z)-1))
          #or = paste(or2, collapse='')  ## minus 1 in max(z) to account for the same position apparently mutated.
          if(paste(or, collapse='') %in% c('1234','3412')){
            novo[[paste(h, refs[ix[1]], alt, sep=':')]] = y2
            dif_ref[[paste(h, refs[ix[1]], alt, sep=':')]] = setdiff(a2, y2)
            dif_alt[[paste(h, refs[ix[1]], alt, sep=':')]] = setdiff(c2, y2)
            bk = mean( c(min(x), max(x),min(z), max(z))[or[2:3]] )
            pairs = rbind(pairs, data.frame(hap=h, ref=refs[ix[1]], d_ref=length(a), alt=alt, d_alt=length(c), d_novo= length(y), dif_ref=paste(min(x), max(x), sep='-'), dif_alt=paste(min(z), max(z), sep='-'), bk=bk))
          }
        }
      }
    }
    if(nrow(pairs) == 0){
      return(NA)
    }else{
      pairs$c_ref = dif_ref
      pairs$c_alt = dif_alt
      pairs$novo = novo
      ##  Filtering -> taking the entries where the minimum number of new mutations is needed to explain the haplotype.
      ix = unlist(tapply(1:nrow(pairs), pairs$hap, function(i) i[pairs$d_novo[i] == min(pairs$d_novo[i])]))
      pairs = pairs[ix,]
      pairs = pairs[ sort(tapply(1:nrow(pairs), as.character(pairs$hap), function(i) i[1])), ] ## Take just one possible recombination event per haplotype.
      return(pairs)
    }
  }
  
  ################
  ##  1. Retrieve data & setting directories. 
  ################
  if(!dir.exists(file.path(w_dir, subset_dir))) stop('Haplotypes directory path needs to be provided or is incorrect')
  sam_tab = process_qlusters.sam_table(file.path(w_dir, subset_dir))  ##  Retrieve directories per samples. 
  if(nrow(sam_tab) ==0 ) stop('SeekDeep qluster output files are not found for any sample.')
  ## Take out replicates known to be bad or to be excluded. 
  if(!is.null(s_names)) sam_tab = subset(sam_tab, ! s_name %in% s_names ) 
  ##  Index for main loop. 
  ix = table(sam_tab$ix) %>% sort(decreasing=TRUE)  ## Number of replicas per sample. 
  sapply(file.path(w_dir, paste(c(summ_dir, file.path(summ_dir, gene), file.path(summ_dir,gene,'pq'), file.path(summ_dir,gene,'graphs')), sep='_')), function(i) if (!dir.exists(i)) dir.create(i)) # Create summary directories. 
  
  ################
  ##  2. Main for loop per sample. 
  ################
  replica_sum = data.frame() ##  Summary about # of reads, # of haplotypes & strain frequency per replica. 
  if(chim == T)  chimeras = data.frame()
  nice_sum = all_pairs_sum = data.frame()
  
  ##  Reconstructed FOR ->  Plots and processing. 
  for(s in 1:length(ix)){ 
    print(names(ix)[s])
    ########
    ##  A.  LOAD DATA from  from each replica for the sample, filter by length & calculate current 
    ########
    freqs = data.frame(); seqs = DNAStringSet() 
    sam_replicas = sam_tab$dir[sam_tab$ix==names(ix)[s]]
    for(replica in sam_replicas){
      tmp = process_qlusters.replica(replica)
      tmp_freqs = tmp$freqs
      
      ## Adding replica info
      replica_sum = rbind(replica_sum, data.frame(id=names(ix)[s], 
                                                  s_name=tmp_freqs$s_name[1], 
                                                  #run=tmp_freqs$run[1],  
                                                  total_rds = sum(tmp_freqs$ClusterSize), 
                                                  tot_hap=nrow(tmp_freqs), 
                                                  filter_rds=NA, filter_hap=NA, stringsAsFactors = F)  )
      
      freqs = rbind(freqs, tmp_freqs)  ## Add frequencies data. 
      seqs = c(seqs, tmp$seqs)  ##Concatenated sequences. 
    }
    ##  Correct replica name.
    freqs$s_name = factor(freqs$s_name, levels=unique(freqs$s_name))
    print(tapply(freqs$ClusterSize, freqs$s_name, sum))
    
    ########
    ## B. SELECTING BEST PAIRS OF REPLICAS by shared haplotypes freq - plots and stats.  
    #### Haplotypes have to be present in at least 2 replicas and not as singletons. 
    ########
    ##  Taking out singletons. 
    seqs = seqs[freqs$ClusterSize > 1]
    freqs = freqs[freqs$ClusterSize > 1, ]
    ##  Having the set of shared haplotypes. 
    rep_num = unique(freqs$s_name) %>% length
    freqs$hap = factor(as.character(seqs)) %>% as.numeric
    shared_hap = acast(subset(freqs, ClusterSize > 1)[,c('s_name', 'hap', 'ClusterSize')], hap ~ s_name, value.var = 'ClusterSize', fill = 0)
    ##  Take out for starters haplotypes that are unique to certain replicate.   rep_num > 1
    if(rep_num > 1){
      aux = rownames(shared_hap)[ apply(shared_hap, 1, function(i) sum(i == 0) <= (ncol(shared_hap) - 2)) ]  ##  Haplotypes present in at least 2 samples. 
      shared_hap = matrix(shared_hap[aux,], ncol = ncol(shared_hap), dimnames = list(aux, colnames(shared_hap)))
      if(length(aux) > 1){
        shared_hap = shared_hap[, apply(shared_hap, 2, sum) > 0] ##  Avoid weird replicas with no overlapping info. 
        shared_hap = shared_hap[order(rowSums(shared_hap), decreasing = T),]  ### Let's no weight the samples nor take the total # of reads. ** Changes between possible pairs of replicates. 
      }
    }  
    ##  Choose best pair of replicates and modify datasets:  shared_hap, freqs & seqs  
    #best_replicas = 1:2; best_cor = 0; best_lm =NA; best_ix = NA
    if(rep_num > 1){
      ########
      ## B.2 CONT.  SELECTING BEST PAIRS OF REPLICAS - plots and stats. 
      ########
      ## Choose the best pair of replicas.
      ##  Best criteria: Best correlation rather than best slope. 
      best_replicas = rep(0, 2); best_cor = 0; best_lm =NA; best_ix = NA
      
      for(i in 1:(ncol(shared_hap) - 1)){
        for(j in (i+1):ncol(shared_hap)){
          shared_ix = which(shared_hap[,i] > 0 & shared_hap[,j] > 0)
          if(length(shared_ix) >= 2){
            shared_tmp = shared_hap[shared_ix, c(i,j)]
            file_graph = file.path(w_dir, summ_dir, gene, 'graphs', paste(names(ix)[s], gene, colnames(shared_hap)[i], colnames(shared_hap)[j], sep='_'))
            pair_stats = process_qlusters.graphs(shared_tmp, gsub('$', '.pdf', file_graph))
            
            ##  Add this info to the all_pairs_sum dataframe. 
            all_pairs_sum = rbind(all_pairs_sum, data.frame(
              id = names(ix)[s], R1 = pair_stats$order[1], R2 = pair_stats$order[2], slope = pair_stats$lm$coefficients, 
              cor = pair_stats$cor   
            ))
            
            ##  Checking if the correlation is better than previusly observed.
            if(pair_stats$cor > best_cor){
              ##  Not sure why I had this expression. -> abs(pair_stats$cor - 1) < abs(best_cor - 1)
              best_replicas = c(i, j)
              best_cor = pair_stats$cor
              best_lm = pair_stats$lm
              best_ix = shared_ix
            }
          }else{
            all_pairs_sum = rbind(all_pairs_sum, data.frame(
              id = names(ix)[s], R1 = colnames(shared_hap)[i], R2 = colnames(shared_hap)[j], slope = NA, cor = NA  
            ))
            if(length(shared_ix) == 1) best_ix = shared_ix
          } ## end of else for 1 or none shared haplotypes. 
        }
      }## End of for pair of replicas
    }else{ ##  Case of a single sample.
      best_replicas = 1; best_ix = 1:nrow(shared_hap)
    }
    ########
    ## B.3 Gathering data from best pair of replicas. 
    ######## 
    ##  If the number of shared haplotypes is at least 1.
    if( length(best_ix) > 1 ){
      ##  Gather best table and stats if more than 1 replicate. 
      shared_hap = matrix(shared_hap[best_ix,best_replicas], nrow = length(best_ix), dimnames = list(rownames(shared_hap)[best_ix], colnames(shared_hap)[best_replicas]))
      rds_rep = colSums(shared_hap) 
      shared_hap = cbind(shared_hap, total = rowSums(shared_hap))
      shared_hap = apply(shared_hap, 2, function(i) i/sum(i))

      ##  Add best slope and confidence interval. 
      ## Slope of the regression.  Try both ways.
      if(class(best_lm) == 'lm'){
        lm_CI = confint(best_lm, names(best_lm$coefficients), level = 0.95)
        sam_slope = best_lm$coefficients
      } else{
        lm_CI = rep(NA, 2); sam_slope = NA
        }
      
      nice_sum = rbind( nice_sum, 
                          data.frame(id = names(ix)[s], R1 = colnames(shared_hap)[1], R2 = colnames(shared_hap)[2], 
                                     R1_rds = rds_rep[1], R2_rds = rds_rep[2], 
                                     hap_num = length(best_ix), sample_freq_cutoff = sum(shared_hap[shared_hap[,'total'] >=cut_off,'total']),
                                     cor = best_cor,  slope = sam_slope, 
                                     CI_min = lm_CI[1], CI_max = lm_CI[2], CI_s1 = lm_CI[1] <= 1 & 1 <= lm_CI[2]   )
      )
      
      ##########
      ## C. SHARED HAPLOTYPES. & COLLAPSING OF NON SHARED HAPS. - CUT_OFF
      ########
      ##  Restrict dataset to only samples that are part of the best replicas. 
      seqs = seqs[freqs$s_name %in% colnames(shared_hap)]  
      freqs = freqs[freqs$s_name %in% colnames(shared_hap), ]
      
      ###  Here is where the cut-off should be applied. - Using shared_hap table. 
      ##  If only those samples that will do to the cut-off in the merge. 
      hap_ix = rownames(shared_hap)[shared_hap[,'total'] >= cut_off] %>% as.numeric
      shared_haps = seqs[sapply(hap_ix, function(i) which(freqs$hap == i)[1])]  ##  Keep only shared haplotypes.
      
      if(ncol(shared_hap)> 2 & length(hap_ix) > 1){
        file_graph = file.path(w_dir, summ_dir, gene, 'graphs', paste(names(ix)[s],"_", gene, '_final.pdf' , sep=''))
        process_qlusters.graphs(shared_hap[as.character(hap_ix),-3], file_graph = file_graph)
      }
      
      ## Order haplotypes per number of reads.
      rds_hap = sapply(1:length(shared_haps), function(i) sum( freqs$ClusterSize[seqs == shared_haps[i]] ))
      shared_haps = shared_haps[order(rds_hap, decreasing = T)]
       
      freqs$hap = NA  #  Add label for shared haplotypes 
      for(i in 1:length(shared_haps)){
        freqs$hap[ seqs == shared_haps[i] ] = i
        #tmp = which(seqs == shared_haps[i])[1]
      }
      ln = width(shared_haps[1])  # Store segment length. 
      non_shared_ix = is.na(freqs$hap)# Non-shared or singletons. 
      
      ## Assign clusters out of the shared set to their closest haplotype in the shared set.  Expensive Step. 
      if(collapse_non_shared == T & sum(rds_rep) < 20000){
        non_shared_ix = which(non_shared_ix)
        for(i in non_shared_ix){
          tmp_seqs = c(seqs[i], shared_haps) %>% subseq(17, ln-15) %>% as.DNAbin 
          t_dist = (tmp_seqs  %>% dist.dna('raw') %>% as.matrix)[-1,1]
          freqs$hap[i] = which(t_dist == min(t_dist))[1]
        }
      }else{
        seqs = seqs[!non_shared_ix]
        freqs = freqs[!non_shared_ix,]
      }
      
      ## Adding data of the reads & haplotypes that passed the filter step per replica.
      replica_sum$filter_rds[replica_sum$id == names(ix)[s]] = tapply(freqs$ClusterSize, freqs$s_name, sum)
      replica_sum$filter_hap[replica_sum$id == names(ix)[s]] = tapply(freqs$hap, freqs$s_name, function(i) length(unique(i)))
      
      ##Dataframe containing all the summary per sample.
      data = data.frame(
        ix = names(ix)[s], #pop_freq = 1,
        rep_num = tapply(freqs$s_name, freqs$hap, function(i) length(unique(i))), rep_frac = 1, 
        read_num = tapply(freqs$ClusterSize, freqs$hap, sum),
        #read_avg_frac = tapply(freqs$clusterFraction, freqs$hap, mean),
        read_adj_frac = 1,  seq = as.character(shared_haps)        )
      
      ############################
      ## D.  Chimeras
      ############################
      
      if(chim == T & nrow(data) >= 3){
        ln = width(seqs)[1]
        sq = DNAStringSet(data$seq) %>% subseq( 15, ln-14); names(sq) = as.character(1:length(sq))
        
        ##Method:  Sort-rel-fre- decreasing / Start comp fparents 1,2 to 3. then go under, for n, parents are all possible pairs i,j < n
        pairs = data.frame()
        
        for(i in 3:nrow(data)){
          ## Potential parents for i.  1:i-1. 
          tmp = find_chimeras(sq[1:i], names(sq)[1:(i-1)])
          if(class(tmp) == "data.frame") pairs = rbind(pairs, tmp)
        }
        
        if(nrow(pairs) > 0 & sum(pairs$d_novo == 0) > 0 ){
          ##  Modify the number of reads. 
          chim_ix = subset(pairs, d_novo == 0 )$hap %>% as.character %>% unique %>% as.numeric
          
          for(h in chim_ix){
            parental = pairs[pairs$hap == h,c('ref', 'alt')] %>% as.numeric
            data$read_num[parental] = data$read_num[parental] + data$read_num[h]/2 %>% round(digits=0)
          }
          ## Restrict the dataset.
          data = data[-chim_ix,]
          # Write to the chimera file. 
          
          chimeras = rbind(chimeras, data.frame(id=names(ix)[s], pairs[pairs$hap %in% chim_ix, 1:9 ])  )
        }
        
      }
      
      ############################
      ## E.  Store data
      ############################
      
      data = data[order(data$read_num, decreasing=T),]
      data$read_adj_frac = data$read_num / sum(data$read_num)
      pq_file = file.path(w_dir, summ_dir, gene, names(ix)[s]) %>% gsub(pattern = '$', replacement = '.csv')
      write.csv(data, pq_file, row.names=F)
    }else{
      ## Worst case scenario... Negative correlation in the only pair or even best pair. 
      message(paste(names(ix)[s], 'sample had incongruent replicas (Correlation below 0), ommited from summary, graph is available'))
    } 
    }
  
  #########
  ###  3. Plots for the  summary. 
  #########
  file_graph = file.path(w_dir, summ_dir, gene, 'graphs')
  
  #ggplot(nice_sum, aes(x = slope, y = cor)) + geom_vline(xintercept = 1, color = 'orange') + geom_point(alpha = 0.7) + theme_bw() + geom_segment(aes(x = CI_min, xend = CI_max, y = cor, yend = cor, color = CI_s1, alpha = 0.5)) + scale_color_manual(values = c('FALSE' = 'red', 'TRUE' = 'black')) + theme(legend.position = 'none') + labs(x = 'Slope - CI', 'Correlation')
  #ggplot(nice_sum, aes(x = slope)) + geom_density(fill = 'skyblue') + geom_point(y = 0, alpha = 0.3) + theme_bw() 
  #ggplot(nice_sum, aes(x = cor)) + geom_density(fill = 'salmon') + geom_point(y = 0, alpha = 0.3) + theme_bw() 
  
  nice_sum$slope2 = round(nice_sum$slope, digits = 3)
  A = ggplot(subset(nice_sum)) + geom_vline(xintercept = 1, color = 'orange') + theme_bw()+
    geom_segment(aes(x = CI_min, xend = CI_max, y = id, yend = id, color = CI_s1)) +
    #geom_segment(aes(x = CI_min, xend = CI_min, y = as.numeric(id) - 0.2, yend = as.numeric(id) + 0.2, color = CI_s1)) +
    #geom_segment(aes(x = CI_max, xend = CI_max, y = as.numeric(id) - 0.2, yend = as.numeric(id) + 0.2, color = CI_s1)) +
    geom_point(aes(x = slope, y = id, color = CI_s1), size = 0.5) + theme(legend.position = 'right') +
    geom_tile(aes(x = -0.5, y = id, fill = cor), color = 'blue') + 
    geom_point(aes(x = -0.5, y = id, color = cor >= 0.7)) +
  #geom_text(aes(x = 0.4, y = id, label = slope2)) +
    scale_fill_gradient(low = 'white', high = 'blue', limits = c(0,1)) + 
    scale_x_continuous(breaks = c(-1.5, 0, 1, 2), labels = c('Cor', '0', '1', '2')) +
    scale_color_manual(values = c('FALSE' = 'red', 'TRUE' = 'black')) + 
    labs(x = 'Slope, 95% CI', y = 'Sample', fill = 'Correlation\n R1 ~ R2', color = '  1 in\nslope CI\n or Cor > 0.7', title = paste('Cluster :', subset_dir))
  
  pdf(paste(file_graph, '/replicas_test_cor_slope_all_reads2.pdf',sep=''), width = 10, height = 6)
  grid.arrange( A %+% nice_sum[1:50,] + theme(legend.position = 'none') , 
                A %+% nice_sum[51:100,] + theme(legend.position = 'none') + labs(title = '', y = NULL), 
                A %+% nice_sum[101:nrow(nice_sum),] + labs(title = '', y = NULL), nrow = 1, widths = c(2,2,3))
  dev.off()
    
  pdf(paste(file_graph,'/replicas_test_cor_slope_all_reads.pdf',sep=''), width = 5, height = 8)
  ggplot(subset(nice_sum)) + geom_vline(xintercept = 1, color = 'orange') + theme_bw()+
    geom_segment(aes(x = CI_min, xend = CI_max, y = id, yend = id, color = CI_s1)) +
    geom_segment(aes(x = CI_min, xend = CI_min, y = as.numeric(id) - 0.2, yend = as.numeric(id) + 0.2, color = CI_s1)) +
    geom_segment(aes(x = CI_max, xend = CI_max, y = as.numeric(id) - 0.2, yend = as.numeric(id) + 0.2, color = CI_s1)) +
    geom_point(aes(x = slope, y = id, color = CI_s1)) + theme(legend.position = 'right') +
    geom_tile(aes(x = -1.5, y = id, fill = cor), color = 'blue') + 
    geom_point(aes(x = -1.5, y = id, color = cor >= 0.7)) +
    scale_fill_gradient(low = 'white', high = 'blue', limits = c(0,1)) + 
    scale_x_continuous(breaks = c(-1.5, 0, 1, 2), labels = c('Cor', '0', '1', '2')) +
    scale_color_manual(values = c('FALSE' = 'red', 'TRUE' = 'black')) + 
    labs(x = 'Slope, 95% CI', y = 'Sample', fill = 'Correlation\n R1 ~ R2', color = '  1 in\nslope CI\n or Cor > 0.7', title = paste('Cluster :', subset_dir))
  dev.off()
  
  ## Num of haplotypes. 
  pdf(paste(file_graph,'/replicas_test_cor_hap_num.pdf',sep=''), width = 5, height = 5)
  ggplot(nice_sum, aes(x = cor, y = hap_num)) + geom_point() + 
    labs(x = 'Correlation', y= 'Number of haplotypes') + scale_y_log10()
  dev.off()
  
  pdf( paste(file_graph,'/replicas_test_slope_correlation.pdf',sep=''), width = 5, height = 5)
  
  ggplot(nice_sum, aes(x = cor)) + geom_density(color='skyblue', fill='skyblue', alpha = 0.4) + 
    geom_point(y=0, alpha=0.1) + labs(title='Replicas correlation') + theme_bw() +
    geom_text(aes(label = id), y = 2, angle = 90, size = 1.5) + labs(x= 'Cor(R1, R2)')
  ggplot(nice_sum, aes(x = slope)) + geom_density(color='salmon', fill='salmon', alpha = 0.4) + 
    geom_point(y=0, alpha=0.1) + labs(title='Regression slope') + theme_bw() +
    geom_text(aes(label = id), y = 1.5, angle = 90, size = 1.5) + labs(x= 'Slope')
  
  ggplot(nice_sum, aes(x = slope, y = cor)) + geom_point(alpha=0.3, color = "red") + labs(title='Regression slope vs correlation') + theme_bw() + theme(legend.position = 'none') +
  geom_text(aes(label = id, size = 1 - cor), angle = 45) + geom_point(alpha=0.3, color = "red") + scale_x_continuous(limits = c(0.2, 1.1))  + scale_y_continuous(limits = c(0.2, 1)) 
  
  dev.off()
  
  ##############
  ## 4. Writing summary tables
  ############
  write.csv(unique(nice_sum), file.path(w_dir, summ_dir, gene, 'best_replicas_summary.csv'), row.names=F)
  write.csv(unique(all_pairs_sum), file.path(w_dir, summ_dir, gene, 'all_pairs_replicas_summary.csv'), row.names=F)
  
  write.csv(unique(replica_sum), file.path(w_dir, summ_dir, gene, '/replicas_summary.csv'), row.names=F)
  if(nrow(chimeras)>0) write.csv(chimeras, file.path(w_dir, summ_dir, gene, '/chimeras.csv'), row.names=F)
  
  weird = setdiff( as.character(unique(replica_sum$id)), as.character(unique(nice_sum$id)) )
  if(length(weird) > 0) write.table(weird, file.path(w_dir, summ_dir, gene, '/incongruent_samples.csv'), row.names=F, col.names = F)
  
}

