##########
## Haplotypes summary -> Generalized version.
## Read the haplotypes summaries for all the samples of a given set and condense them. 
## Creates 4 important structures: total haplotypes per sample, unique haplotype catalog, haplotype sequences and samples summary. 
## Structures connect by haplotype ID. -> Each unique sequence found. 
## Author: Diana Vera Cruz
#####

####  INPUT INFO
##  pq_dir -> directory path for pq summary files. 
## out_dir -> output directory path ->  for each of the raw csv tables & env RData.
## gene -> Gene to be analyzed. (String)
## ref_strains -> file path for reference sequences. (FASTA format)

####  OUTPUT INFO
## data -> data.frame with data from all the samples, one entry per haplotype foun in a given sample and its relative frequency. 
## dsum -> data.frame with stats per sample, including info as number of replicas, number of haplotypes and frequencies of each reference strain. 
## h_ix -> data.frame for the haplotype catalog. IDs for each unique haplotype found in all the samples, and stats for such haplotype.
## seqs -> DNAStringSet with all the sequences of each unique haplotype, named with their correspendant haplotype ID and in the same order as h_ix. 

########
## MAIN FUNCTION
########
hap_summary = function(summ_dir, gene='gB', ref_strains = NULL){  ##  New format. 
  library(Biostrings); library(magrittr); library(ape)
  ############
  ## 0. FUNCTIONS
  ############
  
  ## Function to assign each haplotype to its closest reference strain. 
  #### seqs -> DNAStringSet with sequences to investigate -> SAME LENGTH.
  #### refs -> DNAStringSet with reference sequences.
  ## Returns best reference strain name for each haplotype.
  assign_strain = function(seqs = NULL, rf = NULL){
    library(ape)
    if(is.null(seqs) & is.null(refs)) stop('Incorrect input data')
    
    ln = min(width(seqs)[1], width(rf)[1])
    ix = vector(); i=0; size = 150; d_ix = vector(); dt = matrix(nrow=0, ncol=3)
    while(i < length(seqs)){
      s_dist = c(rf, seqs[(i+1):min(c(i+size, length(seqs)))]) 
      s_dist = s_dist %>% as.DNAbin %>% dist.dna(model = 'raw') %>% as.matrix
      s_dist = s_dist[,1:3] * ln ## Dist to references in number of mismatches.
      ix = c(ix, apply(s_dist, 1, function(i) which(i == min(i))[1] )[-(1:3)] )
      d_ix = c(d_ix, apply(s_dist, 1, function(i) min(i) )[-(1:3)] ) %>% round(digits=0)
      dt = rbind(dt, s_dist[-(1:3),] )
      #ix = c(ix, apply(s_dist[,1:3]*(ln-3), 1, function(i) if(min(i)<= 15) which(i == min(i)) else c(1, 1))[-(1:3)] )  ## Worst case scenario, a sample product with recombination (half an half mm from both replicas)  ## Modified step !!!
      i = i + size
    }
    
    ix = factor(names(rf)[ix], levels=names(rf))
    colnames(dt) = names(rf)
    return(data.frame(strain=ix, dist=d_ix, round(dt, digits=0)))
  }
  
  ## Function to find potential quimeric haplotypes sequences. 
  #### h_seq -> DNAStringSet with sequences to investigate -> SAME LENGTH.
  #### refs -> Names of the sequences in h_seq to use as references. 
  ## Returns potential references pairs that can explain haplotype as chimeric instead of de novo mutations. 
  find_chimeras = function(h_seq, refs){
    d = as.matrix(dist.dna(as.DNAbin(h_seq), 'raw'))*width(h_seq)[1]
    #print(d[refs, refs])
    pairs = data.frame()
    novo = dif_ref = dif_alt = list()
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
  
  #############
  ## 1. Loading data per sample. -> Dataframe data
  #############
  d = file.path(summ_dir, gene, 'pq') %>% list.files(full.names = T)
  data = data.frame() ## Load data per sample.
  for(i in d)   data = rbind(data, read.csv(i, stringsAsFactors = F) )
  data$gene = gene
  data$ix = paste(data$ix, gene, sep='-')
  data$seq = factor(data$seq)
  
  #############
  ## 2. Create a summary vector: dsum  ->  # of replicas &  total reads. 
  #############
  dsum = tapply(data$read_num, data$ix, sum)  
  dsum = data.frame(ix=names(dsum), gene=gene, read_num = unlist(dsum), 
                    rep_num = unlist(tapply(data$rep_num, data$ix, max)), 
                    hap_num = unlist(tapply(data$rep_num, data$ix, length)) )
  dsum$id = gsub(paste('-',gene,sep=""), '', dsum$ix)
  
  #############
  ## 3. Haplotypes (Catalog in h_ix dataframe) 
  #############
  seqs = levels(data$seq) %>% DNAStringSet
  names(seqs) = paste('H', 1:length(seqs), sep='')  ## Tentative name to allow integration.
  data$hap = paste('H', as.numeric(data$seq), sep='') ## adding the id. 
  data = data[, -which(names(data)=='seq')]
  
  h_ix = tapply(data$read_adj_frac, data$hap, max) ##  Creating the other summary --> Haplotypes summary. (similar to h_ix)
  h_ix = data.frame(hap=names(h_ix), max=h_ix, r_max=tapply(data$read_num, data$hap, max), min=tapply(data$read_adj_frac, data$hap, min), 
                    freq=tapply(data$read_adj_frac, data$hap, length) )
  h_ix = h_ix[order(as.numeric(gsub('H', '',h_ix$hap))),]
  
  #############
  ## 4. Loading strains and assigning them to each haplotype & filtering  
  #############
  refs = readDNAStringSet(ref_strains)
  h_ix = data.frame(h_ix, assign_strain(seqs, rf) ) ## Since sequences are not crop to the CDS region.
  cols = sapply(names(refs), grep, names(h_ix), value = T) ## Cols names for distance to each of the reference strains. 
  h_ix$range = sapply(1:nrow(h_ix), function(i) max(h_ix[i,cols]) - min(h_ix[i,cols]) )
  
  ## Add strain column to the data table. 
  data$strain = sapply(data$hap, function(i) h_ix$strain[h_ix$hap==i])
  data = data[order(data$strain, data$ix),] ## Sorting data by strains and frequency.
  
  #############
  ## 5. Contaminants & recombinants 
  #############
  h_ix$filter = NA
  ## Remaining haplotypes. ->  Recombinants not excluded by assumptions about parents existance & relative frequency.
  cols = c('UCD52', 'UCD59','X180.92')
  ideal_haps = subset(h_ix, dist ==0 )[,c('strain', cols)]
  rownames(ideal_haps) = ideal_haps$strain
  tmp = h_ix[,cols] - h_ix$dist - ideal_haps[as.character(h_ix$strain), cols]
  ###  Ideal haps are those that are containing  new mutations from the references. 
  h_ix$filter[which(apply(tmp, 1, sum) == 0 )] = 'New'
  ###  haplotypes with one mutation in the same position or to the same aminoacid that another one.
  h_ix$filter[is.na(h_ix$filter) & h_ix$dist == 1] = 'Potential new'
  ## Also take out samples that are way too far away from the reference. -> Spurious haplotypes. 
  h_ix$filter[h_ix$dist > ceiling(max(h_ix$range)/2)] = 'Spurious'
  
  ### Recombinants -> Only using the canonical sequence as parentals.
  ref_haps = as.character(h_ix$hap[h_ix$dist == 0])
  names(ref_haps) = sapply(ref_haps, function(i) h_ix$strain[h_ix$hap == i])
  
  chim = find_chimeras(seqs[ subset(h_ix, (dist > 1 & is.na(filter)) | dist == 0 )$hap ], ref_haps)
  
  pot_chim  = chim$hap[chim$d_novo == 0]
  h_ix$filter[h_ix$hap %in% pot_chim] = 'chimera'
  
  ### STORING DATA. 
  write.csv(data, paste(summ_dir,'/data_',gene,'.csv',sep=''), row.names=F)
  write.csv(h_ix, paste(summ_dir,'/h_ix_',gene,'.csv',sep=''), row.names=F)
  write.csv(dsum, paste(summ_dir,'/dsum_',gene,'.csv',sep=''), row.names=F)
  writeXStringSet(seqs, paste(summ_dir,'/seqs_',gene,'.fa',sep=''))
}
