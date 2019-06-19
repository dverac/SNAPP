# SNAPP

SNAPP sequencing data pipeline to produce haplotypes from targeted amplicon sequencing

## Overall sequencing data processing. 

Raw sequencing data (Unpaired reads)   -->  Haplotypes per sample (Sequence and frequency)

### Considerations

	• Targeted amplicon sequencing to an specific gene region.
		○ Region can be recovered by merging paired reads by overlapping ends.  
	• Predominantly 2 technical replicates per sample. 
		○ Few samples with more than 2 replicates or only one replicate. 
	• High read depth. 
	• High number of PCR rounds =  High number of PCR errors and increased variability between replicas.
	• Variable initial sample size: Variable viral load at different sampling points in time and fluid. 

The code in this repository does the last step in the complete pipeline, which is haplotypes consolidation.  You can find the complete pipeline description in the wiki section. 

#### Haplotypes consolidation - Merging replicates data for each sample

Per each sample, all its replicates are considered. Shared haplotypes between each possible pair of replicates is obtained after excluding singletons. In general, for samples with more than 2 replicates, the pair with the best correlation between the frequencies of their shared haplotypes are chosen as representatives of each sample. Using the same idea, samples in which the pair of replicates are poorly correlated were flagged since they could lead to inaccurate final frequencies of the shared haplotypes.  

Once the best pair of replicas is chosen, the final relative frequency of each haplotype is the weighted mean of their frequencies, using the read number per replicate as weight. Haplotypes with a frequency over the selected threshold, 0.436% (look haplotypes cut-off  page) are kept and those below the threshold or non shared are collapsed to their closest shared haplotype that satisfies the frequency threshold. Finally, shared haplotypes are tested for chimeric origin. An haplotype is considered a chimera if it is the product of local partial alignments to other two haplotypes present in the same sample  at higher frequencies. 

Check process_qlusters_master.R code for this step. 

#### Database for a set of samples 

Given a collection of samples intended to be analyzed together, we gather all the haplotype information at a given gene for all those samples and condense them in 4 major structures:  
	• data = data.frame that includes all the haplotypes and their observed frequencies in each sample.   
	• dsum = data.frame with summary data for each sample, including total number of haplotypes, relative frequencies given a set of refrence strains and total reads supporting the sample.   
	• h_ix = data.frame for the haplotypes catalog, each unique haplotype found in the collection has a entry with information of number of samples where shouch ahpltoype was found, reference strain, number of nucleotides of difference against reference strain and statstitcis about its relative frequency in the samples.  
	• seqs = DNAStringSet with all the sequences of each unique haplotype, named using the haplotype ID and same order as h_ix.   
	  
Check hap_summary_master.R for this part (Pending...)
