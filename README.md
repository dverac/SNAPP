# SNAPP

SNAPP sequencing data pipeline to produce haplotypes from targeted amplicon sequencing

##Overall sequencing data processing. 

Raw sequencing data (Unpaired reads)   -->  Haplotypes per sample (Sequence and frequency)


### Considerations

	• Targeted amplicon sequencing to an specific gene region.
		○ Region can be recovered by merging paired reads by overlapping ends.  
	• Predominantly 2 technical replicates per sample. 
		○ Few samples with more than 2 replicates or only one replicate. 
	• High read depth. 
	• High number of PCR rounds =  High number of PCR errors and increased variability between replicas.
	• Variable initial sample size: Variable viral load at different sampling points in time and fluid. 
	
	
### Pipeline

#### Paired reads merging - PEAR
For each replicate, their paired reads are merged using PEAR. A software that takes into account the overlapping region as well as the quality scores on the bases located in such region. The merging was performed under default parameters.

`pear -f  sam1_rep1_R1_001.fastq.gz  -r  sam1_rep1_R2_001.fastq.gz -o  sam1_rep1`

#### Fragments quality filtering - Seek Deep Extractor
Fragments are now filter to exclude those with low quality scores and sequences whose ends doesn't match the primers sequences employed for amplification. 
The following parameters were used for all the samples. 

`seekdeep extractor -fastq  sam1_rep1_pear.fastq -id  gene.id  -checkComplement -qualCheck 25 -qualCheckCutOff .75 -variableStart 20 -minLen 350 -maxLen 450 --overWriteDir -dout  sam1_rep1`

#### Haplotypes calculation - Seek Deep Qluster
The filtered fragments from each replicate are now collapsed.  Each unique sequence is associated with their respective number of reads that belong to it. After sorting such clusters by their number of reads a collapsing iteration is performed. Low frequency clusters are compared against more abundant clusters if their sequence only diverge in 1 base and such base is located in a low quality region, the small cluster is added to the major one.  The final clusters are considered afterwards as haplotypes. 

`seekdeep qluster -fastq sam1_rep1_sd_extractor/gene.fastq  -par ./tabs/illumina_lkmer2 --overWriteDir --noHomopolymerWeighting  -dout  sam1/rep1`

#### Haplotypes consolidation - Merging replicates data for each sample

Per each sample, all its replicates are considered. Shared haplotypes between each possible pair of replicates is obtained after excluding singletons. In general, for samples with more than 2 replicates, the pair with the best correlation between the frequencies of their shared haplotypes are chosen as representatives of each sample. Using the same idea, samples in which the pair of replicates are poorly correlated were flagged since they could lead to inaccurate final frequencies of the shared haplotypes.  

Once the best pair of replicas is chosen, the final relative frequency of each haplotype is the weighted mean of their frequencies, using the read number per replicate as weight. Haplotypes with a frequency over the selected threshold, 0.436% (look haplotypes cut-off  page) are kept and those below the threshold or non shared are collapsed to their closest shared haplotype that satisfies the frequency threshold. Finally, shared haplotypes are tested for chimeric origin. An haplotype is considered a chimera if it is the product of local partial alignments to other two haplotypes present in the same sample  at higher frequencies. 
