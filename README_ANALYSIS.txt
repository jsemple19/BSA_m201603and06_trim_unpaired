
To process sequencing data located on the server:

1) Setup:
=========
create a directory in user directory with data of sequencing run e.g 20151202
this will be the base directory for analysis (BASEDIR)
in BASEDIR create a directory called "scripts" and copy scripts into this directory
in BASEDIR create a directory called "genomeSeq" and put in it a copy of the genome file, 
and SNP bed file.

2) Aligning reads to genome with bwa mem:
=========================================
edit bwa_index.sh to have the correct FULL path for BASEDIR (working directory you just 
created) and DATADIR (location on server where fastq.gz files are) (assuming paired end
reads)
edit submitBWAmem_longQ.sh to contain the right number of tasks
e.g. #$ -t 1-8  will work for 8 sequencing samples (each with read1 and read2 files, i.e.
total of 16 fastq files)

submit the job array to the cluster:
qsub submitBWAmem_longQ.sh

This will create sam files in ./samFiles directory

3) Sorting reads with samtools and extracting variant calls with mpileup and bcftools:
======================================================================================
edit script samtools_aln.sh to have the right BASEDIR
edit submitSamtools_longQ.sh to contain the right number of tasks (number of samples 
sequenced)
submit job array to cluster:
qsub submitSamtools_longQ.sh

this will create bam files in bamFiles/ and sort them (.sorted.bam) and the vcf files in 
vcfFiles/ directory

4) Extracting allele counts etc from vcfFiles
=============================================
edit submitVarCollect_longQ.sh to have the correct full path for BASEDIR
edit to contain the right number of tasks (number of samples sequenced)

submit job array to cluster
qsub submitVarCollect_longQ.sh

this will call the custom python script CollectAllSNPs.py
that extracts the read depth (DP) and allele counts (I16) from the info column of the vcf
files, calculates Hawaii allele frequency and writes this data to .txt file in ./finalData

5) Cleaning the data and plotting the results
=============================================
Manually create a file called fileList.txt that contains two columns:
First column contains experiment filenames from finalData/ directory. arranged such that
each two lines represents an experiment (cont then selected). the second column contains 
the name of the experiment to use for filenames and figure titles (e.g. include population 
name and generation number). save this file to finalData/ directory. (suggested format is 
separating info with _ so that it can be used in unix-friendly filenames. for figure 
titles, the _ is replaced with a space.

In R, set the working directory to the scripts directory and run processCounts1.R that calls
functions found in MegaMatePlot.R
These use relative paths, so no need to set BASEDIR.
These will do some basic cleaning of the data (removing loci with <2 total reads, and loci
that whose read depths is >5* median average deviations from the median (extreme outliers
that are likely to be repetitive regions) (this cleaning is very minimal, more stringent
filters can be applied later if necessary).
The number of reads removed by the cleaning is reported in a log_experimentName.txt file.

A pdf file (data_experimentName.pdf) of plots is produced for each experiment with three plots:
a) raw frequencies of Hawaii alleles in control and selected populations
b) frequency difference between control and selected popultaions smoothed with 1000 locus
window
c) significance of difference (-log10P). this is calculated by performing a fisher test at
each locus with the N2 and Hawaii allele counts in control and selected samples. The 
pvalues are adjusted for multiple testing (fdr), and -log10(adjP) is smoothed polynomically 
with a 1000 locus window. (this has the result of reducing the significance of the most 
significant loci, and so is conservative). the fdr threshold is plotted at -log10(0.05)

A pdf file (date_readDepths.pdf) is also produced with histograms of the read depths at different loci for each
sample.

all files are saved in finalData/

When analysis is complete you can delete all directories except for 
scripts/ (so you can recreate analysis) 
finalData/ 
genomeData/ (so you know what you aligned it against)

i.e. you can delete:
/samFiles
/bamFiles
/vcfFiles
/tempData

Even if you wish to be conservatively cautious, you should definitely delete:
samFiles (very big and can be generated from bamFiles) 
and unsorted bamFiles (redundant with sorted bamFiles)
you can do this automatically with use cleanUp.sh
