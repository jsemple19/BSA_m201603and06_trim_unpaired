#! /bin/bash
#last modified 20160701

GENOME_VER="WS230"
BASEDIR=/users/blehner/jsemple/seqResults/m201603and06

# to create env variable for samtools:

export SAMTOOLS_HOME="/software/bl/el6.3/samtools-1.3"
mkdir -p $BASEDIR/bamFiles/$GENOME_VER
mkdir -p $BASEDIR/vcfFiles/$GENOME_VER
GENOMEDIR=/users/blehner/jsemple/seqResults/GenomeBuilds/$GENOME_VER

genomefile=`ls $GENOMEDIR/*.fa`

#index genome
$SAMTOOLS_HOME/samtools faidx $genomefile

# create arrays of file names
samfiles=(`ls $BASEDIR/samFiles/$GENOME_VER/*.sam`)

# ordinal number of samFile input from command line arg
i=$1
let i=i-1
r1=`echo ${samfiles[$i]%_*_aln_pe.sam} | cut -f9 -d'/'`
outbam=`echo $BASEDIR/bamFiles/$GENOME_VER/$r1".bam"`

##convert samfiles to bamfiles (slow):
$SAMTOOLS_HOME/samtools view -b -S -o $outbam ${samfiles[$i]}

##sort the bam file:
outsort=`echo $BASEDIR/bamFiles/$GENOME_VER/$r1".sorted.bam"`	
$SAMTOOLS_HOME/samtools sort $outbam -T $outsort -o $outsort

##count variants at SNP sites from .bed file
#bedfile=`ls $GENOMEDIR/SNVs*.pos`  # 300,000 SNVs from Hawaii genome paper
bedfile=`ls $GENOMEDIR/AnnSNPs*.pos`  # 172,000 SNVs from genome gff annotation

#names of sorted files now end in .bam
#insorted=`echo $outsort".bam"`

#create variant count file
vcffile=`echo $BASEDIR/vcfFiles/$GENOME_VER/$r1"_raw.vcf"`

export BCFTOOLS_HOME="/software/bl/el6.3/bcftools-1.3"

#run samtools mpileup and bcftools
$SAMTOOLS_HOME/samtools mpileup -f $genomefile -l $bedfile -uBg $outsort | $BCFTOOLS_HOME/bcftools view -> $vcffile


