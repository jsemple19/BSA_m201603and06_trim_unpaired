#! /bin/bash
#last modified 20160701
#mpileup under bcftools for multiallelic calling with no prior

GENOME_VER="PRJNA13758.WS250"
BASEDIR=/users/blehner/jsemple/seqResults/m201603and06_trim_unpaired
PAIREDDIR=/users/blehner/jsemple/seqResults/m201603and06_trim_slWin
GENOMEDIR=/users/blehner/jsemple/seqResults/GenomeBuilds/$GENOME_VER

genomefile=`ls $GENOMEDIR/*.fa`

# to create env variable for samtools:
export SAMTOOLS_HOME="/software/bl/el6.3/samtools-1.3.1"

mkdir -p $BASEDIR/bamFiles/$GENOME_VER
mkdir -p $BASEDIR/vcfFiles/$GENOME_VER

#index genome
#$SAMTOOLS_HOME/samtools faidx $genomefile

# create arrays of file names - need to process each read separately
samfiles_up1=(`ls $BASEDIR/samFiles/$GENOME_VER/*_aln_up1.sam`)
samfiles_up2=(`ls $BASEDIR/samFiles/$GENOME_VER/*_aln_up2.sam`)
# ordinal number of samFile input from command line arg
i=$1
let i=i-1
r1=`echo ${samfiles_up1[$i]%*.sam} | cut -f9 -d'/'`
r2=`echo ${samfiles_up2[$i]%*.sam} | cut -f9 -d'/'`
# keep read specifying extension
outbam_up1=`echo $BASEDIR/bamFiles/$GENOME_VER/$r1".bam"`
outbam_up2=`echo $BASEDIR/bamFiles/$GENOME_VER/$r2".bam"`

##convert samfiles to bamfiles (slow):
$SAMTOOLS_HOME/samtools view -b -S -o $outbam_up1 ${samfiles_up1[$i]}
$SAMTOOLS_HOME/samtools view -b -S -o $outbam_up2 ${samfiles_up2[$i]}

##sort the bam files:
outsort_up1=`echo $BASEDIR/bamFiles/$GENOME_VER/$r1".sorted.bam"`
outsort_up2=`echo $BASEDIR/bamFiles/$GENOME_VER/$r2".sorted.bam"`	
$SAMTOOLS_HOME/samtools sort $outbam_up1 -T $outsort_up1 -o $outsort_up1
$SAMTOOLS_HOME/samtools sort $outbam_up2 -T $outsort_up2 -o $outsort_up2

# merge sorted bamfiles of both unpaired and paired reads
r1=`echo ${samfiles_up1[$i]%_*_aln_up1.sam} | cut -f9 -d'/' | perl -i -pe 's/elegans_unpaired/elegans_uppe/g' `
r2=`echo ${samfiles_up2[$i]%_*_aln_up2.sam} | cut -f9 -d'/' | perl -i -pe 's/elegans_unpaired/elegans_uppe/g' `
bamfiles_pe=(`ls $PAIREDDIR/bamFiles/$GENOME_VER/*.sorted.bam`)
r3=`echo ${bamfiles_pe[$i]%*.sorted.bam} | cut -f9 -d'/' | perl -i -pe 's/elegans_paired/elegans_uppe/g' `
outmerge=`echo $BASEDIR/bamFiles/$GENOME_VER/$r1".sorted.bam"`

if [ $r1 = $r2 -a $r2 = $r3 ]
then
	$SAMTOOLS_HOME/samtools merge $outmerge $outsort_up1 $outsort_up2 ${bamfiles_pe[$i]}
else
	echo "not same file name"
fi

##count variants at SNP sites from .bed file
#bedfile=`ls $GENOMEDIR/SNVs*.pos`  # 300,000 SNVs from Hawaii genome paper
#bedfile=`ls $GENOMEDIR/AnnSNPs*.pos` # 172,000 SNVs from genome gff annotation
bedfile=`ls $GENOMEDIR/SNVs_N2-CB4856*.pos`

#names of sorted files now end in .bam
#insorted=`echo $outsort".bam"`

#create variant count file name
vcffile=`echo $BASEDIR/vcfFiles/$GENOME_VER/$r1"_raw.vcf"`

export BCFTOOLS_HOME="/software/bl/el6.3/bcftools-1.3.1"

#run in two steps. first mpileup. then index the outputed bam file. then bcftools multiallelic calling
mpfile=`echo $BASEDIR/vcfFiles/$GENOME_VER/$r1"_mp.bcf"`
$SAMTOOLS_HOME/samtools mpileup -f $genomefile -Q 20 -l $bedfile -BIg $outmerge -o $mpfile

$BCFTOOLS_HOME/bcftools index $mpfile

# multiallelic calling with prior ignored (-P0)
$BCFTOOLS_HOME/bcftools call -R $bedfile -m -P0 -Ov $mpfile -o $vcffile
