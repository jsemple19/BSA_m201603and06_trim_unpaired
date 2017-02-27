#! /bin/bash
#last modified 20160701 to add external ref to GenomeBuilds rather than repeating 
# genomeSeq folder in each dataset and to merge fastq of technical replicates
#files

#PS1="\W>"
#changes prompt to current directory and >
GENOME_VER="WS230"

SOURCEDIR=/users/blehner/sequencing_data/Jennifer_Semple
DATADIR=/users/blehner/sequencing_data/Jennifer_Semple/m201603and06
#BASEDIR=/users/blehner/jsemple/seqResults/m201603and06
GENOMEDIR=/users/blehner/jsemple/seqResults/GenomeBuilds/$GENOME_VER

# to create env variable for bwa:
export BWA_HOME="/software/bl/el6.3/bwa_0.7.12"

#mkdir -p $BASEDIR/samFiles/$GENOME_VER
#mkdir -p $BASEDIR/tempData

Source1=(`ls $SOURCEDIR/2016-03-24/*.fastq.gz`)
Source2=(`ls $SOURCEDIR/2016-06-01/*.fastq.gz`)
#length=16


# i="0"
# while [ $i -lt $length ]; do
# 	f1=`echo ${Source1[$i]%.fastq} | cut -f7 -d'/' | cut -f1-3 -d'_'`
# 	f2=`echo ${Source2[$i]%.fastq} | cut -f7 -d'/' | cut -f1-3 -d'_'`
# 	readNum=`echo ${Source1[$i]%.fastq} | cut -f7 -d'/' | cut -f4 -d'_'`
# 	if ["$f1" != "$f2" ]; then
# 		echo "name mismatch" $f1 $f2
# 	else
# 		outfile=`echo "mrg"$f1"_"$readNum`
#  		zcat ${Source1[$i]} ${Source2[$i]} > $DATADIR/$outfile
#  	fi
#  	i=$[$i+2]
# done

#to run on cluster
i=$1
let i=i-1
f1=`echo ${Source1[$i]%.fastq} | cut -f7 -d'/' | cut -f1-3 -d'_'`
f2=`echo ${Source2[$i]%.fastq} | cut -f7 -d'/' | cut -f1-3 -d'_'`
readNum=`echo ${Source1[$i]%.fastq} | cut -f7 -d'/' | cut -f4 -d'_'`
if [ "$f1" != "$f2" ]; then
	echo "name mismatch" $f1 $f2
else
	outfile=`echo "mrg"$f1"_"$readNum`
	echo "merging" ${Source1[$i]} ${Source2[$i]} "into" $outfile
	zcat ${Source1[$i]} ${Source2[$i]} | gzip > $DATADIR/$outfile
fi

# created merged fastq files as follows:
# cat $SOURCEDIR/2015-12-24/1501_12225_CCTTCA_read1.fastq.gz $SOURCEDIR/2016-02-12/1501_12225_CCTTCA_read1.fastq.gz > $DATADIR/1501mrg_12225_CCTTCA_read1.fastq.gz
# cat $SOURCEDIR/2015-12-24/1501_12225_CCTTCA_read2.fastq.gz $SOURCEDIR/2016-02-12/1501_12225_CCTTCA_read2.fastq.gz > $DATADIR/1501mrg_12225_CCTTCA_read2.fastq.gz
# cat $SOURCEDIR/2015-12-24/1504_12224_AAGGGA_read1.fastq.gz $SOURCEDIR/2016-02-12/1504_12224_AAGGGA_read1.fastq.gz > $DATADIR/1504mrg_12224_AAGGGA_read1.fastq.gz
# cat $SOURCEDIR/2015-12-24/1504_12224_AAGGGA_read2.fastq.gz $SOURCEDIR/2016-02-12/1504_12224_AAGGGA_read2.fastq.gz > $DATADIR/1504mrg_12224_AAGGGA_read2.fastq.gz
# cat $SOURCEDIR/2015-12-24/1572_12229_CCTCGG_read1.fastq.gz $SOURCEDIR/2016-02-12/1572_12229_CCTCGG_read1.fastq.gz > $DATADIR/1572mrg_12229_CCTCGG_read1.fastq.gz
# cat $SOURCEDIR/2015-12-24/1572_12229_CCTCGG_read2.fastq.gz $SOURCEDIR/2016-02-12/1572_12229_CCTCGG_read2.fastq.gz > $DATADIR/1572mrg_12229_CCTCGG_read2.fastq.gz
# 

# then created subsampled files of same depth as new sequencing
#zcat MEmerged*read1* | head -27000000 > MEsubSampled_4998_ACAGTG_read1.fastq
#zcat MEmerged*read2* | head -27000000 > MEsubSampled_4998_ACAGTG_read2.fastq
#zcat 1Fmerged*read1* | head -27000000 > 1FsubSampled_4998_ACAGTG_read1.fastq
#zcat 1Fmerged*read2* | head -27000000 > 1FsubSampled_4998_ACAGTG_read2.fastq


#genomefile=`ls $GENOMEDIR | grep 'c_elegans.*genomic\.fa$'`

#create index (only do once):
#$gunzip $BASEDIR/genomeSeq/$GENOME_VER/$genomefile
#$BWA_HOME/bwa index -p $BASEDIR/genomeSeq/$GENOME_VER/elegans -a bwtsw $BASEDIR/genomeSeq/$GENOME_VER/$genomefile
#$BWA_HOME/bwa index -p $GENOMEDIR/elegans -a bwtsw $GENOMEDIR/$genomefile



