#! /bin/bash
#last modified 20160701 to add external ref to genomeBuilds and merge technical replicate
#files

#PS1="\W>"
#changes prompt to current directory and >
GENOME_VER="WS230"

#SOURCEDIR=/users/blehner/sequencing_data/Jennifer_Semple
DATADIR=/users/blehner/sequencing_data/Jennifer_Semple/m201603and06
BASEDIR=/users/blehner/jsemple/seqResults/m201603and06
GENOMEDIR=/users/blehner/jsemple/seqResults/GenomeBuilds/$GENOME_VER

# to create env variable for bwa:
export BWA_HOME="/software/bl/el6.3/bwa_0.7.12"

mkdir -p $BASEDIR/samFiles/$GENOME_VER
#mkdir -p $BASEDIR/tempData


#genomefile=`ls $GENOMEDIR | grep 'c_elegans.*genomic\.fa$'`

#create index (only do once):
#$gunzip $BASEDIR/genomeSeq/$GENOME_VER/$genomefile
#$BWA_HOME/bwa index -p $BASEDIR/genomeSeq/$GENOME_VER/elegans -a bwtsw $BASEDIR/genomeSeq/$GENOME_VER/$genomefile
#$BWA_HOME/bwa index -p $GENOMEDIR/elegans -a bwtsw $GENOMEDIR/$genomefile


# create arrays of file names
read1_files=(`ls $DATADIR/*read1.fastq.gz`)
read2_files=(`ls $DATADIR/*read2.fastq.gz`)

# number of filenames in a string (separated by spaces) (but chose to use arrays instead)
# length=`echo -n $read1_files | wc -w`

# #length of array
# length=`echo ${#read1_files[@]}`
# let length=length-1
# for i in `seq 0 $length`; do
# 	r1=`echo ${read1_files[$i]%_read1.fastq} | cut -f7 -d'/'`
# 	r2=`echo ${read2_files[$i]%_read2.fastq} | cut -f7 -d'/'`
# 	if [ "$r1" != "$r2" ]; then
# 		echo "name mismatch" $r1 $r2		
# 	else
# 		outfile=`echo "elegans_"$r1"_aln_pe.sam"`
# 		echo $genomefile
# 		$BWA_HOME/bwa mem -t 4 $BASEDIR/genomeSeq/elegans ${read1_files[$i]} ${read2_files[$i]} > $BASEDIR/samFiles/$outfile 		
# 	fi
# done

#to run on cluster:
i=$1
let i=i-1
r1=`echo ${read1_files[$i]%_read1.fastq.gz} | cut -f7 -d'/'`
r2=`echo ${read2_files[$i]%_read2.fastq.gz} | cut -f7 -d'/'`
if [ "$r1" != "$r2" ]; then
	echo "name mismatch" $r1 $r2		
else
	outfile=`echo "elegans_"$r1"_aln_pe.sam"`
	#echo $genomefile
	$BWA_HOME/bwa mem -t 4 $GENOMEDIR/elegans ${read1_files[$i]} ${read2_files[$i]} > $BASEDIR/samFiles/$GENOME_VER/$outfile 		
fi

