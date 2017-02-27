#! /bin/bash
#last modified 20160610

export TRIM_HOME="/users/blehner/jsemple/seqResults/Trimmomatic-0.36"

#reference genome version
#GENOME_VER="WS230"

DATADIR=/users/blehner/sequencing_data/Jennifer_Semple/2016-06-01
NEWDATADIR=/users/blehner/sequencing_data/Jennifer_Semple/2016-06-01_trim
#BASEDIR=/users/blehner/jsemple/seqResults/20160601_trim
#mkdir -p $BASEDIR/samFiles/$GENOME_VER

#genomefile=`ls $BASEDIR/genomeSeq/$GENOME_VER | grep 'c_elegans.*genomic\.fa$'`



read1_files=(`ls $DATADIR/*read1.fastq.gz`)
read2_files=(`ls $DATADIR/*read2.fastq.gz`)


i=$1
let i=i-1
r1=`echo ${read1_files[$i]%_read1.fastq.gz} | cut -f7 -d'/'`
r2=`echo ${read2_files[$i]%_read2.fastq.gz} | cut -f7 -d'/'`
if [ "$r1" != "$r2" ]; then
	echo "name mismatch" $r1 $r2		
else
	outfile1=`echo $NEWDATADIR/"paired_"$r1"_r1.fastq.gz"`
	outfile1u=`echo $NEWDATADIR/"unpaired_"$r1"_r1.fastq.gz"`
	outfile2=`echo $NEWDATADIR/"paired_"$r2"_r2.fastq.gz"`
	outfile2u=`echo $NEWDATADIR/"unpaired_"$r2"_r2.fastq.gz"`
	java -jar $TRIM_HOME/trimmomatic-0.36.jar PE -threads 4 -phred33 ${read1_files[i]} ${read2_files[i]} $outfile1 $outfile1u $outfile2 $outfile2u ILLUMINACLIP:$TRIM_HOME/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36

fi

