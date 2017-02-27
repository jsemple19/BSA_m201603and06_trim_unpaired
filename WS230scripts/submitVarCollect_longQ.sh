#! /bin/bash
#$ -N varCollect
#$ -cwd
#$ -q long-sl65
#$ -l h_rt=160:00:00
#$ -M jennifer.semple@crg.es
#$ -m abe
#$ -o /users/blehner/jsemple/outputs/
#$ -e /users/blehner/jsemple/errors/
#$ -t 1-8

GENOME_VER=WS230
BASEDIR=/users/blehner/jsemple/seqResults/m201603and06
mkdir -p $BASEDIR/finalData/$GENOME_VER
all_files=(`ls $BASEDIR/vcfFiles/$GENOME_VER/*.vcf`)
let i=$SGE_TASK_ID-1
python CollectAllSNPs.py ${all_files[$i]}
