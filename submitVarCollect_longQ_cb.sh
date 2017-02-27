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

MYPYTHON=/software/bl/el6.3/Python-3.5.0

export PATH=/software/bl/el6.3/Python-3.5.0/bin:$PATH
export LD_LIBRARY_PATH=/software/bl/el6.3/Python-3.5.0/lib:$LD_LIBRARY_PATH


GENOME_VER=PRJNA275000.WS250
BASEDIR=/users/blehner/jsemple/seqResults/m201603and06_trim_unpaired
mkdir -p $BASEDIR/finalData/$GENOME_VER
all_files=(`ls $BASEDIR/vcfFiles/$GENOME_VER/*raw.vcf`)
let i=$SGE_TASK_ID-1
$MYPYTHON/python CollectAllSNPs.py ${all_files[$i]}


