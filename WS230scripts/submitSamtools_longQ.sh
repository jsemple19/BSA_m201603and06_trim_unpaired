#! /bin/bash
#$ -N findVars
#$ -cwd
#$ -q long-sl65
#$ -l h_rt=160:00:00
#$ -M jennifer.semple@crg.es
#$ -m abe
#$ -o /users/blehner/jsemple/outputs/
#$ -e /users/blehner/jsemple/errors/
#$ -t 1-8

./samtools_aln.sh $SGE_TASK_ID
