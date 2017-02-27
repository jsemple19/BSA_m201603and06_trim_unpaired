#! /bin/bash
#$ -N trimReads
#$ -cwd
#$ -q long-sl65
#$ -l h_rt=6:00:00 
#$ -l virtual_free=12G 
#$ -M jennifer.semple@crg.es
#$ -m abe
#$ -o /users/blehner/jsemple/outputs/
#$ -e /users/blehner/jsemple/errors/
#$ -t 1-8 

./trimReads.sh $SGE_TASK_ID
