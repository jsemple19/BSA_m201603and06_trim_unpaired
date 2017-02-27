#! /bin/bash
#$ -N concat
#$ -cwd
#$ -q long-sl65
#$ -l h_rt=160:00:00
#$ -o /users/blehner/jsemple/outputs/
#$ -e /users/blehner/jsemple/errors/
#$ -t 1-16

./concatFastq.sh $SGE_TASK_ID
