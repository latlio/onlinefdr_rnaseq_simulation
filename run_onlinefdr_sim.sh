#!/bin/bash
#
#$ -cwd -V -b y
#$ -P vlong
#$ -o tmp.out -e tmp.err
#

R CMD BATCH --no-save --no-restore $1 $SGE_TASK_ID.Rout