#!/bin/bash

#PBS -N masterTbl
#PBS -l nodes=1:ppn=1
#PBS -l mem=200gb
#PBS -l walltime=05:00:00
#PBS -j oe
#PBS -M i.pang@unsw.edu.au
#PBS -m ae

### To check he status of array jobs, this is the command:
# qstat  -t -n1 -u z3371724

module load R/3.6.3

DATE=$(date +"%Y%m%d")

### Create relevant directories

echo "My PBS directory"
echo $PBS_O_WORKDIR

cd $PBS_O_WORKDIR

## First command to run to pick the peaks
COMMAND_TO_RUN="Rscript --vanilla  merge_master_table.R  > \
   $PBS_O_WORKDIR/run_merge_master_table_${DATE}.log 2>&1"

echo $COMMAND_TO_RUN

$COMMAND_TO_RUN

