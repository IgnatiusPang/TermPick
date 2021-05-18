#!/bin/bash

#PBS -N readGRfiles
#PBS -l nodes=1:ppn=1
#PBS -l mem=21gb
#PBS -l walltime=05:00:00
#PBS -j oe
#PBS -M i.pang@unsw.edu.au
#PBS -m ae

### To check he status of array jobs, this is the command:
# qstat  -t -n1 -u z3371724

module load R/3.6.1

DATE=$(date +"%Y%m%d")

### Create relevant directories

echo "My PBS directory"
echo $PBS_O_WORKDIR

cd $PBS_O_WORKDIR

COMMAND_1_TO_RUN="Rscript --vanilla  read_gr_files.R > run_read_gr_files_${DATE}.log 2>&1"

echo $COMMAND_1_TO_RUN

$COMMAND_1_TO_RUN
