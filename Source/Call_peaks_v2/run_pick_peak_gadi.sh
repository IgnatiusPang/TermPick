#!/bin/bash

#PBS -N pickPeak
#PBS -l ncpus=4
#PBS -l mem=128gb
#PBS -l walltime=05:00:00
#PBS -j oe
#PBS -M i.pang@unsw.edu.au
#PBS -m ae
#PBS -P u71
#PBS -l storage=scratch/u71+gdata/u71
#PBS -l jobfs=300GB

### To check he status of array jobs, this is the command:
# qstat  -t -n1 -u z3371724

# https://opus.nci.org.au/display/Help/What+does+exceeded+memory+allocation+mean

module load openmpi #  openmpi/4.0.2(default)

module load gcc/system
module load R/3.6.1

DATE=$(date +"%Y%m%d")

### Create relevant directories

echo "My PBS directory"
echo $PBS_O_WORKDIR

cd $PBS_O_WORKDIR

HPC_CLUSTER="gadi"
NUM_CORES=4
COUNTS=1

for WINDOW_LEN in {10..150..10}
do
  for STD_DEV in {1..10..1}
  do

    echo ARRAY_INDEX  $PBS_ARRAY_INDEX
    echo COUNTS $COUNTS
    echo WINDOW_LEN $WINDOW_LEN
    echo STD_DEV $STD_DEV

    if [[ $COUNTS -eq $PBS_ARRAY_INDEX ]]
    then
      ## First command to run to pick the peaks
      COMMAND_1_TO_RUN="Rscript --vanilla  pick_peak_step.R ${HPC_CLUSTER} ${NUM_CORES} ${WINDOW_LEN} \
        ${STD_DEV} > run_pick_peak_step_win${WINDOW_LEN}_std${STD_DEV}_${DATE}.log 2>&1"

      echo $COMMAND_1_TO_RUN

      $COMMAND_1_TO_RUN

      exit
    fi

    COUNTS=$((COUNTS+1))
  done
done


