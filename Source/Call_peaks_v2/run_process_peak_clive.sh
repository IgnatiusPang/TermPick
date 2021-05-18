#!/bin/bash

#PBS -N procPeak
#PBS -l nodes=1:ppn=1
#PBS -l mem=60gb
#PBS -l walltime=05:00:00
#PBS -j oe
#PBS -M i.pang@unsw.edu.au
#PBS -J 1-110
#PBS -m aej

### To check he status of array jobs, this is the command:
# qstat  -t -n1 -u z3371724

module load R/3.6.3

DATE=$(date +"%Y%m%d")

### Create relevant directories

echo "My PBS directory"
echo $PBS_O_WORKDIR

cd $PBS_O_WORKDIR

HPC_CLUSTER="clive"
COUNTS=1


MY_ARRAY=($(seq 2 1 9))
MY_ARRAY+=($(seq 10 10 30))
printf "%s\n" "${MY_ARRAY[@]}"

for WINDOW_LEN in ${MY_ARRAY[@]}
do
  for STD_DEV in {1..10..1}
  do
    if [[ $COUNTS -eq $PBS_ARRAY_INDEX ]]
    then

                echo ARRAY_INDEX  $PBS_ARRAY_INDEX
                echo COUNTS $COUNTS
                echo WINDOW_LEN $WINDOW_LEN
                echo STD_DEV $STD_DEV

                ## First command to run to pick the peaks
                COMMAND_2_TO_RUN="Rscript --vanilla  process_peak_pick.R  > \
                  $PBS_O_WORKDIR/run_process_peak_step_win${WINDOW_LEN}_std${STD_DEV}_${DATE}.log 2>&1"

                echo $COMMAND_2_TO_RUN

                $COMMAND_2_TO_RUN

        exit 0
    fi

    COUNTS=$((COUNTS+1))
  done
done


