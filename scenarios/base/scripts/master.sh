#!/bin/bash 

# Run batch queue
#sbatch -p csde -A csde runsim.sh 

#Multiple Jobs
#sbatch -p csde -A csde --export=ALL,SIMNO=50,NJOBS=5 runsim.sh 

# Array job
sbatch -p csde -A csde --array=1-5 --export=ALL,SIMNO=80,NJOBS=5 runsim.sh 

