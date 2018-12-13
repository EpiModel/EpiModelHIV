
# Run batch queue
#sbatch -p csde -A csde runsim.sh 

#Multiple Jobs
#sbatch -p csde -A csde --export=SIMNO=1000,NJOBS=4 runsim.sh 

# Array job
sbatch -p csde -A csde --array=1-6 --export=ALL,SIMNO=1000,NJOBS=12 runsim.sh 