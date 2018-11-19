
# Run batch queue
sbatch -p csde -A csde runsim.sh 

# Run checkpoint queue (previously called backfill)
#sbatch -p ckpt -A csde-ckpt runsim.sh

# Pass in environmental variables
sbatch -p csde -A csde --export= runsim.sh 

# Array job
sbatch -p csde -A csde --array=1-5 --export=SIMNO=1000,NJOBS=5, runsim.sh 