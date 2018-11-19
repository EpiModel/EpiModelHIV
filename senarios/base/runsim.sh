

#!/bin/bash 

## Job Name 
#SBATCH --job-name=SHAMP

## Allocation Definition
#SBATCH --account=csde
#SBATCH --partition=csde

## Nodes
#SBATCH --nodes=5

## Tasks per node
#SBATCH --ntasks-per-node=8

## Walltime
#SBATCH --time=6:00:00:00 

## E-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dth2@uw.edu

## Memory per node
#SBATCH --mem=58G

## Specify the working directory
#SBATCH --workdir=/gscratch/csde/deven/SHAMP/scenarios/base


module load r-3.5.0

Rscript sim2.R
