#!/bin/bash 

## Job Name 
#SBATCH --job-name=s50


## Allocation Definition
#SBATCH --account=csde
#SBATCH --partition=csde

## Nodes
#SBATCH --nodes=1

## Tasks per node
#SBATCH --ntasks-per-node=5

## Walltime
#SBATCH --time=24-12 

## E-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dth2@uw.edu

## Memory per node
#SBATCH --mem=52G

## Specify the working directory
#SBATCH --workdir=/gscratch/csde/deven/SHAMP/scenarios/base


. /suppscr/csde/sjenness/spack/share/spack/setup-env.sh
module load gcc-8.2.0-gcc-8.1.0-sh54wqg
module load r-3.5.2-gcc-8.2.0-sby3icq


Rscript sim.R
