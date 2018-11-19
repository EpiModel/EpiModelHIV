#!/bin/bash

### User specs
#PBS -N sim$SIMNO
#PBS -l nodes=1:ppn=16,mem=44gb,feature=16core,walltime=01:00:00
#PBS -o /gscratch/csde/deven/camp/adol/out
#PBS -e /gscratch/csde/deven/camp/adol/out
#PBS -j oe
#PBS -d /gscratch/csde/camp/adol
#PBS -m n
#PBS -M dth2@u.washington.edu

### Standard specs
HYAK_NPE=$(wc -l < $PBS_NODEFILE)
HYAK_NNODES=$(uniq $PBS_NODEFILE | wc -l )
HYAK_TPN=$((HYAK_NPE/HYAK_NNODES))
NODEMEM=`grep MemTotal /proc/meminfo | awk '{print $2}'`
NODEFREE=$((NODEMEM-2097152))
MEMPERTASK=$((NODEFREE/HYAK_TPN))
ulimit -v $MEMPERTASK
export MX_RCACHE=0

### Modules
module load r_3.2.4

### App
Rscript process.R
