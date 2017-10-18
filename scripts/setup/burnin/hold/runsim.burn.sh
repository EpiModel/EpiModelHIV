#!/bin/bash

### User specs
#PBS -N sti-burn
#PBS -l nodes=1:ppn=16,mem=50gb,feature=16core,walltime=02:00:00
#PBS -o /gscratch/csde/sjenness/sti/out
#PBS -e /gscratch/csde/sjenness/sti/out
#PBS -j oe
#PBS -d /gscratch/csde/sjenness/sti
#PBS -m n

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
ALLARGS="${SIMNO} ${PBS_ARRAYID}"

### App
Rscript sim.burn.R ${ALLARGS}
