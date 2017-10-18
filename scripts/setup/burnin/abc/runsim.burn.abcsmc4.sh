#!/bin/bash

### User specs
#PBS -N sti-abc-recalib
#PBS -l nodes=1:ppn=16,mem=50gb,feature=16core,walltime=05:00:00:00
#PBS -o /gscratch/csde/sjenness/stia/out
#PBS -e /gscratch/csde/sjenness/stia/out
#PBS -j oe
#PBS -d /gscratch/csde/sjenness/stia
#PBS -m ae

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
R CMD BATCH --vanilla sim.burn.abcsmc4.R out/sim.burn.abcsmc4.n${NSIM}.p${PACC}.Rout
