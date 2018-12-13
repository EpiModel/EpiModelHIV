#!/bin/bash

### User specs
#PBS -N shamp-abc-calib
#PBS -l nodes=1:ppn=16,mem=50gb,feature=16core,walltime=10:00:00:00
#PBS -o /gscratch/csde/deven/SHAMP/scenarios/base/out
#PBS -e /gscratch/csde/deven/SHAMP/scenarios/base/out
#PBS -j oe
#PBS -d /gscratch/csde/deven/SHAMP/scenarios/base
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
R CMD BATCH --vanilla sim.burn.abc.shamp.R out/sim.burn.abc.shamp.n${NSIM}.p${PACC}.Rout
