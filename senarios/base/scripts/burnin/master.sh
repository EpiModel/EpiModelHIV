#!/bin/bash

# runs abc fit
qsub -q batch runsim.burn.abcsmc2.sh


# runs burnin model
qsub -q batch -t 1-9 -v SIMNO=100,NJOBS=9 runsim.burn.sh


