#!/bin/bash

qsub -q batch -m ae -v NSIM=200,PACC=0.01 runsim.burn.abc.shamp.sh
qsub -q batch -m ae -v NSIM=200,PACC=0.05 runsim.burn.abc.shamp.sh
qsub -q batch -m ae -v NSIM=200,PACC=0.10 runsim.burn.abc.shamp.sh
qsub -q batch -m ae -v NSIM=200,PACC=0.20 runsim.burn.abc.shamp.sh
