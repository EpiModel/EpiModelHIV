#!/bin/bash

qsub -q batch -t 1-10  -v SIMNO=1 runsim.sh
qsub -q batch -t 1-10  -v SIMNO=2 runsim.sh
qsub -q batch -t 1-10  -v SIMNO=3 runsim.sh
qsub -q batch -t 1-10  -v SIMNO=4 runsim.sh
