
library("methods")
library("EpiModelHPC")
library("EpiModelHIV")

args <- commandArgs(trailingOnly = TRUE)
simno <- args[1]
jobno <- args[2]
fsimno <- paste(simno, jobno, sep = ".")
print(fsimno)

load(file = "est/fitsmall.rda")
load(file = "est/data.params.rda")

param <- param_shamp(data.params,temp.adjust=1)
init <- init_shamp()
control <- control_shamp(nsteps = 2600)

netsim_hpc("est/fitsmall.rda", param, init, control, compress = "gzip")
