


library("methods")
library("EpiModelHPC")
library("EpiModelHIV")

args <- commandArgs(trailingOnly = TRUE)
simno <- args[1]
jobno <- args[2]
fsimno <- paste(simno, jobno, sep = ".")
print(fsimno)

load(file = "est/fit.rda")
load(file = "est/data.params.rda")

param <- param_shamp(data.params)
init <- init_shamp()
control <- control_shamp(nsteps = 52)

netsim_hpc("est/fit.rda", param, init, control, compress = "gzip")


####################   SLURM  ##############################
library("methods")
library("EpiModelHPC")
library("EpiModelHIV")

simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")

cat("Array number is ", jobno)
cat("\n fsimno is ", fsimno)

##LOAD FILES
load(file = "est/fit.rda")
load(file = "est/data.params.rda")

# Epidemic model
param <- param_shamp(data.params)
init <- init_shamp()
control <- control_shamp(nsteps = 52)

mod1 <- netsim(est, param, init, control)

print(mod1)

fn <- paste0("sim.", fsimno, ".rda")
save(mod1, file = fn)