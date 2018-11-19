

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

param <- param_shamp(data.params, 
                     VI.foi.scale = 1.8,
                     msm.foi.scale = 3.8,
                     fa.foi.scale = 12,
                     death_stats = FALSE,
                     p.growth = TRUE,
                     p.growth.nsteps = 75,
                     p.growth.size = 10000)

init <- init_shamp(prev.B.f = 0.1,
                   prev.BI.f =0.1,
                   prev.H.f =0.1,
                   prev.HI.f =0.1,
                   prev.W.f = 0.1,
                   prev.B.msf = 0.1,
                   prev.BI.msf =0.1,
                   prev.H.msf =0.1,
                   prev.HI.msf =0.1,
                   prev.W.msf = 0.1,
                   prev.B.msm = 0.1,
                   prev.BI.msm =0.1,
                   prev.H.msm =0.1,
                   prev.HI.msm =0.1,
                   prev.W.msm = 0.1,
                   prev.B.msmf = 0.1,
                   prev.BI.msmf =0.1,
                   prev.H.msmf =0.1,
                   prev.HI.msmf =0.1,
                   prev.W.msmf = 0.1)

control <- control_shamp(nsteps = 52)

mod1 <- netsim(est, param, init, control)

print(mod1)

fn <- paste0("sim.", fsimno, ".rda")
save(mod1, file = fn)