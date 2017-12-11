
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))

## Environmental Arguments
simno <- Sys.getenv("SIMNO")
jobno <- Sys.getenv("PBS_ARRAYID")
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")

## Parameters
load("est/nwstats.rda")

param <- param_msm(nwstats = st,
                   rgc.tprob = 0.398,
                   ugc.tprob = 0.320,
                   rct.tprob = 0.221,
                   uct.tprob = 0.205,
                   rct.asympt.int = 254.1,
                   uct.asympt.int = 254.1,
                   hiv.rgc.rr = 2.78,
                   hiv.ugc.rr = 1.73,
                   hiv.rct.rr = 2.78,
                   hiv.uct.rr = 1.73,
                   hiv.dual.rr = 0.2)
init <- init_msm(nwstats = st)
control <- control_msm(simno = fsimno,
                       nsteps = 3120,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control, verbose = FALSE)
process_simfiles(simno = simno, min.n = njobs, compress = TRUE,
                 outdir = "data/", verbose = FALSE)
