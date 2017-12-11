
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))

## Environmental Arguments
simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("PBS_ARRAYID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")

cov <- as.numeric(Sys.getenv("COV"))
prstiint <- as.numeric(Sys.getenv("PSTIINT"))
rc <- as.numeric(Sys.getenv("RC"))
probtx <- as.numeric(Sys.getenv("PROBTX"))
asymptx <- as.numeric(Sys.getenv("ASYMPTX"))

## Parameters
load("est/nwstats.rda")

param <- param_msm(nwstats = st,

                   rgc.dur.asympt = 35.11851,
                   ugc.dur.asympt = 35.11851,
                   rct.dur.asympt = 44.24538,
                   uct.dur.asympt = 44.24538,
                   hiv.rgc.rr = 2.780673,
                   hiv.ugc.rr = 1.732363,
                   hiv.rct.rr = 2.780673,
                   hiv.uct.rr = 1.732363,

                   prep.coverage = cov,
                   prep.start = 2601,

                   rcomp.prob = rc,
                   rcomp.adh.groups = 2:3,

                   prep.sti.screen.int = prstiint,
                   prep.sti.prob.tx = probtx,

                   gc.asympt.prob.tx = asymptx,
                   ct.asympt.prob.tx = asymptx)

init <- init_msm(st)

control <- control_msm(simno = fsimno,
                       start = 2601,
                       nsteps = 3120,
                       nsims = 16,
                       ncores = 16,
                       initialize.FUN = reinit_msm,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/stimod.burnin.rda", param, init, control,
           compress = FALSE, verbose = FALSE)

process_simfiles(simno = simno, min.n = njobs,
                 outdir = "data/", compress = TRUE)
