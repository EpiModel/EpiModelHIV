
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
library("EpiModelHPC")
sourceDir("source/", FALSE)

## Environmental Arguments
args <- commandArgs(TRUE)
simno <- args[1]
jobno <- args[2]

## Parameters
fsimno <- paste(simno, jobno, sep = ".")
load("est/nwstats.rda")

param <- param_msm(nwstats = st,
                   ai.scale = 1,

                   riskh.start = 1e4,
                   prep.start = 1e4,
                   prep.coverage = 0,

                   rcomp.prob = 0,
                   rcomp.adh.groups = 0:4,
                   rcomp.main.only = FALSE,
                   rcomp.discl.only = FALSE,

                   rgc.tprob = 0.60,
                   ugc.tprob = 0.48,
                   rct.tprob = 0.40,
                   uct.tprob = 0.32,

                   rgc.sympt.prob = 0.16,
                   ugc.sympt.prob = 0.90,
                   rct.sympt.prob = 0.14,
                   uct.sympt.prob = 0.58,

                   rgc.dur.asympt = 300 / 7,
                   ugc.dur.asympt = 240 / 7,
                   gc.dur.tx = 13 / 7,
                   gc.dur.ntx = 185 / 7,

                   rct.dur.asympt = 497 / 7,
                   uct.dur.asympt = 240 / 7,
                   ct.dur.tx = 14 / 7,
                   ct.dur.ntx = 180 / 7,

                   gc.prob.cease = 0,
                   ct.prob.cease = 0,

                   gc.prob.tx = 0.90,
                   ct.prob.tx = 0.85,

                   prep.sti.screen.int = 182,
                   prep.sti.prob.tx = 1,

                   sti.cond.rr = 0.3,

                   hiv.rgc.rr = 2.25,
                   hiv.ugc.rr = 1.50,
                   hiv.rct.rr = 2.25,
                   hiv.uct.rr = 1.50)

init <- init_msm(nwstats = st,
                 prev.B = 0.2,
                 prev.W = 0.2,
                 prev.ugc = 0.1,
                 prev.rgc = 0.1,
                 prev.uct = 0.1,
                 prev.rct = 0.1)

control <- control_msm(simno = fsimno,
                       nsteps = 52 * 25,
                       nsims = 16,
                       ncores = 16,
                       save.int = 1e4,
                       acts.FUN = acts_sti,
                       condoms.FUN = condoms_sti,
                       initialize.FUN = initialize_sti,
                       prep.FUN = prep_sti,
                       prev.FUN = prevalence_sti,
                       riskhist.FUN = riskhist_sti,
                       position.FUN = position_sti,
                       trans.FUN = trans_sti,
                       stitrans.FUN = sti_trans,
                       stirecov.FUN = sti_recov,
                       stitx.FUN = sti_tx,
                       verbose.FUN = verbose_sti,
                       verbose.int = 25,
                       module.order = c("aging.FUN", "deaths.FUN", "births.FUN",
                                        "test.FUN", "tx.FUN", "prep.FUN",
                                        "progress.FUN", "vl.FUN",
                                        "resim_nets.FUN", "disclose.FUN",
                                        "acts.FUN", "condoms.FUN", "riskhist.FUN",
                                        "position.FUN", "trans.FUN", "stitrans.FUN",
                                        "stirecov.FUN", "stitx.FUN", "prev.FUN"))

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = FALSE)
