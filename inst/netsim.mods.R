rm(list = ls())
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))
suppressMessages(library("dplyr"))
sourceDir("R/")


data(est)
data(st)
est
st
param <- param_msm(nwstats = st,
                   ai.scale = 1.11,
                   
                   rsyph.tprob = 0.04668348,
                   usyph.tprob = 0.03598792,
                   
                   hiv.rsyph.rr = 2.98876572, 
                   hiv.usyph.rr = 1.7456618,
                   syph.rhiv.rr = 6.54189295,
                   syph.uhiv.rr = 5.09641658,
                   
                   syph.earlat.rr = 0.5,
                   incu.syph.int = 27,
                   prim.syph.int = 60,
                   seco.syph.int = 120,
                   earlat.syph.int = 365 - 27 - 60 - 120,
                   latelat.syph.int = 9 * 52 * 7,
                   latelatelat.syph.int = 20 * 52 * 7,
                   tert.syph.int = 20 * 52 * 7,
                   syph.tert.prog.prob = 0.15 / (52 * 7 * 20),
                   
                   rgc.tprob = 0.4133300,
                   ugc.tprob = 0.31404720,
                   rct.tprob = 0.1907554,
                   uct.tprob = 0.16394697,
                   
                   hiv.rgc.rr = 2.35,
                   hiv.ugc.rr = 1.35,
                   hiv.rct.rr = 2.35,
                   hiv.uct.rr = 1.35,
                   
                   # adjust prim and seco from 0.1385 each
                   stage.syph.B.prob = c(0.00, 0.20, 0.077, 0.277, 0.22, 0.22, 0.006),
                   stage.syph.W.prob = c(0.00, 0.20, 0.077, 0.277, 0.22, 0.22, 0.006),
                   
                   syph.prim.sympt.prob.tx = 0.35, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008 use 0.45
                   syph.prim.asympt.prob.tx = 0.35,
                   syph.seco.sympt.prob.tx = 0.60, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.seco.asympt.prob.tx = 0.60,
                   syph.earlat.prob.tx = 0.15, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.latelat.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 0.90,
                   syph.tert.asympt.prob.tx = 0.90,
                   
                   hivdx.syph.sympt.tx.rr = 1.45,

                   prep.coverage = 0.4,
                   ept.coverage = 0.4,

                   stianntest.coverage = 0.5,
                   stihighrisktest.coverage = 0.8,

                   prep.start = 5000,
                   stitest.start = 5000,
                   ept.start = 5000,
                   
                   stitest.elig.model = "all",

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months

init <- init_msm(nwstats = st,
                 prev.B = 0.10,
                 prev.W = 0.10,
                 prev.ugc = 0.015,
                 prev.rgc = 0.015,
                 prev.uct = 0.015,
                 prev.rct = 0.015,
                 prev.syph.B = 0.020,
                 prev.syph.W = 0.020)

control <- control_msm(simno = 0.253,
                       nsteps = 2600,
                       nsims = 1,
                       ncores = 1,
                       save.nwstats = TRUE,
                       verbose.int = 1)

sim <- netsim(est, param, init, control)


at <- 1
dat <- initialize_msm(est, param, init, control, s = 1)
# dat <- reinit_msm(sim, param, init, control, s = 1)

at <- at + 1
dat <- aging_msm(dat, at)
dat <- deaths_msm(dat, at)
dat <- births_msm(dat, at)
dat <- test_msm(dat, at)
dat <- test_sti_msm(dat, at)
dat <- tx_msm(dat, at)
dat <- prep_msm(dat, at)
dat <- ept_msm(dat, at)
dat <- progress_msm(dat, at)
dat <- progress_syph_msm(dat, at)
dat <- vl_msm(dat, at)
#dat <- update_vl_msm(dat, at)
dat <- update_aiclass_msm(dat, at)
dat <- update_roleclass_msm(dat, at)
dat <- edges_correct_msm(dat, at)
dat <- simnet_msm(dat, at)
dat <- disclose_msm(dat, at)
dat <- part_msm(dat, at)
dat <- acts_msm(dat, at)
dat <- condoms_msm(dat, at)
dat <- riskhist_msm(dat, at)
dat <- position_msm(dat, at)
dat <- trans_msm(dat, at)
dat <- sti_trans(dat, at)
dat <- sti_recov(dat, at)
dat <- sti_tx(dat, at)
dat <- prevalence_msm(dat, at)
verbose_msm(dat, type = "progress", s = 1, at)

