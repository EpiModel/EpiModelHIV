rm(list=ls())
suppressMessages(library("EpiModelHIV"))
sourceDir("R/")


data(est)
data(st)
est
st
param <- param_msm(nwstats = st,
                   ai.scale = 1.13,
                   syph.tprob = 0.021,
                   rgc.tprob = 0.40,
                   ugc.tprob = 0.30,
                   rct.tprob = 0.24,
                   uct.tprob = 0.2,
                   hiv.syph.rr = 2.00,
                   syph.hiv.rr = 2.40,
                   prep.coverage = 0,
                   ept.coverage = 0,
                   rgc.sympt.prob = 0.16, # Beck
                   ugc.sympt.prob = 0.90, # Beck
                   rct.sympt.prob = 0.14, # Beck
                   uct.sympt.prob = 0.58, # Beck
                   stitest.int = 182) # adjustable for 3 or 6 months
init <- init_msm(nwstats = st, 
                 prev.B = 0.10, 
                 prev.W = 0.10,
                 prev.ugc = 0.013,
                 prev.rgc = 0.013,
                 prev.uct = 0.013,
                 prev.rct = 0.013,
                 prev.syph.B = 0.01,
                 prev.syph.W = 0.01)
control <- control_msm(simno = 0.253, 
                       nsteps = 52*50,
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
dat <- updatenwp_msm(dat, at)
dat <- simnet_msm(dat, at)
dat <- disclose_msm(dat, at)
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

