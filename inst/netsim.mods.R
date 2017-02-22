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
                   ai.scale = 1.12, # 1.11889726, # was 1.13
                   
                   rsyph.tprob = 0.049, # 0.050,
                   usyph.tprob = 0.040, #0.038, 
                   
                   hiv.rsyph.rr = 2.70, 
                   hiv.usyph.rr = 1.70,
                   syph.rhiv.rr = 3.00,
                   syph.uhiv.rr = 2.00,
                   
                   syph.earlat.rr = 0.5, #2/3, 0
                   incu.syph.int = 27,
                   prim.syph.int = 60,
                   seco.syph.int = 120,
                   earlat.syph.int = 365 - 120 - 60 - 27,
                   latelat.syph.int = 9*365,
                   latelatelat.syph.int = 20*365,
                   tert.syph.int = 20*365,
                   immune.syph.int = 5*365,
                   syph.tert.prog.prob = 0.15,
                   
                   rgc.tprob = 0.41, 
                   ugc.tprob = 0.31, 
                   rct.tprob = 0.21, 
                   uct.tprob = 0.15,
                   
                   hiv.rgc.rr = 2.70, #2.780673,
                   hiv.ugc.rr = 1.70, #1.732363,
                   hiv.rct.rr = 2.70, #2.780673,
                   hiv.uct.rr = 1.70, #1.732363,
                   
                   # adjust prim and seco from 0.1385 each
                   stage.syph.B.prob = c(0.00, 0.20, 0.077, 0.277, 0.22, 0.22, 0.006),
                   stage.syph.W.prob = c(0.00, 0.20, 0.077, 0.277, 0.22, 0.22, 0.006),
                   
                   syph.prim.sympt.prob.tx = 0.48201266, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008 use 0.45
                   syph.prim.asympt.prob.tx = 0.00,
                   syph.seco.sympt.prob.tx = 0.67214403, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.seco.asympt.prob.tx = 0.00,
                   syph.earlat.prob.tx = 0.17136638, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.latelat.prob.tx = 0.15437417,
                   syph.tert.sympt.prob.tx = 0.90,
                   syph.tert.asympt.prob.tx = 0.00,
                   
                   prep.coverage = 0,
                   ept.coverage = 0,
                   
                   prep.start = 5000,
                   stitest.start = 5000,
                   
                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months

init <- init_msm(nwstats = st, 
                 prev.B = 0.10, 
                 prev.W = 0.10,
                 prev.ugc = 0.015,
                 prev.rgc = 0.015,
                 prev.uct = 0.015,
                 prev.rct = 0.015,
                 prev.syph.B = 0.02,
                 prev.syph.W = 0.02)

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

