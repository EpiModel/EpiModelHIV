rm(list = ls())
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("dplyr"))
devtools::load_all()

#data(est)
#data(st)
load("est/nwstats.rda")
load("est/fit.rda")
#load("est/stimod.burnin.rda")

param <- param_msm(nwstats = st,

                   ai.scale = 1.04,
                   ai.scale.pospos = 1.04,

                   # STI acquisition
                   rgc.tprob = 0, #0.65,
                   ugc.tprob = 0, #0.55,
                   rct.tprob = 0, #0.29,
                   uct.tprob = 0, #0.23,
                   syph.tprob = 0.06,

                   # HIV acquisition
                   hiv.rgc.rr = 1.75,
                   hiv.ugc.rr = 1.26,
                   hiv.rct.rr = 1.75,
                   hiv.uct.rr = 1.26,
                   hiv.syph.rr = 1.63,

                   syph.incub.sympt.prob = 0,
                   syph.prim.sympt.prob = 0.70,
                   syph.seco.sympt.prob = 0.85,
                   syph.earlat.sympt.prob = 0,
                   syph.latelat.sympt.prob = 0,
                   syph.tert.sympt.prob = 1.0,

                   syph.prim.sympt.prob.tx = 0.80,
                   syph.seco.sympt.prob.tx = 0.80,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   ept.coverage = 0.5,
                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.syph.hivneg.coverage = 0.45,
                   stihighrisktest.gc.hivneg.coverage = 0.1,
                   stihighrisktest.ct.hivneg.coverage = 0.1,
                   stihighrisktest.syph.hivneg.coverage = 0.1,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.61,
                   stianntest.syph.hivpos.coverage = 0.67,
                   stihighrisktest.gc.hivpos.coverage = 0.1,
                   stihighrisktest.ct.hivpos.coverage = 0.1,
                   stihighrisktest.syph.hivpos.coverage = 0.1,

                   prep.start = 7000,
                   stitest.start = 7000,
                   ept.start = 5201,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182,
                   ept.risk.int = 60)

init <- init_msm(nwstats = st,
                 prev.B = 0.10,
                 prev.W = 0.10,
                 prev.ugc = 0.010,
                 prev.rgc = 0.010,
                 prev.uct = 0.010,
                 prev.rct = 0.010,
                 prev.syph.B = 0.010,
                 prev.syph.W = 0.010)

control <- control_msm(nsteps = 2600)

# control <- control_msm(start = 5201,
#                        nsteps = 5210)

#sim <- netsim(est, param, init, control)

# debug(initialize_msm)
# debug(reinit_msm)
# debug(prevalence_msm)

at <- 1
# at <- 5201

dat <- initialize_msm(est, param, init, control, s = 1)
# dat <- reinit_msm(sim, param, init, control, s = 1)

debug(sti_trans_msm)
debug(riskhist_stitest_msm)
debug(sti_test_msm)
debug(sti_recov_msm)
debug(sti_tx_msm)
debug(sti_ept_msm)
debug(part_msm)
# debug(acts_msm)
# debug(hiv_test_msm)
# debug(acts_msm)
# debug(hiv_trans_msm)
# debug(condoms_msm)
# debug(simnet_msm)
# debug(syph_progress_msm)

at <- at + 1
#for (at in at:dat$control$nsteps) {
  dat <- aging_msm(dat, at)
  dat <- deaths_msm(dat, at)
  dat <- births_msm(dat, at)
  dat <- hiv_test_msm(dat, at)
  dat <- sti_test_msm(dat, at)
  dat <- hiv_tx_msm(dat, at)
  dat <- prep_msm(dat, at)
  dat <- hiv_progress_msm(dat, at)
  dat <- syph_progress_msm(dat, at)
  dat <- hiv_vl_msm(dat, at)
  dat <- simnet_msm(dat, at)
  dat <- hiv_disclose_msm(dat, at)
  dat <- part_msm(dat, at)
  dat <- acts_msm(dat, at)
  dat <- condoms_msm(dat, at)
  dat <- riskhist_prep_msm(dat, at)
  dat <- riskhist_stitest_msm(dat, at)
  dat <- position_msm(dat, at)
  dat <- hiv_trans_msm(dat, at)
  dat <- sti_trans_msm(dat, at)
  dat <- sti_recov_msm(dat, at)
  dat <- sti_tx_msm(dat, at)
  dat <- sti_ept_msm(dat, at)
  dat <- prevalence_msm(dat, at)
  cat("\t", at)
#}
