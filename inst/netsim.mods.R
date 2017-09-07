rm(list = ls())
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("dplyr"))
devtools::load_all()

#data(est)
#data(st)
load("est/nwstats.rda")
load("est/fit.rda")

param <- param_msm(nwstats = st,

                   # STI acquisition
                   rgc.tprob = 0.4456,
                   ugc.tprob = 0.3341,
                   rct.tprob = 0.1985,
                   uct.tprob = 0.1787,
                   syph.tprob = 0.17,

                   # HIV acquisition
                   hiv.rgc.rr = 1.80292790,
                   hiv.ugc.rr = 1.1989083,
                   hiv.rct.rr = 1.80292790,
                   hiv.uct.rr = 1.1989083,
                   hiv.syph.rr = 2.00,

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

                   ept.coverage = 0.0,
                   stianntest.gc.hivneg.coverage = 0.1,
                   stianntest.ct.hivneg.coverage = 0.1,
                   stianntest.syph.hivneg.coverage = 0.1,
                   stihighrisktest.gc.hivneg.coverage = 0.0,
                   stihighrisktest.ct.hivneg.coverage = 0.0,
                   stihighrisktest.syph.hivneg.coverage = 0.0,
                   stianntest.gc.hivpos.coverage = 0.1,
                   stianntest.ct.hivpos.coverage = 0.1,
                   stianntest.syph.hivpos.coverage = 0.1,
                   stihighrisktest.gc.hivpos.coverage = 0.0,
                   stihighrisktest.ct.hivpos.coverage = 0.0,
                   stihighrisktest.syph.hivpos.coverage = 0.0,

                   prep.start = 2601,
                   stitest.start = 2601,
                   ept.start = 2601,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182,
                   ept.risk.int = 60)

init <- init_msm(nwstats = st)

control <- control_msm(nsteps = 200)

sim <- netsim(est, param, init, control)

debug(initialize_msm)
at <- 1
dat <- initialize_msm(est, param, init, control, s = 1)
# dat <- reinit_msm(sim, param, init, control, s = 1)

debug(sti_recov_msm)
debug(sti_tx_msm)
debug(hiv_test_msm)
debug(sti_test_msm)

at <- at + 1
for (at in 2:198) {
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
}
