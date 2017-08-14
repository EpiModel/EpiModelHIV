rm(list = ls())
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("dplyr"))
devtools::load_all()

#data(est)
#data(st)
load("est/nwstats.rda")
load("est/fit.rda")

param <- param_msm(nwstats = st,
                   ai.scale = 1.03,

                   syph.earlat.rr = 0.5,
                   incu.syph.int = 27,
                   prim.syph.int = 60,
                   seco.syph.int = 120,
                   earlat.syph.int = 365 - 27 - 60 - 120,
                   latelat.syph.int = 9 * 52 * 7,
                   latelatelat.syph.int = 20 * 52 * 7,
                   tert.syph.int = 20 * 52 * 7,
                   syph.tert.prog.prob = 0.00015625599,

                   # STI acquisition
                   rgc.tprob = 0.4434,
                   ugc.tprob = 0.3343,
                   rct.tprob = 0.2008,
                   uct.tprob = 0.1790,
                   syph.tprob = 0.1424,

                   # HIV acquisition
                   hiv.rgc.rr = 1.80292790,
                   hiv.ugc.rr = 1.1989083,
                   hiv.rct.rr = 1.80292790,
                   hiv.uct.rr = 1.1989083,
                   hiv.syph.rr = 1.62,

                   # HIV transmission
                   hiv.trans.gc.rr = 1,
                   hiv.trans.ct.rr = 1,
                   hiv.trans.syph.rr = 1,

                   syph.prim.sympt.prob.tx = 0.60,
                   syph.seco.sympt.prob.tx = 0.688235,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   syph.prim.asympt.prob.tx = 1.0,
                   syph.seco.asympt.prob.tx = 1.0,
                   syph.earlat.asympt.prob.tx = 1.0,
                   syph.latelat.asympt.prob.tx = 1.0,
                   syph.tert.asympt.prob.tx = 1.0,

                   hivdx.syph.sympt.tx.rr = 1.45,

                   prep.coverage = 0.0,
                   ept.coverage = 0.8,
                   stianntest.coverage = 0.1,
                   stihighrisktest.coverage = 0.1,

                   prep.start = 2600,
                   stitest.start = 2600,
                   ept.start = 5,

                   stitest.elig.model = "sti",

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182,
                   ept.risk.int = 60)

init <- init_msm(nwstats = st,
                 prev.B = 0.08,
                 prev.W = 0.08,
                 prev.ugc = 0.010,
                 prev.rgc = 0.010,
                 prev.uct = 0.010,
                 prev.rct = 0.010,
                 prev.syph.B = 0.010,
                 prev.syph.W = 0.010,
                 stage.syph.B.prob = c(0.00, 0.01, 0.01, 0.10, 0.44, 0.44, 0.00),
                 stage.syph.W.prob = c(0.00, 0.01, 0.01, 0.10, 0.44, 0.44, 0.00))

control <- control_msm(nsteps = 200)

sim <- netsim(est, param, init, control)

at <- 1
dat <- initialize_msm(est, param, init, control, s = 1)
# dat <- reinit_msm(sim, param, init, control, s = 1)

debug(sti_ept_msm)

at <- at + 1
for (at in 2:control$nsteps) {
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
