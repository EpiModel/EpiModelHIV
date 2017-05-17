rm(list = ls())
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("dplyr"))
devtools::load_all()

data(est)
data(st)

param <- param_msm(nwstats = st,
                   ai.scale = 1.11,

                   rsyph.tprob = 0.0675,
                   usyph.tprob = 0.0495,

                   hiv.rsyph.rr = 2.77,
                   hiv.usyph.rr = 2.045,
                   syph.rhiv.rr = 1,
                   syph.uhiv.rr = 1,

                   syph.earlat.rr = 0.5,
                   incu.syph.int = 27,
                   prim.syph.int = 60,
                   seco.syph.int = 120,
                   earlat.syph.int = 365 - 27 - 60 - 120,
                   latelat.syph.int = 9 * 52 * 7,
                   latelatelat.syph.int = 20 * 52 * 7,
                   tert.syph.int = 20 * 52 * 7,
                   syph.tert.prog.prob = 0.00015625599,

                   rgc.tprob = 0.4207,
                   ugc.tprob = 0.3095,
                   rct.tprob = 0.1960,
                   uct.tprob = 0.1654,

                   hiv.rgc.rr = 2.548,
                   hiv.ugc.rr = 1.855,
                   hiv.rct.rr = 2.548,
                   hiv.uct.rr = 1.855,

                   syph.prim.sympt.prob.tx = 0.35,
                   syph.seco.sympt.prob.tx = 0.60,
                   syph.earlat.sympt.prob.tx = 0.15,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 0.90,

                   syph.prim.asympt.prob.tx = 1,
                   syph.seco.asympt.prob.tx = 1,
                   syph.earlat.asympt.prob.tx = 1,
                   syph.latelat.asympt.prob.tx = 1,
                   syph.tert.asympt.prob.tx = 1,

                   hivdx.syph.sympt.tx.rr = 1.45,

                   prep.coverage = 0.0,
                   ept.coverage = 0.0,
                   stianntest.coverage = 0.5,
                   stihighrisktest.coverage = 0.8,

                   prep.start = 5000,
                   stitest.start = 250,
                   ept.start = 5000,

                   stitest.elig.model = "sti",

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
                 prev.syph.B = 0.015,
                 prev.syph.W = 0.015)

control <- control_msm(nsteps = 2600)

sim <- netsim(est, param, init, control)

at <- 1
dat <- initialize_msm(est, param, init, control, s = 1)
# dat <- reinit_msm(sim, param, init, control, s = 1)

# debug(simnet_msm)

# at <- at + 1
for (at in 2:control$nsteps) {
  dat <- aging_msm(dat, at)
  dat <- deaths_msm(dat, at)
  dat <- births_msm(dat, at)
  dat <- hiv_test_msm(dat, at)
  dat <- sti_test_msm(dat, at)
  dat <- hiv_tx_msm(dat, at)
  dat <- prep_msm(dat, at)
  dat <- sti_ept_msm(dat, at)
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
  dat <- riskhist_ept_msm(dat, at)
  dat <- position_msm(dat, at)
  dat <- hiv_trans_msm(dat, at)
  dat <- sti_trans_msm(dat, at)
  dat <- sti_recov_msm(dat, at)
  dat <- sti_tx_msm(dat, at)
  dat <- prevalence_msm(dat, at)
  cat("\t", at)
}
