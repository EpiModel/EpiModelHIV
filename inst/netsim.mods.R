
at <- 1
dat <- initialize.msm(est, param, init, control, s = 1)
# dat <- reinit.msm(sim, param, init, control, s = 1)

at <- at + 1
dat <- aging.msm(dat, at)
dat <- deaths.msm(dat, at)
dat <- births.msm(dat, at)
dat <- test.msm(dat, at)
dat <- tx.msm(dat, at)
dat <- prep.msm(dat, at)
dat <- progress.msm(dat, at)
dat <- update_vl.msm(dat, at)
dat <- update_aiclass.msm(dat, at)
dat <- update_roleclass.msm(dat, at)
dat <- edges_correct.msm(dat, at)
dat <- updatenwp.msm(dat, at)
dat <- simnet.msm(dat, at)
dat <- disclose.msm(dat, at)
dat <- acts.msm(dat, at)
dat <- condoms.msm(dat, at)
dat <- riskhist.msm(dat, at)
dat <- position.msm(dat, at)
dat <- trans.msm(dat, at)
dat <- prevalence.msm(dat, at)
verbose.msm(dat, type = "progress", s = 1, at)

