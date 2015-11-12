
at <- 1
dat <- initialize.mard(est, param, init, control, s = 1)
# dat <- reinit.mard(sim, param, init, control, s = 1)

at <- at + 1
dat <- aging.mard(dat, at)
dat <- deaths.mard(dat, at)
dat <- births.mard(dat, at)
dat <- test.mard(dat, at)
dat <- tx.mard(dat, at)
dat <- prep.mard(dat, at)
dat <- progress.mard(dat, at)
dat <- update_vl.mard(dat, at)
dat <- update_aiclass.mard(dat, at)
dat <- update_roleclass.mard(dat, at)
dat <- edges_correct.mard(dat, at)
dat <- updatenwp.mard(dat, at)
dat <- simnet.mard(dat, at)
dat <- disclose.mard(dat, at)
dat <- acts.mard(dat, at)
dat <- condoms.mard(dat, at)
dat <- riskhist.mard(dat, at)
dat <- position.mard(dat, at)
dat <- trans.mard(dat, at)
dat <- prevalence.mard(dat, at)
verbose.mard(dat, type = "progress", s = 1, at)

