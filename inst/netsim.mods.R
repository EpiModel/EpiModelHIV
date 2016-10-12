rm(list=ls())
suppressMessages(library("EpiModelHIV"))
sourceDir("R/")

data(est)
data(st)
est
st
param <- param_msm(nwstats = st, 
                   ai.scale = 1.323,
                   prep.coverage = 0)
init <- init_msm(nwstats = st, 
                 prev.B = 0.253, 
                 prev.W = 0.253)
control <- control_msm(simno = 0.253, 
                       nsteps = 52,
                       nsims = 5, 
                       ncores = 1,
                       initialize.FUN = initialize_msm, aging.FUN = aging_msm,
                       deaths.FUN = deaths_msm, births.FUN = births_msm, test.FUN = test_msm,
                       tx.FUN = tx_msm, prep.FUN = prep_msm, progress.FUN = progress_msm,
                       progresssyph.FUN = progress_syph_msm,
                       vl.FUN = vl_msm, aiclass.FUN = NULL, roleclass.FUN = NULL,
                       resim_nets.FUN = simnet_msm, disclose.FUN = disclose_msm,
                       acts.FUN = acts_msm, condoms.FUN = condoms_msm,
                       riskhist.FUN = riskhist_msm, position.FUN = position_msm,
                       trans.FUN = trans_msm, stitrans.FUN = sti_trans,
                       stirecov.FUN = sti_recov, stitx.FUN = sti_tx, prev.FUN = prevalence_msm,
                       verbose.FUN = verbose_msm,
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
dat <- tx_msm(dat, at)
dat <- prep_msm(dat, at)
dat <- progress_msm(dat, at)
dat <- progress_syph_msm(dat, at)
dat <- update_vl_msm(dat, at)
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
dat <- prevalence_msm(dat, at)
verbose_msm(dat, type = "progress", s = 1, at)

