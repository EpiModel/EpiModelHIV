rm(list = ls())
suppressMessages(library("EpiModelHIV"))

data(est)
data(st)

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
                       save.nwstats = TRUE,
                       verbose.int = 1)
# sim <- netsim(est, param, init, control)

# debug(stergm_prep)


at <- 1
dat <- initialize_msm(est, param, init, control, s = 1)
# dat <- reinit_msm(sim, param, init, control, s = 1)

# mf <- dat$p[[1]]$model.form
# mf$terms[[4]]


at <- at + 1
dat <- aging_msm(dat, at)       ## <1 ms
dat <- deaths_msm(dat, at)      ## 4 ms
dat <- births_msm(dat, at)      ## 6 ms
dat <- test_msm(dat, at)        ## 2 ms
dat <- tx_msm(dat, at)          ## 3 ms
dat <- prep_msm(dat, at)        ## 2 ms
dat <- progress_msm(dat, at)    ## 2 ms
dat <- vl_msm(dat, at)          ## 3 ms
dat <- simnet_msm(dat, at)      ## 53 ms
dat <- disclose_msm(dat, at)    ## 1 ms
dat <- acts_msm(dat, at)        ## 1 ms
dat <- condoms_msm(dat, at)     ## 2 ms
dat <- riskhist_msm(dat, at)    ## 4 ms
dat <- position_msm(dat, at)    ## 1 ms
dat <- trans_msm(dat, at)       ## 1 ms
dat <- sti_trans(dat, at)       ## 4 ms
dat <- sti_recov(dat, at)       ## 3 ms
dat <- sti_tx(dat, at)          ## 2 ms
dat <- prevalence_msm(dat, at)  ## 1 ms
