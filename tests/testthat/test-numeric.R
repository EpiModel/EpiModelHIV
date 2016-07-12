context("Numeric attributes for C++ code")

rm(list=ls())
devtools::install_github("statnet/EpiModelHIV")
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
                       nsims = 1, 
                       ncores = 1, 
                       save.nwstats = TRUE,
                       verbose.int = 1)
debug(EpiModelHIV:::init_status_msm)
sim <- netsim(est, param, init, control)


test_that("Numeric for ", {
    
    NULL
})