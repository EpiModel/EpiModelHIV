context("Model Runs")

test_that("Burnin model", {

  data(st)
  data(est)

  param <- param_msm(nwstats = st)
  init <- init_msm(nwstats = st)
  control <- control_msm(simno = 1, nsteps = 10, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_is(sim, "netsim")

})

