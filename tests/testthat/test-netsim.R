context("Full netsim Simulation")

data(est)
data(st)

test_that("Testing netsim", {

  param <- param_msm(nwstats = st)
  init <- init_msm(nwstats = st)
  control <- control_msm(nsteps = 5, verbose = FALSE)

  sim <- netsim(est, param, init, control)

  # expect this output on sim
  nm <- c("param", "control", "nwparam", "epi", "stats", "attr", "temp",
          "el", "p")
  expect_identical(names(sim), nm)
  expect_is(sim, "netsim")

})
