context("Model Runs")

test_that("Burnin model", {

  data(st)
  data(est)

  param <- param.mard(nwstats = st)
  init <- init.mard(nwstats = st)
  control <- control.mard(simno = 1, nsteps = 10, verbose = FALSE)

  sim <- netsim(est, param, init, control)
  expect_is(sim, "netsim")

})

test_that("Follow-up model", {

  data(st)
  data(est)

  param <- param.mard(nwstats = st,
                     prep.start = 10)
  init <- init.mard(nwstats = st)
  control <- control.mard(simno = 1, nsteps = 10,
                         save.other = c("attr", "temp", "riskh", "el", "p"),
                         verbose = FALSE)

  sim <- netsim(est, param, init, control)

  param <- param.mard(nwstats = st,
                     prep.start = 10,
                     prep.elig.model = "cdc3",
                     prep.coverage = 0.5,
                     prep.risk.int = 182,
                     prep.class.prob = reallocate_pcp(reall = 0),
                     prep.class.hr = c(1, 0.69, 0.19, 0.05))
  init <- init.mard(nwstats = st)
  control <- control.mard(simno = 1, start = 11, nsteps = 20,
                         verbose = FALSE, initialize.FUN = reinit_msm)

  sim2 <- netsim(sim, param, init, control)
  expect_is(sim2, "netsim")

})
