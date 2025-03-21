# 1. Core Unit Tests ------------------------------------------------------

test_that("Basic parameterization", {
  param <- param_msm(
    netstats          = netstats,
    epistats          = epistats,
    prep.start        = 54,
    riskh.start       = 2
  )

  init <- init_msm(
    prev.ugc = 0.1,
    prev.rct = 0.1,
    prev.rgc = 0.1,
    prev.uct = 0.1
  )

  control <- control_msm(
    nsteps              = 52 * 2, # two years
    nsims               = 1,
    ncores              = 1,
    cumulative.edgelist = TRUE,
    truncate.el.cuml    = 0,
    verbose             = FALSE,
    raw.output          = FALSE,
    save.nwstats        = TRUE
  )

  sim <- netsim(est, param, init, control)
})

test_that("Cumulative edgelist", {
  param <- param_msm(
    netstats          = netstats,
    epistats          = epistats,
    prep.start        = 54,
    riskh.start       = 2
  )

  init <- init_msm(
    prev.ugc = 0.1,
    prev.rct = 0.1,
    prev.rgc = 0.1,
    prev.uct = 0.1
  )

  control <- control_msm(
    nsteps              = 10,
    nsims               = 1,
    ncores              = 1,
    cumulative.edgelist = TRUE,
    truncate.el.cuml    = Inf,
    verbose             = FALSE,
    raw.output          = FALSE
  )
  sim <- netsim(est, param, init, control)

  expect_equal(min(sim$el.cuml$sim1[[3]]$start, na.rm = TRUE), 1)
  expect_equal(max(sim$el.cuml$sim1[[3]]$stop, na.rm = TRUE), 9)
})

test_that("Full prep coverage", {
  param <- param_msm(
    netstats          = netstats,
    epistats          = epistats,
    prep.start        = 54,
    riskh.start       = 2,
    prep.require.lnt  = FALSE,
    prep.discont.int  = rep(0, 3),
    prep.start.prob   = rep(1, 3)
  )

  init <- init_msm(
    prev.ugc = 0.1,
    prev.rct = 0.1,
    prev.rgc = 0.1,
    prev.uct = 0.1
  )

  control <- control_msm(
    nsteps              = 104,
    nsims               = 1,
    ncores              = 1,
    cumulative.edgelist = TRUE,
    truncate.el.cuml    = Inf,
    verbose             = FALSE,
    raw.output          = FALSE
  )

  sim <- netsim(est, param, init, control)

  d_sim <- as.data.frame(sim)
  # >= as eligibility is assessed at each step but stoppage if not eligible
  # occurs on a yearly basis
  expect_true(all(d_sim$prepCurr >= d_sim$prepElig, na.rm = TRUE))
})

test_that("Full diagnostics", {
  param <- param_msm(
    netstats           = netstats,
    epistats           = epistats,
    prep.start         = 54,
    riskh.start        = 2,
    hiv.test.rate = rep(1, 3),
    test.window.int = 0
  )

  init <- init_msm(
    prev.ugc = 0.1,
    prev.rct = 0.1,
    prev.rgc = 0.1,
    prev.uct = 0.1
  )

  control <- control_msm(
    nsteps              = 54,
    nsims               = 1,
    ncores              = 1,
    cumulative.edgelist = TRUE,
    truncate.el.cuml    = Inf,
    verbose             = FALSE,
    raw.output          = FALSE
  )

  sim <- netsim(est, param, init, control)

  d_sim <- as.data.frame(sim)
  # Not exactly equal as transmission occurs after testing in the module order
  expect_true(all(d_sim$i.prev - d_sim$i.prev.dx < 0.05, na.rm = TRUE))


  # test restart
  control <- control_msm(
    start               = 55,
    nsteps              = 56,
    nsims               = 1,
    ncores              = 1,
    cumulative.edgelist = TRUE,
    truncate.el.cuml    = Inf,
    initialize.FUN      = reinit_msm,
    verbose             = FALSE,
    raw.output          = FALSE
  )

  sim2 <- netsim(sim, param, init, control)


})
