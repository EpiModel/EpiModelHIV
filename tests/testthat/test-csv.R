#
# Most inputs are defined by the "testthat/setup.R" file
#

test_that("Default Model Parameters", {
  nsteps  <- year_steps * 2
  control <- control_msm(
    nsteps = year_steps * 2,
    verbose = FALSE
  )

  sim <- netsim(netest, param, init, control)

  control <- control_msm(
    start               = nsteps + 1,
    nsteps              = nsteps * 2 + 1,
    initialize.FUN      = reinit_msm,
    verbose = FALSE
  )

  sim2 <- netsim(sim, param, init, control)
})

test_that("Cumulative edgelist", {
  n_steps <- 10
  control <- control_msm(
    nsteps = n_steps,
    cumulative.edgelist = TRUE,
    truncate.el.cuml = Inf,
    save.cumulative.edgelist = TRUE,
    verbose = FALSE
  )

  sim <- netsim(netest, param, init, control)

  # Check that the cumulative edgelist are there
  #  And that start and stop are the correct values (2, nsteps -1)
  expect_equal(min(sim$cumulative.edgelist$sim1$start, na.rm = TRUE), 2)
  expect_equal(max(sim$cumulative.edgelist$sim1$stop, na.rm = TRUE), n_steps - 1)
})

test_that("Full diagnostics", {
  skip("stochastic, prone to failure")
  param_fd <- param
  param_fd$hiv.test.rate <- rep(1, 3)
  param_fd$test.window.int <- 0

  nsteps  <- year_steps * 2
  control <- control_msm(
    nsteps = year_steps * 2,
    verbose = FALSE
  )

  sim <- netsim(netest, param_fd, init, control)
  d_sim <- as.data.frame(sim)
  # Not exactly equal as transmission occurs after testing in the module order
  expect_true(all(d_sim$i.prev - d_sim$i.prev.dx < 0.05, na.rm = TRUE))
})

test_that("Full prep coverage", {
  skip("stochastic, prone to failure")
  param_pc <- param
  param_pc$prep.discont.int <- rep(0, 3)
  param_pc$prep.start.rate <- rep(1, 3)

  nsteps  <- year_steps * 2
  control <- control_msm(
    nsteps = year_steps * 2,
    verbose = FALSE
  )

  sim <- netsim(netest, param_pc, init, control)
  d_sim <- as.data.frame(sim)

  # >= as eligibility is assessed at each step but stoppage if not eligible
  # occurs on a yearly basis
  expect_true(all(d_sim$prepCurr >= d_sim$prep.indic, na.rm = TRUE))
})
