test_that("Default Model Parameters", {

  params_df <- read.csv(system.file(
    "model_parameters.csv",
    package = "EpiModelHIV"
  )) |>
  dplyr::select(param, value, type)

  # `netsim` Parameters
  param <- param.net(
    data.frame.params   = params_df,
    netstats            = netstats,
    epistats            = epistats,
    prep.start          = 0,
    riskh.start         = 0
  )

  # Initial conditions (default prevalence initialized in epistats)
  # For models without bacterial STIs, these must be initialized here
  # with non-zero values
  init <- init_msm(
    prev.ugc = 0.1,
    prev.rct = 0.1,
    prev.rgc = 0.1,
    prev.uct = 0.1
  )

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
