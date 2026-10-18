# Setup the estimation object before any tests
library(ARTnet)
library(EpiModelHIV)

# Relevant time steps for the simulation
time_unit  <- 7               # number of days in a time step
year_steps <- 364 / time_unit # number of time steps in a year

networks_size   <- 5 * 1e3
control_ergm <- ergm::control.ergm(
  main.method = "Stochastic-Approximation",
  MCMLE.maxit = 500,
  SAN.maxit = 3,
  SAN.nsteps.times = 4,
  MCMC.samplesize = 1e4,
  MCMC.interval = 5e3,
  parallel = 1
)

## Initialize the ARTnet data objects and the networks to be fitted
epistats <- build_epistats(
  geog.lvl = "city",
  geog.cat = "Atlanta",
  init.hiv.prev = c(0.33, 0.137, 0.084),
  race = TRUE,
  time.unit = time_unit
)

netparams <- build_netparams(
  epistats = epistats,
  smooth.main.dur = TRUE
)

netstats <- build_netstats(
  epistats,
  netparams,
  expect.mort = 0.000478213,
  network.size = networks_size
)

nw <- EpiModel::network_initialize(netstats$demog$num)
nw_main <- EpiModel::set_vertex_attribute(
  nw,
  names(netstats$attr),
  netstats$attr
)

nw_casl <- nw_main
nw_inst <- nw_main

## Define and fit the *main* network  model ------------------------------------
model_main <- ~ edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = -1) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("deg.casl", levels = -1) +
  concurrent +
  degrange(from = 3) +
  nodematch("role.class", diff = TRUE, levels = c(1, 2))

# Target Stats
netstats_main <- c(
  edges                = netstats$main$edges,
  nodematch_age.grp    = netstats$main$nodematch_age.grp,
  nodefactor_age.grp   = netstats$main$nodefactor_age.grp[-1],
  nodematch_race       = netstats$main$nodematch_race_diffF,
  nodefactor_race      = netstats$main$nodefactor_race[-1],
  nodefactor_deg.casl  = netstats$main$nodefactor_deg.casl[-1],
  concurrent           = netstats$main$concurrent,
  degrange             = 0,
  nodematch_role.class = c(0, 0)
) |> unname()

# Fit model
fit_main <- EpiModel::netest(
  nw_main,
  formation = model_main,
  target.stats = netstats_main,
  coef.diss = netstats$main$diss.byage,
  set.control.ergm = control_ergm,
  verbose = FALSE
) |> EpiModel::trim_netest()

# Keep only the necessary objects
rm(model_main, netstats_main)

## Define and fit the *casual* network  model ----------------------------------
model_casl <- ~ edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = -5) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("deg.main", levels = -3) +
  concurrent +
  degrange(from = 4) +
  nodematch("role.class", diff = TRUE, levels = c(1, 2))

# Target Stats
netstats_casl <- c(
  edges                = netstats$casl$edges,
  nodematch_age.grp    = netstats$casl$nodematch_age.grp,
  nodefactor_age.grp   = netstats$casl$nodefactor_age.grp[-5],
  nodematch_race       = netstats$casl$nodematch_race_diffF,
  nodefactor_race      = netstats$casl$nodefactor_race[-1],
  nodefactor_deg.main  = netstats$casl$nodefactor_deg.main[-3],
  concurrent           = netstats$casl$concurrent,
  degrange             = 0,
  nodematch_role.class = c(0, 0)
) |> unname()

# Fit model
fit_casl <- EpiModel::netest(
  nw_casl,
  formation = model_casl,
  target.stats = netstats_casl,
  coef.diss = netstats$casl$diss.byage,
  set.control.ergm = control_ergm,
  verbose = FALSE
) |> trim_netest()

# Keep only the necessary objects
rm(model_casl, netstats_casl)

## Define and fit the *one-off* network  model ---------------------------------

# Formula
model_ooff <- ~ edges +
  nodematch("age.grp", diff = FALSE) +
  nodefactor("age.grp", levels = -1) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("risk.grp", levels = -5) +
  nodefactor("deg.tot", levels = -1) +
  nodematch("role.class", diff = TRUE, levels = c(1, 2))

# Target Stats
netstats_ooff <- c(
  edges                = netstats$inst$edges,
  nodematch_age.grp    = sum(netstats$inst$nodematch_age.grp),
  nodefactor_age.grp   = netstats$inst$nodefactor_age.grp[-1],
  nodematch_race       = netstats$inst$nodematch_race_diffF,
  nodefactor_race      = netstats$inst$nodefactor_race[-1],
  nodefactor_risk.grp  = netstats$inst$nodefactor_risk.grp[-5],
  nodefactor_deg.tot   = netstats$inst$nodefactor_deg.tot[-1],
  nodematch_role.class = c(0, 0)
) |> unname()

# Fit model
fit_ooff <- EpiModel::netest(
  nw_inst,
  formation = model_ooff,
  target.stats = netstats_ooff,
  coef.diss = dissolution_coefs(~ offset(edges), 1),
  set.control.ergm = control_ergm,
  verbose = FALSE
) |> trim_netest()

# Keep only the necessary objects
rm(model_ooff, netstats_ooff)

# format and trim the data -----------------------------------------------------
netest <- list(fit_main = fit_main, fit_casl = fit_casl, fit_ooff = fit_ooff)
# Reduce the size of netstats and epistats
netstats <- ARTnet::trim_netstats(netstats)
epistats <- ARTnet::trim_epistats(epistats)

# Code to run after all the tests
withr::defer(
  {
    rm(list = ls()) # placeholder, not useful
  },
  teardown_env()
)
