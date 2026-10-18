# nolint start
# This file contains the Partner Identification and Services module
#
#
# In order to use this module you need to:
#
# - activate it:
# by adding `partident.FUN = partident_msm` to the list of
# modules in `control_msm`. This module should be added at the end, just
# before the `prev.FUN` and `verbose.FUN` modules.
#
# - Setup the attributes it needs:
# add the following elements to the `get_default_attrs` list:
#
# ```r
#     part.scrnd = -Inf,
#     part.ident = -Inf,
#     part.ident.counter = 0,
# ```
#
# - Setup the required parameters:
# here are their CSV format you could add to the `model_parameters.csv` file
#
# ```csv
# part.ident.start,Inf,numeric,,,,,Time step to initiate the partner identification process. Default of `Inf` indicates no partner identification.
# part.index.window.int,0,numeric,,,,,Number of prior time steps for which a newly diagnosed index patient is eligible for partner notification. Default of 0 indicates index eligible only if newly diagnosed at current time step.
# part.index.degree,1,numeric,,,,,The minimum number of partnerships an index patient may have for partner elicitation. Default of `1` indicating all partnerships under consideration.
# part.index.prob,0.666,numeric,,,,,Probability that an index patient will initiate partner services.
# part.ident.main.window.int,52,numeric,,,,,The number of time steps in the past that a main partnership qualifies for partner identification.
# part.ident.casl.window.int,52,numeric,,,,,The number of time steps in the past that a casual partnership qualifies for partner identification.
# part.ident.ooff.window.int,52,numeric,,,,,The number of time steps in the past that a one-off partnership qualifies for partner identification.
# part.ident.main.prob,0.157,numeric,,,,,Probability that a main partner is identified through partner identification.
# part.ident.casl.prob,0.08518719,numeric,,,,,Probability that a casual partner is identified through partner identification.
# part.ident.ooff.prob,0.04448852,numeric,,,,,Probability that a one-time partner is identified through partner identification.
# part.hiv.test.rate_1,0.394,numeric,,,,,Mean probability per time step of HIV testing for a Black MSM who has been identified as the partner of an index HIV+ MSM.
# part.hiv.test.rate_2,0.394,numeric,,,,,Mean probability per time step of HIV testing for a Hispanic MSM who has been identified as the partner of an index HIV+ MSM.
# part.hiv.test.rate_3,0.394,numeric,,,,,Mean probability per time step of HIV testing for a White MSM who has been identified as the partner of an index HIV+ MSM.
# part.prep.start.prob_1,0,numeric,,,,,Probability of a Black MSM identified through partner identification starting PrEP given current indications at the current time step.
# part.prep.start.prob_2,0,numeric,,,,,Probability of a Hispanic MSM identified through partner identification starting PrEP given current indications at the current time step.
# part.prep.start.prob_3,0,numeric,,,,,Probability of a White MSM identified through partner identification starting PrEP given current indications at the current time step.
# part.tx.init.rate_1,0.387,numeric,,,,,Probability that a Black MSM who has been identified as the partner of an incident HIV+ MSM will initiate treatment during the current time step.
# part.tx.init.rate_2,0.387,numeric,,,,,Probability that a Hispanic MSM who has been identified as the partner of an incident HIV+ MSM will initiate treatment during the current time step.
# part.tx.init.rate_3,0.387,numeric,,,,,Probability that a White MSM who has been identified as the partner of an incident HIV+ MSM will initiate treatment during the current time step.
# part.tx.reinit.rate_1,0,numeric,,,,,"Probability per time step that a Black MSM who has been identified through partner identification, stopped treatment will restart treatment."
# part.tx.reinit.rate_2,0,numeric,,,,,"Probability per time step that a Hispanic MSM who has been identified through partner identification, stopped treatment will restart treatment."
# part.tx.reinit.rate_3,0,numeric,11,,,,"Probability per time step that a White MSM who has been identified through partner identification, stopped treatment will restart treatment."
# ```
# nolint end


#' @title Partner Identification Module
#'
#' @description Module function for identifying partners of incident HIV+ MSM.
#' This modules also runs the partner services associated. These are independant
#' of the other systems. (e.g. getting screened through partner services does
#' not prevent you from being screened through the normal process on the same
#' step)
#'
#' @inheritParams aging_msm
#'
#' @details
#' This partner identification module handles the elicitation of partners of
#' incident HIV+ MSM for the purposes of HIV testing ([`hivtest_msm`]).
#'
#' Attributes modified:
#'   - `part.ident`: time of identification, NA if not
#'   - `part.ident.counter`: number of time identified
#' Epi trackers:
#'   - elig_indexes
#'   - found_indexes
#'   - elig_partners
#'   - found_partners
#'
#' @return
#' This function returns the `dat` object with at updated `part.ident` attribute detailing time of
#' identification.
#'
#'
#' @export
#'
partident_msm <- function(dat, at) {
  if (at < get_param(dat, "part.ident.start")) {
    return(dat)
  }

  # Input
  ## Attributes
  diag.status <- get_attr(dat, "diag.status")
  diag.time   <- get_attr(dat, "diag.time")

  # Parameters
  part.index.window.int    <- get_param(dat, "part.index.window.int")
  part.index.degree        <- get_param(dat, "part.index.degree")
  part.index.prob          <- get_param(dat, "part.index.prob")

  # Concatenate Partnership Parameters
  part.ident.int <- c(
    get_param(dat, "part.ident.main.window.int"),
    get_param(dat, "part.ident.casl.window.int"),
    get_param(dat, "part.ident.ooff.window.int")
  )
  part.ident.prob <- c(
    get_param(dat, "part.ident.main.prob"),
    get_param(dat, "part.ident.casl.prob"),
    get_param(dat, "part.ident.ooff.prob")
  )

  ## Process

  # Who was diagnosed last step. Identification occurs before testing
  previous_step <- at - 1

  # Identify individuals with recent HIV+ diagnoses
  recent_diag <- previous_step - diag.time

  # set the epi to 0, they will be overwritten if needed. (prevent NAs in early
  # return cases)
  dat <- set_epi(dat, "elig_indexes", at, 0)
  dat <- set_epi(dat, "found_indexes", at, 0)
  dat <- set_epi(dat, "elig_partners", at, 0)
  dat <- set_epi(dat, "found_partners", at, 0)

  hivpos.at <- which(diag.status == 1 & recent_diag <= part.index.window.int)
  if (length(hivpos.at) == 0) {
    return(dat)
  }

  # Partnerships involving an HIV positive as index
  hivpos_pel <- get_partners(dat, hivpos.at, only.active.nodes = TRUE)

  # Remove partnerships of index patients with mean degree < part.index.degree
  keep_ids <- which(table(hivpos_pel$index) >= part.index.degree)
  keep_ids <- names(keep_ids)
  hivpos_uid <- intersect(keep_ids, hivpos_pel$index)

  # Each index is found with probability `part.index.prob`
  found_uid <- hivpos_uid[runif(length(hivpos_uid)) < part.index.prob]
  hivpos_pel <- hivpos_pel[hivpos_pel$index %in% found_uid, ]

  # Set trackers for the indexes
  dat <- set_epi(dat, "elig_indexes", at, length(hivpos_uid))
  dat <- set_epi(dat, "found_indexes", at, length(found_uid))

  if (nrow(hivpos_pel) == 0) {
    return(dat)
  }

  # Subset plist.ident by partners who are in a qualified relationship,
  # defined by relationship time ("stop" is NA for ongoing relationships)
  relationship_age <- ifelse(is.na(hivpos_pel$stop), 0, at - hivpos_pel$stop)
  age_thresh <- part.ident.int[hivpos_pel$network]
  relationship_elig <- relationship_age <= age_thresh
  hivpos_pel <- hivpos_pel[relationship_elig, ]

  if (nrow(hivpos_pel) == 0) {
    return(dat)
  }

  # Set trackers for the eligible partners
  elig_partners <- unique(hivpos_pel$partner)
  dat <- set_epi(dat, "elig_partners", at, length(elig_partners))

  if (length(elig_partners) == 0) {
    return(dat)
  }

  ## Attempt to identify eligible partners
  part.ident.prob.type <- part.ident.prob[hivpos_pel$network]
  found_relationship <- runif(nrow(hivpos_pel)) < part.ident.prob.type
  hivpos_pel <- hivpos_pel[found_relationship, ]

  identified_ids <- get_posit_ids(dat, unique(hivpos_pel$partner))

  dat <- set_epi(dat, "found_partners", at, length(identified_ids))

  # Update attributes
  dat <- set_attr(dat, "part.ident", at, posit_ids = identified_ids)
  counter <- get_attr(dat, "part.ident.counter", posit_ids = identified_ids)
  dat <- set_attr(dat, "part.ident.counter", counter + 1, posit_ids = identified_ids)

  dat <- partident_hivtest(dat)
  dat <- partident_hivtx(dat)
  dat <- partident_prep(dat)

  return(dat)
}

partident_hivtest <- function(dat) {
  # Attributes
  status      <- get_attr(dat, "status")
  diag.status <- get_attr(dat, "diag.status")
  part.ident  <- get_attr(dat, "part.ident")
  race        <- get_attr(dat, "race")
  # Parameters
  part.hiv.test.rate <- get_param(dat, "part.hiv.test.rate")
  at <- get_current_timestep(dat)

  ## Process
  elig_ids <- which(part.ident == at & diag.status == 0)
  prob.screen <- part.hiv.test.rate[race[elig_ids]]
  screened_ids <- elig_ids[runif(length(elig_ids)) < prob.screen]

  ## Outputs
  dat <- set_attr(dat, "part.scrnd", at, posit_ids = screened_ids)
  dat <- set_pos_hiv_test_attrs(dat, screened_ids[status[screened_ids] == 1])
  dat <- set_neg_hiv_test_attrs(dat, screened_ids[status[screened_ids] == 0])

  return(dat)
}

partident_hivtx <- function(dat) {
  # Attributes
  diag.status      <- get_attr(dat, "diag.status")
  tx.status        <- get_attr(dat, "tx.status")
  cuml.time.on.tx  <- get_attr(dat, "cuml.time.on.tx")
  part.ident       <- get_attr(dat, "part.ident")
  race             <- get_attr(dat, "race")
  # Parameters
  part.tx.init.rate   <- get_param(dat, "part.tx.init.rate")
  part.tx.reinit.rate <- get_param(dat, "part.tx.reinit.rate")
  at <- get_current_timestep(dat)

  tx_elig <- part.ident == at & tx.status == 0 & diag.status == 1

  # Initiation
  elig_ids <- which(tx_elig & cuml.time.on.tx == 0)

  elig_rates <- part.tx.init.rate[race[elig_ids]]
  tx_start_ids <- elig_ids[runif(length(elig_ids)) < elig_rates]

  # Restarting
  elig_ids <- which(tx_elig & cuml.time.on.tx > 0)
  elig_rates <- part.tx.reinit.rate[race[elig_ids]]
  tx_restart_ids <- elig_ids[runif(length(elig_ids)) < elig_rates]

  ## Outputs
  dat <- set_attr(dat, "part.tx.init.time", at, posit_ids = tx_start_ids)
  dat <- set_start_hiv_tx_attrs(dat, tx_start_ids)
  dat <- set_attr(dat, "part.tx.reinit.time", at, posit_ids = tx_restart_ids)
  dat <- set_restart_hiv_tx_attrs(dat, tx_restart_ids)

  return(dat)
}

partident_prep <- function(dat) {
  # Attributes
  part.scrnd   <- get_attr(dat, "part.scrnd")
  prep.indic   <- get_attr(dat, "prep.indic")
  prep         <- get_attr(dat, "prep")
  race         <- get_attr(dat, "race")

  # Parameters
  part.prep.start.rate <- get_param(dat, "part.prep.start.rate")
  at <- get_current_timestep(dat)

  ## Process
  elig_ids <- which(part.scrnd == at & prep == 0 & prep.indic == 1)
  elig_rates <- part.prep.start.rate[race[elig_ids]]
  prep_start_ids <- elig_ids[runif(length(elig_ids)) < elig_rates]

  ## Outputs
  dat <- set_epi(dat, "part.prep.start", at, length(prep_start_ids))
  dat <- set_start_prep_attrs(dat, prep_start_ids)
  dat <- set_attr(dat, "prep.start.part", at, posit_ids = prep_start_ids)

  return(dat)
}
