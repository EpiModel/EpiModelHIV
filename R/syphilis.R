#' @title Initialize The Syphilis Epidemic
#'
#' @inheritParams aging_msm
#'
#' @details
#' The Syphilis epidemic is currently initialize within the higher sexual risk
#' group. I.e., the nodes with one offs partners. This is to emulate the fact
#' that syphilis is mostly present in a "core groupe" of highly sexually active
#' individuals.
#'
#' @return
#' Returns the `dat` object after the Syphilis epidemic has been initialized.
#'
syph_init <- function(dat) {
  n_new <- length(get_attr(dat, "active"))
  prev.syph <- get_init(dat, "prev.syph")

  # init syph only on higher risk grps
  high_risk_ids <- which(get_attr(dat, "risk.grp") >= 3)

  syph_ids <- sample(high_risk_ids, size = round(prev.syph * n_new), FALSE)
  dat <- set_inf_syph_attrs(dat, syph_ids)
  syph.test.last <- 1 - rgeom(n_new, get_param(dat, "syph.screen.hivneg.rate"))
  dat <- set_attr(dat, "syph.test.last", syph.test.last)
  return(dat)
}

#' @title Handles the Transmission of Syphilis
#'
#' @details
#' Within this model, Syphilis is only infectious during the primary and
#' secondary stage (stages number 2 and 3).
#'
#' @inheritParams aging_msm
#'
#' @return
#' Returns the `dat` object with updated attributes for Syphilis infections.
#'
syphilis_transmission <- function(dat) {
  # Setup ----------------------------------------------------------------------
  at <- get_current_timestep(dat)

  ## Attributes
  race <- get_attr(dat, "race")
  syph.inf <- get_attr(dat, "syph.inf")
  syph.stage <- get_attr(dat, "syph.stage")

  ## Parameters
  syph.prob <- get_param(dat, "syph.prob")
  syph.cond.eff.or <- get_param(dat, "syph.cond.eff.or")
  cond.fail.or <- get_param(dat, "cond.fail.or")

  # Process --------------------------------------------------------------------

  # Pull act list
  al <- dat[["temp"]][["al"]]
  ## ins variable coding
  # ins = 1 : p1 is insertive partner; ins = 0 : p1 is receptive partner

  # syphilis is  only infectious during the primary & secondary stage
  infectious <- syph.stage %in% c(2, 3)

  # acts where p1 can transmit to p2
  p1i_acts <- which(infectious[al[, "p1"]])           # p1 is infectious
  p1i_acts <- p1i_acts[syph.inf[al[p1i_acts, "p2"]] == 0] # p2 not infected
  # acts where p2 can transmit to p1
  p2i_acts <- which(infectious[al[, "p2"]])           # p2 is infectious
  p2i_acts <- p2i_acts[syph.inf[al[p2i_acts, "p1"]] == 0] # p1 not infected

  # Store the acts so that infection always is from p1 to p2
  all_acts <- c(p1i_acts, p2i_acts)
  elig_acts <- matrix(
    c(
      al[p1i_acts, "p1"], al[p2i_acts, "p2"], # first col: infected partner
      al[p1i_acts, "p2"], al[p2i_acts, "p1"], # second col: susceptible partner
      al[all_acts, "uai"], # third col: condomless act flag
      # 4th column: insertive partner
      ifelse(al[all_acts, "ins"] == 1, al[all_acts, "p1"], al[all_acts, "p2"])
    ),
    ncol = 4
  )
  colnames(elig_acts) <- c("infectious", "susceptible", "uai", "ins_partner")

  # probability of infection ---------------------------------------------------
  log_odds <- prob_to_log_odds(rep(syph.prob, nrow(elig_acts)))

  # Condom modifier: UAI == 0 -> a condom is used
  cond_acts <- which(elig_acts[, "uai"] == 0)
  race_ins_partner <- race[elig_acts[cond_acts, "ins_partner"]]
  condom_or <- 1 - (syph.cond.eff.or - cond.fail.or[race_ins_partner])
  log_odds[cond_acts] <- log_odds[cond_acts] + log(condom_or)

  # NOTE: we could add some transmission in the early latent and latent stages.
  # It requires adding the stages 4 and 5 to the `infectious` variables above.
  # And setting an odds ratio for these stages below.

  inf_prob <- log_odds_to_prob(log_odds)

  # Stochastic transmission
  inf_acts <- which(runif(nrow(elig_acts)) < inf_prob)
  new_inf_ids <- elig_acts[inf_acts, "susceptible"]
  new_inf_ids <- unique(new_inf_ids) # deduplicate (multiple acts/partnerships)

  # Output ---------------------------------------------------------------------
  # Update attr
  dat <- set_inf_syph_attrs(dat, new_inf_ids)
  # Summary stats
  dat <- set_epi(dat, "incid.syph", at, length(new_inf_ids))

  return(dat)
}

#' @title Handles the Natural History of Syphilis
#'
#' @details
#' Syphilis natural history is modeled here in 6 stages:
#' 1. exposed: not infectious, not detectable, no visible symptoms
#' 2. primary: infectious, detectable, possibly visible symptoms
#' 3. secondary: infectious, detectable, possibly visible symptoms
#' 4. early latent: not infectious, detectable, no visible symptoms
#' 5. latent: not infectious, detectable, no visible symptoms
#' 6. tertiary: not infectious, detectable, always visible symptoms
#'
#' The transition from one stage to the next is stochastic with probability
#' derived from the mean interval duration of each stage, in the parameter:
#' `syph.prog.int`
#'
#' Syphilis visible symptoms are only possible during the primary and secondary
#' stage for now. This is governed by the `syph.sympt.prob` parameter vector.
#'
#' Recovery from Syphilis in our model is never spontaneous. It always requires
#' treatment. Upon treatment, the infection immediately stops.
#'
#' NOTE: At this point, there is no protective period after recovery.
#'
#' @inheritParams aging_msm
#'
#' @return
#' Returns the `dat` object with updated attributes for Syphilis infections.
#'
syphilis_progression <- function(dat) {
  # Setup ----------------------------------------------------------------------
  at <- get_current_timestep(dat)
  ## Attributes
  syph.stage <- get_attr(dat, "syph.stage") # 1->6: incub then up to tertiary
  syph.stage.last <- get_attr(dat, "syph.stage.last")
  syph.sympt <- get_attr(dat, "syph.sympt")
  syph.tx <- get_attr(dat, "syph.tx")
  ## Parameters
  syph.prog.int <- get_param(dat, "syph.prog.int") # mean duration of each stage
  syph.sympt.prob <- get_param(dat, "syph.sympt.prob") # probability of visible symptom at each stage

  # Process --------------------------------------------------------------------

  # Progression to next stage --------------------------------------------------

  # Stochastic Progression
  syph.prog.rate <- 1 / syph.prog.int
  # with prob syph.prog.rate[stage] except for tertiary
  ids_elig <- which(syph.stage < 6)
  prog_rate <- syph.prog.rate[syph.stage[ids_elig]]
  prog_ids <- ids_elig[runif(length(ids_elig)) < prog_rate]
  syph.stage[prog_ids] <- syph.stage[prog_ids] + 1
  prog_stage <- syph.stage[prog_ids]
  syph.stage.last[prog_ids] <- at

  # assign visible symptoms
  sympt_prob <- syph.sympt.prob[prog_stage]
  syph.sympt[prog_ids] <- as.numeric(runif(length(prog_ids)) < sympt_prob)

  # Syphilis Recovery ----------------------------------------------------
  # through TX
  ids_recov_tx <- which(syph.tx == 1)
  # NOTE: currently, recovery is only possible via treatment
  ids_recov <- ids_recov_tx

  # Output ---------------------------------------------------------------------

  # Update attr
  dat <- set_attr(dat, "syph.stage", syph.stage)
  dat <- set_attr(dat, "syph.stage.last", syph.stage.last)
  dat <- set_attr(dat, "syph.sympt", syph.sympt)
  dat <- set_recov_syph_attrs(dat, ids_recov)

  return(dat)
}

#' @title Handles the Screening for Syphilis
#'
#' @details
#' Screening for Syphilis can occur at every timestep for all
#' nodes not diagnosed with the disease. Different probablity are used for nodes
#' with an HIV positive diagnostic, for nodes on PrEP and for nodes not
#' diagnosed with HIV and not on PrEP.
#'
#' Additionally, the PrEP starters will test on the step they initiate PrEP.
#' (PrEP occurs before STIs)
#'
#' @inheritParams aging_msm
#'
#' @return
#' Returns the `dat` object with updated testing attributes for Syphilis
#'
syphilis_testing <- function(dat) {
  # Setup ----------------------------------------------------------------------
  at <- get_current_timestep(dat)
  ## Attributes
  syph.inf        <- get_attr(dat, "syph.inf")
  syph.dx         <- get_attr(dat, "syph.dx")
  syph.stage      <- get_attr(dat, "syph.stage")
  prep.any        <- get_attr(dat, "prep.any")
  prep            <- get_attr(dat, "prep")
  prep.start.last <- get_attr(dat, "prep.start.last")
  hiv.dx          <- get_attr(dat, "diag.status")

  prep.otc <- get_attr(dat, "prep.otc")
  prep.otc.start.last <- get_attr(dat, "prep.otc.start.last")

  ## Parameters
  screen.hivneg.rate <- get_param(dat, "syph.screen.hivneg.rate")
  screen.hivpos.rate <- get_param(dat, "syph.screen.hivpos.rate")
  screen.prep.rate   <- get_param(dat, "sti.screen.prep.rate")

  screen.prep.otc.rate <- get_param(dat, "sti.screen.prep.otc.rate")
  prep.otc.always.sti.tst <- get_param(dat, "prep.otc.always.sti.tst")

  # Process --------------------------------------------------------------------
  screen_ids <- vector(mode = "list", 4)
  names(screen_ids) <- c("nprep", "prep", "hivpos", "prep.otc")

  # HIV Negative - not on PrEP
  elig_ids <- which(syph.dx == 0 & hiv.dx == 0 & prep.any == 0)
  screen_ids$nprep <- elig_ids[runif(length(elig_ids)) < screen.hivneg.rate]

  # HIV Negative - on PrEP
  elig_ids <- which(syph.dx == 0 & prep == 1)
  screen_ids$prep <- elig_ids[runif(length(elig_ids)) < screen.prep.rate]
  # unconditionnaly add those who started PrEP this step
  prep_start_ids <- elig_ids[prep.start.last[elig_ids] == at]
  screen_ids$prep <- unique(c(screen_ids$prep, prep_start_ids))

  # HIV Positive
  elig_ids <- which(syph.dx == 0 & hiv.dx == 1)
  screen_ids$hivpos <- elig_ids[runif(length(elig_ids)) < screen.hivpos.rate]

  # HIV Negative - on PrEP OTC
  elig_ids <- which(syph.dx == 0 & prep.otc == 1)
  screen_ids$prep.otc <- elig_ids[runif(length(elig_ids)) < screen.prep.otc.rate]
  if (prep.otc.always.sti.tst == 1) {
    prep_otc_start_ids <- elig_ids[prep.otc.start.last[elig_ids] == at]
    screen_ids$prep.otc <- unique(c(screen_ids$prep.otc, prep_otc_start_ids))
  }

  ## Outputs -------------------------------------------------------------------
  # Simple combine as they are exclusive
  final_screen_ids <- unlist(screen_ids)
  dat <- set_attr(dat, "syph.test.last", at, final_screen_ids)
  # diagnosed are: screened && infected at the correct site
  dx_ids <- final_screen_ids[syph.inf[final_screen_ids] == 1]
  dx_ids <- dx_ids[syph.stage[dx_ids] != 1] # no dx possible in `exposed`
  dat <- set_attr(dat, "syph.dx", 1, dx_ids)
  dat <- set_attr(dat, "syph.dx.last", at, dx_ids)

  if (get_control(dat, "debug")) {
    for (grp in names(screen_ids)) {
      dat <- set_epi(dat, paste0("dbg_syph_tst_", grp), at,
                     length(screen_ids[[grp]]))
    }
    dat <- set_epi(dat, "dbg_syph_dx", at, length(dx_ids))
  }

  return(dat)
}

#' @title Handles the Treatment Attribution for Syphilis
#'
#' @details
#' Treatment attribution for Syphilis can happend via two
#' non-exclusive mecanismes.
#'
#' Symptomatic treatment is a recurring process. A symptomatic individual have a
#' weekly probability of getting treated `syph.sympt.tx.prob`.
#'
#' Asymptomatic treatment is a one time process. When a node is diagnosed
#' with Syphilis it is automatically treated. Therefore, upon diagnostic a node
#' is treated on the same step as the testing part happens before treatment.
#'
#' @inheritParams aging_msm
#'
#' @return
#' Returns the `dat` object with updated treatment attributes for Syphilis
#'
syphilis_treatment <- function(dat) {
  # Setup ----------------------------------------------------------------------
  at <- get_current_timestep(dat)
  ## Attributes
  syph.sympt    <- get_attr(dat, "syph.sympt")
  syph.dx.last  <- get_attr(dat, "syph.dx.last")
  race          <- get_attr(dat, "race")
  ## Parameters
  syph.sympt.tx.prob  <- get_param(dat, "syph.sympt.tx.prob")

  # Process --------------------------------------------------------------------

  ## Symptomatic treatment eligible
  elig_ids <- which(syph.sympt == 1)
  rates <- syph.sympt.tx.prob[race[elig_ids]]
  sympt_tx_ids <- elig_ids[runif(length(elig_ids)) < rates]

  ## Asymptomatic diagnosed this step -> automatic treatment
  dx_ids <- which(syph.dx.last == at)

  ## Combine all of the treated (potential overlap)
  tx_ids <- unique(c(sympt_tx_ids, dx_ids))

  # Output ---------------------------------------------------------------------
  dat <- set_attr(dat, "syph.tx",        1, tx_ids)
  dat <- set_attr(dat, "syph.dx",        1, tx_ids)
  dat <- set_attr(dat, "syph.test.last", at, tx_ids)
  dat <- set_attr(dat, "syph.dx.last",   at, tx_ids)

  return(dat)
}

#' Utility to Set The Default Attributes For Newly Syphilis Infected
#'
#' @param ids Positional IDs of the newly infected
#'
#' @inheritParams aging_msm
#'
#' @noRd
set_inf_syph_attrs <- function(dat, ids) {
  at <- get_current_timestep(dat)
  dat <- set_attr(dat, "syph.inf", 1, ids)
  dat <- set_attr(dat, "syph.stage.last", at, ids)
  dat <- set_attr(dat, "syph.inf.last", at, ids)
  inf.count <- get_attr(dat, "syph.inf.count", ids)
  dat <- set_attr(dat, "syph.inf.count", inf.count + 1, ids)
  dat <- set_attr(dat, "syph.stage", 1, ids)
  return(dat)
}

#' Utility to Reset The Default Attributes For Newly Syphilis Recovered
#'
#' @param ids Positional IDs of the newly infected
#'
#' @inheritParams aging_msm
#'
#' @noRd
set_recov_syph_attrs <- function(dat, ids) {
  dat <- set_attr(dat, "syph.inf", 0, ids)
  dat <- set_attr(dat, "syph.sympt", 0, ids)
  dat <- set_attr(dat, "syph.dx", 0, ids)
  dat <- set_attr(dat, "syph.tx", 0, ids)
  dat <- set_attr(dat, "syph.stage", NA, ids)
  return(dat)
}
