#' @title Renal Function Module
#'
#' @description Module function managing the renal function of MSMs in the model.
#'
#' @inheritParams aging_msm
#'
#' @export
#'
renal_msm <- function(dat, at) {
  dat <- gfr_decline(dat)

  # Attributes
  gfr <- get_attr(dat, "gfr")
  gfr.prep.baseline <- get_attr(dat, "gfr.prep.baseline")
  prep.any <- get_attr(dat, "prep.any")
  age <- get_attr(dat, "age")

  # Parameters
  gfr.decline.prep.rate <- get_param(dat, "gfr.decline.prep.rate")
  gfr.decline.age.gt50.or <- get_param(dat, "gfr.decline.age.gt50.or")
  gfr.decline.gfr.lt90.or <- get_param(dat, "gfr.decline.gfr.lt90.or")
  gfr.prep.recov.rate <- get_param(dat, "gfr.prep.recov.rate")

  # prep decline
  #   prob base on age > 50 & GFR < 90 (2 HR) + base prob
  #   use ~2 per 100 pyar for a start as baseline
  drop_elig_ids <- which(gfr >= 60 & prep.any == 1)
  log_rates <- rep(prob_to_log_odds(gfr.decline.prep.rate),
                   length(drop_elig_ids))
  yo50_eids <- which(age[drop_elig_ids] > 50)
  log_rates[yo50_eids] <- log_rates[yo50_eids] + log(gfr.decline.age.gt50.or)
  gfr90_eids <- which(gfr[drop_elig_ids] < 90)
  log_rates[gfr90_eids] <- log_rates[gfr90_eids] + log(gfr.decline.gfr.lt90.or)
  rates <- log_odds_to_prob(log_rates)
  drop_ids <- drop_elig_ids[runif(length(rates)) < rates]

  gfr.prep.baseline[drop_ids] <- gfr[drop_ids]
  gfr[drop_ids] <- 50

  # recovery
  #   per week prob so 75% recov in 8 weeks
  recov_elig_ids <- which(prep.any == 0 & !is.na(gfr.prep.baseline))
  recov_ids <- recov_elig_ids[runif(length(recov_elig_ids)) <
                                gfr.prep.recov.rate]

  # reset recovs to their baseline
  gfr[recov_ids] <- gfr.prep.baseline[recov_ids]
  gfr.prep.baseline[recov_ids] <- NA

  # Outputs --------------------------------------------------------------------
  dat <- set_attr(dat, "gfr", gfr)
  dat <- set_attr(dat, "gfr.prep.baseline", gfr.prep.baseline)

  if (get_control(dat, "debug")) {

    dat <- set_epi(dat, "dbg_gfr_drop_elig", at, length(drop_elig_ids))
    dat <- set_epi(dat, "dbg_gfr_drop", at, length(drop_ids))
    # above 50 yo
    yo50_pwar <- length(yo50_eids) # / 52 # for pyar
    yo50_drop <- length(intersect(drop_elig_ids[yo50_eids], drop_ids))
    dat <- set_epi(dat, "dbg_gfr_drop_yo50", at, yo50_drop)
    dat <- set_epi(dat, "dbg_yo50_pwar", at, yo50_pwar)
    # below 50 yo
    yo30_pwar <- length(drop_elig_ids) - yo50_pwar
    yo30_drop <- length(drop_ids) - yo50_drop
    dat <- set_epi(dat, "dbg_gfr_drop_yo30", at, yo30_drop)
    dat <- set_epi(dat, "dbg_yo30_pwar", at, yo30_pwar)

    gfr90_pwar <- length(gfr90_eids) # / 52 # for pyar
    gfr90_drop <- length(intersect(drop_elig_ids[gfr90_eids], drop_ids))
    dat <- set_epi(dat, "dbg_gfr_drop_gfr90", at, gfr90_drop)
    dat <- set_epi(dat, "dbg_gfr90_pwar", at, gfr90_pwar)
    # below 50 gfr
    gfr100_pwar <- length(drop_elig_ids) - gfr90_pwar
    gfr100_drop <- length(drop_ids) - gfr90_drop
    dat <- set_epi(dat, "dbg_gfr_drop_gfr100", at, gfr100_drop)
    dat <- set_epi(dat, "dbg_gfr100_pwar", at, gfr100_pwar)

    # Recov
    dat <- set_epi(dat, "dbg_gfr_recov_elig", at, length(recov_elig_ids))
    dat <- set_epi(dat, "dbg_gfr_recov", at, length(recov_ids))
  }

  return(dat)
}

  # NOTE: Decline rates based upon age
  #   3 possible states: 100 (>90), 75 (>60), 50 (<60)
  # - < 50yo
  # - > 50yo
gfr_decline <- function(dat) {
  # Attributes
  gfr <- get_attr(dat, "gfr")
  age <- get_attr(dat, "age")
  gfr.prep.baseline <- get_attr(dat, "gfr.prep.baseline")

  # Parameters:
  gfr.90.decline.rate <- get_param(dat, "gfr.90.decline.rate")
  gfr.60.decline.rate <- get_param(dat, "gfr.60.decline.rate")

  age_grps <- cut(age, c(15, 50, 65), labels = FALSE, right = FALSE)

  # Individual > 60 becoming < 60
  elig_ids <- which(gfr > 60 & gfr < 90)
  rates <- gfr.60.decline.rate[age_grps[elig_ids]]
  decline_60_ids <- elig_ids[runif(length(elig_ids)) < rates]
  gfr[decline_60_ids] <- 50

  # Individual > 90 becoming < 90 (75)
  elig_ids <- which(gfr > 90)
  rates <- gfr.90.decline.rate[age_grps[elig_ids]]
  decline_90_ids <- elig_ids[runif(length(elig_ids)) < rates]
  gfr[decline_90_ids] <- 75

  # Output
  dat <- set_attr(dat, "gfr", gfr)

  if (get_control(dat, "debug")) {
    # GFR per age grp
    at <- get_current_timestep(dat)
    age <- get_attr(dat, "age")
    pop_ids <- which(is.na(gfr.prep.baseline))
    age_breaks <- c(15, 30, 50, 65)
    age_grps <- cut(age[pop_ids], age_breaks, right = FALSE)
    gfr_grp_90 <- tapply(gfr[pop_ids], age_grps[pop_ids], \(x) mean(x < 90))
    gfr_grp_60 <- tapply(gfr[pop_ids], age_grps[pop_ids], \(x) mean(x < 60))

    for (i in seq_along(gfr_grp_60)) {
      dat <- set_epi(
        dat,
        paste0("dbg_gfr60_prop_", names(gfr_grp_60)[i]),
        at, gfr_grp_60[i]
      )
      dat <- set_epi(
        dat,
        paste0("dbg_gfr90_prop_", names(gfr_grp_90)[i]),
        at, gfr_grp_90[i]
      )
    }
  }

  return(dat)
}