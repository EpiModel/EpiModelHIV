#' @title PrEP OTC Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging_msm
#'
#' @export
#'
prep_otc_msm <- function(dat, at) {
  if (get_param(dat, "prep.otc.adhr.dist.realloc")) {
    dat <- prep_otc_realloc_adhr(dat)
  }
  dat <- prep_otc_indications(dat)
  dat <- prep_otc_gfr_stop(dat)
  dat <- prep_otc_stop(dat)
  dat <- prep_otc_start(dat)
  dat <- prep_otc_switch(dat)
  if (get_control(dat, "debug")) {
    dat <- dbg_prep_otc(dat)
  }
  return(dat)
}

prep_otc_realloc_adhr <- function(dat) {
  prep.otc.class <- get_attr(dat, "prep.otc.class")
  dist <- get_param(dat, "prep.otc.adhr.dist")
  prep.class <- sample(3, length(prep.otc.class), TRUE, dist)
  dat <- set_attr(dat, "prep.otc.class", prep.class)
  dat <- set_param(dat, "prep.otc.adhr.dist.realloc", FALSE)
  return(dat)
}

prep_otc_start <- function(dat) {
  # Attributes
  race <- get_attr(dat, "race")
  prep.otc.indic <- get_attr(dat, "prep.otc.indic")
  prep.any <- get_attr(dat, "prep.any")
  # Parameters
  prep.otc.start.rate <- get_param(dat, "prep.otc.start.rate")

  elig_ids <- which(prep.otc.indic == 1 & prep.any == 0)
  if (get_param(dat, "prep.otc.hard.indications") == 1) {
    gfr <- get_attr(dat, "gfr")
    elig_ids <- elig_ids[gfr[elig_ids] >= 60]
  }
  rates <- prep.otc.start.rate[race[elig_ids]]
  start_ids <- elig_ids[runif(length(elig_ids)) < rates]

  # Output
  dat <- set_start_prep_otc_attrs(dat, start_ids)
  return(dat)
}

prep_otc_stop <- function(dat) {
  at <- get_current_timestep(dat)
  # Attributes
  race <- get_attr(dat, "race")
  prep.otc.indic <- get_attr(dat, "prep.otc.indic")
  prep.otc <- get_attr(dat, "prep.otc")
  prep.otc.risk.last <- get_attr(dat, "prep.otc.risk.last")
  # Parameters
  prep.otc.risk.reassess.int <- get_param(dat, "prep.otc.risk.reassess.int")
  prep.otc.discont.int <- get_param(dat, "prep.otc.discont.int")
  prep.otc.discont.rate <- 1 - 2^(-1 / prep.otc.discont.int)

  # Indication lapse - is only checked every `prep.risk.reassess.int`
  #                  - HIV test on previous step
  risk_assess <- prep.otc == 1 &
    (at - prep.otc.risk.last >= prep.otc.risk.reassess.int)
  stop_elig <- risk_assess & prep.otc.indic == 0
  stop_elig_ids <- which(stop_elig)

  # Random Discontinuation
  elig_ids <- which(prep.otc == 1 & !stop_elig)
  rates <- prep.otc.discont.rate[race[elig_ids]]
  stop_rand_ids <- elig_ids[runif(length(elig_ids)) < rates]

  # Combine - they are exclusives
  stop_ids <- c(stop_elig_ids, stop_rand_ids)

  # Output
  # Epis
  dat <- set_epi(dat, "prep.otc.stop.rng", at, length(stop_rand_ids))
  dat <- set_epi(dat, "prep.otc.stop.inelig", at, length(stop_elig_ids))

  prep.otc.start.last <- get_attr(dat, "prep.otc.start.last")
  dat <- set_epi(
    dat,
    "prep.otc.dur.mean",
    at,
    mean(at - prep.otc.start.last[stop_ids])
  )
  dat <- set_epi(
    dat,
    "prep.otc.dur.mean.rng",
    at,
    mean(at - prep.otc.start.last[stop_rand_ids])
  )
  dat <- set_epi(
    dat,
    "prep.otc.dur.mean.inelig",
    at,
    mean(at - prep.otc.start.last[stop_elig_ids])
  )

  # Attrs
  dat <- set_attr(dat, "prep.otc.risk.last", at, posit_ids = which(risk_assess))
  dat <- set_stop_prep_otc_attrs(dat, stop_ids)

  return(dat)
}

# Indications:
#   - not on prep
#   - not diagnosed with HIV
#   - sexually active (any sex 6M)
prep_otc_indications <- function(dat) {
  if (get_param(dat, "prep.otc.hard.indications") == 1) {
    # Use same indications as normal PrEP
    indic <- get_attr(dat, "prep.indic")
    dat <- set_attr(dat, "prep.otc.indic", indic)
    return(dat)
  }

  # Else, relaxed indications
  #
  # Attributes
  diag.status <- get_attr(dat, "diag.status")
  act.last <- get_attr(dat, "act.last")

  # Parameters
  at <- get_current_timestep(dat)
  time.unit <- get_param(dat, "time.unit")
  d6m <- (364 / 2) / time.unit

  prep.otc.indic <- at - act.last <= d6m & diag.status == 0

  dat <- set_attr(dat, "prep.otc.indic", as.numeric(prep.otc.indic))

  return(dat)
}

prep_otc_switch <- function(dat) {
  # Attributes
  prep <- get_attr(dat, "prep")
  prep.indic <- get_attr(dat, "prep.indic")
  prep.otc <- get_attr(dat, "prep.otc")
  prep.otc.indic <- get_attr(dat, "prep.otc.indic")
  # Params
  prep.std.switch.otc.prob <- get_param(dat, "prep.std.switch.otc.prob")
  prep.otc.switch.std.prob <- get_param(dat, "prep.otc.switch.std.prob")

  # Switch from OTC to STD
  switch_prob <- prep.otc.switch.std.prob
  otc_elig_ids <- which(prep.otc == 1 & prep.indic == 1)
  otc_switch_ids <- otc_elig_ids[runif(length(otc_elig_ids)) < switch_prob]

  # Switch from STD to OTC
  switch_prob <- prep.std.switch.otc.prob
  std_elig_ids <- which(prep == 1 & prep.otc.indic == 1)
  std_switch_ids <- std_elig_ids[runif(length(std_elig_ids)) < switch_prob]

  # Switch OTC -> STD
  dat <- set_stop_prep_otc_attrs(dat, otc_switch_ids)
  dat <- set_start_prep_attrs(dat, otc_switch_ids)

  # Switch STD -> OTC
  dat <- set_stop_prep_attrs(dat, std_switch_ids)
  dat <- set_start_prep_otc_attrs(dat, std_switch_ids)

  if (get_control(dat, "debug")) {
    at <- get_current_timestep(dat)
    dat <- set_epi(dat, "dbg_switch_otc_std_elig", at, length(otc_elig_ids))
    dat <- set_epi(dat, "dbg_switch_otc_std", at, length(otc_switch_ids))
    dat <- set_epi(dat, "dbg_switch_std_otc_elig", at, length(std_elig_ids))
    dat <- set_epi(dat, "dbg_switch_std_otc", at, length(std_switch_ids))
  }

  return(dat)
}

set_start_prep_otc_attrs <- function(dat, ids) {
  at <- get_current_timestep(dat)
  dat <- set_attr(dat, "prep.otc.start.last", at, posit_ids = ids)
  dat <- set_attr(dat, "prep.otc.risk.last", at, posit_ids = ids)
  dat <- set_attr(dat, "prep.otc", 1, posit_ids = ids)
  dat <- set_attr(dat, "prep.any", 1, posit_ids = ids)
  counter <- get_attr(dat, "prep.otc.start.count", posit_ids = ids)
  dat <- set_attr(dat, "prep.otc.start.count", counter + 1, posit_ids = ids)

  pclass <- get_attr(dat, "prep.otc.class", posit_ids = ids)
  dat <- set_attr(dat, "prep.class", pclass, posit_ids = ids)

  return(dat)
}

set_stop_prep_otc_attrs <- function(dat, ids) {
  at <- get_current_timestep(dat)
  dat <- set_attr(dat, "prep.otc", 0, posit_ids = ids)
  dat <- set_attr(dat, "prep.any", 0, posit_ids = ids)
  dat <- set_attr(dat, "prep.otc.risk.last", -Inf, posit_ids = ids)
  dat <- set_attr(dat, "prep.otc.start.last", -Inf, posit_ids = ids)
  dat <- set_attr(dat, "prep.otc.stop.last", at, posit_ids = ids)
  dat <- set_attr(dat, "prep.class", -1, posit_ids = ids)
  return(dat)
}

dbg_prep_otc <- function(dat) {
  at <- get_current_timestep(dat)
  dat <- with(get_attr_list(dat), {
    dat <- set_epi(dat, "dbg_prep_otc_cur", at, sum(prep.otc == 1))
    dat <- set_epi(dat, "dbg_prep_otc_elig", at, sum(prep.otc.indic == 1))
    dat <- set_epi(
      dat,
      "dbg_prep_otc_starts",
      at,
      sum(prep.otc.start.last == at)
    )
    dat <- set_epi(dat, "dbg_prep_otc_stops", at, sum(prep.otc.stop.last == at))

    any_otc_ids <- which(prep.otc.start.count > 0)
    dat <- set_epi(dat, "dbg_prep_otc_any", at, length(any_otc_ids))
    dat <- set_epi(
      dat,
      "dbg_prep_otc_eps_mean",
      at,
      mean(prep.otc.start.count[any_otc_ids])
    )
    dat <- set_epi(
      dat,
      "dbg_prep_otc_eps_sd",
      at,
      sd(prep.otc.start.count[any_otc_ids])
    )
    dat <- set_epi(
      dat,
      "dbg_prep_otc_eps_q1",
      at,
      quantile(prep.otc.start.count[any_otc_ids], 0.25)
    )
    dat <- set_epi(
      dat,
      "dbg_prep_otc_eps_q2",
      at,
      quantile(prep.otc.start.count[any_otc_ids], 0.50)
    )
    dat <- set_epi(
      dat,
      "dbg_prep_otc_eps_q3",
      at,
      quantile(prep.otc.start.count[any_otc_ids], 0.75)
    )
    dat <- set_epi(
      dat,
      "dbg_prep_otc_eps_max",
      at,
      max(prep.otc.start.count[any_otc_ids])
    )
    dat
  })
  return(dat)
}

# Manages the stoppage rules for PrEP OTC users due to GFR decline
# NOTE: share the gfr.prep.ref
prep_otc_gfr_stop <- function(dat) {
  at <- get_current_timestep(dat)
  # Attributes
  prep <- get_attr(dat, "prep.otc")
  age <- get_attr(dat, "age")
  gfr <- get_attr(dat, "gfr")
  gfr.prep.ref <- get_attr(dat, "gfr.prep.ref")
  prep.gfr.risk.last <- get_attr(dat, "prep.gfr.risk.last")
  # Parameters
  low.int <- get_param(dat, "prep.otc.gfr.low.risk.int")
  high.int <- get_param(dat, "prep.otc.gfr.high.risk.int")

  hig_risk_ids <- which(prep == 1 & (age >= 50 | gfr.prep.ref < 90))
  low_risk_ids <- which(prep == 1 & (age < 50 & gfr.prep.ref >= 90))

  if (get_param(dat, "prep.otc.gfr.risk.rng") == 1) {
    assess_ids <- c(
      hig_risk_ids[runif(length(hig_risk_ids)) < 1 / high.int],
      low_risk_ids[runif(length(low_risk_ids)) < 1 / low.int]
    )
  } else {
    assess_ids <- c(
      hig_risk_ids[at - prep.gfr.risk.last[hig_risk_ids] >= high.int],
      low_risk_ids[at - prep.gfr.risk.last[low_risk_ids] >= low.int]
    )
  }

  stop_ids <- assess_ids[gfr[assess_ids] < 60]

  gfr.prep.ref[assess_ids] <- gfr[assess_ids]
  dat <- set_attr(dat, "gfr.prep.ref", gfr.prep.ref)
  dat <- set_attr(dat, "prep.gfr.risk.last", at, posit_ids = assess_ids)
  dat <- set_stop_prep_otc_attrs(dat, stop_ids)

  if (get_control(dat, "debug")) {
    at <- get_current_timestep(dat)
    dat <- set_epi(dat, "dbg_prep_otc_gfr_stop", at, length(stop_ids))
  }

  return(dat)
}
