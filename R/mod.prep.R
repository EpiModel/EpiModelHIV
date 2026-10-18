
#' @title PrEP Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging_msm
#'
#' @export
#'
prep_msm <- function(dat, at) {
  dat <- prep_indications(dat)
  dat <- prep_gfr_stop(dat)
  dat <- prep_stop(dat)
  dat <- prep_start(dat)
  if (get_control(dat, "debug")) {
    dat <- dbg_prep(dat)
  }
  return(dat)
}

prep_start <- function(dat) {
  # Attributes
  race <- get_attr(dat, "race")
  prep.indic <- get_attr(dat, "prep.indic")
  prep.any <- get_attr(dat, "prep.any")
  gfr <- get_attr(dat, "gfr")
  # Parameters
  prep.start.rate <- get_param(dat, "prep.start.rate")

  elig_ids <- which(prep.indic == 1 & prep.any == 0 & gfr >= 60)
  rates <- prep.start.rate[race[elig_ids]]
  start_ids <- elig_ids[runif(length(elig_ids)) < rates]

  # Output
  dat <- set_start_prep_attrs(dat, start_ids)
  return(dat)
}

prep_stop <- function(dat) {
  at <- get_current_timestep(dat)
  # Attributes
  race <- get_attr(dat, "race")
  prep.indic <- get_attr(dat, "prep.indic")
  prep <- get_attr(dat, "prep")
  prep.risk.last <- get_attr(dat, "prep.risk.last")
  # Parameters
  prep.risk.reassess.int <- get_param(dat, "prep.risk.reassess.int")
  prep.discont.int <- get_param(dat, "prep.discont.int")
  prep.discont.rate <- 1 - 2^(-1 / prep.discont.int)

  # Indication lapse - is only checked every `prep.risk.reassess.int`
  #                  - HIV test on previous step
  risk_assess <- prep == 1 & (
    at - prep.risk.last >= prep.risk.reassess.int
  )
  stop_elig <- risk_assess & prep.indic == 0
  stop_elig_ids <- which(stop_elig)

  # Random Discontinuation
  elig_ids <- which(prep == 1 & !stop_elig)
  rates <- prep.discont.rate[race[elig_ids]]
  stop_rand_ids <- elig_ids[runif(length(elig_ids)) < rates]

  # Combine - they are exclusives
  stop_ids <- c(stop_elig_ids, stop_rand_ids)

  # Output
  # Epis
  dat <- set_epi(dat, "prep.stop.rng", at, length(stop_rand_ids))
  dat <- set_epi(dat, "prep.stop.inelig", at, length(stop_elig_ids))

  prep.start.last <- get_attr(dat, "prep.start.last")
  dat <- set_epi(dat, "prep.dur.mean", at,
                 mean(at - prep.start.last[stop_ids]))
  dat <- set_epi(dat, "prep.dur.mean.rng", at,
                 mean(at - prep.start.last[stop_rand_ids]))
  dat <- set_epi(dat, "prep.dur.mean.inelig", at,
                 mean(at - prep.start.last[stop_elig_ids]))

  # Attrs
  dat <- set_attr(dat, "prep.risk.last", at, posit_ids = which(risk_assess))
  dat <- set_stop_prep_attrs(dat, stop_ids)

  return(dat)
}

# Manages the stoppage rules for PrEP users due to GFR decline
prep_gfr_stop <- function(dat) {
  at <- get_current_timestep(dat)
  # Attributes
  prep <- get_attr(dat, "prep")
  age <- get_attr(dat, "age")
  gfr <- get_attr(dat, "gfr")
  gfr.prep.ref <- get_attr(dat, "gfr.prep.ref")
  prep.gfr.risk.last <- get_attr(dat, "prep.gfr.risk.last")
  # Parameters
  time.unit <- get_param(dat, "time.unit")

  low.int <- (364 / time.unit)
  high.int <- (364 / time.unit) / 2

  assess_ids <- which(
    prep == 1 & (
      # high risk, test every 6 months
      ((age >= 50 | gfr.prep.ref < 90) & at - prep.gfr.risk.last >= high.int)  |
        # low risk, test every 12 months
        ((age < 50 & gfr.prep.ref >= 90) & at - prep.gfr.risk.last >= low.int)
    )
  )

  stop_ids <- assess_ids[gfr[assess_ids] < 60]

  gfr.prep.ref[assess_ids] <- gfr[assess_ids]
  dat <- set_attr(dat, "gfr.prep.ref", gfr.prep.ref)
  dat <- set_attr(dat, "prep.gfr.risk.last", at, posit_ids = assess_ids)
  dat <- set_stop_prep_attrs(dat, stop_ids)

  if (get_control(dat, "debug")) {
    at <- get_current_timestep(dat)
    dat <- set_epi(dat, "dbg_prep_gfr_stop", at, length(stop_ids))
  }

  return(dat)
}

# Indications:
# - any sex 6 month + >=1 HIV+ non UVL partner
# - any sex 6 month + any STI
# - any UAI 6 month + >=1 unknown status partner
# - no gfr.lt60 in past 6 month
#
# Should we deal with unkown?
# Should we deal with PrEP outside of these "hard" indications?
#     "Discuss PrEP, prescribe if resquested"
prep_indications <- function(dat) {
  # Attributes
  diag.status <- get_attr(dat, "diag.status")
  prep.indic <- get_attr(dat, "prep.indic")
  act.last <- get_attr(dat, "act.last")
  uai.last <- get_attr(dat, "uai.last")
  usupp.partner.last <- get_attr(dat, "usupp.partner.last")
  unknown.partner.last <- get_attr(dat, "unknown.partner.last")
  vl.last.usupp <- get_attr(dat, "vl.last.usupp")
  last.neg.test <- get_attr(dat, "last.neg.test")
  syph.dx.last <- get_attr(dat, "syph.dx.last")
  gono.dx.last <- get_attr(dat, "gono.dx.last")
  chla.dx.last <- get_attr(dat, "chla.dx.last")

  # Params
  at <- get_current_timestep(dat)
  time.unit <- get_param(dat, "time.unit")

  # Act List
  # columns: p1, p2, ptype, uai, pid
  al <- unique(dat$temp$al[, c("p1", "p2", "uai")])

  # Compute last acts and last unprotected act
  any_act_ids <- unique(c(al[, c("p1", "p2")]))
  any_uai_ids <- unique(c(al[al[, "uai"] == 1, c("p1", "p2")]))
  act.last[any_act_ids] <- at
  uai.last[any_uai_ids] <- at

  # 6 months in time steps
  d6m <- (364 / 2) / time.unit

  # any partner HIV+ not suppressed past 6 months?
  usupp_6m <- diag.status == 1 & at - vl.last.usupp <= d6m
  usupp_6m_part_ids <- unique(c(
    al[, "p1"][usupp_6m[al[, "p2"]]],
    al[, "p2"][usupp_6m[al[, "p1"]]]
  ))
  usupp.partner.last[usupp_6m_part_ids] <- at

  # any partner of unknown HIV status? (last neg test > 6 month)
  unknown <- diag.status == 0 & at - last.neg.test > d6m
  unknown_parts_ids <- unique(c(
    al[, "p1"][unknown[al[, "p2"]]],
    al[, "p2"][unknown[al[, "p1"]]]
  ))
  unknown.partner.last[unknown_parts_ids] <- at

  # any STI diag in past 6 month
  sti_diag_6m <-
    at - syph.dx.last <= d6m |
    at - gono.dx.last <= d6m |
    at - chla.dx.last <= d6m

  # Finally calc the actual PrEP eligibility
  # nolint start
  prep.indic <- diag.status == 0 & (
    at - uai.last <= d6m & at - unknown.partner.last <= d6m |
    at - act.last <= d6m & (sti_diag_6m | at - usupp.partner.last <= d6m))

  # nolint end

  # Output
  dat <- set_attr(dat, "act.last", act.last)
  dat <- set_attr(dat, "uai.last", uai.last)
  dat <- set_attr(dat, "usupp.partner.last", usupp.partner.last)
  dat <- set_attr(dat, "unknown.partner.last", unknown.partner.last)
  dat <- set_attr(dat, "prep.indic", as.numeric(prep.indic))

  return(dat)
}

set_start_prep_attrs <- function(dat, ids) {
  at <- get_current_timestep(dat)
  dat <- set_attr(dat, "prep.start.last", at, posit_ids = ids)
  dat <- set_attr(dat, "prep.risk.last", at, posit_ids = ids)
  dat <- set_attr(dat, "prep", 1, posit_ids = ids)
  dat <- set_attr(dat, "prep.any", 1, posit_ids = ids)
  counter <- get_attr(dat, "prep.start.count", posit_ids = ids)
  dat <- set_attr(dat, "prep.start.count", counter + 1, posit_ids = ids)

  pclass <- get_attr(dat, "prep.std.class", posit_ids = ids)
  dat <- set_attr(dat, "prep.class", pclass, posit_ids = ids)

  gfr <- get_attr(dat, "gfr", posit_ids = ids)
  dat <- set_attr(dat, "gfr.prep.ref", gfr, posit_ids = ids)

  return(dat)
}

set_stop_prep_attrs <- function(dat, ids) {
  at <- get_current_timestep(dat)
  dat <- set_attr(dat, "prep", 0, posit_ids = ids)
  dat <- set_attr(dat, "prep.any", 0, posit_ids = ids)
  dat <- set_attr(dat, "prep.risk.last", -Inf, posit_ids = ids)
  dat <- set_attr(dat, "prep.start.last", -Inf, posit_ids = ids)
  dat <- set_attr(dat, "prep.stop.last", at, posit_ids = ids)
  dat <- set_attr(dat, "prep.class", -1, posit_ids = ids)
  return(dat)
}

dbg_prep <- function(dat) {
  at <- get_current_timestep(dat)
  dat <- with(get_attr_list(dat), {
    dat <- set_epi(dat, "dbg_prep_cur", at, sum(prep == 1))
    dat <- set_epi(dat, "dbg_prep_elig", at, sum(prep.indic == 1))
    dat <- set_epi(dat, "dbg_prep_starts", at, sum(prep.start.last == at))
    dat <- set_epi(dat, "dbg_prep_stops", at, sum(prep.stop.last == at))

    any_prep_ids <- which(prep.start.count > 0)
    dat <- set_epi(dat, "dbg_prep_any", at, length(any_prep_ids))
    dat <- set_epi(dat, "dbg_prep_eps_mean", at,
                   mean(prep.start.count[any_prep_ids]))
    dat <- set_epi(dat, "dbg_prep_eps_sd", at,
                   sd(prep.start.count[any_prep_ids]))
    dat <- set_epi(dat, "dbg_prep_eps_q1", at,
                   quantile(prep.start.count[any_prep_ids], 0.25))
    dat <- set_epi(dat, "dbg_prep_eps_q2", at,
                   quantile(prep.start.count[any_prep_ids], 0.50))
    dat <- set_epi(dat, "dbg_prep_eps_q3", at,
                   quantile(prep.start.count[any_prep_ids], 0.75))
    dat <- set_epi(dat, "dbg_prep_eps_max", at,
                   max(prep.start.count[any_prep_ids]))
    dat
  })
  return(dat)
}