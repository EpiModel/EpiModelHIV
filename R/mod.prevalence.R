
#' @title Prevalence Calculations within Time Steps
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Summary statistic calculations are of two broad forms: prevalence and
#' incidence. This function establishes the summary statistic vectors for both
#' prevalence and incidence at time 1, and then calculates the prevalence
#' statistics for times 2 onward. Incidence statistics (e.g., number of new
#' infections or deaths) are calculated within the modules as they depend on
#' vectors that are not stored external to the module.
#'
#' @return
#' This function returns the `dat` object with an updated summary of current
#' attributes stored in `dat$epi`.
#'
#' @keywords module msm
#'
#' @export
#'
prevalence_msm <- function(dat, at) {

  # Attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  diag.status <- get_attr(dat, "diag.status")
  race <- get_attr(dat, "race")
  age <- get_attr(dat, "age")

  prep.indic <- get_attr(dat, "prep.indic")
  prep <- get_attr(dat, "prep")
  prep.any <- get_attr(dat, "prep.any")
  prep.otc <- get_attr(dat, "prep.otc")
  prep.otc.indic <- get_attr(dat, "prep.otc.indic")

  gfr <- get_attr(dat, "gfr")

  gono.rect <- get_attr(dat, "gono.rect")
  gono.uret <- get_attr(dat, "gono.uret")
  chla.rect <- get_attr(dat, "chla.rect")
  chla.uret <- get_attr(dat, "chla.uret")
  syph <- get_attr(dat, "syph.inf")

  # Parameter
  time.unit <- get_param(dat, "time.unit")

  # Pop Size / Demog
  num <- sum(active == 1, na.rm = TRUE)
  dat <- set_epi(dat, "num", at, num)
  dat <- set_epi(dat, "num.B", at, sum(race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "num.H", at, sum(race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "num.W", at, sum(race == 3, na.rm = TRUE))
  dat <- set_epi(dat, "s.num", at, sum(status == 0, na.rm = TRUE))
  dat <- set_epi(dat, "i.num", at, sum(status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "i.num.B", at, sum(status == 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "i.num.H", at, sum(status == 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "i.num.W", at, sum(status == 1 & race == 3, na.rm = TRUE))
  dat <- set_epi(dat, "age.mean", at, mean(age, na.rm = TRUE))

  i.num <- sum(status == 1, na.rm = TRUE)
  # Prev and Incid
  dat <- set_epi(dat, "i.prev", at,  i.num / num)
  dat <- set_epi(dat, "i.prev.B", at,
                 sum(race == 1 & status == 1, na.rm = TRUE) /
                   sum(race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "i.prev.H", at,
                 sum(race == 2 & status == 1, na.rm = TRUE) /
                   sum(race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "i.prev.W", at,
                 sum(race == 3 & status == 1, na.rm = TRUE) /
                   sum(race == 3, na.rm = TRUE))

  dat <- set_epi(dat, "i.prev.dx", at,
                 sum(diag.status == 1, na.rm = TRUE) / num)
  dat <- set_epi(dat, "i.prev.dx.B", at,
                 sum(race == 1 & diag.status == 1, na.rm = TRUE) /
                   sum(race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "i.prev.dx.H", at,
                 sum(race == 2 & diag.status == 1, na.rm = TRUE) /
                   sum(race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "i.prev.dx.W", at,
                 sum(race == 3 & diag.status == 1, na.rm = TRUE) /
                   sum(race == 3, na.rm = TRUE))

  incid   <- get_epi(dat, "incid", at)
  incid.B <- get_epi(dat, "incid.B", at)
  incid.H <- get_epi(dat, "incid.H", at)
  incid.W <- get_epi(dat, "incid.W", at)
  dat <- set_epi(dat, "ir100", at,
                 (incid / sum(status == 0, incid, na.rm = TRUE)) *
                   100 * (364 / time.unit))
  dat <- set_epi(dat, "ir100.B", at,
                 (incid.B / sum(status == 0 & race == 1, incid.B, na.rm = TRUE)) *
                   100 * (364 / time.unit))
  dat <- set_epi(dat, "ir100.H", at,
                 (incid.H / sum(status == 0 & race == 2, incid.H, na.rm = TRUE)) *
                   100 * (364 / time.unit))
  dat <- set_epi(dat, "ir100.W", at,
                 (incid.W / sum(status == 0 & race == 3, incid.W, na.rm = TRUE)) *
                   100 * (364 / time.unit))

  # PrEP
  dat <- set_epi(dat, "prep.indic", at, sum(prep.indic == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.indic.B", at,
                 sum(prep.indic == 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.indic.H", at,
                 sum(prep.indic == 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "prep.indic.W", at,
                 sum(prep.indic == 1 & race == 3, na.rm = TRUE))

  dat <- set_epi(dat, "prepCurr", at, sum(prep == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prepCurr.B", at,
                 sum(prep == 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prepCurr.H", at,
                 sum(prep == 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "prepCurr.W", at,
                 sum(prep == 1 & race == 3, na.rm = TRUE))

  # PrEP OTC
  dat <- set_epi(dat, "prep.otc.indic", at, sum(prep.otc.indic == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.otc.indic.B", at,
                 sum(prep.otc.indic == 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.otc.indic.H", at,
                 sum(prep.otc.indic == 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "prep.otc.indic.W", at,
                 sum(prep.otc.indic == 1 & race == 3, na.rm = TRUE))

  dat <- set_epi(dat, "prep.otcCurr", at, sum(prep.otc == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.otcCurr.B", at,
                 sum(prep.otc == 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.otcCurr.H", at,
                 sum(prep.otc == 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "prep.otcCurr.W", at,
                 sum(prep.otc == 1 & race == 3, na.rm = TRUE))

  dat <- set_epi(dat, "prep.otc.std.indic", at, sum(
      prep.indic == 1 & prep.otc == 1, na.rm = TRUE))


  # On PrEP with GFR<60
  dat <- set_epi(dat, "prep_any_gfr_lt60", at,
                 sum(prep.any == 1 & gfr < 60) / sum(prep.any == 1))
  dat <- set_epi(dat, "prep_std_gfr_lt60", at,
                 sum(prep == 1 & gfr < 60) / sum(prep == 1))
  dat <- set_epi(dat, "prep_otc_gfr_lt60", at,
                 sum(prep.otc == 1 & gfr < 60) / sum(prep.otc == 1))

  # STIs
  incid.gono.rect <- get_epi(dat, "incid.gono.rect", at)
  incid.gono.uret <- get_epi(dat, "incid.gono.uret", at)
  incid.chla.rect <- get_epi(dat, "incid.chla.rect", at)
  incid.chla.uret <- get_epi(dat, "incid.chla.uret", at)
  incid.syph <- get_epi(dat, "incid.syph", at)

  prev.gono <- sum((gono.rect == 1 | gono.uret == 1), na.rm = TRUE) / num
  prev.chla <- sum((chla.rect == 1 | chla.uret == 1), na.rm = TRUE) / num
  prev.syph <- sum((syph == 1), na.rm = TRUE) / num

  ir100.gono.rect <- (incid.gono.rect / sum(gono.rect == 0, incid.gono.rect, na.rm = TRUE)) *
    100 * (364 / time.unit)
  ir100.gono.uret <- (incid.gono.uret / sum(gono.uret == 0, incid.gono.uret, na.rm = TRUE)) *
    100 * (364 / time.unit)
  ir100.gono <- ir100.gono.rect + ir100.gono.uret

  ir100.chla.rect <- (incid.chla.rect / sum(chla.rect == 0, incid.chla.rect, na.rm = TRUE)) *
    100 * (364 / time.unit)
  ir100.chla.uret <- (incid.chla.uret / sum(chla.uret == 0, incid.chla.uret, na.rm = TRUE)) *
    100 * (364 / time.unit)
  ir100.chla <- ir100.chla.rect + ir100.chla.uret

  ir100.syph <- (incid.syph / sum(syph == 0, incid.syph, na.rm = TRUE)) *
    100 * (364 / time.unit)

  ir100.sti <- ir100.gono + ir100.chla + ir100.syph

  dat <- set_epi(dat, "prev.gono", at, prev.gono)
  dat <- set_epi(dat, "prev.chla", at, prev.chla)
  dat <- set_epi(dat, "prev.syph", at, prev.syph)
  dat <- set_epi(dat, "ir100.gono", at, ir100.gono)
  dat <- set_epi(dat, "ir100.chla", at, ir100.chla)
  dat <- set_epi(dat, "ir100.syph", at, ir100.syph)
  dat <- set_epi(dat, "ir100.sti", at, ir100.sti)

  # Person time on PrEP with HIV+ (in person-year)
  ptime <- sum(status & prep.any) / 52
  dat <- set_epi(dat, "hiv.prep.py", at, ptime )

  return(dat)
}
