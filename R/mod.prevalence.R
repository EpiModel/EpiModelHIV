
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

  prepElig  <- get_attr(dat, "prepElig")
  prepStat  <- get_attr(dat, "prepStat")
  prepClass <- get_attr(dat, "prepClass")

  rGC <- get_attr(dat, "rGC")
  uGC <- get_attr(dat, "uGC")
  rCT <- get_attr(dat, "rCT")
  uCT <- get_attr(dat, "uCT")

  # EDP Attributes
  edp.prepStat <- get_attr(dat, "edp.prepStat")
  edp.prepElig <- get_attr(dat, "edp.prepElig")

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
    sum(race == 1 & status == 1, na.rm = TRUE) / sum(race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "i.prev.H", at,
    sum(race == 2 & status == 1, na.rm = TRUE) / sum(race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "i.prev.W", at,
    sum(race == 3 & status == 1, na.rm = TRUE) / sum(race == 3, na.rm = TRUE))

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
  dat <- set_epi(dat, "ir100", at, (
      incid / sum(status == 0, incid, na.rm = TRUE)) * 100 * (364 / time.unit))
  dat <- set_epi(dat, "ir100.B", at,
    (incid.B / sum(status == 0 & race == 1, incid.B, na.rm = TRUE)) *
    100 * (364 / time.unit))
  dat <- set_epi(dat, "ir100.H", at,
    (incid.H / sum(status == 0 & race == 2, incid.H, na.rm = TRUE)) *
    100 * (364 / time.unit))
  dat <- set_epi(dat, "ir100.W", at,
    (incid.W / sum(status == 0 & race == 3, incid.W, na.rm = TRUE)) *
      100 * (364 / time.unit))

  dat <- set_epi(dat, "prepElig", at, sum(prepElig == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prepElig.B", at,
    sum(prepElig == 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prepElig.H", at,
    sum(prepElig == 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "prepElig.W", at,
    sum(prepElig == 1 & race == 3, na.rm = TRUE))

  dat <- set_epi(dat, "edp.prepElig", at, sum(edp.prepElig == 1, na.rm = TRUE))
  dat <- set_epi(dat, "edp.prepElig.B", at,
                 sum(edp.prepElig == 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "edp.prepElig.H", at,
                 sum(edp.prepElig == 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "edp.prepElig.W", at,
                 sum(edp.prepElig == 1 & race == 3, na.rm = TRUE))

  dat <- set_epi(dat, "prepCurr", at, sum(prepStat == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prepCurr.B", at,
    sum(prepStat == 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prepCurr.H", at,
    sum(prepStat == 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "prepCurr.W", at,
    sum(prepStat == 1 & race == 3, na.rm = TRUE))

  dat <- set_epi(dat, "edp.prepCurr", at, sum(edp.prepStat == 1, na.rm = TRUE))
  dat <- set_epi(dat, "edp.prepCurr.B", at,
                 sum(edp.prepStat == 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "edp.prepCurr.H", at,
                 sum(edp.prepStat == 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "edp.prepCurr.W", at,
                 sum(edp.prepStat == 1 & race == 3, na.rm = TRUE))

  dat <- set_epi(dat, "prepCurr.low", at, sum(prepStat == 1 & prepClass == 1))
  dat <- set_epi(dat, "prepCurr.med", at, sum(prepStat == 1 & prepClass == 2))
  dat <- set_epi(dat, "prepCurr.hi", at, sum(prepStat == 1 & prepClass == 3))

  # STIs
  incid.rgc <- get_epi(dat, "incid.rgc", at)
  incid.ugc <- get_epi(dat, "incid.ugc", at)
  incid.rct <- get_epi(dat, "incid.rct", at)
  incid.uct <- get_epi(dat, "incid.uct", at)

  prev.gc <- sum((rGC == 1 | uGC == 1), na.rm = TRUE) / num
  prev.ct <- sum((rCT == 1 | uCT == 1), na.rm = TRUE) / num

  ir100.rgc <- (incid.rgc / sum(rGC == 0, incid.rgc, na.rm = TRUE)) *
    100 * (364 / time.unit)
  ir100.ugc <- (incid.ugc / sum(uGC == 0, incid.ugc, na.rm = TRUE)) *
    100 * (364 / time.unit)
  ir100.gc <- ir100.rgc + ir100.ugc

  ir100.rct <- (incid.rct / sum(rCT == 0, incid.rct, na.rm = TRUE)) *
    100 * (364 / time.unit)
  ir100.uct <- (incid.uct / sum(uCT == 0, incid.uct, na.rm = TRUE)) *
    100 * (364 / time.unit)
  ir100.ct <- ir100.rct + ir100.uct

  ir100.sti <- ir100.gc + ir100.ct

  dat <- set_epi(dat, "prev.gc", at, prev.gc)
  dat <- set_epi(dat, "prev.ct", at, prev.ct)
  dat <- set_epi(dat, "ir100.gc", at, ir100.gc)
  dat <- set_epi(dat, "ir100.ct", at, ir100.ct)
  dat <- set_epi(dat, "ir100.sti", at, ir100.sti)

  return(dat)
}
