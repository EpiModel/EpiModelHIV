
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

  prepElig <- get_attr(dat, "prepElig")
  prepStat <- get_attr(dat, "prepStat")

  rGC <- get_attr(dat, "rGC")
  uGC <- get_attr(dat, "uGC")
  rCT <- get_attr(dat, "rCT")
  uCT <- get_attr(dat, "uCT")

  sexStigmaClass <- get_attr(dat, "sexStigmaClass")
  mdd.status <- get_attr(dat, "mdd.status")
  mdd.diag.status <- get_attr(dat, "mdd.diag.status")
  tx.status <- get_attr(dat, "tx.status")



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

  dat <- set_epi(dat, "prepCurr", at, sum(prepStat == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prepCurr.B", at,
    sum(prepStat == 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prepCurr.H", at,
    sum(prepStat == 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "prepCurr.W", at,
    sum(prepStat == 1 & race == 3, na.rm = TRUE))

  dat <- set_epi(dat, "prepCov", at, sum(prepStat == 1, na.rm = TRUE) /
                   sum(prepElig == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prepCov.mdd1", at,
                 sum(prepStat == 1 & mdd.status == 1, na.rm = TRUE) /
                   sum(prepElig == 1 & mdd.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prepCov.mdd0", at,
                 sum(prepStat == 1 & mdd.status == 0, na.rm = TRUE) /
                   sum(prepElig == 1 & mdd.status == 0, na.rm = TRUE))

  dat <- set_epi(dat, "prepCov.stigma1", at,
                 sum(prepStat == 1 & sexStigmaClass == 1, na.rm = TRUE) /
                   sum(prepElig == 1 & sexStigmaClass == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prepCov.stigma2", at,
                 sum(prepStat == 1 & sexStigmaClass == 2, na.rm = TRUE) /
                   sum(prepElig == 1 & sexStigmaClass == 2, na.rm = TRUE))
  dat <- set_epi(dat, "prepCov.stigma3", at,
                 sum(prepStat == 1 & sexStigmaClass == 3, na.rm = TRUE) /
                   sum(prepElig == 1 & sexStigmaClass == 3, na.rm = TRUE))
  dat <- set_epi(dat, "prepCov.stigma4", at,
                 sum(prepStat == 1 & sexStigmaClass == 4, na.rm = TRUE) /
                   sum(prepElig == 1 & sexStigmaClass == 4, na.rm = TRUE))


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


  #stigma
  dat <- set_epi(dat, "stigma.1.num", at, sum(sexStigmaClass == 1, na.rm = TRUE))
  dat <- set_epi(dat, "stigma.2.num", at, sum(sexStigmaClass == 2, na.rm = TRUE))
  dat <- set_epi(dat, "stigma.3.num", at, sum(sexStigmaClass == 3, na.rm = TRUE))
  dat <- set_epi(dat, "stigma.4.num", at, sum(sexStigmaClass == 4, na.rm = TRUE))

  s.num <- sum(status == 0, na.rm = TRUE)

  dat <- set_epi(dat, "stigma.1.prpall", at, sum(sexStigmaClass == 1, na.rm = TRUE)/ num)
  dat <- set_epi(dat, "stigma.2.prpall", at, sum(sexStigmaClass == 2, na.rm = TRUE)/ num)
  dat <- set_epi(dat, "stigma.3.prpall", at, sum(sexStigmaClass == 3, na.rm = TRUE)/ num)
  dat <- set_epi(dat, "stigma.4.prpall", at, sum(sexStigmaClass == 4, na.rm = TRUE)/ num)

  dat <- set_epi(dat, "stigma.1.prphiv1", at,
                 sum(sexStigmaClass == 1 & status == 1, na.rm = TRUE)/ i.num)
  dat <- set_epi(dat, "stigma.2.prphiv1", at,
                 sum(sexStigmaClass == 2& status == 1, na.rm = TRUE)/ i.num)
  dat <- set_epi(dat, "stigma.3.prphiv1", at,
                 sum(sexStigmaClass == 3 & status == 1, na.rm = TRUE)/ i.num)
  dat <- set_epi(dat, "stigma.4.prphiv1", at,
                 sum(sexStigmaClass == 4 & status == 1, na.rm = TRUE)/ i.num)

  dat <- set_epi(dat, "stigma.1.prphiv0", at,
                 sum(sexStigmaClass == 1 & status == 0, na.rm = TRUE)/ s.num)
  dat <- set_epi(dat, "stigma.2.prphiv0", at,
                 sum(sexStigmaClass == 2& status == 0, na.rm = TRUE)/ s.num)
  dat <- set_epi(dat, "stigma.3.prphiv0", at,
                 sum(sexStigmaClass == 3 & status == 0, na.rm = TRUE)/ s.num)
  dat <- set_epi(dat, "stigma.4.prphiv0", at,
                 sum(sexStigmaClass == 4 & status == 0, na.rm = TRUE)/ s.num)

  #mdd
  dat <- set_epi(dat, "mdd.num", at, sum(mdd.status == 1, na.rm = T))
  dat <- set_epi(dat, "mdd.prpall", at, sum(mdd.status == 1, na.rm = T) / num)

  dat <- set_epi(dat, "mdd.prphiv1", at,
                 sum(mdd.status == 1 & status == 1, na.rm = TRUE)/ i.num)
  dat <- set_epi(dat, "mdd.prphiv0", at,
                 sum(mdd.status == 1 & status == 0, na.rm = TRUE)/ s.num)

  dat <- set_epi(
    dat, "mdd.diag.prphiv1", at,
    sum(mdd.diag.status == 1 & status == 1, na.rm = TRUE) /
      sum(mdd.status == 1 & status == 1, na.rm = TRUE)
  )
  dat <- set_epi(
    dat, "mdd.diag.prphiv0", at,
    sum(mdd.diag.status == 1 & status == 0, na.rm = TRUE) /
      sum(mdd.status == 1 & status == 0, na.rm = TRUE)
  )


  dat <- set_epi(dat, "mdd.prpstigma1", at,
                 sum(mdd.status == 1 & sexStigmaClass == 1, na.rm = TRUE)/
                   sum(sexStigmaClass == 1, na.rm = TRUE))
  dat <- set_epi(dat, "mdd.prpstigma2", at,
                 sum(mdd.status == 1 & sexStigmaClass == 2, na.rm = TRUE)/
                   sum(sexStigmaClass == 2, na.rm = TRUE))
  dat <- set_epi(dat, "mdd.prpstigma3", at,
                 sum(mdd.status == 1 & sexStigmaClass == 3, na.rm = TRUE)/
                   sum(sexStigmaClass == 3, na.rm = TRUE))
  dat <- set_epi(dat, "mdd.prpstigma4", at,
                 sum(mdd.status == 1 & sexStigmaClass == 4, na.rm = TRUE)/
                   sum(sexStigmaClass == 4, na.rm = TRUE))


  dat <- set_epi(dat, "mdd.prprace1", at,
                 sum(mdd.status == 1 & race == 1, na.rm = TRUE)/
                   sum(race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "mdd.prprace2", at,
                 sum(mdd.status == 1 & race == 2, na.rm = TRUE)/
                   sum(race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "mdd.prprace3", at,
                 sum(mdd.status == 1 & race == 3, na.rm = TRUE)/
                   sum(race == 3, na.rm = TRUE))


  #ir - mdd & sexstigmaclass
  incid.mdd1 <- get_epi(dat, "incid.mdd1", at)
  incid.mdd0 <- get_epi(dat, "incid.mdd0", at)
  incid.stigma1 <- get_epi(dat, "incid.stigma1", at)
  incid.stigma2 <- get_epi(dat, "incid.stigma2", at)
  incid.stigma3 <- get_epi(dat, "incid.stigma3", at)
  incid.stigma4 <- get_epi(dat, "incid.stigma4", at)

  dat <- set_epi(dat, "mdd.ir100.1", at,
                 (incid.mdd1 /
                    sum(status == 0 & mdd.status == 1, na.rm = T)) * 100)
  dat <- set_epi(dat, "mdd.ir100.0", at,
                 (incid.mdd1 /
                    sum(status == 0 & mdd.status == 0, na.rm = T)) * 100)

  dat <- set_epi(dat, "stigma.ir100.1", at,
                 (incid.stigma1 /
                    sum(status == 0 & sexStigmaClass == 1, na.rm = T)) * 100)
  dat <- set_epi(dat, "stigma.ir100.2", at,
                 (incid.stigma2 /
                    sum(status == 0 & sexStigmaClass == 2, na.rm = T)) * 100)
  dat <- set_epi(dat, "stigma.ir100.3", at,
                 (incid.stigma3 /
                    sum(status == 0 & sexStigmaClass == 3, na.rm = T)) * 100)
  dat <- set_epi(dat, "stigma.ir100.4", at,
                 (incid.stigma4 /
                    sum(status == 0 & sexStigmaClass == 4, na.rm = T)) * 100)


  #diag and ART coverage
  dat <- set_epi(dat, "hivdiag", at, sum(status == 1 & diag.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "hivdiagCov", at,
                 sum(status == 1 & diag.status == 1, na.rm = TRUE) /
                   sum(status == 1, na.rm = TRUE))

  dat <- set_epi(dat, "hivdiag.mdd1", at,
                 sum(status == 1 & diag.status == 1 & mdd.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "hivdiagCov.mdd1", at,
                 sum(status == 1 & diag.status == 1 & mdd.status == 1, na.rm = TRUE) /
                   sum(status == 1 & mdd.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "hivdiag.mdd0", at,
                 sum(status == 1 & diag.status == 1 & mdd.status == 0, na.rm = TRUE))
  dat <- set_epi(dat, "hivdiagCov.mdd0", at,
                 sum(status == 1 & diag.status == 1 & mdd.status == 0, na.rm = TRUE) /
                   sum(status == 1 & mdd.status == 0, na.rm = TRUE))

  dat <- set_epi(dat, "hivdiag.stigma1", at,
                 sum(status == 1 & diag.status == 1 & sexStigmaClass == 1, na.rm = TRUE))
  dat <- set_epi(dat, "hivdiagCov.stigma1", at,
                 sum(status == 1 & diag.status == 1 & sexStigmaClass == 1, na.rm = TRUE) /
                   sum(status == 1 & sexStigmaClass == 1, na.rm = TRUE))
  dat <- set_epi(dat, "hivdiag.stigma2", at,
                 sum(status == 1 & diag.status == 1 & sexStigmaClass == 2, na.rm = TRUE))
  dat <- set_epi(dat, "hivdiagCov.stigma2", at,
                 sum(status == 1 & diag.status == 1 & sexStigmaClass == 2, na.rm = TRUE) /
                   sum(status == 1 & sexStigmaClass == 2, na.rm = TRUE))
  dat <- set_epi(dat, "hivdiag.stigma3", at,
                 sum(status == 1 & diag.status == 1 & sexStigmaClass == 3, na.rm = TRUE))
  dat <- set_epi(dat, "hivdiagCov.stigma3", at,
                 sum(status == 1 & diag.status == 1 & sexStigmaClass == 3, na.rm = TRUE) /
                   sum(status == 1 & sexStigmaClass == 3, na.rm = TRUE))
  dat <- set_epi(dat, "hivdiag.stigma4", at,
                 sum(status == 1 & diag.status == 1 & sexStigmaClass == 4, na.rm = TRUE))
  dat <- set_epi(dat, "hivdiagCov.stigma4", at,
                 sum(status == 1 & diag.status == 1 & sexStigmaClass == 4, na.rm = TRUE) /
                   sum(status == 1 & sexStigmaClass == 4, na.rm = TRUE))

  dat <- set_epi(dat, "artCurr", at, sum(tx.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "artCov", at, sum(tx.status == 1, na.rm = TRUE) /
                   sum(status == 1 & diag.status == 1, na.rm = TRUE))

  dat <- set_epi(dat, "artCurr.mdd1", at, sum(tx.status == 1 & mdd.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "artCov.mdd1", at, sum(tx.status == 1 & mdd.status == 1, na.rm = TRUE) /
                   sum(status == 1 & diag.status == 1 & mdd.status == 1, na.rm = TRUE))
  dat <- set_epi(dat, "artCurr.mdd0", at, sum(tx.status == 1 & mdd.status == 0, na.rm = TRUE))
  dat <- set_epi(dat, "artCov.mdd0", at, sum(tx.status == 1 & mdd.status == 0, na.rm = TRUE) /
                   sum(status == 1 & diag.status == 1 & mdd.status == 0, na.rm = TRUE))

  dat <- set_epi(dat, "artCurr.stigma1", at,
                 sum(tx.status == 1 & sexStigmaClass == 1, na.rm = TRUE))
  dat <- set_epi(dat, "artCov.stigma1", at,
                 sum(tx.status == 1 & sexStigmaClass == 1, na.rm = TRUE) /
                   sum(status == 1 & diag.status == 1 & sexStigmaClass == 1, na.rm = TRUE))
  dat <- set_epi(dat, "artCurr.stigma2", at,
                 sum(tx.status == 1 & sexStigmaClass == 2, na.rm = TRUE))
  dat <- set_epi(dat, "artCov.stigma2", at,
                 sum(tx.status == 1 & sexStigmaClass == 2, na.rm = TRUE) /
                   sum(status == 1 & diag.status == 1 & sexStigmaClass == 2, na.rm = TRUE))
  dat <- set_epi(dat, "artCurr.stigma3", at,
                 sum(tx.status == 1 & sexStigmaClass == 3, na.rm = TRUE))
  dat <- set_epi(dat, "artCov.stigma3", at,
                 sum(tx.status == 1 & sexStigmaClass == 3, na.rm = TRUE) /
                   sum(status == 1 & diag.status == 1 & sexStigmaClass == 3, na.rm = TRUE))
  dat <- set_epi(dat, "artCurr.stigma4", at,
                 sum(tx.status == 1 & sexStigmaClass == 4, na.rm = TRUE))
  dat <- set_epi(dat, "artCov.stigma4", at,
                 sum(tx.status == 1 & sexStigmaClass == 4, na.rm = TRUE) /
                   sum(status == 1 & diag.status == 1 & sexStigmaClass == 4, na.rm = TRUE))

  v.supp <- get_epi(dat, "vSupp", at)
  v.supp.prp <- v.supp / sum(status == 1 & diag.status == 1, na.rm = TRUE)
  dat <- set_epi(dat, "vSuppCov", at, v.supp.prp)


  return(dat)
}
