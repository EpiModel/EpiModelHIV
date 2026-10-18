
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
#' This function returns the \code{dat} object with an updated summary of current
#' attributes stored in \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
prevalence_msm <- function(dat, at) {

  active <- get_attr(dat, "active")
  dat <- set_epi(dat, "num", at, sum(active == 1, na.rm = TRUE))

  # active <- get_attr(dat, "active")
  # status <- get_attr(dat, "status")
  # diag.status <- get_attr(dat, "diag.status")
  # diag.time <- get_attr(dat, "diag.time")
  # inf.time <- get_attr(dat, "inf.time")
  # race <- get_attr(dat, "race")
  # age <- get_attr(dat, "age")
  # tx.init.time <- get_attr(dat, "tx.init.time")
  # vl <- get_attr(dat, "vl")
  # last.neg.test <- get_attr(dat, "last.neg.test")

  # prepElig <- get_attr(dat, "prepElig")
  # prepStat <- get_attr(dat, "prepStat")
  # prepClass <- get_attr(dat, "prepClass")

  rGC <- get_attr(dat, "rGC")
  uGC <- get_attr(dat, "uGC")
  rCT <- get_attr(dat, "rCT")
  uCT <- get_attr(dat, "uCT")

  # # Pop Size / Demog
  # dat <- set_epi(dat, "num", at, sum(active == 1, na.rm = TRUE))
  # dat <- set_epi(dat, "num.B", at, sum(race == 1, na.rm = TRUE))
  # dat <- set_epi(dat, "num.H", at, sum(race == 2, na.rm = TRUE))
  # dat <- set_epi(dat, "num.W", at, sum(race == 3, na.rm = TRUE))
  # dat <- set_epi(dat, "age.mean", at, mean(age, na.rm = TRUE))
  # dat <- set_epi(dat, "s.num", at, sum(status == 0, na.rm = TRUE))
  # dat <- set_epi(dat, "i.num", at, sum(status == 1, na.rm = TRUE))
  # dat <- set_epi(dat, "i.num.B", at, sum(status == 1 & race == 1, na.rm = TRUE))
  # dat <- set_epi(dat, "i.num.H", at, sum(status == 1 & race == 2, na.rm = TRUE))
  # dat <- set_epi(dat, "i.num.W", at, sum(status == 1 & race == 3, na.rm = TRUE))

  # dat <- set_epi(dat, "i.num.dx", at, sum(diag.status == 1, na.rm = TRUE))

  # # Prev / Incid
  # dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  # dat$epi$i.prev.B[at] <- sum(race == 1 & status == 1, na.rm = TRUE) /
  #   sum(race == 1, na.rm = TRUE)
  # dat$epi$i.prev.H[at] <- sum(race == 2 & status == 1, na.rm = TRUE) /
  #   sum(race == 2, na.rm = TRUE)
  # dat$epi$i.prev.W[at] <- sum(race == 3 & status == 1, na.rm = TRUE) /
  #   sum(race == 3, na.rm = TRUE)

  # dat$epi$i.prev.dx[at] <- sum(diag.status == 1, na.rm = TRUE) / dat$epi$num[at]
  # dat$epi$i.prev.dx.B[at] <- sum(race == 1 & diag.status == 1, na.rm = TRUE) /
  #   sum(race == 1, na.rm = TRUE)
  # dat$epi$i.prev.dx.H[at] <- sum(race == 2 & diag.status == 1, na.rm = TRUE) /
  #   sum(race == 2, na.rm = TRUE)
  # dat$epi$i.prev.dx.W[at] <- sum(race == 3 & diag.status == 1, na.rm = TRUE) /
  #   sum(race == 3, na.rm = TRUE)

  # dat$epi$ir100[at] <- (dat$epi$incid[at] / sum(status == 0, dat$epi$incid[at], na.rm = TRUE)) * 5200
  # dat$epi$ir100.B[at] <- (dat$epi$incid.B[at] / sum(status == 0 & race == 1,
  #                                                   dat$epi$incid.B[at], na.rm = TRUE)) * 5200
  # dat$epi$ir100.H[at] <- (dat$epi$incid.H[at] / sum(status == 0 & race == 2,
  #                                                   dat$epi$incid.H[at], na.rm = TRUE)) * 5200
  # dat$epi$ir100.W[at] <- (dat$epi$incid.W[at] / sum(status == 0 & race == 3,
  #                                                   dat$epi$incid.W[at], na.rm = TRUE)) * 5200


  # # Care continuum stats (primary)
  # dat$epi$cc.dx[at] <- sum(diag.status == 1 & inf.time >= 2, na.rm = TRUE) /
  #                      sum(status == 1 & inf.time >= 2, na.rm = TRUE)
  # dat$epi$cc.dx.B[at] <- sum(diag.status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE) /
  #                        sum(status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE)
  # dat$epi$cc.dx.H[at] <- sum(diag.status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE) /
  #                        sum(status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)
  # dat$epi$cc.dx.W[at] <- sum(diag.status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE) /
  #                        sum(status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE)

  # dat$epi$cc.linked1m[at] <- sum(dat$attr$tx.init.time - dat$attr$diag.time <= 4 & dat$attr$diag.time >= 2, na.rm = TRUE) /
  #                            sum(dat$attr$diag.status == 1 & dat$attr$diag.time >= 2, na.rm = TRUE)
  # dat$epi$cc.linked1m.B[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 1, na.rm = TRUE) /
  #                              sum(diag.status == 1 & diag.time >= 2 & race == 1, na.rm = TRUE)
  # dat$epi$cc.linked1m.H[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 2, na.rm = TRUE) /
  #                              sum(diag.status == 1 & diag.time >= 2 & race == 2, na.rm = TRUE)
  # dat$epi$cc.linked1m.W[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 3, na.rm = TRUE) /
  #                              sum(diag.status == 1 & diag.time >= 2 & race == 3, na.rm = TRUE)

  # dat$epi$cc.linked1m.int[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380, na.rm = TRUE) /
  #                                sum(dat$attr$diag.status == 1 & diag.time >= 3380, na.rm = TRUE)
  # dat$epi$cc.linked1m.int.B[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 1, na.rm = TRUE) /
  #                                  sum(diag.status == 1 & diag.time >= 3380 & race == 1, na.rm = TRUE)
  # dat$epi$cc.linked1m.int.H[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 2, na.rm = TRUE) /
  #                                  sum(diag.status == 1 & diag.time >= 3380 & race == 2, na.rm = TRUE)
  # dat$epi$cc.linked1m.int.W[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 3, na.rm = TRUE) /
  #                                  sum(diag.status == 1 & diag.time >= 3380 & race == 3, na.rm = TRUE)

  # dat$epi$cc.vsupp[at] <- sum(vl <= log10(200) & diag.status == 1 & diag.time >= 2, na.rm = TRUE) /
  #                         sum(diag.status == 1 & diag.time >= 2, na.rm = TRUE)
  # dat$epi$cc.vsupp.B[at] <- sum(vl <= log10(200) & diag.status == 1 &
  #                               diag.time >= 2 & race == 1, na.rm = TRUE) /
  #                           sum(diag.status == 1 & diag.time >= 2 & race == 1, na.rm = TRUE)
  # dat$epi$cc.vsupp.H[at] <- sum(vl <= log10(200) & diag.status == 1 &
  #                               diag.time >= 2 & race == 2, na.rm = TRUE) /
  #                           sum(diag.status == 1 & diag.time >= 2 & race == 2, na.rm = TRUE)
  # dat$epi$cc.vsupp.W[at] <- sum(vl <= log10(200) & diag.status == 1 &
  #                               diag.time >= 2 & race == 3, na.rm = TRUE) /
  #                           sum(diag.status == 1 & diag.time >= 2 & race == 3, na.rm = TRUE)

  # dat$epi$cc.vsupp.all[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2, na.rm = TRUE) /
  #                             sum(status == 1 & inf.time >= 2, na.rm = TRUE)
  # dat$epi$cc.vsupp.all.B[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE) /
  #                               sum(status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE)
  # dat$epi$cc.vsupp.all.H[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE) /
  #                               sum(status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)
  # dat$epi$cc.vsupp.all.W[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE) /
  #                               sum(status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE)

  # dat$epi$cc.HIV.mr[at] <- (dat$epi$dep.HIV[at]/dat$epi$i.num[at])*52

  # # Care continuum stats (secondary)
  # dat$epi$cc.test.int[at] <- mean(at - last.neg.test[diag.status == 0], na.rm = TRUE)
  # dat$epi$cc.test.int.B[at] <- mean(at - last.neg.test[diag.status == 0 & race == 1], na.rm = TRUE)
  # dat$epi$cc.test.int.H[at] <- mean(at - last.neg.test[diag.status == 0 & race == 2], na.rm = TRUE)
  # dat$epi$cc.test.int.W[at] <- mean(at - last.neg.test[diag.status == 0 & race == 3], na.rm = TRUE)

  # dat$epi$cc.dx.delay[at] <- mean(diag.time[diag.time >= 2] - inf.time[diag.time >= 2], na.rm = TRUE)
  # dat$epi$cc.dx.delay.B[at] <- mean(diag.time[diag.time >= 2 & race == 1] -
  #                                     inf.time[diag.time >= 2 & race == 1], na.rm = TRUE)
  # dat$epi$cc.dx.delay.H[at] <- mean(diag.time[diag.time >= 2 & race == 2] -
  #                                     inf.time[diag.time >= 2 & race == 2], na.rm = TRUE)
  # dat$epi$cc.dx.delay.W[at] <- mean(diag.time[diag.time >= 2 & race == 3] -
  #                                     inf.time[diag.time >= 2 & race == 3], na.rm = TRUE)

  # dat$epi$cc.dx.delay.int[at] <- mean(diag.time[diag.time >= 3380] -
  #                                       inf.time[diag.time >= 3380], na.rm = TRUE)
  # dat$epi$cc.dx.delay.int.B[at] <- mean(diag.time[diag.time >= 3380 & race == 1] -
  #                                     inf.time[diag.time >= 3380 & race == 1], na.rm = TRUE)
  # dat$epi$cc.dx.delay.int.H[at] <- mean(diag.time[diag.time >= 3380 & race == 2] -
  #                                     inf.time[diag.time >= 3380 & race == 2], na.rm = TRUE)
  # dat$epi$cc.dx.delay.int.W[at] <- mean(diag.time[diag.time >= 3380 & race == 3] -
  #                                     inf.time[diag.time >= 3380 & race == 3], na.rm = TRUE)

  # # HIV stage
  # # dat$epi$hstage.acute[at] <- sum(stage %in% 1:2 & diag.time >= 2, na.rm = TRUE) /
  # #                             sum(status == 1 & diag.time >= 2, na.rm = TRUE)
  # # dat$epi$hstage.chronic[at] <- sum(stage == 3 & diag.time >= 2, na.rm = TRUE) /
  # #                               sum(status == 1 & diag.time >= 2, na.rm = TRUE)
  # # dat$epi$hstage.aids[at] <- sum(stage == 4 & diag.time >= 2, na.rm = TRUE) /
  # #                            sum(status == 1 & diag.time >= 2, na.rm = TRUE)

  # dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  # dat$epi$prepElig.B[at] <- sum(prepElig == 1 & race == 1, na.rm = TRUE)
  # dat$epi$prepElig.H[at] <- sum(prepElig == 1 & race == 2, na.rm = TRUE)
  # dat$epi$prepElig.W[at] <- sum(prepElig == 1 & race == 3, na.rm = TRUE)

  # dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  # dat$epi$prepCurr.B[at] <- sum(prepStat == 1 & race == 1, na.rm = TRUE)
  # dat$epi$prepCurr.H[at] <- sum(prepStat == 1 & race == 2, na.rm = TRUE)
  # dat$epi$prepCurr.W[at] <- sum(prepStat == 1 & race == 3, na.rm = TRUE)

  # dat$epi$prepCurr.hadr[at] <- sum(prepStat == 1 & prepClass == 3, na.rm = TRUE)

  # STIs
  dat$epi$prev.gc[at] <- sum((rGC == 1 | uGC == 1), na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ct[at] <- sum((rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at]
  ir100.rgc <- (dat$epi$incid.rgc[at]/sum(rGC == 0, dat$epi$incid.rgc[at], na.rm = TRUE))*5200
  ir100.ugc <- (dat$epi$incid.ugc[at]/sum(uGC == 0, dat$epi$incid.ugc[at], na.rm = TRUE))*5200
  dat$epi$ir100.gc[at] <- ir100.rgc + ir100.ugc
  ir100.rct <- (dat$epi$incid.rct[at]/sum(rCT == 0, dat$epi$incid.rct[at], na.rm = TRUE))*5200
  ir100.uct <- (dat$epi$incid.uct[at]/sum(uCT == 0, dat$epi$incid.uct[at], na.rm = TRUE))*5200
  dat$epi$ir100.ct[at] <- ir100.rct + ir100.uct
  dat$epi$ir100.sti[at] <- dat$epi$ir100.gc[at] + dat$epi$ir100.ct[at]

  # # Partner Demographics
  # identified <- which(at - dat$attr$part.ident <= 52)
  # dat$epi$part.identified[at] <- length(identified)
  # screened <- which(at - dat$attr$part.scrnd <= 52)
  # dat$epi$part.screened[at] <- length(screened)
  # prep.start.part <- which(at - dat$attr$prep.start.part <= 52)
  # dat$epi$prep.part[at] <- length(prep.start.part)
  # part.tx.init.time <- get_attr(dat, "part.tx.init.time")
  # part.tx.reinit.time <- get_attr(dat, "part.tx.reinit.time")
  # dat$epi$part.init[at] <- sum(part.tx.init.time - diag.time <= 52 & diag.time >= 2, na.rm = TRUE)
  # dat$epi$part.reinit[at] <- sum(part.tx.reinit.time - diag.time <= 52 & diag.time >= 2, na.rm = TRUE)

  epi_trackers <- get_param(dat, "epi_trackers")
  for (nm in names(epi_trackers)) {
    dat <- set_epi(dat, nm, at, epi_trackers[[nm]](dat, at))
  }

  return(dat)
}
