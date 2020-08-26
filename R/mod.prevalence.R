
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

  active <- dat$attr$active
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  diag.stage <- dat$attr$diag.stage
  diag.time <- dat$attr$diag.time
  aids.time <- dat$attr$aids.time
  inf.time <- dat$attr$inf.time
  race <- dat$attr$race
  age <- dat$attr$age
  tx.init.time <- dat$attr$tx.init.time
  vl <- dat$attr$vl
  vl.last.usupp <- dat$attr$vl.last.usupp
  last.neg.test <- dat$attr$last.neg.test
  stage <- dat$attr$stage

  prepElig <- dat$attr$prepElig
  prepElig.la <- dat$attr$prepElig.la

  prepStat <- dat$attr$prepStat
  prepStat.la <- dat$attr$prepStat.la

  prepClass <- dat$attr$prepClass

  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT

  # Pop Size / Demog
  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$num.B[at] <- sum(race == 1, na.rm = TRUE)
  dat$epi$num.H[at] <- sum(race == 2, na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == 3, na.rm = TRUE)
  dat$epi$age.mean[at] <- mean(age, na.rm = TRUE)
  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$i.num.B[at] <- sum(status == 1 & race == 1, na.rm = TRUE)
  dat$epi$i.num.H[at] <- sum(status == 1 & race == 2, na.rm = TRUE)
  dat$epi$i.num.W[at] <- sum(status == 1 & race == 3, na.rm = TRUE)

  dat$epi$i.num.dx[at] <- sum(diag.status == 1, na.rm = TRUE)

  # Prev / Incid
  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.B[at] <- sum(race == 1 & status == 1, na.rm = TRUE) / sum(race == 1, na.rm = TRUE)
  dat$epi$i.prev.H[at] <- sum(race == 2 & status == 1, na.rm = TRUE) / sum(race == 2, na.rm = TRUE)
  dat$epi$i.prev.W[at] <- sum(race == 3 & status == 1, na.rm = TRUE) / sum(race == 3, na.rm = TRUE)

  dat$epi$i.prev.dx[at] <- sum(diag.status == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$i.prev.dx.B[at] <- sum(race == 1 & diag.status == 1, na.rm = TRUE) / sum(race == 1, na.rm = TRUE)
  dat$epi$i.prev.dx.H[at] <- sum(race == 2 & diag.status == 1, na.rm = TRUE) / sum(race == 2, na.rm = TRUE)
  dat$epi$i.prev.dx.W[at] <- sum(race == 3 & diag.status == 1, na.rm = TRUE) / sum(race == 3, na.rm = TRUE)

  dat$epi$ir100[at] <- (dat$epi$incid[at] / sum(status == 0, dat$epi$incid[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.B[at] <- (dat$epi$incid.B[at] / sum(status == 0 & race == 1, dat$epi$incid.B[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.H[at] <- (dat$epi$incid.H[at] / sum(status == 0 & race == 2, dat$epi$incid.H[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.W[at] <- (dat$epi$incid.W[at] / sum(status == 0 & race == 3, dat$epi$incid.W[at], na.rm = TRUE)) * 5200


  # Care continuum stats (primary)
  dat$epi$cc.dx[at] <- sum(diag.status == 1 & inf.time >= 2, na.rm = TRUE) /
                       sum(status == 1 & inf.time >= 2, na.rm = TRUE)
  dat$epi$cc.dx.B[at] <- sum(diag.status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE) /
                         sum(status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE)
  dat$epi$cc.dx.H[at] <- sum(diag.status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE) /
                         sum(status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)
  dat$epi$cc.dx.W[at] <- sum(diag.status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE) /
                         sum(status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.dx.aids[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                  aids.time - diag.time <= 52, na.rm = TRUE) /
                            sum(diag.status == 1 & inf.time >= 2, na.rm = TRUE)
  dat$epi$cc.dx.aids.B[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                    aids.time - diag.time <= 52 & race == 1, na.rm = TRUE) /
                              sum(diag.status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE)
  dat$epi$cc.dx.aids.H[at] <- sum(diag.status == 1 & stage == 4 &  inf.time >= 2 &
                                    aids.time - diag.time <= 52 & race == 2, na.rm = TRUE) /
                              sum(diag.status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)
  dat$epi$cc.dx.aids.W[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                    aids.time - diag.time <= 52 & race == 3, na.rm = TRUE) /
                              sum(diag.status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.linked1m[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2, na.rm = TRUE) /
                             sum(dat$attr$diag.status == 1 & diag.time >= 2, na.rm = TRUE)
  dat$epi$cc.linked1m.B[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 1, na.rm = TRUE) /
                               sum(diag.status == 1 & diag.time >= 2 & race == 1, na.rm = TRUE)
  dat$epi$cc.linked1m.H[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 2, na.rm = TRUE) /
                               sum(diag.status == 1 & diag.time >= 2 & race == 2, na.rm = TRUE)
  dat$epi$cc.linked1m.W[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 3, na.rm = TRUE) /
                               sum(diag.status == 1 & diag.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.linked1m.int[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380, na.rm = TRUE) /
                                 sum(dat$attr$diag.status == 1 & diag.time >= 3380, na.rm = TRUE)
  dat$epi$cc.linked1m.int.B[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 1, na.rm = TRUE) /
                                   sum(diag.status == 1 & diag.time >= 3380 & race == 1, na.rm = TRUE)
  dat$epi$cc.linked1m.int.H[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 2, na.rm = TRUE) /
                                   sum(diag.status == 1 & diag.time >= 3380 & race == 2, na.rm = TRUE)
  dat$epi$cc.linked1m.int.W[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 3, na.rm = TRUE) /
                                   sum(diag.status == 1 & diag.time >= 3380 & race == 3, na.rm = TRUE)

  dat$epi$cc.vsupp[at] <- sum(vl <= log10(200) & diag.status == 1 & diag.time >= 2, na.rm = TRUE) /
                          sum(diag.status == 1 & diag.time >= 2, na.rm = TRUE)
  dat$epi$cc.vsupp.B[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                diag.time >= 2 & race == 1, na.rm = TRUE) /
                            sum(diag.status == 1 & diag.time >= 2 & race == 1, na.rm = TRUE)
  dat$epi$cc.vsupp.H[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                diag.time >= 2 & race == 2, na.rm = TRUE) /
                            sum(diag.status == 1 & diag.time >= 2 & race == 2, na.rm = TRUE)
  dat$epi$cc.vsupp.W[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                diag.time >= 2 & race == 3, na.rm = TRUE) /
                            sum(diag.status == 1 & diag.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.vsupp.all[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2, na.rm = TRUE) /
                              sum(status == 1 & inf.time >= 2, na.rm = TRUE)
  dat$epi$cc.vsupp.all.B[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE) /
                                sum(status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE)
  dat$epi$cc.vsupp.all.H[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE) /
                                sum(status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)
  dat$epi$cc.vsupp.all.W[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE) /
                                sum(status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.vsupp.dur1y[at] <- 1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
                                         diag.status == 1, na.rm = TRUE) /
                                     sum(diag.status == 1 & diag.time >= 2, na.rm = TRUE))
  dat$epi$cc.vsupp.dur1y.B[at] <- 1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
                                         diag.status == 1 & race == 1, na.rm = TRUE) /
                                     sum(diag.status == 1 & diag.time >= 2 & race == 1, na.rm = TRUE))
  dat$epi$cc.vsupp.dur1y.H[at] <- 1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
                                           diag.status == 1 & race == 2, na.rm = TRUE) /
                                       sum(diag.status == 1 & diag.time >= 2 & race == 2, na.rm = TRUE))
  dat$epi$cc.vsupp.dur1y.W[at] <- 1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
                                           diag.status == 1 & race == 3, na.rm = TRUE) /
                                       sum(diag.status == 1 & diag.time >= 2 & race == 3, na.rm = TRUE))

  dat$epi$cc.HIV.mr[at] <- (dat$epi$dep.HIV[at]/dat$epi$i.num[at])*52

  # Care continuum stats (secondary)
  dat$epi$cc.test.int[at] <- mean(at - last.neg.test[diag.status == 0], na.rm = TRUE)
  dat$epi$cc.test.int.B[at] <- mean(at - last.neg.test[diag.status == 0 & race == 1], na.rm = TRUE)
  dat$epi$cc.test.int.H[at] <- mean(at - last.neg.test[diag.status == 0 & race == 2], na.rm = TRUE)
  dat$epi$cc.test.int.W[at] <- mean(at - last.neg.test[diag.status == 0 & race == 3], na.rm = TRUE)

  dat$epi$cc.dx.delay[at] <- mean(diag.time[diag.time >= 2] - inf.time[diag.time >= 2], na.rm = TRUE)
  dat$epi$cc.dx.delay.B[at] <- mean(diag.time[diag.time >= 2 & race == 1] -
                                      inf.time[diag.time >= 2 & race == 1], na.rm = TRUE)
  dat$epi$cc.dx.delay.H[at] <- mean(diag.time[diag.time >= 2 & race == 2] -
                                      inf.time[diag.time >= 2 & race == 2], na.rm = TRUE)
  dat$epi$cc.dx.delay.W[at] <- mean(diag.time[diag.time >= 2 & race == 3] -
                                      inf.time[diag.time >= 2 & race == 3], na.rm = TRUE)

  dat$epi$cc.dx.delay.int[at] <- mean(diag.time[diag.time >= 3380] - inf.time[diag.time >= 3380], na.rm = TRUE)
  dat$epi$cc.dx.delay.int.B[at] <- mean(diag.time[diag.time >= 3380 & race == 1] -
                                      inf.time[diag.time >= 3380 & race == 1], na.rm = TRUE)
  dat$epi$cc.dx.delay.int.H[at] <- mean(diag.time[diag.time >= 3380 & race == 2] -
                                      inf.time[diag.time >= 3380 & race == 2], na.rm = TRUE)
  dat$epi$cc.dx.delay.int.W[at] <- mean(diag.time[diag.time >= 3380 & race == 3] -
                                      inf.time[diag.time >= 3380 & race == 3], na.rm = TRUE)

  # dat$epi$cc.tx.any1y[at] <- sum((at - dat$attr$tx.period.last <= 52), na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.tx.any1y.B[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 1, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & race == 1, na.rm = TRUE)
  # dat$epi$cc.tx.any1y.H[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 2, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & race == 2, na.rm = TRUE)
  # dat$epi$cc.tx.any1y.W[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 3, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & race == 3, na.rm = TRUE)

  # dat$epi$cc.dx.delay[at] <- mean(dat$attr$diag.time - dat$attr$inf.time, na.rm = TRUE)
  # dat$epi$cc.testpy[at] <- 1-sum((at - dat$attr$last.neg.test) > 52 & status == 0,
  #     is.na(dat$attr$last.neg.test) & status == 0, na.rm = TRUE) /
  #   sum(status == 0)
  # dat$epi$cc.linked[at] <- sum(dat$attr$cuml.time.on.tx > 0, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.tx[at] <- sum(dat$attr$tx.status == 1, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.tx.ret3m[at] <- sum((at - dat$attr$tx.period.last) <= 52 &
  #       (dat$attr$tx.period.last - dat$attr$tx.period.first) > 13, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt1[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 1, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 1, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt2[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 2, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 2, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt3[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 3, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 3, na.rm = TRUE)


  # HIV screening outcomes
  dat$epi$mean.neg.tests[at] <- mean(dat$attr$num.neg.tests[diag.status == 0], na.rm = TRUE)
  dat$epi$mean.neg.tests.B[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & race == 1], na.rm = TRUE)
  dat$epi$mean.neg.tests.H[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & race == 2], na.rm = TRUE)
  dat$epi$mean.neg.tests.W[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & race == 3], na.rm = TRUE)

  dat$epi$test.past.year[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0, na.rm = TRUE) /
    sum(diag.status == 0, na.rm = TRUE)
  dat$epi$test.past.year.B[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & race == 1, na.rm = TRUE) /
    sum(diag.status == 0 & race == 1, na.rm = TRUE)
  dat$epi$test.past.year.H[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & race == 2, na.rm = TRUE) /
    sum(diag.status == 0 & race == 2, na.rm = TRUE)
  dat$epi$test.past.year.W[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & race == 3, na.rm = TRUE) /
    sum(diag.status == 0 & race == 3, na.rm = TRUE)

  # HIV stage
  dat$epi$hstage.acute[at] <- sum(stage %in% 1:2 & diag.time >= 2, na.rm = TRUE) /
                              sum(status == 1 & diag.time >= 2, na.rm = TRUE)
  dat$epi$hstage.chronic[at] <- sum(stage == 3 & diag.time >= 2, na.rm = TRUE) /
                                sum(status == 1 & diag.time >= 2, na.rm = TRUE)
  dat$epi$hstage.aids[at] <- sum(stage == 4 & diag.time >= 2, na.rm = TRUE) /
                             sum(status == 1 & diag.time >= 2, na.rm = TRUE)

  dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  dat$epi$prepElig.B[at] <- sum(prepElig == 1 & race == 1, na.rm = TRUE)
  dat$epi$prepElig.H[at] <- sum(prepElig == 1 & race == 2, na.rm = TRUE)
  dat$epi$prepElig.W[at] <- sum(prepElig == 1 & race == 3, na.rm = TRUE)

  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prep.laCurr[at] <- sum(prepStat.la == 1, na.rm = TRUE)

  dat$epi$prepCurr.B[at] <- sum(prepStat == 1 & race == 1, na.rm = TRUE)
  dat$epi$prepCurr.H[at] <- sum(prepStat == 1 & race == 2, na.rm = TRUE)
  dat$epi$prepCurr.W[at] <- sum(prepStat == 1 & race == 3, na.rm = TRUE)

  dat$epi$prepCurr.hadr[at] <- sum(prepStat == 1 & prepClass == 3, na.rm = TRUE)

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

  dat$epi$i.prev.prep[at] <- sum(status == 1 & (prepStat == 1 | prepStat.la == 1), na.rm = TRUE) /
                             sum(prepStat == 1 | prepStat.la == 1, na.rm = TRUE)
  dat$epi$i.prev.oprep[at] <- sum(status == 1 & (prepStat == 1), na.rm = TRUE) /
                              sum(prepStat == 1, na.rm = TRUE)
  dat$epi$i.prev.iprep[at] <- sum(status == 1 & (prepStat.la == 1), na.rm = TRUE) /
                              sum(prepStat.la == 1, na.rm = TRUE)

  dat$epi$ir100.prep[at] <- (dat$epi$incid.prep[at] /
                             sum(status == 0 & (prepStat == 1),
                                 na.rm = TRUE)) * 5200

  dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepElig.la[at] <- sum(prepElig.la == 1, na.rm = TRUE)
  dat$epi$prepCurr.la[at] <- sum(prepStat.la == 1, na.rm = TRUE)
  dat$epi$prepCurr.8w.la[at] <- sum(at - dat$attr$prepTimeLastInj <= 8, na.rm = TRUE)

  dat$epi$waning_laprep[at] <- epi_waning_laprep(dat, at, 1:3)
  dat$epi$waning_laprep.B[at] <- epi_waning_laprep(dat, at, 1)
  dat$epi$waning_laprep.H[at] <- epi_waning_laprep(dat, at, 2)
  dat$epi$waning_laprep.W[at] <- epi_waning_laprep(dat, at, 3)
  return(dat)
}

#' @export
#' @rdname prevalence_msm
prevalence_het <- function(dat, at) {

  status <- dat$attr$status
  male <- dat$attr$male
  age <- dat$attr$age

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  # Initialize vectors
  if (at == 1) {
    dat$epi$i.num <- rNA
    dat$epi$num <- rNA

    dat$epi$i.num.male <- rNA
    dat$epi$i.num.feml <- rNA
    dat$epi$i.prev.male <- rNA
    dat$epi$i.prev.feml <- rNA

    dat$epi$num.male <- rNA
    dat$epi$num.feml <- rNA
    dat$epi$meanAge <- rNA
    dat$epi$propMale <- rNA

    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA
  }

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$num[at] <- length(status)

  dat$epi$i.num.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.num.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE)
  dat$epi$i.prev.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE) /
    sum(male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE) /
    sum(male == 0, na.rm = TRUE)

  dat$epi$num.male[at] <- sum(male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(male == 0, na.rm = TRUE)
  dat$epi$meanAge[at] <- mean(age, na.rm = TRUE)
  dat$epi$propMale[at] <- mean(male, na.rm = TRUE)

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}

epi_waning_laprep <- function(dat, at, r_ind) {
  with(dat$attr, {
    cond <- at - prepTimeLastInj <= 52 & at - prepTimeLastInj > 8
    pop <- race %in% r_ind & prepStat.la != 1

    sum(cond & pop, na.rm = TRUE)
  })
}
