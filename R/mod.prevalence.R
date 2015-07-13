
#' @title Prevalence Calculations within Time Steps
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging.mard
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
#' @keywords module
#' @export
#'
prevalence.mard <- function(dat, at) {

  # Variables
  active <- dat$attr$active
  race <- dat$attr$race
  status <- dat$attr$status
  age <- dat$attr$age
  diag.status <- dat$attr$diag.status
  tx.status <- dat$attr$tx.status
  time.on.tx <- dat$attr$cum.time.on.tx
  vl.full.supp <- dat$param$vl.full.supp
  vl.part.supp <- dat$param$vl.part.supp
  vl <- dat$attr$vl
  tt.traj <- dat$attr$tt.traj
  stage <- dat$attr$stage
  prepStat <- dat$attr$prepStat
  riskg <- dat$attr$riskg

  prevfull <- dat$control$prevfull
  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  if (at == 1) {
    dat$epi$num <- dat$epi$num.B <- dat$epi$num.W <- rNA
    dat$epi$i.num <- dat$epi$i.num.B <- dat$epi$i.num.W <- rNA
    dat$epi$i.prev <- dat$epi$i.prev.B <- dat$epi$i.prev.W <- rNA
    dat$epi$nBirths <- dat$epi$nBirths.B <- dat$epi$nBirths.W <- rNA
    dat$epi$dth.gen <- dat$epi$dth.dis <- rNA
    dat$epi$dth.gen.B <- dat$epi$dth.gen.W <- rNA
    dat$epi$dth.dis.B <- dat$epi$dth.dis.W <- rNA
    dat$epi$incid <- dat$epi$incid.B <- dat$epi$incid.W <- rNA
    dat$epi$incid.acte <- dat$epi$incid.chrn <- dat$epi$incid.aids <- rNA
    dat$epi$tx.naive <- dat$epi$tx.naive.B <- dat$epi$tx.naive.W <- rNA
    dat$epi$tx.full.supp <- dat$epi$tx.full.supp.B <- dat$epi$tx.full.supp.W <- rNA
    dat$epi$acute <- dat$epi$acute.B <- dat$epi$acute.W <- rNA
    dat$epi$chronic <- dat$epi$chronic.B <- dat$epi$chronic.W <- rNA
    dat$epi$aids <- dat$epi$aids.B <- dat$epi$aids.W <- rNA
    dat$epi$i.prev.rg1 <- dat$epi$i.prev.rg2 <- dat$epi$i.prev.rg3 <-
      dat$epi$i.prev.rg4 <- dat$epi$i.prev.rg5 <- rNA

    if (prevfull == TRUE) {
      dat$epi$num.yng <- dat$epi$num.old <- rNA
      dat$epi$num.B.yng <- dat$epi$num.B.old <- rNA
      dat$epi$num.W.yng <- dat$epi$num.W.old <- rNA
      dat$epi$i.num.yng <- dat$epi$i.num.old <- rNA
      dat$epi$i.num.B.yng <- dat$epi$i.num.B.old <- rNA
      dat$epi$i.num.W.yng <- dat$epi$i.num.W.old <- rNA
      dat$epi$i.prev.yng <- dat$epi$i.prev.old <- rNA
      dat$epi$i.prev.B.yng <- dat$epi$i.prev.B.old <- rNA
      dat$epi$i.prev.W.yng <- dat$epi$i.prev.W.old <- rNA
      dat$epi$nBirths.yng <- dat$epi$nBirths.old <- rNA
      dat$epi$nBirths.B.yng <- dat$epi$nBirths.B.old <- rNA
      dat$epi$nBirths.W.yng <- dat$epi$nBirths.w.old <- rNA
      dat$epi$dth.gen.yng <- dat$epi$dth.gen.old <- rNA
      dat$epi$dth.gen.B.yng <- dat$epi$dth.gen.B.old <- rNA
      dat$epi$dth.gen.W.yng <- dat$epi$dth.gen.W.old <- rNA
      dat$epi$dth.dis.yng <- dat$epi$dth.dis.old <- rNA
      dat$epi$dth.dis.B.yng <- dat$epi$dth.dis.B.old <- rNA
      dat$epi$dth.dis.W.yng <- dat$epi$dth.dis.W.old <- rNA
      dat$epi$incid.yng <- dat$epi$incid.old <- rNA
      dat$epi$incid.B.yng <- dat$epi$incid.B.old <-  rNA
      dat$epi$incid.W.yng <- dat$epi$incid.W.old <-  rNA
      dat$epi$undiag <- dat$epi$undiag.B <- dat$epi$undiag.W <- rNA
      dat$epi$undiag.yng <- dat$epi$undiag.old <- rNA
      dat$epi$undiag.B.yng <- dat$epi$undiag.B.old <- rNA
      dat$epi$undiag.W.yng <- dat$epi$undiag.W.old <- rNA
      dat$epi$diag <- dat$epi$diag.B <- dat$epi$diag.W <- rNA
      dat$epi$diag.yng <- dat$epi$diag.old <- rNA
      dat$epi$diag.B.yng <- dat$epi$diag.B.old <- rNA
      dat$epi$diag.W.yng <- dat$epi$diag.W.old <- rNA
      dat$epi$tx.naive.yng <- dat$epi$tx.naive.old <- rNA
      dat$epi$tx.naive.B.yng <- dat$epi$tx.naive.B.old <- rNA
      dat$epi$tx.naive.W.yng <- dat$epi$tx.naive.W.old <- rNA
      dat$epi$tx.full.supp.yng <- dat$epi$tx.full.supp.old <- rNA
      dat$epi$tx.full.supp.B.yng <- dat$epi$tx.full.supp.B.old <- rNA
      dat$epi$tx.full.supp.W.yng <- dat$epi$tx.full.supp.W.old <- rNA
      dat$epi$tx.part.supp <- dat$epi$tx.part.supp.B <-
        dat$epi$tx.part.supp.W <- rNA
      dat$epi$tx.part.supp.yng <- dat$epi$tx.part.supp.old <- rNA
      dat$epi$tx.part.supp.B.yng <- dat$epi$tx.part.supp.B.old <- rNA
      dat$epi$tx.part.supp.W.yng <- dat$epi$tx.part.supp.W.old <- rNA
      dat$epi$tx.influx.full <- dat$epi$tx.influx.full.B <-
        dat$epi$tx.influx.full.W <- rNA
      dat$epi$tx.influx.full.yng <- dat$epi$tx.influx.full.old <- rNA
      dat$epi$tx.influx.full.B.yng <- dat$epi$tx.influx.full.B.old <- rNA
      dat$epi$tx.influx.full.W.yng <- dat$epi$tx.influx.full.W.old <- rNA
      dat$epi$tx.influx.part <- dat$epi$tx.influx.part.B <-
        dat$epi$tx.influx.part.W <- rNA
      dat$epi$tx.influx.part.yng <- dat$epi$tx.influx.part.old <- rNA
      dat$epi$tx.influx.part.B.yng <- dat$epi$tx.influx.part.B.old <- rNA
      dat$epi$tx.influx.part.W.yng <- dat$epi$tx.influx.part.W.old <- rNA
      dat$epi$off.tx <- dat$epi$off.tx.B <- dat$epi$off.tx.W <- rNA
      dat$epi$off.tx.yng <- dat$epi$off.tx.old <- rNA
      dat$epi$off.tx.B.yng <- dat$epi$off.tx.B.old <- rNA
      dat$epi$off.tx.W.yng <- dat$epi$off.tx.W.old <- rNA
      dat$epi$acute.yng <- dat$epi$acute.old <- rNA
      dat$epi$acute.B.yng <- dat$epi$acute.B.old <- rNA
      dat$epi$acute.W.yng <- dat$epi$acute.W.old <- rNA
      dat$epi$chronic.yng <- dat$epi$chronic.old <- rNA
      dat$epi$chronic.B.yng <- dat$epi$chronic.B.old <- rNA
      dat$epi$chronic.W.yng <- dat$epi$chronic.W.old <- rNA
      dat$epi$aids.yng <- dat$epi$aids.old <- rNA
      dat$epi$aids.B.yng <- dat$epi$aids.B.old <- rNA
      dat$epi$aids.W.yng <- dat$epi$aids.W.old <- rNA
    }

    dat$epi$prepCurr <- dat$epi$prepEver <-
      dat$epi$prepCov <- dat$epi$prepStart <- rNA
    dat$epi$incid.prep0 <- dat$epi$incid.prep1 <-
      dat$epi$i.num.prep0 <- dat$epi$i.num.prep1 <-
      dat$epi$i.prev.prep0 <- dat$epi$i.prev.prep1 <- rNA
  }


  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$num.B[at] <- sum(active == 1 & race == "B", na.rm = TRUE)
  dat$epi$num.W[at] <- sum(active == 1 & race == "W", na.rm = TRUE)
  dat$epi$i.num[at] <- sum(active == 1 & status == 1, na.rm = TRUE)
  dat$epi$i.num.B[at] <- sum(active == 1 & status == 1 & race == "B", na.rm = TRUE)
  dat$epi$i.num.W[at] <- sum(active == 1 & status == 1 & race == "W", na.rm = TRUE)
  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.B[at] <- dat$epi$i.num.B[at] / dat$epi$num.B[at]
  dat$epi$i.prev.W[at] <- dat$epi$i.num.W[at] / dat$epi$num.W[at]
  dat$epi$i.prev.rg1 <- sum(active == 1 & status == 1 & riskg == 1, na.rm = TRUE)/
                        sum(active == 1 & riskg == 1, na.rm = TRUE)
  dat$epi$i.prev.rg2 <- sum(active == 1 & status == 1 & riskg == 2, na.rm = TRUE)/
                        sum(active == 1 & riskg == 2, na.rm = TRUE)
  dat$epi$i.prev.rg3 <- sum(active == 1 & status == 1 & riskg == 3, na.rm = TRUE)/
                        sum(active == 1 & riskg == 3, na.rm = TRUE)
  dat$epi$i.prev.rg4 <- sum(active == 1 & status == 1 & riskg == 4, na.rm = TRUE)/
                        sum(active == 1 & riskg == 4, na.rm = TRUE)
  dat$epi$i.prev.rg5 <- sum(active == 1 & status == 1 & riskg == 5, na.rm = TRUE)/
                        sum(active == 1 & riskg == 5, na.rm = TRUE)

  dat$epi$tx.naive[at] <- sum(active == 1 & time.on.tx == 0, na.rm = TRUE)
  dat$epi$tx.naive.B[at] <- sum(active == 1 & time.on.tx == 0 & race == "B",
                                na.rm = TRUE)
  dat$epi$tx.naive.W[at] <- sum(active == 1 & time.on.tx == 0 & race == "W",
                                na.rm = TRUE)
  dat$epi$tx.full.supp[at] <- sum(active == 1 & tx.status == 1 &
                                    vl %in% vl.full.supp, na.rm = TRUE)
  dat$epi$tx.full.supp.B[at] <- sum(active == 1 & tx.status == 1 &
                                      vl %in% vl.full.supp & race == "B",
                                    na.rm = TRUE)
  dat$epi$tx.full.supp.W[at] <- sum(active == 1 & tx.status == 1 &
                                      vl %in% vl.full.supp & race == "W",
                                    na.rm = TRUE)
  dat$epi$acute[at] <- sum(active == 1 & stage %in% c("AR", "AF"), na.rm = TRUE)
  dat$epi$acute.B[at] <- sum(active == 1 & stage %in% c("AR", "AF") & race == "B",
                             na.rm = TRUE)
  dat$epi$acute.W[at] <- sum(active == 1 & stage %in% c("AR", "AF") & race == "W",
                             na.rm = TRUE)
  dat$epi$chronic[at] <- sum(active == 1 & stage == "C", na.rm = TRUE)
  dat$epi$chronic.B[at] <- sum(active == 1 & stage == "C" & race == "B",
                               na.rm = TRUE)
  dat$epi$chronic.W[at] <- sum(active == 1 & stage == "C" & race == "W",
                               na.rm = TRUE)
  dat$epi$aids[at] <- sum(active == 1 & stage == "D", na.rm = TRUE)
  dat$epi$aids.B[at] <- sum(active == 1 & stage == "D" & race == "B",
                            na.rm = TRUE)
  dat$epi$aids.W[at] <- sum(active == 1 & stage == "D" & race == "W",
                            na.rm = TRUE)


  if (prevfull == TRUE) {
    dat$epi$num.yng[at] <- sum(active == 1 & age < 30, na.rm = TRUE)
    dat$epi$num.old[at] <- sum(active == 1 & age >= 30, na.rm = TRUE)
    dat$epi$num.B.yng[at] <- sum(active == 1 & race == "B" & age < 30, na.rm = TRUE)
    dat$epi$num.B.old[at] <- sum(active == 1 & race == "B" & age >= 30, na.rm = TRUE)
    dat$epi$num.W.yng[at] <- sum(active == 1 & race == "W" & age < 30, na.rm = TRUE)
    dat$epi$num.W.old[at] <- sum(active == 1 & race == "W" & age >= 30, na.rm = TRUE)
    dat$epi$i.num.yng[at] <- sum(active == 1 & status == 1 & age < 30, na.rm = TRUE)
    dat$epi$i.num.old[at] <- sum(active == 1 & status == 1 & age >= 30, na.rm = TRUE)
    dat$epi$i.num.B.yng[at] <- sum(active == 1 & status == 1 & race == "B" &
                                     age < 30, na.rm = TRUE)
    dat$epi$i.num.B.old[at] <- sum(active == 1 & status == 1 & race == "B" &
                                     age >= 30, na.rm = TRUE)
    dat$epi$i.num.W.yng[at] <- sum(active == 1 & status == 1 & race == "W" &
                                     age < 30, na.rm = TRUE)
    dat$epi$i.num.W.old[at] <- sum(active == 1 & status == 1 & race == "W" &
                                     age >= 30, na.rm = TRUE)
    dat$epi$i.prev.yng[at] <- dat$epi$i.num.yng[at] / dat$epi$num.yng[at]
    dat$epi$i.prev.old[at] <- dat$epi$i.num.old[at] / dat$epi$num.old[at]
    dat$epi$i.prev.B.yng[at] <- dat$epi$i.num.B.yng[at] / dat$epi$num.B.yng[at]
    dat$epi$i.prev.B.old[at] <- dat$epi$i.num.B.old[at] / dat$epi$num.B.old[at]
    dat$epi$i.prev.W.yng[at] <- dat$epi$i.num.W.yng[at] / dat$epi$num.W.yng[at]
    dat$epi$i.prev.W.old[at] <- dat$epi$i.num.W.old[at] / dat$epi$num.W.old[at]
    dat$epi$undiag[at] <- sum(active == 1 & status == 1 & diag.status == 0,
                              na.rm = TRUE)
    dat$epi$undiag.B[at] <- sum(active == 1 & status == 1 & diag.status == 0 &
                                  race == "B", na.rm = TRUE)
    dat$epi$undiag.W[at] <- sum(active == 1 & status == 1 & diag.status == 0 &
                                  race == "W", na.rm = TRUE)
    dat$epi$undiag.yng[at] <- sum(active == 1 & status == 1 & diag.status == 0 &
                                    age < 30, na.rm = TRUE)
    dat$epi$undiag.old[at] <- sum(active == 1 & status == 1 & diag.status == 0 &
                                    age >= 30, na.rm = TRUE)
    dat$epi$undiag.B.yng[at] <- sum(active == 1 & status == 1 & diag.status == 0 &
                                      race == "B" & age < 30, na.rm = TRUE)
    dat$epi$undiag.B.old[at] <- sum(active == 1 & status == 1 & diag.status == 0 &
                                      race == "B" & age >= 30, na.rm = TRUE)
    dat$epi$undiag.W.yng[at] <- sum(active == 1 & status == 1 & diag.status == 0 &
                                      race == "W" & age < 30, na.rm = TRUE)
    dat$epi$undiag.W.old[at] <- sum(active == 1 & status == 1 & diag.status == 0 &
                                      race == "W" & age >= 30, na.rm = TRUE)
    dat$epi$diag[at] <- sum(active == 1 & status == 1 & diag.status == 1,
                            na.rm = TRUE)
    dat$epi$diag.B[at] <- sum(active == 1 & status == 1 & diag.status == 1 &
                                race == "B", na.rm = TRUE)
    dat$epi$diag.W[at] <- sum(active == 1 & status == 1 & diag.status == 1 &
                                race == "W", na.rm = TRUE)
    dat$epi$diag.yng[at] <- sum(active == 1 & status == 1 & diag.status == 1 &
                                  age < 30, na.rm = TRUE)
    dat$epi$diag.old[at] <- sum(active == 1 & status == 1 & diag.status == 1 &
                                  age >= 30, na.rm = TRUE)
    dat$epi$diag.B.yng[at] <- sum(active == 1 & status == 1 & diag.status == 1 &
                                    race == "B" & age < 30, na.rm = TRUE)
    dat$epi$diag.B.old[at] <- sum(active == 1 & status == 1 & diag.status == 1 &
                                    race == "B" & age >= 30, na.rm = TRUE)
    dat$epi$diag.W.yng[at] <- sum(active == 1 & status == 1 & diag.status == 1 &
                                    race == "W" & age < 30, na.rm = TRUE)
    dat$epi$diag.W.old[at] <- sum(active == 1 & status == 1 & diag.status == 1 &
                                    race == "W" & age >= 30, na.rm = TRUE)
    dat$epi$tx.naive.yng[at] <- sum(active == 1 & time.on.tx == 0 & age < 30,
                                    na.rm = TRUE)
    dat$epi$tx.naive.old[at] <- sum(active == 1 & time.on.tx == 0 & age >= 30,
                                    na.rm = TRUE)
    dat$epi$tx.naive.B.yng[at] <- sum(active == 1 & time.on.tx == 0 & race == "B" &
                                        age < 30, na.rm = TRUE)
    dat$epi$tx.naive.B.old[at] <- sum(active == 1 & time.on.tx == 0 & race == "B" &
                                        age >= 30, na.rm = TRUE)
    dat$epi$tx.naive.W.yng[at] <- sum(active == 1 & time.on.tx == 0 & race == "W" &
                                        age < 30, na.rm = TRUE)
    dat$epi$tx.naive.W.old[at] <- sum(active == 1 & time.on.tx == 0 & race == "W" &
                                        age >= 30, na.rm = TRUE)
    dat$epi$tx.full.supp.yng[at] <- sum(active == 1 & tx.status == 1 &
                                          vl %in% vl.full.supp & age < 30,
                                        na.rm = TRUE)
    dat$epi$tx.full.supp.old[at] <- sum(active == 1 & tx.status == 1 &
                                          vl %in% vl.full.supp & age >= 30,
                                        na.rm = TRUE)
    dat$epi$tx.full.supp.B.yng[at] <- sum(active == 1 & tx.status == 1 &
                                            vl %in% vl.full.supp & race == "B" &
                                            age < 30, na.rm = TRUE)
    dat$epi$tx.full.supp.B.old[at] <- sum(active == 1 & tx.status == 1 &
                                            vl %in% vl.full.supp & race == "B" &
                                            age >= 30, na.rm = TRUE)
    dat$epi$tx.part.supp[at] <- sum(active == 1 & tx.status == 1 &
                                      vl %in% vl.part.supp, na.rm = TRUE)
    dat$epi$tx.part.supp.B[at] <- sum(active == 1 & tx.status == 1 &
                                        vl %in% vl.part.supp & race == "B",
                                      na.rm = TRUE)
    dat$epi$tx.part.supp.W[at] <- sum(active == 1 & tx.status == 1 &
                                        vl %in% vl.part.supp & race == "W",
                                      na.rm = TRUE)
    dat$epi$tx.part.supp.yng[at] <- sum(active == 1 & tx.status == 1 &
                                          vl %in% vl.part.supp & age < 30,
                                        na.rm = TRUE)
    dat$epi$tx.part.supp.old[at] <- sum(active == 1 & tx.status == 1 &
                                          vl %in% vl.part.supp & age >= 30,
                                        na.rm = TRUE)
    dat$epi$tx.part.supp.B.yng[at] <- sum(active == 1 & tx.status == 1 &
                                            vl %in% vl.part.supp & race == "B" &
                                            age < 30, na.rm = TRUE)
    dat$epi$tx.part.supp.B.old[at] <- sum(active == 1 & tx.status == 1 &
                                            vl %in% vl.part.supp & race == "B" &
                                            age >= 30, na.rm = TRUE)
    dat$epi$tx.part.supp.W.yng[at] <- sum(active == 1 & tx.status == 1 &
                                            vl %in% vl.part.supp & race == "W" &
                                            age < 30, na.rm = TRUE)
    dat$epi$tx.part.supp.W.old[at] <- sum(active == 1 & tx.status == 1 &
                                            vl %in% vl.part.supp & race == "W" &
                                            age >= 30, na.rm = TRUE)
    dat$epi$tx.influx.full[at] <- sum(active == 1 & tx.status == 1 &
                                        !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                        tt.traj == "YF", na.rm = TRUE)
    dat$epi$tx.influx.full.B[at] <- sum(active == 1 & tx.status == 1 &
                                          !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                          tt.traj == "YF" & race == "B",
                                        na.rm = TRUE)
    dat$epi$tx.influx.full.W[at] <- sum(active == 1 & tx.status == 1 &
                                          !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                          tt.traj == "YF" & race == "W",
                                        na.rm = TRUE)
    dat$epi$tx.influx.full.yng[at] <- sum(active == 1 & tx.status == 1 &
                                            !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                            tt.traj == "YF" & age < 30,
                                          na.rm = TRUE)
    dat$epi$tx.influx.full.old[at] <- sum(active == 1 & tx.status == 1 &
                                            !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                            tt.traj == "YF" & age >= 30,
                                          na.rm = TRUE)
    dat$epi$tx.influx.full.B.yng[at] <- sum(active == 1 & tx.status == 1 &
                                              !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                              tt.traj == "YF" & race == "B" &
                                              age < 30, na.rm = TRUE)
    dat$epi$tx.influx.full.B.old[at] <- sum(active == 1 & tx.status == 1 &
                                              !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                              tt.traj == "YF" & race == "B" &
                                              age >= 30, na.rm = TRUE)
    dat$epi$tx.influx.full.W.yng[at] <- sum(active == 1 & tx.status == 1 &
                                              !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                              tt.traj == "YF" & race == "W" &
                                              age < 30, na.rm = TRUE)
    dat$epi$tx.influx.full.W.old[at] <- sum(active == 1 & tx.status == 1 &
                                              !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                              tt.traj == "YF" & race == "W" &
                                              age >= 30, na.rm = TRUE)
    dat$epi$tx.influx.part[at] <- sum(active == 1 & tx.status == 1 &
                                        !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                        tt.traj == "YP", na.rm = TRUE)
    dat$epi$tx.influx.part.B[at] <- sum(active == 1 & tx.status == 1 &
                                          !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                          tt.traj == "YP" & race == "B",
                                        na.rm = TRUE)
    dat$epi$tx.influx.part.W[at] <- sum(active == 1 & tx.status == 1 &
                                          !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                          tt.traj == "YP" & race == "W",
                                        na.rm = TRUE)
    dat$epi$tx.influx.part.yng[at] <- sum(active == 1 & tx.status == 1 &
                                            !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                            tt.traj == "YP" & age < 30,
                                          na.rm = TRUE)
    dat$epi$tx.influx.part.old[at] <- sum(active == 1 & tx.status == 1 &
                                            !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                            tt.traj == "YP" & age >= 30,
                                          na.rm = TRUE)
    dat$epi$tx.influx.part.B.yng[at] <- sum(active == 1 & tx.status == 1 &
                                              !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                              tt.traj == "YP" & race == "B" &
                                              age < 30, na.rm = TRUE)
    dat$epi$tx.influx.part.B.old[at] <- sum(active == 1 & tx.status == 1 &
                                              !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                              tt.traj == "YP" & race == "B" &
                                              age >= 30, na.rm = TRUE)
    dat$epi$tx.influx.part.W.yng[at] <- sum(active == 1 & tx.status == 1 &
                                              !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                              tt.traj == "YP" & race == "W" &
                                              age < 30, na.rm = TRUE)
    dat$epi$tx.influx.part.W.old[at] <- sum(active == 1 & tx.status == 1 &
                                              !(vl %in% c(vl.full.supp, vl.part.supp)) &
                                              tt.traj == "YP" & race == "W" &
                                              age >= 30, na.rm = TRUE)

    dat$epi$off.tx[at] <- sum(active == 1 & tx.status == 0 & time.on.tx > 0,
                              na.rm = TRUE)
    dat$epi$off.tx.B[at] <- sum(active == 1 & tx.status == 0 & time.on.tx > 0 &
                                  race == "B", na.rm = TRUE)
    dat$epi$off.tx.W[at] <- sum(active == 1 & tx.status == 0 & time.on.tx > 0 &
                                  race == "W", na.rm = TRUE)
    dat$epi$off.tx.yng[at] <- sum(active == 1 & tx.status == 0 & time.on.tx > 0 &
                                    age < 30, na.rm = TRUE)
    dat$epi$off.tx.old[at] <- sum(active == 1 & tx.status == 0 & time.on.tx > 0 &
                                    age >= 30, na.rm = TRUE)
    dat$epi$off.tx.B.yng[at] <- sum(active == 1 & tx.status == 0 & time.on.tx > 0 &
                                      race == "B" & age < 30, na.rm = TRUE)
    dat$epi$off.tx.B.old[at] <- sum(active == 1 & tx.status == 0 & time.on.tx > 0 &
                                      race == "B" & age >= 30, na.rm = TRUE)
    dat$epi$off.tx.W.yng[at] <- sum(active == 1 & tx.status == 0 & time.on.tx > 0 &
                                      race == "W" & age < 30, na.rm = TRUE)
    dat$epi$off.tx.W.old[at] <- sum(active == 1 & tx.status == 0 & time.on.tx > 0 &
                                      race == "W" & age >= 30, na.rm = TRUE)
    dat$epi$acute.yng[at] <- sum(active == 1 & stage %in% c("AR", "AF") & age < 30,
                                 na.rm = TRUE)
    dat$epi$acute.old[at] <- sum(active == 1 & stage %in% c("AR", "AF") & age >= 30,
                                 na.rm = TRUE)
    dat$epi$acute.B.yng[at] <- sum(active == 1 & stage %in% c("AR", "AF") &
                                     race == "B" & age < 30, na.rm = TRUE)
    dat$epi$acute.B.old[at] <- sum(active == 1 & stage %in% c("AR", "AF") &
                                     race == "B" & age >= 30, na.rm = TRUE)
    dat$epi$acute.W.yng[at] <- sum(active == 1 & stage %in% c("AR", "AF") &
                                     race == "W" & age < 30, na.rm = TRUE)
    dat$epi$acute.W.old[at] <- sum(active == 1 & stage %in% c("AR", "AF") &
                                     race == "W" & age >= 30, na.rm = TRUE)
    dat$epi$chronic.yng[at] <- sum(active == 1 & stage == "C" & age < 30,
                                   na.rm = TRUE)
    dat$epi$chronic.old[at] <- sum(active == 1 & stage == "C" & age >= 30,
                                   na.rm = TRUE)
    dat$epi$chronic.B.yng[at] <- sum(active == 1 & stage == "C" & race == "B" &
                                       age < 30, na.rm = TRUE)
    dat$epi$chronic.B.old[at] <- sum(active == 1 & stage == "C" & race == "B" &
                                       age >= 30, na.rm = TRUE)
    dat$epi$chronic.W.yng[at] <- sum(active == 1 & stage == "C" & race == "W" &
                                       age < 30, na.rm = TRUE)
    dat$epi$chronic.W.old[at] <- sum(active == 1 & stage == "C" & race == "W" &
                                       age >= 30, na.rm = TRUE)
    dat$epi$aids.yng[at] <- sum(active == 1 & stage == "D" & age < 30,
                                na.rm = TRUE)
    dat$epi$aids.old[at] <- sum(active == 1 & stage == "D" & age >= 30,
                                na.rm = TRUE)
    dat$epi$aids.B.yng[at] <- sum(active == 1 & stage == "D" & race == "B" &
                                    age < 30, na.rm = TRUE)
    dat$epi$aids.B.old[at] <- sum(active == 1 & stage == "D" & race == "B" &
                                    age >= 30, na.rm = TRUE)
    dat$epi$aids.W.yng[at] <- sum(active == 1 & stage == "D" & race == "W" &
                                    age < 30, na.rm = TRUE)
    dat$epi$aids.W.old[at] <- sum(active == 1 & stage == "D" & race == "W" &
                                    age >= 30, na.rm = TRUE)

  }


  dat$epi$prepCurr[at] <- sum(active == 1 & prepStat == 1, na.rm = TRUE)
  dat$epi$prepEver[at] <- sum(active == 1 & dat$attr$prepEver == 1, na.rm = TRUE)
  dat$epi$i.num.prep0[at] <- sum(active == 1 & (is.na(prepStat) | prepStat == 0) &
                                   status == 1, na.rm = TRUE)
  dat$epi$i.num.prep1[at] <- sum(active == 1 & prepStat == 1 & status == 1,
                                 na.rm = TRUE)
  dat$epi$i.prev.prep0[at] <- dat$epi$i.num.prep0[at] /
                              sum(active == 1 & (is.na(prepStat) | prepStat == 0),
                                  na.rm = TRUE)
  if (at == 1) {
    dat$epi$i.prev.prep1[1] <- 0
  } else {
    dat$epi$i.prev.prep1[at] <- dat$epi$i.num.prep1[at] /
                                sum(active == 1 & prepStat == 1, na.rm = TRUE)
  }

  return(dat)
}
