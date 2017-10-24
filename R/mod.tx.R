
#' @title Treatment Module
#'
#' @description Module function for anti-retroviral treatment initiation and
#'              adherence over time.
#'
#' @inheritParams aging_camplc
#'
#' @details
#' Persons enter into the simulation with one of four ART "patterns": never
#' tested, tested but never treated, treated and achieving partial HIV viral
#' suppression, or treated with full viral suppression (these types are stored
#' as individual-level attributes in \code{tt.traj}). This module initiates ART
#' for treatment naive persons in the latter two types, and then cycles them on
#' and off treatment conditional on empirical race-specific adherence rates. ART
#' initiation, non-adherence, and restarting are all stochastically simulated
#' based on binomial statistical models.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{tx.status},
#' \code{tx.init.time}, \code{cum.time.on.tx}, \code{cum.time.off.tx} attributes.
#'
#' @keywords module msm
#'
#' @export
#'
tx_msm <- function(dat, at) {

  ## Variables

  # Attributes
  race <- dat$attr$race
  status <- dat$attr$status
  tx.status <- dat$attr$tx.status
  diag.status <- dat$attr$diag.status
  tt.traj <- dat$attr$tt.traj
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  stage <- dat$attr$stage
  asmm <- dat$attr$asmm

  # Parameters
  tx.init.B.prob.msm <- dat$param$tx.init.B.prob
  tx.init.W.prob.msm <- dat$param$tx.init.W.prob
  tx.halt.B.prob.msm <- dat$param$tx.halt.B.prob
  tx.halt.W.prob.msm <- dat$param$tx.halt.W.prob
  tx.reinit.B.prob.msm <- dat$param$tx.reinit.B.prob
  tx.reinit.W.prob.msm <- dat$param$tx.reinit.W.prob
  
  tx.init.B.prob.asmm <- dat$param$tx.init.B.prob.asmm
  tx.init.W.prob.asmm <- dat$param$tx.init.W.prob.asmm
  tx.halt.B.prob.asmm <- dat$param$tx.halt.B.prob.asmm
  tx.halt.W.prob.asmm <- dat$param$tx.halt.W.prob.asmm
  tx.reinit.B.prob.asmm <- dat$param$tx.reinit.B.prob.asmm
  tx.reinit.W.prob.asmm <- dat$param$tx.reinit.W.prob.asmm


  ## Initiation
  tx.init.elig.B.msm <- which(race == "B" & status == 1 & asmm == 0 &
                          tx.status == 0 & diag.status == 1 &
                          tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                          stage != 4)
  tx.init.B.msm <- tx.init.elig.B.msm[rbinom(length(tx.init.elig.B.msm), 1,
                                     tx.init.B.prob.msm) == 1]

  tx.init.elig.W.msm <- which(race == "W" & status == 1 & asmm == 0 &
                          tx.status == 0 & diag.status == 1 &
                          tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                          stage != 4)
  tx.init.W.msm <- tx.init.elig.W.msm[rbinom(length(tx.init.elig.W.msm), 1,
                                     tx.init.W.prob.msm) == 1]

  
  tx.init.elig.B.asmm <- which(race == "B" & status == 1 & asmm == 1 &
                            tx.status == 0 & diag.status == 1 &
                            tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                            stage != 4)
  tx.init.B.asmm <- tx.init.elig.B.asmm[rbinom(length(tx.init.elig.B.asmm), 1,
                                     tx.init.B.prob.asmm) == 1]
  
  tx.init.elig.W.asmm <- which(race == "W" & status == 1 & asmm == 0 &
                            tx.status == 0 & diag.status == 1 &
                            tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                            stage != 4)
  tx.init.W.asmm <- tx.init.elig.W.asmm[rbinom(length(tx.init.elig.W.asmm), 1,
                                     tx.init.W.prob.asmm) == 1]
  
  tx.init.msm <- c(tx.init.B.msm, tx.init.W.msm)
  tx.init.asmm <- c(tx.init.B.asmm, tx.init.W.asmm)
  tx.init <- c(tx.init.B.msm, tx.init.W.msm, tx.init.B.asmm, tx.init.W.asmm)

  dat$attr$tx.status[tx.init] <- 1
  dat$attr$tx.init.time[tx.init] <- at


  ## Halting
  tx.halt.elig.B.msm <- which(race == "B" & tx.status == 1 & asmm ==0)
  tx.halt.B.msm <- tx.halt.elig.B.msm[rbinom(length(tx.halt.elig.B.msm), 1,
                                     tx.halt.B.prob.msm) == 1]

  tx.halt.elig.W.msm <- which(race == "W" & tx.status == 1 & asmm ==0)
  tx.halt.W.msm <- tx.halt.elig.W.msm[rbinom(length(tx.halt.elig.W.msm),
                                     1, tx.halt.W.prob.msm) == 1]
  
  tx.halt.elig.B.asmm <- which(race == "B" & tx.status == 1 & asmm ==1)
  tx.halt.B.asmm <- tx.halt.elig.B.asmm[rbinom(length(tx.halt.elig.B.asmm), 1,
                                     tx.halt.B.prob.asmm) == 1]
  
  tx.halt.elig.W.asmm <- which(race == "W" & tx.status == 1 & asmm ==0)
  tx.halt.W.asmm <- tx.halt.elig.W.asmm[rbinom(length(tx.halt.elig.W.asmm),
                                     1, tx.halt.W.prob.asmm) == 1]
  
  tx.halt.msm <- c(tx.halt.B.msm, tx.halt.W.msm)
  tx.halt.asmm <- c(tx.halt.B.asmm, tx.halt.W.asmm)
  tx.halt <- c(tx.halt.B.msm, tx.halt.W.msm, tx.halt.B.asmm, tx.halt.W.asmm)
  
  dat$attr$tx.status[tx.halt] <- 0


  ## Restarting
  tx.reinit.elig.B.msm <- which(race == "B" & tx.status == 0 &
                            cum.time.on.tx > 0 & stage != 4)
  tx.reinit.B.msm <- tx.reinit.elig.B.msm[rbinom(length(tx.reinit.elig.B.msm),
                                         1, tx.reinit.B.prob.msm) == 1]

  tx.reinit.elig.W.msm <- which(race == "W" & tx.status == 0 &
                            cum.time.on.tx > 0 & stage != 4)
  tx.reinit.W.msm <- tx.reinit.elig.W.msm[rbinom(length(tx.reinit.elig.W.msm),
                                         1, tx.reinit.W.prob.msm) == 1]

  
  tx.reinit.elig.B.asmm <- which(race == "B" & tx.status == 0 &
                              cum.time.on.tx > 0 & stage != 4)
  tx.reinit.B.asmm <- tx.reinit.elig.B.asmm[rbinom(length(tx.reinit.elig.B.asmm),
                                         1, tx.reinit.B.prob.asmm) == 1]
  
  tx.reinit.elig.W.asmm <- which(race == "W" & tx.status == 0 &
                              cum.time.on.tx > 0 & stage != 4)
  tx.reinit.W.asmm <- tx.reinit.elig.W.asmm[rbinom(length(tx.reinit.elig.W.asmm),
                                         1, tx.reinit.W.prob.asmm) == 1]
  
  
  tx.reinit.msm <- c(tx.reinit.B.msm, tx.reinit.W.msm)
  tx.reinit.asmm <- c(tx.reinit.B.msm, tx.reinit.W.msm)
  tx.reinit <- c(tx.reinit.B.msm, tx.reinit.W.msm, tx.reinit.B.msm, tx.reinit.W.msm)
  
  dat$attr$tx.status[tx.reinit] <- 1


  ## Other output
  dat$attr$cum.time.on.tx <- dat$attr$cum.time.on.tx +
                             ((dat$attr$tx.status == 1) %in% TRUE)
  dat$attr$cum.time.off.tx <- dat$attr$cum.time.off.tx +
                              ((dat$attr$tx.status == 0) %in% TRUE)

  ## Summary statistics
  dat$epi$tx.init.inc.msm[at] <- length(tx.init.msm)
  dat$epi$tx.halt.inc.msm[at] <- length(tx.halt.msm)
  dat$epi$tx.resm.inc.msm[at] <- length(tx.reinit.msm)
  
  dat$epi$tx.init.inc.asmm[at] <- length(tx.init.asmm)
  dat$epi$tx.halt.inc.asmm[at] <- length(tx.halt.asmm)
  dat$epi$tx.resm.inc.asmm[at] <- length(tx.reinit.asmm)
  
  dat$epi$tx.init.inc[at] <- length(tx.init)
  dat$epi$tx.halt.inc[at] <- length(tx.halt)
  dat$epi$tx.resm.inc[at] <- length(tx.reinit)

  return(dat)
}

