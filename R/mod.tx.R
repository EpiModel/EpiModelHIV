
#' @title Treatment Module
#'
#' @description Module function for anti-retroviral treatment initiation and
#'              adherence over time.
#'
#' @inheritParams aging.mard
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
#' @keywords module
#' @export
#'
tx.mard <- function(dat, at) {

  ## Variables

  # Attributes
  active <- dat$attr$active
  race <- dat$attr$race
  status <- dat$attr$status
  tx.status <- dat$attr$tx.status
  diag.status <- dat$attr$diag.status
  tt.traj <- dat$attr$tt.traj
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  stage <- dat$attr$stage

  # Parameters
  tx.init.B.prob <- dat$param$tx.init.B.prob
  tx.init.W.prob <- dat$param$tx.init.W.prob
  tx.halt.B.prob <- dat$param$tx.halt.B.prob
  tx.halt.W.prob <- dat$param$tx.halt.W.prob
  tx.reinit.B.prob <- dat$param$tx.reinit.B.prob
  tx.reinit.W.prob <- dat$param$tx.reinit.W.prob


  ## Initiation
  tx.init.elig.B <- which(active == 1 & race == "B" & status == 1 &
                          tx.status == 0 & diag.status == 1 &
                          tt.traj %in% c("YP", "YF") & cum.time.on.tx == 0 &
                          stage != "D")
  tx.init.B <- tx.init.elig.B[rbinom(length(tx.init.elig.B), 1,
                                     tx.init.B.prob) == 1]

  tx.init.elig.W <- which(active == 1 & race == "W" & status == 1 &
                          tx.status == 0 & diag.status == 1 &
                          tt.traj %in% c("YP", "YF") & cum.time.on.tx == 0 &
                          stage != "D")
  tx.init.W <- tx.init.elig.W[rbinom(length(tx.init.elig.W), 1,
                                     tx.init.W.prob) == 1]

  tx.init <- c(tx.init.B, tx.init.W)

  dat$attr$tx.status[tx.init] <- 1
  dat$attr$tx.init.time[tx.init] <- at


  ## Halting
  tx.halt.elig.B <- which(active == 1 & race == "B" & tx.status == 1)
  tx.halt.B <- tx.halt.elig.B[rbinom(length(tx.halt.elig.B), 1,
                                     tx.halt.B.prob) == 1]

  tx.halt.elig.W <- which(active == 1 & race == "W" & tx.status == 1)
  tx.halt.W <- tx.halt.elig.W[rbinom(length(tx.halt.elig.W),
                                     1, tx.halt.W.prob) == 1]
  tx.halt <- c(tx.halt.B, tx.halt.W)
  dat$attr$tx.status[tx.halt] <- 0


  ## Restarting
  tx.reinit.elig.B <- which(active == 1 & race == "B" & tx.status == 0 &
                            cum.time.on.tx > 0 & stage != "D")
  tx.reinit.B <- tx.reinit.elig.B[rbinom(length(tx.reinit.elig.B),
                                         1, tx.reinit.B.prob) == 1]

  tx.reinit.elig.W <- which(active == 1 & race == "W" & tx.status == 0 &
                            cum.time.on.tx > 0 & stage != "D")
  tx.reinit.W <- tx.reinit.elig.W[rbinom(length(tx.reinit.elig.W),
                                         1, tx.reinit.W.prob) == 1]

  tx.reinit <- c(tx.reinit.B, tx.reinit.W)
  dat$attr$tx.status[tx.reinit] <- 1


  ## Other output
  idsAct <- which(active == 1)
  dat$attr$cum.time.on.tx[idsAct] <- dat$attr$cum.time.on.tx[idsAct] +
                                     ((dat$attr$tx.status[idsAct] == 1) %in% TRUE)
  dat$attr$cum.time.off.tx[idsAct] <- dat$attr$cum.time.off.tx[idsAct] +
                                     ((dat$attr$tx.status[idsAct] == 0) %in% TRUE)

  ## Summary statistics
  dat$epi$tx.init.inc[at] <- length(tx.init)
  dat$epi$tx.halt.inc[at] <- length(tx.halt)
  dat$epi$tx.resm.inc[at] <- length(tx.reinit)

  return(dat)
}
