
#' @title Treatment Module
#'
#' @description Module function for anti-retroviral treatment initiation and
#'              adherence over time.
#'
#' @inheritParams aging_msm
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

  # Parameters
  tx.init.B.prob <- dat$param$tx.init.B.prob
  tx.init.W.prob <- dat$param$tx.init.W.prob
  tx.halt.B.prob <- dat$param$tx.halt.B.prob
  tx.halt.W.prob <- dat$param$tx.halt.W.prob
  tx.reinit.B.prob <- dat$param$tx.reinit.B.prob
  tx.reinit.W.prob <- dat$param$tx.reinit.W.prob


  ## Initiation
  tx.init.elig.B <- which(race == "B" & status == 1 &
                          tx.status == 0 & diag.status == 1 &
                          tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                          stage != 4)
  tx.init.B <- tx.init.elig.B[rbinom(length(tx.init.elig.B), 1,
                                     tx.init.B.prob) == 1]

  tx.init.elig.W <- which(race == "W" & status == 1 &
                          tx.status == 0 & diag.status == 1 &
                          tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                          stage != 4)
  tx.init.W <- tx.init.elig.W[rbinom(length(tx.init.elig.W), 1,
                                     tx.init.W.prob) == 1]

  tx.init <- c(tx.init.B, tx.init.W)

  dat$attr$tx.status[tx.init] <- 1
  dat$attr$tx.init.time[tx.init] <- at


  ## Halting
  tx.halt.elig.B <- which(race == "B" & tx.status == 1)
  tx.halt.B <- tx.halt.elig.B[rbinom(length(tx.halt.elig.B), 1,
                                     tx.halt.B.prob) == 1]

  tx.halt.elig.W <- which(race == "W" & tx.status == 1)
  tx.halt.W <- tx.halt.elig.W[rbinom(length(tx.halt.elig.W),
                                     1, tx.halt.W.prob) == 1]
  tx.halt <- c(tx.halt.B, tx.halt.W)
  dat$attr$tx.status[tx.halt] <- 0


  ## Restarting
  tx.reinit.elig.B <- which(race == "B" & tx.status == 0 &
                            cum.time.on.tx > 0 & stage != 4)
  tx.reinit.B <- tx.reinit.elig.B[rbinom(length(tx.reinit.elig.B),
                                         1, tx.reinit.B.prob) == 1]

  tx.reinit.elig.W <- which(race == "W" & tx.status == 0 &
                            cum.time.on.tx > 0 & stage != 4)
  tx.reinit.W <- tx.reinit.elig.W[rbinom(length(tx.reinit.elig.W),
                                         1, tx.reinit.W.prob) == 1]

  tx.reinit <- c(tx.reinit.B, tx.reinit.W)
  dat$attr$tx.status[tx.reinit] <- 1


  ## Other output
  dat$attr$cum.time.on.tx <- dat$attr$cum.time.on.tx +
                             ((dat$attr$tx.status == 1) %in% TRUE)
  dat$attr$cum.time.off.tx <- dat$attr$cum.time.off.tx +
                              ((dat$attr$tx.status == 0) %in% TRUE)

  ## Summary statistics
  dat$epi$tx.init.inc[at] <- length(tx.init)
  dat$epi$tx.halt.inc[at] <- length(tx.halt)
  dat$epi$tx.resm.inc[at] <- length(tx.reinit)

  return(dat)
}


#' @title HIV Anti-Retroviral Treatment Module
#'
#' @description Module function for simulating HIV therapy after diagnosis,
#'              including adherence and non-adherence to ART.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
tx_het <- function(dat, at) {

  # Variables ---------------------------------------------------------------
  dxStat <- dat$attr$dxStat
  txStat <- dat$attr$txStat
  txStartTime <- dat$attr$txStartTime
  txStops <- dat$attr$txStops
  txTimeOn <- dat$attr$txTimeOn
  txTimeOff <- dat$attr$txTimeOff
  txCD4start <- dat$attr$txCD4start

  cd4Count <- dat$attr$cd4Count
  tx.elig.cd4 <- dat$param$tx.elig.cd4
  tx.coverage <- dat$param$tx.coverage

  txType <- dat$attr$txType
  tx.adhere.full <- dat$param$tx.adhere.full
  tx.adhere.part <- dat$param$tx.adhere.part


  # Start tx for tx naive ---------------------------------------------------

  ## Calculate tx coverage
  allElig <- which((cd4Count < tx.elig.cd4 | !is.na(txStartTime)))
  txCov <- sum(!is.na(txStartTime[allElig]))/length(allElig)
  if (is.nan(txCov)) {
    txCov <- 0
  }

  idsElig <- which(dxStat == 1 & txStat == 0 &
                   is.na(txStartTime) & cd4Count < tx.elig.cd4)
  nElig <- length(idsElig)
  idsTx <- NULL


  ## Treatment coverage
  nStart <- max(0, min(nElig, round((tx.coverage - txCov) * length(allElig))))
  if (nStart > 0) {
    idsTx <- ssample(idsElig, nStart)
  }


  ## Treatment type assignment
  if (length(idsTx) > 0) {
    needtxType <- which(is.na(txType[idsTx]))
    if (length(needtxType) > 0) {
      txType[idsTx[needtxType]] <- rbinom(length(needtxType), 1, tx.adhere.full)
    }
    if (tx.adhere.part == 0) {
      idsTx <- intersect(idsTx, which(txType == 1))
    }
  }

  if (length(idsTx) > 0) {
    txStat[idsTx] <- 1
    txStartTime[idsTx] <- at
    txStops[idsTx] <- 0
    txTimeOn[idsTx] <- 0
    txTimeOff[idsTx] <- 0
    txCD4start[idsTx] <- cd4Count[idsTx]
  }


  # Stop tx -----------------------------------------------------------------
  idsStop <- NULL
  idsEligStop <- which(dat$attr$txStat == 1 & txType == 0)
  nEligStop <- length(idsEligStop)
  if (nEligStop > 0) {
    vecStop <- which(rbinom(nEligStop, 1, (1 - tx.adhere.part)) == 1)
    if (length(vecStop) > 0) {
      idsStop <- idsEligStop[vecStop]
      txStat[idsStop] <- 0
      txStops[idsStop] <- txStops[idsStop] + 1
    }
  }


  # Restart tx --------------------------------------------------------------
  idsRest <- NULL
  idsEligRest <- which(dat$attr$txStat == 0 & txStops > 0)
  nEligRest <- length(idsEligRest)
  if (nEligRest > 0) {
    vecRes <- which(rbinom(nEligRest, 1, tx.adhere.part) == 1)
    if (length(vecRes) > 0) {
      idsRest <- idsEligRest[vecRes]
      txStat[idsRest] <- 1
      dat$attr$vlSlope[idsRest] <- NA
    }
  }


  # Output ------------------------------------------------------------------
  idsOnTx <- which(txStat == 1)
  idsOffTx <- which(txStat == 0 & !is.na(txStartTime))
  txTimeOn[idsOnTx] <- txTimeOn[idsOnTx] + 1
  txTimeOff[idsOffTx] <- txTimeOff[idsOffTx] + 1

  dat$attr$txStat <- txStat
  dat$attr$txStartTime <- txStartTime
  dat$attr$txStops <- txStops
  dat$attr$txTimeOn <- txTimeOn
  dat$attr$txTimeOff <- txTimeOff
  dat$attr$txType <- txType
  dat$attr$txCD4start <- txCD4start

  return(dat)
}

