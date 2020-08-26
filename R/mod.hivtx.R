
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
#' \code{tx.init.time}, \code{cuml.time.on.tx}, \code{cuml.time.off.tx} attributes.
#'
#' @keywords module msm
#'
#' @export
#'
hivtx_msm <- function(dat, at) {

  # Attributes
  race <- dat$attr$race
  status <- dat$attr$status
  tx.status <- dat$attr$tx.status
  diag.status <- dat$attr$diag.status
  cuml.time.on.tx <- dat$attr$cuml.time.on.tx
  cuml.time.off.tx <- dat$attr$cuml.time.off.tx
  tx.period.first <- dat$attr$tx.period.first
  tx.period.last <- dat$attr$tx.period.last
  tx.init.time <- dat$attr$tx.init.time
  tt.traj <- dat$attr$tt.traj

  # Parameters
  tx.init.prob <- dat$param$tx.init.prob
  tx.halt.part.prob <- dat$param$tx.halt.part.prob
  tx.reinit.part.prob <- dat$param$tx.reinit.part.prob
  tx.halt.full.rr <- dat$param$tx.halt.full.rr
  tx.halt.dur.rr <- dat$param$tx.halt.dur.rr
  tx.reinit.full.rr <- dat$param$tx.reinit.full.rr
  tx.reinit.dur.rr <- dat$param$tx.reinit.full.rr

  ## Initiation
  tx.init.elig <- which(status == 1 &
                        tx.status == 0 &
                        diag.status == 1 &
                        cuml.time.on.tx == 0)
  rates <- tx.init.prob[race[tx.init.elig]]
  tx.init <- tx.init.elig[rbinom(length(tx.init.elig), 1, rates) == 1]

  ## Halting
  tx.halt.part.elig <- which(tx.status == 1 & tt.traj == 1)
  rates.part <- tx.halt.part.prob[race[tx.halt.part.elig]]
  tx.halt.part <- tx.halt.part.elig[rbinom(length(tx.halt.part.elig), 1, rates.part) == 1]

  tx.halt.full.elig <- which(tx.status == 1 & tt.traj == 2)
  rates.full <- tx.halt.part.prob[race[tx.halt.full.elig]] * tx.halt.full.rr[race[tx.halt.full.elig]]
  tx.halt.full <- tx.halt.full.elig[rbinom(length(tx.halt.full.elig), 1, rates.full) == 1]

  tx.halt.dur.elig <- which(tx.status == 1 & tt.traj == 3)
  rates.dur <- tx.halt.part.prob[race[tx.halt.dur.elig]] * tx.halt.dur.rr[race[tx.halt.dur.elig]]
  tx.halt.dur <- tx.halt.dur.elig[rbinom(length(tx.halt.dur.elig), 1, rates.dur) == 1]

  tx.halt <- c(tx.halt.part, tx.halt.full, tx.halt.dur)

  ## Restarting
  tx.reinit.part.elig <- which(tx.status == 0 & tt.traj == 1 &
                               cuml.time.on.tx > 0)
  rates.part <- tx.reinit.part.prob[race[tx.reinit.part.elig]]
  tx.reinit.part <- tx.reinit.part.elig[rbinom(length(tx.reinit.part.elig), 1, rates.part) == 1]

  tx.reinit.full.elig <- which(tx.status == 0 & tt.traj == 2 &
                               cuml.time.on.tx > 0)
  rates.full <- tx.reinit.part.prob[race[tx.reinit.full.elig]] * tx.reinit.full.rr[race[tx.reinit.full.elig]]
  tx.reinit.full <- tx.reinit.full.elig[rbinom(length(tx.reinit.full.elig), 1, rates.full) == 1]

  tx.reinit.dur.elig <- which(tx.status == 0 & tt.traj == 3 &
                              cuml.time.on.tx > 0)
  rates.dur <- tx.reinit.part.prob[race[tx.reinit.dur.elig]] * tx.reinit.dur.rr[race[tx.reinit.dur.elig]]
  tx.reinit.dur <- tx.reinit.dur.elig[rbinom(length(tx.reinit.dur.elig), 1, rates.dur) == 1]

  tx.reinit <- c(tx.reinit.part, tx.reinit.full, tx.reinit.dur)

  ## Update Attributes
  tx.status[tx.init] <- 1
  tx.status[tx.halt] <- 0
  tx.status[tx.reinit] <- 1

  cuml.time.on.tx[which(tx.status == 1)] <- cuml.time.on.tx[which(tx.status == 1)] + 1
  cuml.time.off.tx[which(tx.status == 0)] <- cuml.time.off.tx[which(tx.status == 0)] + 1

  tx.init.time[tx.init] <- at

  idsSetPeriod <- union(tx.init, tx.reinit)
  tx.period.first[idsSetPeriod] <- at
  tx.period.last[idsSetPeriod] <- at

  idsContPeriod <- setdiff(which(tx.status == 1), idsSetPeriod)
  tx.period.last[idsContPeriod] <- at

  dat$attr$tt.traj <- tt.traj
  dat$attr$tx.status <- tx.status
  dat$attr$cuml.time.on.tx <- cuml.time.on.tx
  dat$attr$cuml.time.off.tx <- cuml.time.off.tx
  dat$attr$tx.period.first <- tx.period.first
  dat$attr$tx.period.last <- tx.period.last
  dat$attr$tx.init.time <- tx.init.time

  dat$epi$mean.tx.on[at] <- mean(cuml.time.on.tx, na.rm = TRUE)
  dat$epi$mean.tx.off[at] <- mean(cuml.time.off.tx, na.rm = TRUE)

  dat$epi$mean.tx.on.part[at] <- mean(cuml.time.on.tx[tt.traj == 1], na.rm = TRUE)
  dat$epi$mean.tx.off.part[at] <- mean(cuml.time.off.tx[tt.traj == 1], na.rm = TRUE)

  return(dat)
}


#' @export
#' @rdname hivtx_msm
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
