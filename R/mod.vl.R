
#' @title Viral Load Module
#'
#' @description Module function for updating HIV viral load.
#'
#' @inheritParams aging_msm
#'
#' @details
#' HIV viral load varies over time as a function of time since infection and ART
#' history. In the absence of ART, VL rises during the acute rising stage and
#' falls during the acute falling stage, until it reaches a set-point value in
#' chronic stage infection. VL again rises during AIDS stage disease until the
#' point of death.
#'
#' For persons who have ever initated treatment (\code{tt.traj} is \code{3} or
#' \code{4}), VL changes depending on current ART use in that time step.
#' Current use is associated with a reduction in VL, with the rates of decline
#' and nadirs dependent on partial or full suppression levels. Current
#' non-adherence is associated with an equal level of increase to VL. All persons
#' who have reached AIDS, regardless of how they arrived, have a similar rate of
#' VL increase.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{vl} attribute.
#'
#' @keywords module msm
#'
#' @export
#'
vl_msm <- function(dat, at) {

  ## Variables

  # Attributes
  inf.time.bp <- at - dat$attr$inf.time
  # cum.time.off.tx <- dat$attr$cum.time.off.tx
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  status <- dat$attr$status
  tt.traj <- dat$attr$tt.traj
  stage <- dat$attr$stage
  vl <- dat$attr$vl
  tx.status <- dat$attr$tx.status

  # Parameters
  vlard <- dat$param$vl.acute.rise.int
  vlap <- dat$param$vl.acute.peak
  vlafd <- dat$param$vl.acute.fall.int
  vlsp <- dat$param$vl.set.point
  vldo <- dat$param$vl.aids.onset
  vldd <- dat$param$vl.aids.int
  vlf  <- dat$param$vl.fatal
  vl.full.supp <- dat$param$vl.full.supp
  full.supp.down.slope <- dat$param$full.supp.down.slope
  vl.part.supp <- dat$param$vl.part.supp
  part.supp.down.slope <- dat$param$part.supp.down.slope
  full.supp.up.slope <- dat$param$full.supp.up.slope
  part.supp.up.slope <- dat$param$part.supp.up.slope
  # max.time.off.tx.part <- dat$param$max.time.off.tx.part
  # max.time.on.tx.part <- dat$param$max.time.on.tx.part

  # Calculations
  vlds <- (vlf - vlsp) / vldd
  # part.tx.score <-  (cum.time.off.tx / max.time.off.tx.part) +
  #                   (cum.time.on.tx / max.time.on.tx.part)


  ## Process

  # 1. tx-naive men
  target <- which(status == 1 & cum.time.on.tx == 0)
  inf.time.bp.tn <- inf.time.bp[target]
  new.vl <- (inf.time.bp.tn <= vlard) * (vlap * inf.time.bp.tn / vlard) +
            (inf.time.bp.tn > vlard) * (inf.time.bp.tn <= vlard + vlafd) *
               ((vlsp - vlap) * (inf.time.bp.tn - vlard) / vlafd + vlap) +
            (inf.time.bp.tn > vlard + vlafd) * (inf.time.bp.tn <= vldo) * (vlsp) +
            (inf.time.bp.tn > vldo) * (vlsp + (inf.time.bp.tn - vldo) * vlds)
  vl[target] <- new.vl

  # 2. men on tx, tt.traj=full, not yet escaped
  target <- which(tx.status == 1 & tt.traj == 4 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmax(current.vl - full.supp.down.slope, vl.full.supp)
  vl[target] <- new.vl

  # 3. men on tx, tt.traj=part, not yet escaped
  target <- which(tx.status == 1 & tt.traj == 3 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmax(current.vl - part.supp.down.slope, vl.part.supp)
  vl[target] <- new.vl

  # 4. men off tx, not naive, tt.traj=full, not yet escaped
  target <- which(tx.status == 0 & tt.traj == 4 &
                  cum.time.on.tx > 0 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmin(current.vl + full.supp.up.slope, vlsp)
  vl[target] <- new.vl

  # 5. men off tx, not naive, tt.traj=part, not yet escaped
  target <- which(tx.status == 0 & tt.traj == 3 &
                  cum.time.on.tx > 0 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmin(current.vl + part.supp.up.slope, vlsp)
  vl[target] <- new.vl

  # 6. men on tx, tt.traj=full, escaped
  # Doesn't exist.

  # 7. men on tx, tt.traj=part, escaped
  target <- which(tx.status == 1 &
                  tt.traj == 3 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vlds
  vl[target] <- new.vl

  # 8. men off tx, tt.traj=full, and escaped
  target <- which(tx.status == 0 & tt.traj == 4 &
                  cum.time.on.tx > 0 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vlds
  vl[target] <- new.vl

  # 9. men off tx, tt.traj=part, and escaped
  target <- which(tx.status == 0 & tt.traj == 3 &
                  cum.time.on.tx > 0 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vlds
  vl[target] <- new.vl


  ## Output
  dat$attr$vl <- vl

  return(dat)
}



#' @title Viral Load Module
#'
#' @description Module function for simulating progression of HIV viral load in
#'              natural disease dynamics and in the presence of ART.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
vl_het <- function(dat, at) {

  ## Common variables
  status <- dat$attr$status
  infTime <- dat$attr$infTime


  # Assign base VL ----------------------------------------------------------
  if (is.null(dat$attr$vlLevel)) {
    dat$attr$vlLevel <- rep(NA, length(status))
    dat$attr$vlSlope <- rep(NA, length(status))
  }
  vlLevel <- dat$attr$vlLevel

  idsEligAsn <- which(status == 1 & is.na(vlLevel))
  if (length(idsEligAsn) > 0) {
    vlLevel[idsEligAsn] <- expected_vl(male = dat$attr$male[idsEligAsn],
                                       age = dat$attr$age[idsEligAsn],
                                       ageInf = dat$attr$ageInf[idsEligAsn],
                                       param = dat$param)
  }


  # Update natural VL -------------------------------------------------------
  txStartTime <- dat$attr$txStartTime
  idsEligUpd <- which(status == 1 &
                      infTime < at & is.na(txStartTime))

  if (length(idsEligUpd) > 0) {
    vlLevel[idsEligUpd] <- expected_vl(male = dat$attr$male[idsEligUpd],
                                       age = dat$attr$age[idsEligUpd],
                                       ageInf = dat$attr$ageInf[idsEligUpd],
                                       param = dat$param)
  }

  # VL decline with ART -----------------------------------------------------
  txStat <- dat$attr$txStat
  idsEligTx <- which(status == 1 & infTime < at & txStat == 1)
  if (length(idsEligTx) > 0) {
    tx.vlsupp.time <- dat$param$tx.vlsupp.time
    tx.vlsupp.level <- dat$param$tx.vlsupp.level

    vlSlope <- dat$attr$vlSlope
    needSlope <- intersect(idsEligTx, which(is.na(vlSlope)))

    vl.slope <- vlSlope
    if (length(needSlope) > 0) {
      vl.diff <- pmin(tx.vlsupp.level - vlLevel[needSlope], 0)
      vl.slope[needSlope] <- vl.diff / tx.vlsupp.time
      dat$attr$vlSlope[needSlope] <- vl.slope[needSlope]
    }

    vlLevel[idsEligTx] <- pmax(vlLevel[idsEligTx] + vl.slope[idsEligTx], tx.vlsupp.level)
  }


  # VL rebound post ART -----------------------------------------------------
  idsEligNoTx <- which(status == 1 &
                       txStat == 0 & !is.na(txStartTime))
  if (length(idsEligNoTx) > 0) {
    tx.vlsupp.time <- dat$param$tx.vlsupp.time

    expVl <- expected_vl(male = dat$attr$male[idsEligNoTx],
                         age = dat$attr$age[idsEligNoTx],
                         ageInf = dat$attr$ageInf[idsEligNoTx],
                         param = dat$param)

    vl.slope <- dat$attr$vlSlope

    vlLevel[idsEligNoTx] <- pmin(vlLevel[idsEligNoTx] - vl.slope[idsEligNoTx], expVl)
  }

  dat$attr$vlLevel <- vlLevel

  return(dat)
}


expected_vl <- function(male, age, ageInf, param) {

  timeInf <- (age - ageInf) * (365 / param$time.unit)

  slope1 <- param$vl.acute.peak / param$vl.acute.topeak
  slope2 <- (param$vl.setpoint - param$vl.acute.peak) /
    (param$vl.acute.toset - param$vl.acute.topeak)

  sl3denom <- expected_cd4(method = "timeto",
                           cd4Count1 = 200, cd4Count2 = 25,
                           male = male, age = age, ageInf = ageInf,
                           time.unit = param$time.unit)
  slope3 <- (param$vl.aidsmax - param$vl.setpoint) / sl3denom

  setptTime <- param$vl.acute.topeak + param$vl.acute.toset
  aidsTime <- expected_cd4(method = "timeto", cd4Count1 = 200,
                           male = male, age = age, ageInf = ageInf,
                           time.unit = param$time.unit)

  gp <- 1 * (timeInf <= param$vl.acute.topeak) +
    2 * (timeInf > param$vl.acute.topeak & timeInf <= setptTime) +
    3 * (timeInf > setptTime & timeInf <= aidsTime) +
    4 * (timeInf > aidsTime)

  vlLevel <- rep(NA, length(timeInf))
  vlLevel[gp == 1] <- timeInf[gp == 1] * slope1
  vlLevel[gp == 2] <- pmax(param$vl.setpoint,
                           param$vl.acute.peak +
                             (timeInf[gp == 2] - param$vl.acute.topeak) * slope2)
  vlLevel[gp == 3] <- param$vl.setpoint
  vlLevel[gp == 4] <- pmin(param$vl.aidsmax,
                           param$vl.setpoint +
                             (timeInf[gp == 4] - aidsTime[gp == 4]) * slope3[gp == 4])

  return(vlLevel)
}
