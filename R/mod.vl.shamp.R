
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
#' @keywords module shamp
#'
#' @export
#'
vl_shamp <- function(dat, at) {

  ## Variables

  # Attributes
  inf.time.bp <- at - dat$attr$inf.time
  #cum.time.off.tx <- dat$attr$cum.time.off.tx
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

  # 1. tx-naive persons
  target <- which(status == 1 & cum.time.on.tx == 0)
  inf.time.bp.tn <- inf.time.bp[target]
  new.vl <- (inf.time.bp.tn <= vlard) * (vlap * inf.time.bp.tn / vlard) +
            (inf.time.bp.tn > vlard) * (inf.time.bp.tn <= vlard + vlafd) *
               ((vlsp - vlap) * (inf.time.bp.tn - vlard) / vlafd + vlap) +
            (inf.time.bp.tn > vlard + vlafd) * (inf.time.bp.tn <= vldo) * (vlsp) +
            (inf.time.bp.tn > vldo) * (vlsp + (inf.time.bp.tn - vldo) * vlds)
  vl[target] <- new.vl

  # 2. persons on tx, tt.traj=full, not yet escaped
  target <- which(tx.status == 1 & tt.traj == 4 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmax(current.vl - full.supp.down.slope, vl.full.supp)
  vl[target] <- new.vl

  # 3. persons on tx, tt.traj=part, not yet escaped
  target <- which(tx.status == 1 & tt.traj == 3 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmax(current.vl - part.supp.down.slope, vl.part.supp)
  vl[target] <- new.vl

  # 4. persons off tx, not naive, tt.traj=full, not yet escaped
  target <- which(tx.status == 0 & tt.traj == 4 &
                  cum.time.on.tx > 0 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmin(current.vl + full.supp.up.slope, vlsp)
  vl[target] <- new.vl

  # 5. persons off tx, not naive, tt.traj=part, not yet escaped
  target <- which(tx.status == 0 & tt.traj == 3 &
                  cum.time.on.tx > 0 & stage != 4)
  current.vl <- vl[target]
  new.vl <- pmin(current.vl + part.supp.up.slope, vlsp)
  vl[target] <- new.vl

  # 6. persons on tx, tt.traj=full, escaped
  # Doesn't exist.

  # 7. persons on tx, tt.traj=part, escaped
  target <- which(tx.status == 1 &
                  tt.traj == 3 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vlds
  vl[target] <- new.vl

  # 8. persons off tx, tt.traj=full, and escaped
  target <- which(tx.status == 0 & tt.traj == 4 &
                  cum.time.on.tx > 0 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vlds
  vl[target] <- new.vl

  # 9. persons off tx, tt.traj=part, and escaped
  target <- which(tx.status == 0 & tt.traj == 3 &
                  cum.time.on.tx > 0 & stage == 4)
  current.vl <- vl[target]
  new.vl <- current.vl + vlds
  vl[target] <- new.vl


  ## Output
  dat$attr$vl <- vl

  return(dat)
}

