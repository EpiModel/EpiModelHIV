
#' @title Disease Progression Module
#'
#' @description Module function for HIV disease progression through acute,
#'              chronic and AIDS stages.
#'
#' @inheritParams aging_msm
#'
#' @details
#' HIV disease is divided into four stages: acute rising, acute falling, chronic
#' and AIDS. Acute rising is the time from infection to peak viremia, while
#' acute falling is the time from peak viremia to chronic stage infection with
#' an established set-point HIV viral load.
#'
#' The time spent in chronic stage infection, and thus the time from infection
#' to AIDS, depends on ART history. For ART-naive persons, time to AIDS is
#' established by the \code{vl.aids.onset.int} parameter. For persons ever on ART
#' who fall into the partially suppressed category (the \code{tt.traj} attribute
#' is \code{3}), time to AIDS depends on the sum of two ratios: time on
#' treatment over maximum time on treatment plus time off treatment over maximum
#' time off treatment. For persons ever on ART who fall into the fully
#' suppressed category (\code{tt.traj=4}), time to AIDS depends on whether the
#' cumulative time off treatment exceeds a time threshold specified in the
#' \code{max.time.off.tx.full} parameter.
#'
#' @return
#' This function returns the \code{dat} object after updating the disease stage
#' of infected individuals.
#'
#' @keywords module msm
#'
#' @export
#'
hiv_progress_msm <- function(dat, at) {

  ## Variables

  # Attributes
  status <- dat$attr$status
  time.since.inf <- at - dat$attr$inf.time
  time.since.diag <- at - dat$attr$diag.time
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  cum.time.off.tx <- dat$attr$cum.time.off.tx
  stage <- dat$attr$stage
  stage.time <- dat$attr$stage.time
  stage.time.ar.ndx <- dat$attr$stage.time.ar.ndx
  stage.time.ar.dx <- dat$attr$stage.time.ar.dx
  stage.time.af.ndx <- dat$attr$stage.time.af.ndx
  stage.time.af.dx <- dat$attr$stage.time.af.dx
  stage.time.early.chronic.ndx <- dat$attr$stage.time.early.chronic.ndx
  stage.time.early.chronic.dx.yrone <- dat$attr$stage.time.early.chronic.dx.yrone
  stage.time.early.chronic.dx.yrstwotolate <- dat$attr$stage.time.early.chronic.dx.yrstwotolate
  stage.time.early.chronic.art <- dat$attr$stage.time.early.chronic.art
  stage.time.late.chronic.ndx <- dat$attr$stage.time.late.chronic.ndx
  stage.time.late.chronic.dx <- dat$attr$stage.time.late.chronic.dx
  stage.time.late.chronic.art <- dat$attr$stage.time.late.chronic.art
  stage.time.aids.ndx <- dat$attr$stage.time.aids.ndx
  stage.time.aids.dx <- dat$attr$stage.time.aids.dx
  stage.time.aids.art <- dat$attr$stage.time.aids.art
  tt.traj <- dat$attr$tt.traj
  tx.status <- dat$attr$tx.status
  diag.status <- dat$attr$diag.status

  # Parameters
  vl.acute.rise.int <- dat$param$vl.acute.rise.int
  vl.acute.fall.int <- dat$param$vl.acute.fall.int
  vl.aids.onset.int <- dat$param$vl.aids.onset.int
  max.time.off.tx.part <- dat$param$max.time.off.tx.part
  max.time.on.tx.part <- dat$param$max.time.on.tx.part
  max.time.off.tx.full <- dat$param$max.time.off.tx.full

  # Seventy percent value for early/late chronic split
  early.chronic.int <- floor(0.7 * dat$param$max.time.off.tx.full.int)

  ## Process
  # Current stage
  AR.ndx <- which(stage == 1 & diag.status == 0)
  AF.ndx <- which(stage == 2 & diag.status == 0)
  early.chronic.ndx <- which(stage == 3 & diag.status == 0 & cum.time.off.tx <= early.chronic.int)
  late.chronic.ndx <- which(stage == 3 & diag.status == 0 & cum.time.off.tx > early.chronic.int)
  aids.ndx <- which(stage == 4 & diag.status == 0)

  AR.dx <- which(stage == 1 & diag.status == 1)
  AF.dx <- which(stage == 2 & diag.status == 1)
  early.chronic.dx.yrone <- which(stage == 3 & diag.status == 1 &
                                    time.since.diag <= 52 &
                                    cum.time.off.tx <= early.chronic.int &
                                    tx.status == 0)
  early.chronic.dx.yrstwotolate <- which(stage == 3 & diag.status == 1 &
                                           time.since.diag > 52 &
                                           cum.time.off.tx <= early.chronic.int &
                                           tx.status == 0)
  late.chronic.dx <- which(stage == 3 & diag.status == 1 &
                             cum.time.off.tx > early.chronic.int &
                             tx.status == 0)
  aids.dx <- which(stage == 4 & diag.status == 1 & tx.status == 0)

  early.chronic.art <- which(stage == 3 & diag.status == 1 &
                               cum.time.off.tx <= early.chronic.int &
                               tx.status == 1)
  late.chronic.art <- which(stage == 3 & diag.status == 1 &
                              cum.time.off.tx > early.chronic.int &
                              tx.status == 1)
  aids.art <- which(stage == 4 & diag.status == 1 & tx.status == 1)

  # Population numbers (person-time contributed at each time step)
  dat$epi$stage.time.ar.ndx[at] <- length(AR.ndx)
  dat$epi$stage.time.ar.dx[at] <- length(AR.dx)
  dat$epi$stage.time.af.ndx[at] <- length(AF.ndx)
  dat$epi$stage.time.af.dx[at] <- length(AF.dx)
  dat$epi$stage.time.early.chronic.ndx[at] <- length(early.chronic.ndx)
  dat$epi$stage.time.early.chronic.dx.yrone[at] <- length(early.chronic.dx.yrone)
  dat$epi$stage.time.early.chronic.dx.yrstwotolate[at] <- length(early.chronic.dx.yrstwotolate)
  dat$epi$stage.time.early.chronic.art[at] <- length(early.chronic.art)
  dat$epi$stage.time.late.chronic.ndx[at] <- length(late.chronic.ndx)
  dat$epi$stage.time.late.chronic.dx[at] <- length(late.chronic.dx)
  dat$epi$stage.time.late.chronic.art[at] <- length(late.chronic.art)
  dat$epi$stage.time.aids.ndx[at] <- length(aids.ndx)
  dat$epi$stage.time.aids.dx[at] <- length(aids.dx)
  dat$epi$stage.time.aids.art[at] <- length(aids.art)
  dat$epi$time.hivneg[at] <- length(which(status == 0))

  # Increment time step
  stage.time[status == 1] <- stage.time[status == 1] + 1
  stage.time.ar.ndx[AR.ndx] <- stage.time.ar.ndx[AR.ndx] + 1
  stage.time.ar.dx[AR.dx] <- stage.time.ar.dx[AR.dx] + 1
  stage.time.af.ndx[AF.ndx] <- stage.time.af.ndx[AF.ndx] + 1
  stage.time.af.dx[AF.dx] <- stage.time.af.dx[AF.dx] + 1
  stage.time.early.chronic.ndx[early.chronic.ndx] <- stage.time.early.chronic.ndx[early.chronic.ndx] + 1
  stage.time.early.chronic.dx.yrone[early.chronic.dx.yrone] <- stage.time.early.chronic.dx.yrone[early.chronic.dx.yrone] + 1
  stage.time.early.chronic.dx.yrstwotolate[early.chronic.dx.yrstwotolate] <-
    stage.time.early.chronic.dx.yrstwotolate[early.chronic.dx.yrstwotolate] + 1
  stage.time.early.chronic.art[early.chronic.art] <- stage.time.early.chronic.art[early.chronic.art] + 1
  stage.time.late.chronic.ndx[late.chronic.ndx] <- stage.time.late.chronic.ndx[late.chronic.ndx] + 1
  stage.time.late.chronic.dx[late.chronic.dx] <- stage.time.late.chronic.dx[late.chronic.dx] + 1
  stage.time.late.chronic.art[late.chronic.art] <- stage.time.late.chronic.art[late.chronic.art] + 1
  stage.time.aids.ndx[aids.ndx] <- stage.time.aids.ndx[aids.ndx] + 1
  stage.time.aids.dx[aids.dx] <- stage.time.aids.dx[aids.dx] + 1
  stage.time.aids.art[aids.art] <- stage.time.aids.art[aids.art] + 1

  # Change stage to Acute Falling
  toAF.ndx <- which(time.since.inf == (vl.acute.rise.int + 1) & diag.status == 0)
  toAF.dx <- which(time.since.inf == (vl.acute.rise.int + 1) & diag.status == 1)
  toAF <- which(time.since.inf == (vl.acute.rise.int + 1))
  stage[toAF.ndx] <- 2
  stage[toAF.dx] <- 2
  stage.time[toAF] <- 0

  # Change stage to Chronic
  toC <- which(time.since.inf == (vl.acute.rise.int + vl.acute.fall.int + 1))
  toC.ndx <- which(time.since.inf == (vl.acute.rise.int + vl.acute.fall.int + 1) & diag.status == 0)
  toC.dx <- which(time.since.inf == (vl.acute.rise.int + vl.acute.fall.int + 1) & diag.status == 1)

  stage[toC.ndx] <- 3
  stage[toC.dx] <- 3
  stage.time[toC] <- 0

  # Change stage to AIDS
  aids.tx.naive.ndx <- which(status == 1 & cum.time.on.tx == 0 &
                         (time.since.inf >= vl.aids.onset.int) & stage != 4 &
                           diag.status == 0)
  aids.tx.naive.dx <- which(status == 1 & cum.time.on.tx == 0 &
                                 (time.since.inf >= vl.aids.onset.int) & stage != 4
                                  & diag.status == 1 & tx.status == 0)
  aids.tx.naive.art <- which(status == 1 & cum.time.on.tx == 0 &
                                 (time.since.inf >= vl.aids.onset.int) & stage != 4
                                  & diag.status == 1 & tx.status == 1)

  part.tx.score <- (cum.time.off.tx / max.time.off.tx.part) +
                   (cum.time.on.tx / max.time.on.tx.part)

  aids.part.escape.ndx <- which(cum.time.on.tx > 0 & tt.traj == 3 &
                                stage == 3 & part.tx.score >= 1 & stage != 4 &
                                diag.status == 0)
  aids.part.escape.dx <- which(cum.time.on.tx > 0 & tt.traj == 3 &
                                stage == 3 & part.tx.score >= 1 & stage != 4 &
                                diag.status == 1 & tx.status == 0)
  aids.part.escape.art <- which(cum.time.on.tx > 0 & tt.traj == 3 &
                                stage == 3 & part.tx.score >= 1 & stage != 4 &
                                diag.status == 1 & tx.status == 1)

  aids.off.tx.full.escape.ndx <- which(tx.status == 0 & tt.traj == 4 &
                                        cum.time.on.tx > 0 & cum.time.off.tx >= max.time.off.tx.full &
                                        stage != 4 & diag.status == 0)
  aids.off.tx.full.escape.dx <- which(tx.status == 0 & tt.traj == 4 &
                                        cum.time.on.tx > 0 & cum.time.off.tx >= max.time.off.tx.full &
                                        stage != 4 & diag.status == 1)
  aids.off.tx.full.escape.art <- which(tx.status == 1 & tt.traj == 4 &
                                        cum.time.on.tx > 0 & cum.time.off.tx >= max.time.off.tx.full &
                                        stage != 4 & diag.status == 1)

  isAIDS <- c(aids.tx.naive.ndx, aids.tx.naive.dx, aids.tx.naive.art,
              aids.part.escape.ndx, aids.part.escape.dx, aids.part.escape.art,
              aids.off.tx.full.escape.ndx, aids.off.tx.full.escape.dx)
  isAIDS.ndx <- c(aids.tx.naive.ndx, aids.part.escape.ndx,
                  aids.off.tx.full.escape.ndx)
  isAIDS.dx <- c(aids.tx.naive.dx, aids.part.escape.dx,
                 aids.off.tx.full.escape.dx)
  isAIDS.art <- c(aids.tx.naive.art, aids.part.escape.art, aids.off.tx.full.escape.art)

  stage[isAIDS.ndx] <- 4
  stage[isAIDS.dx] <- 4
  stage[isAIDS.art] <- 4
  stage.time[isAIDS] <- 0

  ## Output
  # Individual attribute: time in stage
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time
  dat$attr$stage.time.ar.ndx <- stage.time.ar.ndx
  dat$attr$stage.time.ar.dx <- stage.time.ar.dx
  dat$attr$stage.time.af.ndx <- stage.time.af.ndx
  dat$attr$stage.time.af.dx <- stage.time.af.dx
  dat$attr$stage.time.early.chronic.dx.yrone <- stage.time.early.chronic.dx.yrone
  dat$attr$stage.time.early.chronic.dx.yrstwotolate <- stage.time.early.chronic.dx.yrstwotolate
  dat$attr$stage.time.early.chronic.art <- stage.time.early.chronic.art
  dat$attr$stage.time.late.chronic.ndx <- stage.time.late.chronic.ndx
  dat$attr$stage.time.late.chronic.dx <- stage.time.late.chronic.dx
  dat$attr$stage.time.late.chronic.art <- stage.time.late.chronic.art
  dat$attr$stage.time.aids.ndx <- stage.time.aids.ndx
  dat$attr$stage.time.aids.dx <- stage.time.aids.dx
  dat$attr$stage.time.aids.art <- stage.time.aids.art

  if (at < dat$param$prep.start) {
      dat$attr$time.off.prep[dat$attr$prepStat == 0] <- dat$attr$time.off.prep[dat$attr$prepStat == 0] + 1
  }

  return(dat)
}

#' @title Disease Progression Module
#'
#' @description Module function for Syphilis disease progression through
#'              multiple stages.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Syphilis disease is divided into multiple stages: incubating, primary,
#' secondary, early latent, late latent, tertiary, and remission.
#'
#' The time spent in chronic stage infection, and thus the time from infection
#' to AIDS, depends on ART history. For ART-naive persons, time to AIDS is
#' established by the \code{vl.aids.onset.int} parameter. For persons ever on ART
#' who fall into the partially suppressed category (the \code{tt.traj} attribute
#' is \code{3}), time to AIDS depends on the sum of two ratios: time on
#' treatment over maximum time on treatment plus time off treatment over maximum
#' time off treatment.
#' For persons ever on ART who fall into the fully suppressed category
#' (\code{tt.traj=4}), time to AIDS depends on whether the cumulative time
#' off treatment exceeds a time threshold specified in the
#' \code{max.time.off.tx.full} parameter.
#'
#' @return
#' This function returns the \code{dat} object after updating the disease stage
#' of infected individuals.
#'
#' @keywords module msm syphilis
#'
#' @export
#'
syph_progress_msm <- function(dat, at) {

  ## Variables

  # Attributes
  syphilis <- dat$attr$syphilis
  stage.syph <- dat$attr$stage.syph
  stage.time.syph <- dat$attr$stage.time.syph
  syph.sympt <- dat$attr$syph.sympt
  syph.incub.tx <- dat$attr$syph.incub.tx
  syph.prim.tx <- dat$attr$syph.prim.tx
  syph.seco.tx <- dat$attr$syph.seco.tx
  syph.earlat.tx <- dat$attr$syph.earlat.tx
  syph.latelat.tx <- dat$attr$syph.latelat.tx

  # Parameters
  incu.syph.int <- dat$param$incu.syph.int
  prim.syph.int <- dat$param$prim.syph.int
  seco.syph.int <- dat$param$seco.syph.int
  earlat.syph.int <- dat$param$earlat.syph.int

  syph.prim.sympt.prob <- dat$param$syph.prim.sympt.prob
  syph.seco.sympt.prob <- dat$param$syph.seco.sympt.prob
  syph.earlat.sympt.prob <- dat$param$syph.earlat.sympt.prob
  syph.latelat.sympt.prob <- dat$param$syph.latelat.sympt.prob
  syph.tert.sympt.prob <- dat$param$syph.tert.sympt.prob
  syph.tert.prog.prob <- dat$param$syph.tert.prog.prob

  ## Process

  # Increment time unit
  stage.time.syph[which(syphilis == 1)] <- stage.time.syph[which(syphilis == 1)] + 1

  # Change stage to Primary and assign symptoms
  toPrim <- which(stage.time.syph == (incu.syph.int + 1) &
                  stage.syph == 1 &
                  syphilis == 1)
  stage.syph[toPrim] <- 2
  stage.time.syph[toPrim] <- 0
  syph.incub.tx[toPrim] <- NA
  syph.sympt[toPrim] <- rbinom(length(toPrim), 1, syph.prim.sympt.prob)

  # Change stage to Secondary and assign symptoms
  toSeco <- which(stage.time.syph == (prim.syph.int + 1) &
                  stage.syph == 2 &
                  syphilis == 1)
  stage.syph[toSeco] <- 3
  stage.time.syph[toSeco] <- 0
  syph.prim.tx[toSeco] <- 0
  syph.sympt[toSeco] <- NA
  syph.sympt[toSeco] <- rbinom(length(toSeco), 1, syph.seco.sympt.prob)

  # Change stage to Early Latent and assign symptoms
  toEarLat <- which(stage.time.syph == (seco.syph.int + 1) &
                    stage.syph == 3 &
                    syphilis == 1)
  stage.syph[toEarLat] <- 4
  stage.time.syph[toEarLat] <- 0
  syph.seco.tx[toEarLat] <- NA
  syph.sympt[toEarLat] <- NA
  syph.sympt[toEarLat] <- rbinom(length(toEarLat), 1, syph.earlat.sympt.prob)

  # Change stage to Late Latent and assign symptoms
  toLateLat <- which(stage.time.syph == (earlat.syph.int + 1) &
                     stage.syph == 4 &
                     syphilis == 1)
  stage.syph[toLateLat] <- 5
  stage.time.syph[toLateLat] <- 0
  syph.earlat.tx[toLateLat] <- NA
  syph.sympt[toLateLat] <- NA
  syph.sympt[toLateLat] <- rbinom(length(toLateLat), 1, syph.latelat.sympt.prob)

  # Change stage to tertiary for fraction of those in late late latent
  toTert <- which(stage.syph == 5 &
                  syphilis == 1)
  toTert <- which(rbinom(length(toTert), 1, syph.tert.prog.prob) == 1)
  stage.syph[toTert] <- 6
  stage.time.syph[toTert] <- 0
  syph.latelat.tx[toTert] <- NA
  syph.sympt[toTert] <- NA
  syph.sympt[toTert] <- rbinom(length(toTert), 1, syph.tert.sympt.prob)

  ## Output
  dat$attr$syph.incub.tx <- syph.incub.tx
  dat$attr$syph.prim.tx <- syph.prim.tx
  dat$attr$syph.seco.tx <- syph.seco.tx
  dat$attr$syph.earlat.tx <- syph.earlat.tx
  dat$attr$syph.latelat.tx <- syph.latelat.tx
  dat$attr$stage.syph <- stage.syph
  dat$attr$stage.time.syph <- stage.time.syph
  dat$attr$syph.sympt <- syph.sympt

  return(dat)
}
