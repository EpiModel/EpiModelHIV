
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
#' established by the \code{vl.aids.onset} parameter. For persons ever on ART
#' who fall into the partially suppressed category (the \code{tt.traj} attribute
#' is \code{3}), time to AIDS depends on the sum of two ratios: time on
#' treatment over maximum time on treatment plus time off treatment over maximum
#' time off treatment. For persons ever on ART who fall into the fully
#' suppressed cateogry (\code{tt.traj=4}), time to AIDS depends on whether the
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
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  cum.time.off.tx <- dat$attr$cum.time.off.tx
  stage <- dat$attr$stage
  stage.time <- dat$attr$stage.time
  stage.time.ar.ndx <- dat$attr$stage.time.ar.ndx
  stage.time.ar.dx <- dat$attr$stage.time.ar.dx
  stage.time.ar.art <- dat$attr$stage.time.ar.art
  stage.time.af.ndx <- dat$attr$stage.time.af.ndx
  stage.time.af.dx <- dat$attr$stage.time.af.dx
  stage.time.af.art <- dat$attr$stage.time.af.art
  stage.time.chronic.ndx <- dat$attr$stage.time.chronic.ndx
  stage.time.chronic.dx <- dat$attr$stage.time.chronic.dx
  stage.time.chronic.art <- dat$attr$stage.time.chronic.art
  stage.time.aids.ndx <- dat$attr$stage.time.aids.ndx
  stage.time.aids.dx <- dat$attr$stage.time.aids.dx
  stage.time.aids.art <- dat$attr$stage.time.aids.art
  tt.traj <- dat$attr$tt.traj
  tx.status <- dat$attr$tx.status
  diag.status <- dat$attr$diag.status

  # Parameters
  vl.acute.rise.int <- dat$param$vl.acute.rise.int
  vl.acute.fall.int <- dat$param$vl.acute.fall.int
  vl.aids.onset <- dat$param$vl.aids.onset
  max.time.off.tx.part <- dat$param$max.time.off.tx.part
  max.time.on.tx.part <- dat$param$max.time.on.tx.part

  max.time.off.tx.full <- dat$param$max.time.off.tx.full

  ## Process
  # Current stage
  AR.ndx <- which(stage == 1 & diag.status == 0)
  AF.ndx <- which(stage == 2 & diag.status == 0)
  chronic.ndx <- which(stage == 3 & diag.status == 0)
  aids.ndx <- which(stage == 4 & diag.status == 0)

  AR.dx <- which(stage == 1 & diag.status == 1 & tx.status == 0)
  AF.dx <- which(stage == 2 & diag.status == 1 & tx.status == 0)
  chronic.dx <- which(stage == 3 & diag.status == 1 & tx.status == 0)
  aids.dx <- which(stage == 4 & diag.status == 1 & tx.status == 0)

  AR.art <- which(stage == 1 & diag.status == 1 & tx.status == 1)
  AF.art <- which(stage == 2 & diag.status == 1 & tx.status == 1)
  chronic.art <- which(stage == 3 & diag.status == 1 & tx.status == 1)
  aids.art <- which(stage == 4 & diag.status == 1 & tx.status == 1)

  # Population numbers (person-time contributed at each time step)
  dat$epi$stage.time.ar.ndx[at] <- length(AR.ndx)
  dat$epi$stage.time.ar.dx[at] <- length(AR.dx)
  dat$epi$stage.time.ar.art[at] <- length(AR.art)
  dat$epi$stage.time.af.ndx[at] <- length(AF.ndx)
  dat$epi$stage.time.af.dx[at] <- length(AF.dx)
  dat$epi$stage.time.af.art[at] <- length(AF.art)
  dat$epi$stage.time.chronic.ndx[at] <- length(chronic.ndx)
  dat$epi$stage.time.chronic.dx[at] <- length(chronic.dx)
  dat$epi$stage.time.chronic.art[at] <- length(chronic.art)
  dat$epi$stage.time.aids.ndx[at] <- length(aids.ndx)
  dat$epi$stage.time.aids.dx[at] <- length(aids.dx)
  dat$epi$stage.time.aids.art[at] <- length(aids.art)
  dat$epi$time.hivneg[at] <- length(which(status == 0))

  # Increment day
  stage.time[status == 1] <- stage.time[status == 1] + 1
  stage.time.ar.ndx[AR.ndx] <- stage.time.ar.ndx[AR.ndx] + 1
  stage.time.ar.dx[AR.dx] <- stage.time.ar.dx[AR.dx] + 1
  stage.time.ar.art[AR.art] <- stage.time.ar.art[AR.art] + 1
  stage.time.af.ndx[AF.ndx] <- stage.time.af.ndx[AF.ndx] + 1
  stage.time.af.dx[AF.dx] <- stage.time.af.art[AF.dx] + 1
  stage.time.af.art[AF.art] <- stage.time.af.art[AF.art] + 1
  stage.time.chronic.ndx[chronic.ndx] <- stage.time.chronic.ndx[chronic.ndx] + 1
  stage.time.chronic.dx[chronic.dx] <- stage.time.chronic.dx[chronic.dx] + 1
  stage.time.chronic.art[chronic.art] <- stage.time.chronic.art[chronic.art] + 1
  stage.time.aids.ndx[aids.ndx] <- stage.time.aids.ndx[aids.ndx] + 1
  stage.time.aids.dx[aids.dx] <- stage.time.aids.dx[aids.dx] + 1
  stage.time.aids.art[aids.art] <- stage.time.aids.art[aids.art] + 1

  # Change stage to Acute Falling
  toAF.ndx <- which(time.since.inf == (vl.acute.rise.int + 1) & diag.status == 0)
  toAF.dx <- which(time.since.inf == (vl.acute.rise.int + 1) & diag.status == 1 & tx.status == 0)
  toAF.art <- which(time.since.inf == (vl.acute.rise.int + 1) & diag.status == 1 & tx.status == 1)
  toAF <- which(time.since.inf == (vl.acute.rise.int + 1))
  stage[toAF.ndx] <- 2
  stage[toAF.dx] <- 2
  stage[toAF.art] <- 2
  stage.time[toAF] <- 1
  stage.time.af.ndx[toAF.ndx] <- 1
  stage.time.af.dx[toAF.dx] <- 1
  stage.time.af.art[toAF.art] <- 1

  # Change stage to Chronic
  toC <- which(time.since.inf == (vl.acute.rise.int + vl.acute.fall.int + 1))
  toC.ndx <- which(time.since.inf == (vl.acute.rise.int + vl.acute.fall.int + 1) & diag.status == 0)
  toC.dx <- which(time.since.inf == (vl.acute.rise.int + vl.acute.fall.int + 1) & diag.status == 1 & tx.status == 0)
  toC.art <- which(time.since.inf == (vl.acute.rise.int + vl.acute.fall.int + 1) & diag.status == 1 & tx.status == 1)

  stage[toC.ndx] <- 3
  stage[toC.dx] <- 3
  stage[toC.art] <- 3
  stage.time[toC] <- 1
  stage.time.chronic.ndx[toC.ndx] <- 1
  stage.time.chronic.dx[toC.dx] <- 1
  stage.time.chronic.art[toC.art] <- 1

  # Change stage to AIDS
  aids.tx.naive.ndx <- which(status == 1 & cum.time.on.tx == 0 &
                         (time.since.inf >= vl.aids.onset) & stage != 4 & diag.status == 0)
  aids.tx.naive.dx <- which(status == 1 & cum.time.on.tx == 0 &
                                 (time.since.inf >= vl.aids.onset) & stage != 4
                                  & diag.status == 1 & tx.status == 0)
  aids.tx.naive.art <- which(status == 1 & cum.time.on.tx == 0 &
                                 (time.since.inf >= vl.aids.onset) & stage != 4
                                  & diag.status == 1 & tx.status == 1)

  part.tx.score <- (cum.time.off.tx / max.time.off.tx.part) +
                   (cum.time.on.tx / max.time.on.tx.part)

  aids.part.escape.ndx <- which(cum.time.on.tx > 0 & tt.traj == 3 &
                                stage == 3 & part.tx.score >= 1 & stage != 4)
  aids.part.escape.dx <- which(cum.time.on.tx > 0 & tt.traj == 3 &
                                stage == 3 & part.tx.score >= 1 & stage != 4 &
                                diag.status == 1 & tx.status == 0)
  aids.part.escape.art <- which(cum.time.on.tx > 0 & tt.traj == 3 &
                                stage == 3 & part.tx.score >= 1 & stage != 4 &
                                diag.status == 1 & tx.status == 1)

  aids.off.tx.full.escape.ndx <- which(tx.status == 0 & tt.traj == 4 &
                                        cum.time.on.tx > 0 & cum.time.off.tx >= max.time.off.tx.full &
                                        stage != 4)
  aids.off.tx.full.escape.dx <- which(tx.status == 0 & tt.traj == 4 &
                                        cum.time.on.tx > 0 & cum.time.off.tx >= max.time.off.tx.full &
                                        stage != 4 & diag.status == 1 & tx.status == 0)
  aids.off.tx.full.escape.art <- which(tx.status == 0 & tt.traj == 4 & cum.time.on.tx > 0 &
                                        cum.time.off.tx >= max.time.off.tx.full &
                                        stage != 4 & diag.status == 1 & tx.status == 1)

  isAIDS <- c(aids.tx.naive.ndx, aids.tx.naive.dx, aids.tx.naive.art,
              aids.part.escape.ndx, aids.part.escape.dx, aids.part.escape.art,
              aids.off.tx.full.escape.ndx, aids.off.tx.full.escape.dx,
              aids.off.tx.full.escape.art)
  isAIDS.ndx <- c(aids.tx.naive.ndx, aids.part.escape.ndx,
                  aids.off.tx.full.escape.ndx)
  isAIDS.dx <- c(aids.tx.naive.dx, aids.part.escape.dx,
                 aids.off.tx.full.escape.dx)
  isAIDS.art <- c(aids.tx.naive.art, aids.part.escape.art,
                  aids.off.tx.full.escape.art)

  stage[isAIDS.ndx] <- 4
  stage[isAIDS.dx] <- 4
  stage[isAIDS.art] <- 4
  stage.time[isAIDS] <- 1
  stage.time.aids.ndx[isAIDS.ndx] <- 1
  stage.time.aids.dx[isAIDS.dx] <- 1
  stage.time.aids.art[isAIDS.art] <- 1

  ## Output
  # Individual attribute: time in stage
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time
  dat$attr$stage.time.ar.ndx <- stage.time.ar.ndx
  dat$attr$stage.time.ar.dx <- stage.time.ar.dx
  dat$attr$stage.time.ar.art <- stage.time.ar.art
  dat$attr$stage.time.af.ndx <- stage.time.af.ndx
  dat$attr$stage.time.af.dx <- stage.time.af.dx
  dat$attr$stage.time.af.art <- stage.time.af.art
  dat$attr$stage.time.chronic.ndx <- stage.time.chronic.ndx
  dat$attr$stage.time.chronic.dx <- stage.time.chronic.dx
  dat$attr$stage.time.chronic.art <- stage.time.chronic.art
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
#' established by the \code{vl.aids.onset} parameter. For persons ever on ART
#' who fall into the partially suppressed category (the \code{tt.traj} attribute
#' is \code{3}), time to AIDS depends on the sum of two ratios: time on
#' treatment over maximum time on treatment plus time off treatment over maximum
#' time off treatment.
#' For persons ever on ART who fall into the fully suppressed cateogry
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

  # Parameters
  incu.syph.int <- dat$param$incu.syph.int
  prim.syph.int <- dat$param$prim.syph.int
  seco.syph.int <- dat$param$seco.syph.int
  earlat.syph.int <- dat$param$earlat.syph.int
  latelat.syph.int <- dat$param$latelat.syph.int

  syph.prim.sympt.prob <- dat$param$syph.prim.sympt.prob
  syph.seco.sympt.prob <- dat$param$syph.seco.sympt.prob
  syph.earlat.sympt.prob <- dat$param$syph.earlat.sympt.prob
  syph.latelat.sympt.prob <- dat$param$syph.latelat.sympt.prob
  syph.tert.sympt.prob <- dat$param$syph.tert.sympt.prob
  syph.tert.prog.prob <- dat$param$syph.tert.prog.prob

  stage.prim.sympt <- dat$attr$stage.prim.sympt
  stage.seco.sympt <- dat$attr$stage.seco.sympt
  stage.earlat.sympt <- dat$attr$stage.earlat.sympt
  stage.latelat.sympt <- dat$attr$stage.latelat.sympt
  stage.latelatelat.sympt <- dat$attr$stage.latelatelat.sympt
  stage.tert.sympt <- dat$attr$stage.tert.sympt

  ## Process

  # Increment time unit
  stage.time.syph[syphilis == 1] <- stage.time.syph[syphilis == 1] + 1

  # Change stage to Primary and assign symptoms
  toPrim <- which(stage.time.syph == (incu.syph.int + 1) &
                  stage.syph == 1 &
                  syphilis == 1)
  stage.syph[toPrim] <- 2
  stage.time.syph[toPrim] <- 0
  stage.prim.sympt[toPrim] <- rbinom(length(toPrim), 1, syph.prim.sympt.prob)

  # Change stage to Secondary and assign symptoms
  toSeco <- which(stage.time.syph == (prim.syph.int + 1) &
                  stage.syph == 2 &
                  syphilis == 1)
  stage.syph[toSeco] <- 3
  stage.time.syph[toSeco] <- 0
  stage.seco.sympt[toSeco] <- rbinom(length(toSeco), 1, syph.seco.sympt.prob)
  stage.prim.sympt[toSeco] <- NA

  # Change stage to Early Latent and assign symptoms
  toEarLat <- which(stage.time.syph == (seco.syph.int + 1) &
                    stage.syph == 3 &
                    syphilis == 1)
  stage.syph[toEarLat] <- 4
  stage.time.syph[toEarLat] <- 0
  stage.earlat.sympt[toEarLat] <- rbinom(length(toEarLat), 1, syph.earlat.sympt.prob)
  stage.seco.sympt[toEarLat] <- NA

  # Change stage to Late Latent and assign symptoms
  toLateLat <- which(stage.time.syph == (earlat.syph.int + 1) &
                     stage.syph == 4 &
                     syphilis == 1)
  stage.syph[toLateLat] <- 5
  stage.time.syph[toLateLat] <- 0
  stage.latelat.sympt[toLateLat] <- rbinom(length(toLateLat), 1, syph.latelat.sympt.prob)
  stage.earlat.sympt[toLateLat] <- NA

  # Change stage to late late latent (functions the same way as late latent)
  tolatelate <- which(stage.time.syph == (latelat.syph.int + 1) &
                      stage.syph == 5 &
                      syphilis == 1)
  stage.syph[tolatelate] <- 6
  stage.time.syph[tolatelate] <- 0
  stage.latelatelat.sympt[tolatelate] <- rbinom(length(tolatelate), 1, syph.latelat.sympt.prob)
  stage.latelat.sympt[tolatelate] <- NA

  # Change stage to tertiary for fraction of those in late late latent
  toTert <- which(stage.time.syph >= 1 &
                  stage.syph == 6 &
                  syphilis == 1)
  toTert <- which(rbinom(length(toTert), 1, syph.tert.prog.prob) == 1)
  stage.syph[toTert] <- 7
  stage.time.syph[toTert] <- 0
  stage.tert.sympt[toTert] <- rbinom(length(toTert), 1, syph.tert.sympt.prob)
  stage.latelatelat.sympt[toTert] <- NA


  ## Output
  dat$attr$stage.syph <- stage.syph
  dat$attr$stage.time.syph <- stage.time.syph
  dat$attr$stage.prim.sympt <- stage.prim.sympt
  dat$attr$stage.seco.sympt <- stage.seco.sympt
  dat$attr$stage.earlat.sympt <- stage.earlat.sympt
  dat$attr$stage.latelat.sympt <- stage.latelat.sympt
  dat$attr$stage.latelatelat.sympt <- stage.latelatelat.sympt
  dat$attr$stage.tert.sympt <- stage.tert.sympt

    return(dat)
}
