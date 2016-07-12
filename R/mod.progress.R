
#' @title Disease Progression Module
#'
#' @description Module function for HIV disease progression through acute, chronic
#'              and AIDS stages.
#'
#' @inheritParams aging_msm
#'
#' @details
#' HIV disease is divided into four stages: acute rising, acute falling, chronic
#' and AIDS. Acute rising is the time from infection to peak viremia, while
#' acute falling is the time from peak viremia to chronic stage infection with
#' an established set-point HIV viral load.
#'
#' The time spent in chronic stage infection, and thus the time from infection to
#' AIDS, depends on ART history. For ART-naive persons, time to AIDS is established
#' by the \code{vl.aids.onset} parameter. For persons ever on ART who fall into
#' the partially suppressed category (the \code{tt.traj} attribute is \code{3}),
#' time to AIDS depends on the sum of two ratios: time on treatment over maximum
#' time on treatment plus time off treatment over maximum time off treatment.
#' For persons ever on ART who fall into the fully suppressed cateogry
#' (\code{tt.traj=4}), time to AIDS depends on whether the cumulative time
#' off treatment exceeds a time threshold specified in the \code{max.time.off.tx.full}
#' parameter.
#'
#' @return
#' This function returns the \code{dat} object after updating the disease stage
#' of infected individuals.
#'
#' @keywords module msm
#' 
#' @export
#'
progress_msm <- function(dat, at) {

  ## Variables

  # Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  time.since.inf <- at - dat$attr$inf.time
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  cum.time.off.tx <- dat$attr$cum.time.off.tx
  stage <- dat$attr$stage
  stage.time <- dat$attr$stage.time
  tt.traj <- dat$attr$tt.traj
  tx.status <- dat$attr$tx.status

  # Parameters
  vl.acute.rise.int <- dat$param$vl.acute.rise.int
  vl.acute.fall.int <- dat$param$vl.acute.fall.int
  vl.aids.onset <- dat$param$vl.aids.onset
  max.time.off.tx.part <- dat$param$max.time.off.tx.part
  max.time.on.tx.part <- dat$param$max.time.on.tx.part

  max.time.off.tx.full <- dat$param$max.time.off.tx.full


  ## Process

  # Increment day
  stage.time[active == 1] <- stage.time[active == 1] + 1

  # Change stage to Acute Falling
  toAF <- which(active == 1 & time.since.inf == (vl.acute.rise.int + 1))
  stage[toAF] <- 2
  stage.time[toAF] <- 1

  # Change stage to Chronic
  toC <- which(active == 1 & time.since.inf == (vl.acute.rise.int +
                                                vl.acute.fall.int + 1))
  stage[toC] <- 3
  stage.time[toC] <- 1

  # Change stage to AIDS
  aids.tx.naive <- which(active == 1 & status == 1 & cum.time.on.tx == 0 &
                         (time.since.inf >= vl.aids.onset) & stage != 4)

  part.tx.score <- (cum.time.off.tx / max.time.off.tx.part) +
                   (cum.time.on.tx / max.time.on.tx.part)

  aids.part.escape <- which(active == 1 & cum.time.on.tx > 0 & tt.traj == 3 &
                            stage == 3 & part.tx.score >= 1 & stage != 4)

  aids.off.tx.full.escape <- which(active == 1 & tx.status == 0 & tt.traj == 4 &
                                   cum.time.on.tx > 0 &
                                   cum.time.off.tx >= max.time.off.tx.full &
                                   stage != 4)

  isAIDS <- c(aids.tx.naive, aids.part.escape, aids.off.tx.full.escape)
  stage[isAIDS] <- 4
  stage.time[isAIDS] <- 1


  ## Output
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time

  return(dat)
}
