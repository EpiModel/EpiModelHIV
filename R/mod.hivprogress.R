
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
#' the partially suppressed category (the \code{tt.traj} attribute is \code{1}),
#' time to AIDS depends on the sum of two ratios: time on treatment over maximum
#' time on treatment plus time off treatment over maximum time off treatment.
#' For persons ever on ART who fall into the fully suppressed cateogry
#' (\code{tt.traj=2}), time to AIDS depends on whether the cumulative time
#' off treatment exceeds a time threshold specified in the \code{max.time.off.tx.full.int}
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
hivprogress_msm <- function(dat, at) {

  ## Input
  # Attributes
  active           <- get_attr(dat, "active")
  status           <- get_attr(dat, "status")
  inf.time         <- get_attr(dat, "inf.time")
  cuml.time.on.tx  <- get_attr(dat, "cuml.time.on.tx")
  cuml.time.off.tx <- get_attr(dat, "cuml.time.off.tx")
  stage            <- get_attr(dat, "stage")
  stage.time       <- get_attr(dat, "stage.time")
  aids.time        <- get_attr(dat, "aids.time")
  tt.traj          <- get_attr(dat, "tt.traj")
  tx.status        <- get_attr(dat, "tx.status")

  time.since.inf <- at - inf.time

  # Parameters
  vl.acute.rise.int           <- get_param(dat, "vl.acute.rise.int")
  vl.acute.fall.int           <- get_param(dat, "vl.acute.fall.int")
  vl.aids.onset.int           <- get_param(dat, "vl.aids.onset.int")
  max.time.off.tx.partial.int <- get_param(dat, "max.time.off.tx.partial.int")
  max.time.on.tx.partial.int  <- get_param(dat, "max.time.on.tx.partial.int")
  max.time.off.tx.full.int    <- get_param(dat, "max.time.off.tx.full.int")

  ## Process

  # Increment day
  stage.time[active] <- stage.time[active] + 1

  # Change stage to Acute Falling
  toAF <- which(
    active &
    stage == 1 &
    time.since.inf >= (vl.acute.rise.int + 1)
  )
  stage[toAF] <- 2
  stage.time[toAF] <- 1

  # Change stage to Chronic
  toC <- which(
    active &
    stage == 2 &
    time.since.inf >= (vl.acute.rise.int + vl.acute.fall.int + 1)
  )
  stage[toC] <- 3
  stage.time[toC] <- 1

  # Change stage to AIDS
  aids.tx.naive <- which(
    active &
    status == 1 &
    cuml.time.on.tx == 0 &
    time.since.inf >= vl.aids.onset.int &
    stage != 4
  )

  part.tx.score <- (cuml.time.off.tx / max.time.off.tx.partial.int) +
                   (cuml.time.on.tx / max.time.on.tx.partial.int)

  aids.part.escape <- which(
    active &
    cuml.time.on.tx > 0 &
    tt.traj == 1 &
    stage == 3 &
    part.tx.score >= 1 &
    stage != 4
  )

  aids.off.tx.full.escape <- which(
    active &
    tx.status &
    tt.traj %in% 2:3 &
    cuml.time.on.tx > 0 &
    cuml.time.off.tx >= max.time.off.tx.full.int &
    stage != 4
  )

  isAIDS <- c(aids.tx.naive, aids.part.escape, aids.off.tx.full.escape)
  stage[isAIDS] <- 4
  stage.time[isAIDS] <- 1
  aids.time[isAIDS] <- at

  ## Output
  # Set Attributes
  dat <- set_attr(dat, "stage", stage)
  dat <- set_attr(dat, "stage.time", stage.time)
  dat <- set_attr(dat, "aids.time", aids.time)

  # Set Epi Trackers
  dat <- set_epi(dat, "new.aids.tot", at, length(isAIDS))
  dat <- set_epi(dat, "new.aids.part", at, length(aids.part.escape))
  dat <- set_epi(dat, "new.aids.full", at, length(aids.off.tx.full.escape))

  return(dat)
}
