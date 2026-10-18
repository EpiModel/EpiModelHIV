
#' @title Treatment Module
#'
#' @description Module function for antiretroviral treatment initiation and
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
#' \code{tx.init.time}, \code{cuml.time.on.tx}, \code{cuml.time.off.tx}
#' attributes.
#'
#' @keywords module msm
#'
#' @export
#'
hivtx_msm <- function(dat, at) {

  ## Input
  # Attributes
  race             <- get_attr(dat, "race")
  status           <- get_attr(dat, "status")
  tx.status        <- get_attr(dat, "tx.status")
  diag.status      <- get_attr(dat, "diag.status")
  cuml.time.on.tx  <- get_attr(dat, "cuml.time.on.tx")
  cuml.time.off.tx <- get_attr(dat, "cuml.time.off.tx")
  tx.init.time     <- get_attr(dat, "tx.init.time")
  tt.traj          <- get_attr(dat, "tt.traj")

  part.ident       <- get_attr(dat, "part.ident")

  # Parameters
  tx.init.prob           <- get_param(dat, "tx.init.prob")
  tx.halt.partial.prob   <- get_param(dat, "tx.halt.partial.prob")
  tx.halt.full.rr        <- get_param(dat, "tx.halt.full.rr")
  tx.halt.durable.rr     <- get_param(dat, "tx.halt.durable.rr")
  tx.reinit.partial.prob <- get_param(dat, "tx.reinit.partial.prob")
  tx.reinit.full.rr      <- get_param(dat, "tx.reinit.full.rr")
  tx.reinit.durable.rr   <- get_param(dat, "tx.reinit.full.rr")

  part.tx.init.prob      <- get_param(dat, "part.tx.init.prob")
  part.tx.halt.prob      <- get_param(dat, "part.tx.halt.prob")
  part.tx.reinit.prob    <- get_param(dat, "part.tx.reinit.prob")

  ## Process
  # Initiation (partner services)
  part.tx.init.elig <- which(
    status == 1 &
    tx.status == 0 &
    diag.status == 1 &
    cuml.time.on.tx == 0 &
    part.ident == at
  )

  rates.init.part <- part.tx.init.prob[race[part.tx.init.elig]]
  part.tx.init <- part.tx.init.elig[
    runif(length(part.tx.init.elig)) < rates.init.part
  ]

  ## Initiation (general)
  tx.init.elig <- which(
    status == 1 &
    tx.status == 0 &
    diag.status == 1 &
    cuml.time.on.tx == 0
  )

  tx.init.elig <- setdiff(tx.init.elig, part.tx.init.elig)
  rates.init.gen <- tx.init.prob[race[tx.init.elig]]
  tx.init <- tx.init.elig[
    runif(length(tx.init.elig)) < rates.init.gen
  ]

  tx.init.all <- c(tx.init, part.tx.init)

  ## Halting (partner services)

  tx.halt.part.elig <- numeric()
  tx.halt.part <- numeric()

  ## Halting (general)
  tx.halt.partial.elig <- which(tx.status == 1 & tt.traj == 1)
  tx.halt.partial.elig <- setdiff(tx.halt.partial.elig, tx.halt.part.elig)
  rates.halt.partial <- tx.halt.partial.prob[race[tx.halt.partial.elig]]
  tx.halt.partial <- tx.halt.partial.elig[
    runif(length(tx.halt.partial.elig)) < rates.halt.partial
  ]

  tx.halt.full.elig <- which(tx.status == 1 & tt.traj == 2)
  tx.halt.full.elig <- setdiff(tx.halt.full.elig, tx.halt.part.elig)
  rates.halt.full <- tx.halt.partial.prob[race[tx.halt.full.elig]] *
                     tx.halt.full.rr[race[tx.halt.full.elig]]
  tx.halt.full <- tx.halt.full.elig[
    runif(length(tx.halt.full.elig)) < rates.halt.full
  ]

  tx.halt.durable.elig <- which(tx.status == 1 & tt.traj == 3)
  tx.halt.durable.elig <- setdiff(tx.halt.durable.elig, tx.halt.part.elig)
  rates.halt.durable <- tx.halt.partial.prob[race[tx.halt.durable.elig]] *
                        tx.halt.durable.rr[race[tx.halt.durable.elig]]
  tx.halt.durable <- tx.halt.durable.elig[
    runif(length(tx.halt.durable.elig)) < rates.halt.durable
  ]

  tx.halt <- c(tx.halt.part, tx.halt.partial, tx.halt.full, tx.halt.durable)

  ## Restarting (partner services)
  part.tx.reinit.elig <- which(
    tx.status == 0 &
    part.ident == at &
    cuml.time.on.tx > 0
  )
  rates.reinit.part <- part.tx.reinit.prob[race[part.tx.reinit.elig]]
  tx.reinit.part <- part.tx.reinit.elig[
    runif(length(part.tx.reinit.elig)) < rates.reinit.part
  ]

  ## Restarting (general)
  tx.reinit.partial.elig <- which(
    tx.status == 0 &
    tt.traj == 1 &
    cuml.time.on.tx > 0
  )
  tx.reinit.partial.elig <- setdiff(tx.reinit.partial.elig, part.tx.reinit.elig)
  rates.reinit.partial <- tx.reinit.partial.prob[race[tx.reinit.partial.elig]]
  tx.reinit.partial <- tx.reinit.partial.elig[
    runif(length(tx.reinit.partial.elig)) < rates.reinit.partial
  ]

  tx.reinit.full.elig <- which(
    tx.status == 0 &
    tt.traj == 2 &
    cuml.time.on.tx > 0
  )
  tx.reinit.full.elig <- setdiff(tx.reinit.full.elig, part.tx.reinit.elig)
  rates.reinit.full <- tx.reinit.partial.prob[race[tx.reinit.full.elig]] *
                       tx.reinit.full.rr[race[tx.reinit.full.elig]]
  tx.reinit.full <- tx.reinit.full.elig[
    runif(length(tx.reinit.full.elig)) < rates.reinit.full
  ]

  tx.reinit.durable.elig <- which(
    tx.status == 0 &
    tt.traj == 3 &
    cuml.time.on.tx > 0
  )
  tx.reinit.durable.elig <- setdiff(
    tx.reinit.durable.elig,
    tx.reinit.partial.elig
  )
  rates.reinit.durable <- tx.reinit.partial.prob[race[tx.reinit.durable.elig]] *
                          tx.reinit.durable.rr[race[tx.reinit.durable.elig]]
  tx.reinit.durable <- tx.reinit.durable.elig[
    runif(length(tx.reinit.durable.elig)) < rates.reinit.durable
  ]

  tx.reinit <- c(
    tx.reinit.part,
    tx.reinit.partial,
    tx.reinit.full,
    tx.reinit.durable
  )

  ## Update Attributes
  tx.status[tx.init.all] <- 1
  tx.status[tx.halt] <- 0
  tx.status[tx.reinit] <- 1

  cuml.time.on.tx[which(tx.status == 1)] <- cuml.time.on.tx[which(tx.status == 1)] + 1
  cuml.time.off.tx[which(tx.status == 0)] <- cuml.time.off.tx[which(tx.status == 0)] + 1

  tx.init.time[tx.init.all] <- at

  dat <- set_attr(dat, "tt.traj", tt.traj)
  dat <- set_attr(dat, "tx.status", tx.status)
  dat <- set_attr(dat, "cuml.time.on.tx", cuml.time.on.tx)
  dat <- set_attr(dat, "cuml.time.off.tx", cuml.time.off.tx)
  dat <- set_attr(dat, "tx.init.time", tx.init.time)

  dat <- set_attr(dat, "part.tx.init.time", at, posit_ids = part.tx.init)
  dat <- set_attr(dat, "part.tx.reinit.time", at, posit_ids = tx.reinit.part)

  dat <- set_epi(dat, "mean.tx.on", at, mean(cuml.time.on.tx, na.rm = TRUE))
  dat <- set_epi(dat, "mean.tx.off", at, mean(cuml.time.off.tx, na.rm = TRUE))
  dat <- set_epi(dat, "mean.tx.on.part", at,
                 mean(cuml.time.on.tx[tt.traj == 1], na.rm = TRUE))
  dat <- set_epi(dat, "mean.tx.off.part", at,
                 mean(cuml.time.off.tx[tt.traj == 1], na.rm = TRUE))

  return(dat)
}
