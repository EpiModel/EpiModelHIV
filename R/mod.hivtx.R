
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
#' as individual-level attributes in `tt.traj`). This module initiates ART
#' for treatment naive persons in the latter two types, and then cycles them on
#' and off treatment conditional on empirical race-specific adherence rates. ART
#' initiation, non-adherence, and restarting are all stochastically simulated
#' based on binomial statistical models.
#'
#' @return
#' This function returns the `dat` object with updated `tx.status`, `tx.init.time`,
#' `cuml.time.on.tx`, `cuml.time.off.tx` attributes.
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
  tx.init.rate           <- get_param(dat, "tx.init.rate")
  tx.halt.partial.rate   <- get_param(dat, "tx.halt.partial.rate")
  tx.halt.full.or        <- get_param(dat, "tx.halt.full.or")
  tx.halt.durable.or     <- get_param(dat, "tx.halt.durable.or")
  tx.reinit.partial.rate <- get_param(dat, "tx.reinit.partial.rate")
  tx.reinit.full.or      <- get_param(dat, "tx.reinit.full.or")
  tx.reinit.durable.or   <- get_param(dat, "tx.reinit.durable.or")

  part.tx.init.rate      <- get_param(dat, "part.tx.init.rate")
  part.tx.reinit.rate    <- get_param(dat, "part.tx.reinit.rate")

  ## Process
  # Initiation (partner services)
  part.tx.init.elig <- which(
    status == 1 &
      tx.status == 0 &
      diag.status == 1 &
      cuml.time.on.tx == 0 &
      part.ident == at
  )

  rates.init.part <- part.tx.init.rate[race[part.tx.init.elig]]
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
  rates.init.gen <- tx.init.rate[race[tx.init.elig]]
  tx.init <- tx.init.elig[
    runif(length(tx.init.elig)) < rates.init.gen
  ]

  tx.init.all <- c(tx.init, part.tx.init)

  ## Halting (general)
  tx.halt.partial.elig <- which(tx.status == 1 & tt.traj == 1)
  rates.halt.partial <- tx.halt.partial.rate[race[tx.halt.partial.elig]]
  tx.halt.partial <- tx.halt.partial.elig[
    runif(length(tx.halt.partial.elig)) < rates.halt.partial
  ]

  tx.halt.full.elig <- which(tx.status == 1 & tt.traj == 2)

  rates.halt.full <- tx.halt.partial.rate[race[tx.halt.full.elig]]
  halt.full.lo <- log(rates.halt.full / (1 - rates.halt.full))
  halt.full.lo <- halt.full.lo + log(tx.halt.full.or[race[tx.halt.full.elig]])
  rates.halt.full <- plogis(halt.full.lo)

  tx.halt.full <- tx.halt.full.elig[
    runif(length(tx.halt.full.elig)) < rates.halt.full
  ]

  tx.halt.durable.elig <- which(tx.status == 1 & tt.traj == 3)

  rates.halt.durable <- tx.halt.partial.rate[race[tx.halt.durable.elig]]
  halt.durable.lo <- log(rates.halt.durable / (1 - rates.halt.durable))
  halt.durable.lo <- halt.durable.lo +
    log(tx.halt.durable.or[race[tx.halt.durable.elig]])
  rates.halt.durable <- plogis(halt.durable.lo)

  tx.halt.durable <- tx.halt.durable.elig[
    runif(length(tx.halt.durable.elig)) < rates.halt.durable
  ]

  tx.halt <- c(tx.halt.partial, tx.halt.full, tx.halt.durable)

  ## Restarting (partner services)
  part.tx.reinit.elig <- which(
    tx.status == 0 &
      part.ident == at &
      cuml.time.on.tx > 0
  )
  rates.reinit.part <- part.tx.reinit.rate[race[part.tx.reinit.elig]]
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
  rates.reinit.partial <- tx.reinit.partial.rate[race[tx.reinit.partial.elig]]
  tx.reinit.partial <- tx.reinit.partial.elig[
    runif(length(tx.reinit.partial.elig)) < rates.reinit.partial
  ]

  tx.reinit.full.elig <- which(
    tx.status == 0 &
      tt.traj == 2 &
      cuml.time.on.tx > 0
  )
  tx.reinit.full.elig <- setdiff(tx.reinit.full.elig, part.tx.reinit.elig)

  rates.reinit.full <- tx.reinit.partial.rate[race[tx.reinit.full.elig]]
  reinit.full.lo <- log(rates.reinit.full / (1 - rates.reinit.full))
  reinit.full.lo <- reinit.full.lo +
    log(tx.reinit.full.or[race[tx.reinit.full.elig]])
  rates.reinit.full <- plogis(reinit.full.lo)

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

  rates.reinit.durable <- tx.reinit.partial.rate[race[tx.reinit.durable.elig]]
  reinit.durable.lo <- log(rates.reinit.durable / (1 - rates.reinit.durable))
  reinit.durable.lo <- reinit.durable.lo +
    log(tx.reinit.durable.or[race[tx.reinit.durable.elig]])
  rates.reinit.durable <- plogis(reinit.durable.lo)

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

  dat <- set_attr(dat, "tx.status", tx.status)
  dat <- set_attr(dat, "cuml.time.on.tx", cuml.time.on.tx)
  dat <- set_attr(dat, "cuml.time.off.tx", cuml.time.off.tx)
  dat <- set_attr(dat, "tx.init.time", tx.init.time)

  dat <- set_attr(dat, "part.tx.init.time", at, posit_ids = part.tx.init)
  dat <- set_attr(dat, "part.tx.reinit.time", at, posit_ids = tx.reinit.part)

  dat <- set_epi(dat, "mean.tx.on", at, mean(cuml.time.on.tx, na.rm = TRUE))
  dat <- set_epi(dat, "mean.tx.off", at, mean(cuml.time.off.tx, na.rm = TRUE))
  # dat <- set_epi(dat, "mean.tx.on.part", at,
  #                mean(cuml.time.on.tx[tt.traj == 1], na.rm = TRUE))
  # dat <- set_epi(dat, "mean.tx.off.part", at,
  #                mean(cuml.time.off.tx[tt.traj == 1], na.rm = TRUE))

  return(dat)
}
