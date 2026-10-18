
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

  part.ident         <- get_attr(dat, "part.ident")
  mdeCurr            <- get_attr(dat, "mdeCurr")
  mdd.tx.status      <- get_attr(dat, "mdd.tx.status")
  mdd.currtx.tottime <- get_attr(dat, "mdd.currtx.tottime")


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

  txinit.rate.mult       <- get_param(dat, "txinit.rate.mult")
  txhalt.rate.mult       <- get_param(dat, "txhalt.rate.mult")
  minadeq.txdur          <- get_param(dat, "minadeq.txdur")



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
  #add depression effects
  if(at >= get_param(dat, "mhefx.start")) {
    mdeCurr.elig <- mdeCurr[part.tx.init.elig]

    #remove mdd impact for those on mdd tx, who have received at least
    #minimally adequate tx
    minadeqmddtx <- which(mdd.tx.status[part.tx.init.elig] == 1 &
                            mdd.currtx.tottime[part.tx.init.elig] >= minadeq.txdur)
    mdeCurr.elig[minadeqmddtx] <- 0

    txinit.mult <- txinit.rate.mult[mdeCurr.elig + 1]
    rates.init.part <- exp(log(rates.init.part) + log(txinit.mult))
  }

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
  #add depression effects
  if(at >= get_param(dat, "mhefx.start")) {
    mdeCurr.elig <- mdeCurr[tx.init.elig]

    #remove mdd impact for those on mdd tx, who have received at least
    #minimally adequate tx
    minadeqmddtx <- which(mdd.tx.status[tx.init.elig] == 1 &
                            mdd.currtx.tottime[tx.init.elig] >= minadeq.txdur)
    mdeCurr.elig[minadeqmddtx] <- 0

    txinit.mult <- txinit.rate.mult[mdeCurr.elig + 1]
    rates.init.gen <- exp(log(rates.init.gen) + log(txinit.mult))
  }

  tx.init <- tx.init.elig[
    runif(length(tx.init.elig)) < rates.init.gen
  ]


  tx.init.all <- c(tx.init, part.tx.init)



  ## Halting (general)
  tx.halt.partial.elig <- which(tx.status == 1 & tt.traj == 1)
  rates.halt.partial <- tx.halt.partial.rate[race[tx.halt.partial.elig]]
  #add depression effects
  if(at >= get_param(dat, "mhefx.start")) {
    mdeCurr.elig <- mdeCurr[tx.halt.partial.elig]

    #remove mdd impact for those on mdd tx, who have received at least
    #minimally adequate tx
    minadeqmddtx <- which(mdd.tx.status[tx.halt.partial.elig] == 1 &
                            mdd.currtx.tottime[tx.halt.partial.elig] >= minadeq.txdur)
    mdeCurr.elig[minadeqmddtx] <- 0

    txhalt.mult <- txhalt.rate.mult[mdeCurr.elig + 1]
    rates.halt.partial <- plogis(qlogis(rates.halt.partial) + log(txhalt.mult))
  }

  tx.halt.partial <- tx.halt.partial.elig[
    runif(length(tx.halt.partial.elig)) < rates.halt.partial
  ]



  tx.halt.full.elig <- which(tx.status == 1 & tt.traj == 2)
  rates.halt.full <- tx.halt.partial.rate[race[tx.halt.full.elig]]
  halt.full.lo <- log(rates.halt.full / (1 - rates.halt.full))
  halt.full.lo <- halt.full.lo + log(tx.halt.full.or[race[tx.halt.full.elig]])
  rates.halt.full <- plogis(halt.full.lo)
  #add depression effects
  if(at >= get_param(dat, "mhefx.start")) {
    mdeCurr.elig <- mdeCurr[tx.halt.full.elig]

    #remove mdd impact for those on mdd tx, who have received at least
    #minimally adequate tx
    minadeqmddtx <- which(mdd.tx.status[tx.halt.full.elig] == 1 &
                            mdd.currtx.tottime[tx.halt.full.elig] >= minadeq.txdur)
    mdeCurr.elig[minadeqmddtx] <- 0

    txhalt.mult <- txhalt.rate.mult[mdeCurr.elig + 1]
    rates.halt.full <- plogis(qlogis(rates.halt.full) + log(txhalt.mult))
  }

  tx.halt.full <- tx.halt.full.elig[
    runif(length(tx.halt.full.elig)) < rates.halt.full
  ]



  tx.halt.durable.elig <- which(tx.status == 1 & tt.traj == 3)
  rates.halt.durable <- tx.halt.partial.rate[race[tx.halt.durable.elig]]
  halt.durable.lo <- log(rates.halt.durable / (1 - rates.halt.durable))
  halt.durable.lo <- halt.durable.lo + log(tx.halt.durable.or[race[tx.halt.durable.elig]])
  rates.halt.durable <- plogis(halt.durable.lo)
  #add depression effects
  if(at >= get_param(dat, "mhefx.start")) {
    mdeCurr.elig <- mdeCurr[tx.halt.durable.elig]

    #remove mdd impact for those on mdd tx, who have received at least
    #minimally adequate tx
    minadeqmddtx <- which(mdd.tx.status[tx.halt.durable.elig] == 1 &
                            mdd.currtx.tottime[tx.halt.durable.elig] >= minadeq.txdur)
    mdeCurr.elig[minadeqmddtx] <- 0

    txhalt.mult <- txhalt.rate.mult[mdeCurr.elig + 1]
    rates.halt.durable <- plogis(qlogis(rates.halt.durable) + log(txhalt.mult))
  }

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
  #add depression effects
  if(at >= get_param(dat, "mhefx.start")) {
    mdeCurr.elig <- mdeCurr[part.tx.reinit.elig]

    #remove mdd impact for those on mdd tx, who have received at least
    #minimally adequate tx
    minadeqmddtx <- which(mdd.tx.status[part.tx.reinit.elig] == 1 &
                            mdd.currtx.tottime[part.tx.reinit.elig] >= minadeq.txdur)
    mdeCurr.elig[minadeqmddtx] <- 0

    txinit.mult <- txinit.rate.mult[mdeCurr.elig + 1]
    rates.reinit.part <- exp(log(rates.reinit.part) + log(txinit.mult))
  }

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
  #add depression effects
  if(at >= get_param(dat, "mhefx.start")) {
    mdeCurr.elig <- mdeCurr[tx.reinit.partial.elig]

    #remove mdd impact for those on mdd tx, who have received at least
    #minimally adequate tx
    minadeqmddtx <- which(mdd.tx.status[tx.reinit.partial.elig] == 1 &
                            mdd.currtx.tottime[tx.reinit.partial.elig] >= minadeq.txdur)
    mdeCurr.elig[minadeqmddtx] <- 0

    txinit.mult <- txinit.rate.mult[mdeCurr.elig + 1]
    rates.reinit.partial <- plogis(qlogis(rates.reinit.partial) + log(txinit.mult))
  }

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
  #add depression effects
  if(at >= get_param(dat, "mhefx.start")) {
    mdeCurr.elig <- mdeCurr[tx.reinit.full.elig]

    #remove mdd impact for those on mdd tx, who have received at least
    #minimally adequate tx
    minadeqmddtx <- which(mdd.tx.status[tx.reinit.full.elig] == 1 &
                            mdd.currtx.tottime[tx.reinit.full.elig] >= minadeq.txdur)
    mdeCurr.elig[minadeqmddtx] <- 0

    txinit.mult <- txinit.rate.mult[mdeCurr.elig + 1]
    rates.reinit.full <- plogis(qlogis(rates.reinit.full) + log(txinit.mult))
  }

  tx.reinit.full <- tx.reinit.full.elig[
    runif(length(tx.reinit.full.elig)) < rates.reinit.full
  ]



  tx.reinit.durable.elig <- which(
    tx.status == 0 &
    tt.traj == 3 &
    cuml.time.on.tx > 0
  )
  tx.reinit.durable.elig <- setdiff(tx.reinit.durable.elig, tx.reinit.partial.elig)

  rates.reinit.durable <- tx.reinit.partial.rate[race[tx.reinit.durable.elig]]
  reinit.durable.lo <- log(rates.reinit.durable / (1 - rates.reinit.durable))
  reinit.durable.lo <- reinit.durable.lo +
    log(tx.reinit.durable.or[race[tx.reinit.durable.elig]])
  rates.reinit.durable <- plogis(reinit.durable.lo)
  #add depression effects
  if(at >= get_param(dat, "mhefx.start")) {
    mdeCurr.elig <- mdeCurr[tx.reinit.durable.elig]

    #remove mdd impact for those on mdd tx, who have received at least
    #minimally adequate tx
    minadeqmddtx <- which(mdd.tx.status[tx.reinit.durable.elig] == 1 &
                            mdd.currtx.tottime[tx.reinit.durable.elig] >= minadeq.txdur)
    mdeCurr.elig[minadeqmddtx] <- 0

    txinit.mult <- txinit.rate.mult[mdeCurr.elig + 1]
    rates.reinit.durable <- plogis(qlogis(rates.reinit.durable) + log(txinit.mult))
  }

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

  #epi trackers
  dat <- set_epi(dat, "mean.tx.on", at, mean(cuml.time.on.tx, na.rm = TRUE))
  dat <- set_epi(dat, "mean.tx.off", at, mean(cuml.time.off.tx, na.rm = TRUE))

  dat <- set_epi(dat, "all.init.tx", at, length(c(tx.init, part.tx.init)))
  dat <- set_epi(dat, "elig.gen.start.tx", at, length(tx.init.elig))
  dat <- set_epi(dat, "gen.start.tx", at, length(tx.init))
  dat <- set_epi(dat, "elig.part.start.tx", at, length(part.tx.init.elig))
  dat <- set_epi(dat, "part.start.tx", at, length(part.tx.init))

  dat <- set_epi(dat, "all.reinit.tx", at, length(tx.reinit))
  dat <- set_epi(dat, "gen.elig.for.reinit", at, length(c(tx.reinit.partial.elig,
                                                          tx.reinit.full.elig,
                                                          tx.reinit.durable.elig)))
  dat <- set_epi(dat, "gen.reinit.tx", at, length(c(tx.reinit.partial,
                                                    tx.reinit.full,
                                                    tx.reinit.durable)))
  dat <- set_epi(dat, "part.elig.for.reinit", at, length(part.tx.reinit.elig))
  dat <- set_epi(dat, "part.reinit.tx", at, length(tx.reinit.part))

  dat <- set_epi(dat, "txStop", at, length(tx.halt))



  return(dat)
}
