
#' @title STI Transmission Module
#'
#' @description Stochastically simulates GC/CT transmission given the current
#'              state of the edgelist.
#'
#' @inheritParams aging_msm
#'
#' @export
#'
stitrans_msm <- function(dat, at) {

  ## Inputs
  # Attributes
  race <- get_attr(dat, "race")

  # Current infection state
  rGC <- get_attr(dat, "rGC")
  uGC <- get_attr(dat, "uGC")
  rCT <- get_attr(dat, "rCT")
  uCT <- get_attr(dat, "uCT")

  # n Times infected
  rGC.timesInf <- get_attr(dat, "rGC.timesInf")
  uGC.timesInf <- get_attr(dat, "uGC.timesInf")
  rCT.timesInf <- get_attr(dat, "rCT.timesInf")
  uCT.timesInf <- get_attr(dat, "uCT.timesInf")

  # Infection time
  rGC.infTime <- get_attr(dat, "rGC.infTime")
  uGC.infTime <- get_attr(dat, "uGC.infTime")
  rCT.infTime <- get_attr(dat, "rCT.infTime")
  uCT.infTime <- get_attr(dat, "uCT.infTime")

  # Infection symptoms (non-varying)
  rGC.sympt <- get_attr(dat, "rGC.sympt")
  uGC.sympt <- get_attr(dat, "uGC.sympt")
  rCT.sympt <- get_attr(dat, "rCT.sympt")
  uCT.sympt <- get_attr(dat, "uCT.sympt")

  # Parameters
  # Acquisition probabilities given contact with infected man
  rgc.prob <- get_param(dat, "rgc.prob")
  ugc.prob <- get_param(dat, "ugc.prob")
  rct.prob <- get_param(dat, "rct.prob")
  uct.prob <- get_param(dat, "uct.prob")

  # Probability of symptoms given infection
  rgc.sympt.prob <- get_param(dat, "rgc.sympt.prob")
  ugc.sympt.prob <- get_param(dat, "ugc.sympt.prob")
  rct.sympt.prob <- get_param(dat, "rct.sympt.prob")
  uct.sympt.prob <- get_param(dat, "uct.sympt.prob")

  # Relative risk of infection given condom use during act
  sti.cond.eff.or    <- get_param(dat, "sti.cond.eff.or")
  sti.cond.fail.or   <- get_param(dat, "sti.cond.fail.or")

  # Pull act list
  al <- dat[["temp"]][["al"]]

  ## ins variable coding
  # ins = 0 : p2 is insertive
  # ins = 1 : p1 is insertive
  # ins = 2 : both p1 and p2 are insertive


  # Rectal GC -----------------------------------------------------------

  # Requires: uGC in insertive man, and no rGC in receptive man
  p1Inf_rgc <- which(
    uGC[al[, "p1"]] == 1 &
    uGC.infTime[al[, "p1"]] < at &
    rGC[al[, "p2"]] == 0 &
    al[, "ins"] %in% c(1, 2)
  )
  p2Inf_rgc <- which(
    uGC[al[, "p2"]] == 1 &
    uGC.infTime[al[, "p2"]] < at &
    rGC[al[, "p1"]] == 0 &
    al[, "ins"] %in% c(0, 2)
  )
  allActs_rgc <- c(p1Inf_rgc, p2Inf_rgc)

  # UAI modifier
  uai_rgc <- al[allActs_rgc, "uai"]
  tprob_rgc <- rep(rgc.prob, length(allActs_rgc))

  # Transform to log odds
  tlo_rgc <- log(tprob_rgc / (1 - tprob_rgc))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_rgc, "p1"]], race[al[p2Inf_rgc, "p2"]])
  condom.or <- rep(NA, length(races))
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.or[ids.race] <- 1 - (sti.cond.eff.or - sti.cond.fail.or[i])
  }

  tlo_rgc[uai_rgc == 0] <- tlo_rgc[uai_rgc == 0] + log(condom.or[uai_rgc == 0])

  # Back-transform to probability
  tprob_rgc <- plogis(tlo_rgc)

  # Stochastic transmission
  trans_rgc <- runif(length(allActs_rgc)) < tprob_rgc

  # Determine the infected partner
  idsInf_rgc <- numeric()
  if (sum(trans_rgc) > 0) {
    transAL_rgc <- al[allActs_rgc[trans_rgc == 1], , drop = FALSE]
    idsInf_rgc <- union(
      intersect(al[p1Inf_rgc, "p2"], transAL_rgc[, "p2"]),
      intersect(al[p2Inf_rgc, "p1"], transAL_rgc[, "p1"])
    )
    stopifnot(all(rGC[idsInf_rgc] == 0))
  }

  # Update attributes
  rGC[idsInf_rgc] <- 1
  rGC.infTime[idsInf_rgc] <- at
  rGC.sympt[idsInf_rgc] <- runif(length(idsInf_rgc)) < rgc.sympt.prob
  rGC.timesInf[idsInf_rgc] <- rGC.timesInf[idsInf_rgc] + 1


  # Urethral GC ---------------------------------------------------------

  # Requires: rGC in receptive man, and no uGC in insertive man
  p1Inf_ugc <- which(
    rGC[al[, "p1"]] == 1 &
    rGC.infTime[al[, "p1"]] < at &
    uGC[al[, "p2"]] == 0 &
    al[, "ins"] %in% c(0, 2)
  )
  p2Inf_ugc <- which(
    rGC[al[, "p2"]] == 1 &
    rGC.infTime[al[, "p2"]] < at &
    uGC[al[, "p1"]] == 0 &
    al[, "ins"] %in% c(1, 2)
  )
  allActs_ugc <- c(p1Inf_ugc, p2Inf_ugc)

  # UAI modifier
  uai_ugc <- al[allActs_ugc, "uai"]
  tprob_ugc <- rep(ugc.prob, length(allActs_ugc))

  # Transform to log odds
  tlo_ugc <- log(tprob_ugc / (1 - tprob_ugc))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_ugc, "p2"]], race[al[p2Inf_ugc, "p1"]])
  condom.or <- rep(NA, length(races))
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.or[ids.race] <- 1 - (sti.cond.eff.or - sti.cond.fail.or[i])
  }

  tlo_ugc[uai_ugc == 0] <- tlo_ugc[uai_ugc == 0] + log(condom.or[uai_ugc == 0])

  # Back-transform to probability
  tprob_ugc <- plogis(tlo_ugc)

  # Stochastic transmission
  trans_ugc <- runif(length(allActs_ugc)) < tprob_ugc

  # Determine the newly infected partner
  idsInf_ugc <- numeric()
  if (sum(trans_ugc) > 0) {
    transAL_ugc <- al[allActs_ugc[trans_ugc == 1],  , drop = FALSE]
    idsInf_ugc <- union(
      intersect(al[p1Inf_ugc, "p2"], transAL_ugc[, "p2"]),
      intersect(al[p2Inf_ugc, "p1"], transAL_ugc[, "p1"])
    )
    stopifnot(all(uGC[idsInf_ugc] == 0))
  }

  # Update attributes
  uGC[idsInf_ugc] <- 1
  uGC.infTime[idsInf_ugc] <- at
  uGC.sympt[idsInf_ugc] <- runif(length(idsInf_ugc)) < ugc.sympt.prob
  uGC.timesInf[idsInf_ugc] <- uGC.timesInf[idsInf_ugc] + 1


  # Rectal CT -----------------------------------------------------------

  # Requires: uCT in insertive man, and no rCT in receptive man
  p1Inf_rct <- which(
    uCT[al[, "p1"]] == 1 &
      uCT.infTime[al[, "p1"]] < at &
      rCT[al[, "p2"]] == 0 &
      al[, "ins"] %in% c(1, 2)
  )
  p2Inf_rct <- which(
    uCT[al[, "p2"]] == 1 &
    uCT.infTime[al[, "p2"]] < at &
    rCT[al[, "p1"]] == 0 &
    al[, "ins"] %in% c(0, 2)
  )
  allActs_rct <- c(p1Inf_rct, p2Inf_rct)

  # UAI modifier
  uai_rct <- al[allActs_rct, "uai"]
  tprob_rct <- rep(rct.prob, length(allActs_rct))

  # Transform to log odds
  tlo_rct <- log(tprob_rct / (1 - tprob_rct))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_rct, "p1"]], race[al[p2Inf_rct, "p2"]])
  condom.or <- rep(NA, length(races))
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.or[ids.race] <- 1 - (sti.cond.eff.or - sti.cond.fail.or[i])
  }

  tlo_rct[uai_rct == 0] <- tlo_rct[uai_rct == 0] + log(condom.or[uai_rct == 0])

  # Back-transform to probability
  tprob_rct <- plogis(tlo_rct)

  # Stochastic transmission
  trans_rct <- runif(length(allActs_rct)) < tprob_rct

  # Determine the newly infected partner
  idsInf_rct <- numeric()
  if (sum(trans_rct) > 0) {
    transAL_rct <- al[allActs_rct[trans_rct == 1],  , drop = FALSE]
    idsInf_rct <- union(
      intersect(al[p1Inf_rct, "p2"], transAL_rct[, "p2"]),
      intersect(al[p2Inf_rct, "p1"], transAL_rct[, "p1"])
    )
    stopifnot(all(rCT[idsInf_rct] == 0))
  }

  # Update attributes
  rCT[idsInf_rct] <- 1
  rCT.infTime[idsInf_rct] <- at
  rCT.sympt[idsInf_rct] <- runif(length(idsInf_rct)) < rct.sympt.prob
  rCT.timesInf[idsInf_rct] <- rCT.timesInf[idsInf_rct] + 1


  # Urethral CT ---------------------------------------------------------

  # Requires: rCT in receptive man, and no uCT in insertive man
  p1Inf_uct <- which(
    rCT[al[, "p1"]] == 1 &
    rCT.infTime[al[, "p1"]] < at &
    uCT[al[, "p2"]] == 0 &
    al[, "ins"] %in% c(0, 2)
  )
  p2Inf_uct <- which(
    rCT[al[, "p2"]] == 1 &
    rCT.infTime[al[, "p2"]] < at &
    uCT[al[, "p1"]] == 0 &
    al[, "ins"] %in% c(1, 2)
  )
  allActs_uct <- c(p1Inf_uct, p2Inf_uct)

  # UAI modifier
  uai_uct <- al[allActs_uct, "uai"]
  tprob_uct <- rep(uct.prob, length(allActs_uct))

  # Transform to log odds
  tlo_uct <- log(tprob_uct / (1 - tprob_uct))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_uct, "p2"]], race[al[p2Inf_uct, "p1"]])
  condom.or <- rep(NA, length(races))
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.or[ids.race] <- 1 - (sti.cond.eff.or - sti.cond.fail.or[i])
  }

  tlo_uct[uai_uct == 0] <- tlo_uct[uai_uct == 0] + log(condom.or[uai_uct == 0])

  # Back-transform to probability
  tprob_uct <- plogis(tlo_uct)

  # Stochastic transmission
  trans_uct <- runif(length(allActs_uct)) < tprob_uct

  # Determine the newly infected partner
  idsInf_uct <- numeric()
  if (sum(trans_uct) > 0) {
    transAL_uct <- al[allActs_uct[trans_uct == 1],  , drop = FALSE]
    idsInf_uct <- union(
      intersect(al[p1Inf_uct, "p2"], transAL_uct[, "p2"]),
      intersect(al[p2Inf_uct, "p1"], transAL_uct[, "p1"])
    )
    stopifnot(all(uCT[idsInf_uct] == 0))
  }

  # Update attributes
  uCT[idsInf_uct] <- 1
  uCT.infTime[idsInf_uct] <- at
  uCT.sympt[idsInf_uct] <- runif(length(idsInf_uct)) < uct.sympt.prob
  uCT.timesInf[idsInf_uct] <- uCT.timesInf[idsInf_uct] + 1


  # Output --------------------------------------------------------------

  # attributes
  dat <- set_attr(dat, "rGC", rGC)
  dat <- set_attr(dat, "uGC", uGC)
  dat <- set_attr(dat, "rCT", rCT)
  dat <- set_attr(dat, "uCT", uCT)

  dat <- set_attr(dat, "rGC.infTime", rGC.infTime)
  dat <- set_attr(dat, "uGC.infTime", uGC.infTime)
  dat <- set_attr(dat, "rCT.infTime", rCT.infTime)
  dat <- set_attr(dat, "uCT.infTime", uCT.infTime)

  dat <- set_attr(dat, "rGC.timesInf", rGC.timesInf)
  dat <- set_attr(dat, "uGC.timesInf", uGC.timesInf)
  dat <- set_attr(dat, "rCT.timesInf", rCT.timesInf)
  dat <- set_attr(dat, "uCT.timesInf", uCT.timesInf)

  dat <- set_attr(dat, "rGC.sympt", rGC.sympt)
  dat <- set_attr(dat, "uGC.sympt", uGC.sympt)
  dat <- set_attr(dat, "rCT.sympt", rCT.sympt)
  dat <- set_attr(dat, "uCT.sympt", uCT.sympt)


  # Summary stats
  dat <- set_epi(dat, "incid.gc", at, length(union(idsInf_rgc, idsInf_ugc)))
  dat <- set_epi(dat, "incid.rgc", at, length(idsInf_rgc))
  dat <- set_epi(dat, "incid.ugc", at, length(idsInf_ugc))

  dat <- set_epi(dat, "incid.ct", at, length(union(idsInf_rct, idsInf_uct)))
  dat <- set_epi(dat, "incid.rct", at, length(idsInf_rct))
  dat <- set_epi(dat, "incid.uct", at, length(idsInf_uct))

  # Check all infected have all STI attributes
  stopifnot(
    all(!is.na(rGC.infTime[rGC == 1])),
    all(!is.na(rGC.sympt[rGC == 1])),
    all(!is.na(uGC.infTime[uGC == 1])),
    all(!is.na(uGC.sympt[uGC == 1])),
    all(!is.na(rCT.infTime[rCT == 1])),
    all(!is.na(rCT.sympt[rCT == 1])),
    all(!is.na(uCT.infTime[uCT == 1])),
    all(!is.na(uCT.sympt[uCT == 1]))
  )

  return(dat)
}
