
#' @title PrEP Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging_msm
#'
#' @export
#'
prep_msm <- function(dat, at) {

  # Function Selection ------------------------------------------------------

  if (at >= get_param(dat, "riskh.start")) {
    dat <- riskhist_msm(dat, at)
  } else {
    return(dat)
  }

  if (at < get_param(dat, "prep.start")) {
    return(dat)
  }

  ## Input
  # Core Attributes
  active      <- get_attr(dat, "active")
  diag.status <- get_attr(dat, "diag.status")
  lnt         <- get_attr(dat, "last.neg.test")
  part.scrnd  <- get_attr(dat, "part.scrnd")
  race        <- get_attr(dat, "race")

  # Refresh PrEP duration and PrEP discontinuation duration attributes
  if (at == 3901) {
    n = length(active)
    dat <- set_attr(dat, "prepDuration", rep(NA, n))
    dat <- set_attr(dat, "prepStopDuration", rep(NA, n))
    dat <- set_attr(dat, "prepStopTime", rep(NA, n))
  }

  # PrEP Attributes
  prepElig           <- get_attr(dat, "prepElig")
  prepStat           <- get_attr(dat, "prepStat")
  prepClass          <- get_attr(dat, "prepClass")
  prepLastRisk       <- get_attr(dat, "prepLastRisk")
  prepStartTime      <- get_attr(dat, "prepStartTime")
  prep.start.counter <- get_attr(dat, "prep.start.counter")
  prep.stop.counter  <- get_attr(dat, "prep.stop.counter")
  stopper.speed      <- get_attr(dat, "stopper.speed")
  prepStopTime       <- get_attr(dat, "prepStopTime")
  prepDuration       <- get_attr(dat, "prepDuration")
  prepStopDuration   <- get_attr(dat, "prepStopDuration")

  # Re-assign stopper.speed for persistence scenarios
  stopper.reassign <- get_param(dat, "stopper.reassign")
  if (at == 3901 & stopper.reassign) {
    stopper.speed <- make_stopper_speed(dat, race)
    dat <- set_attr(dat, "stopper.speed", stopper.speed)
  }

  # Indications
  ind1 <- get_attr(dat, "prep.ind.uai.mono")
  ind2 <- get_attr(dat, "prep.ind.uai.conc")
  ind3 <- get_attr(dat, "prep.ind.sti")

  # Parameters
  prep.start.prob           <- get_param(dat, "prep.start.prob")
  part.prep.start.prob      <- get_param(dat, "part.prep.start.prob")
  prep.adhr.dist            <- get_param(dat, "prep.adhr.dist")
  prep.require.lnt          <- get_param(dat, "prep.require.lnt")
  prep.risk.reassess.int    <- get_param(dat, "prep.risk.reassess.int")
  prep.risk.int             <- get_param(dat, "prep.risk.int")

  twind <- at - prep.risk.int

  # New PrEP parameters----------------------------------------------------
  prep.reinit.start         <- get_param(dat, "prep.reinit.start")
  prep.reinit.slow.rate     <- get_param(dat, "prep.reinit.slow.rate")
  prep.reinit.fast.rate     <- get_param(dat, "prep.reinit.fast.rate")
  prep.fast.discont.rate    <- get_param(dat, "prep.fast.discont.rate")
  prep.slow.discont.rate    <- get_param(dat, "prep.slow.discont.rate")

  # No indications in window-----------------------------------------------
  idsNoIndic <- which(
    (ind1 < twind | is.na(ind1)) &
      (ind2 < twind | is.na(ind2)) &
      (ind3 < twind | is.na(ind3))
  )
  base.cond.no <- which(!active | diag.status == 1)
  idsNoIndic <- union(idsNoIndic, base.cond.no)

  # Indications in window
  idsIndic <- which(ind1 >= twind | ind2 >= twind | ind3 >= twind)
  base.cond.yes <- which(active & diag.status == 0)
  idsIndic <- intersect(idsIndic, base.cond.yes)

  # Set eligibility to 1 if indications
  prepElig[idsIndic] <- 1

  # Set eligibility to 0 if no indications
  prepElig[idsNoIndic] <- 0

  ## Stoppage ------------------------------------------------------------------

  # Indication lapse
  # Rules = None, instant, yearly (CDC guidelines)
  # generalize for when risk reassessment should happen in general language
  idsRiskAssess <- which(
    active &
      prepStat == 1 &
      lnt == at &
      (at - prepLastRisk) >= prep.risk.reassess.int
  )
  prepLastRisk[idsRiskAssess] <- at
  idsStpInd <- intersect(idsNoIndic, idsRiskAssess)

  # Random discontinuation
  idsEligStpRand.fast <- which(active & prepStat == 1 & stopper.speed == 1)
  prob.fast.discont <- prep.fast.discont.rate[race[idsEligStpRand.fast]]
  vecStpRand.fast <- runif(length(idsEligStpRand.fast)) < prob.fast.discont
  idsStpRand.fast <- idsEligStpRand.fast[vecStpRand.fast]

  idsEligStpRand.slow <- which(active & prepStat == 1 & stopper.speed == 0)
  prob.slow.discont <- prep.slow.discont.rate[race[idsEligStpRand.slow]]
  vecStpRand.slow <- runif(length(idsEligStpRand.slow)) < prob.slow.discont
  idsStpRand.slow <- idsEligStpRand.slow[vecStpRand.slow]

  idsStpRand <- c(idsStpRand.fast, idsStpRand.slow)

  # Diagnosis
  idsStpDx <- which(active & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(!active & prepStat == 1)

  # Reset PrEP status
  idsStp <- c(idsStpInd, idsStpDx, idsStpDth)

  # Update attributes for stoppers
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA

  # Update attributes for stoppers due to random discontinuation
  prepStat[idsStpRand] <- 0
  prepLastRisk[idsStpRand] <- NA
  prepDuration[idsStpRand] <- at - prepStartTime[idsStpRand]
  prepStopTime[idsStpRand] <- at

  prepStartTime[idsStp] <- NA
  prepStopTime[idsStp] <- at
  prepStartTime[idsStpRand] <- NA

  # Update number of random discontinuations
  prep.stop.counter[idsStpRand] <- ifelse(
    is.na(prep.stop.counter[idsStpRand]),
    yes = 1,
    no = prep.stop.counter[idsStpRand] + 1
  )


  ## Initiation ----------------------------------------------------------------

  ## Eligibility ##

  # Indications to start: identified partners
  idsStartPart <- numeric()
  idsEligStartPart <- numeric()

  idsEligStartPart <- which(
    part.scrnd == at &
      prepStat == 0 &
      prepElig == 1
  )

  prob.start.part <- part.prep.start.prob[race[idsEligStartPart]]
  vecStartPart <- runif(length(idsEligStartPart)) < prob.start.part
  idsStartPart <- idsEligStartPart[vecStartPart]

  #NOTE: for validation only
  dat <- set_epi(dat, "prepStartPart", at, length(idsStartPart))

  # Update partner attributes
  dat <- set_attr(dat, "prep.start.part", at, posit_ids = idsStartPart)

  # Indications to start: general prep initiation
  if (prep.require.lnt) {
    idsEligStart.new         <- which(prepElig & prepStat == 0 & lnt == at & is.na(prep.start.counter))
    idsEligStart.reinit.fast <- which(prepElig & prepStat == 0 & lnt == at & prep.start.counter > 0 & stopper.speed == 1)
    idsEligStart.reinit.slow <- which(prepElig & prepStat == 0 & lnt == at & prep.start.counter > 0 & stopper.speed == 0)
  } else {
    idsEligStart.new         <- which(prepElig & prepStat == 0 & is.na(prep.start.counter))
    idsEligStart.reinit.fast <- which(prepElig & prepStat == 0 & prep.start.counter > 0 & stopper.speed == 1)
    idsEligStart.reinit.slow <- which(prepElig & prepStat == 0 & prep.start.counter > 0 & stopper.speed == 0)
    idsEligStart <- which(prepElig & prepStat == 0)
  }

  idsEligStart.new <- setdiff(idsEligStart.new, idsEligStartPart)
  idsEligStart.reinit.fast <- setdiff(idsEligStart.reinit.fast, idsEligStartPart)
  idsEligStart.reinit.slow <- setdiff(idsEligStart.reinit.slow, idsEligStartPart)

  prob.start <- prep.start.prob[race[idsEligStart.new]]
  vecStart <- runif(length(idsEligStart.new)) < prob.start
  idsStart.new <- idsEligStart.new[vecStart]

  if (at < prep.reinit.start) {
    prob.start.fast <- prep.start.prob[race[idsEligStart.reinit.fast]]
    vecStart.reinit.fast <- runif(length(idsEligStart.reinit.fast)) < prob.start.fast
    idsStart.reinit.fast <- idsEligStart.reinit.fast[vecStart.reinit.fast]

    prob.start.slow <- prep.start.prob[race[idsEligStart.reinit.slow]]
    vecStart.reinit.slow <- runif(length(idsEligStart.reinit.slow)) < prob.start.slow
    idsStart.reinit.slow <- idsEligStart.reinit.slow[vecStart.reinit.slow]
  } else {
    prob.start.reinit.fast <- prep.reinit.fast.rate[race[idsEligStart.reinit.fast]]
    vecStart.reinit.fast <- runif(length(idsEligStart.reinit.fast)) < prob.start.reinit.fast
    idsStart.reinit.fast <- idsEligStart.reinit.fast[vecStart.reinit.fast]

    prob.start.reinit.slow <- prep.reinit.slow.rate[race[idsEligStart.reinit.slow]]
    vecStart.reinit.slow <- runif(length(idsEligStart.reinit.slow)) < prob.start.reinit.slow
    idsStart.reinit.slow <- idsEligStart.reinit.slow[vecStart.reinit.slow]
  }

  ## Attribute Update ##
  idsPrepStart <- c(idsStartPart, idsStart.new, idsStart.reinit.fast, idsStart.reinit.slow)
  idsReinitiate <- c(idsStart.reinit.fast, idsStart.reinit.slow)

  # Set attributes for starters
  if (length(idsPrepStart) > 0) {
    prepStat[idsPrepStart] <- 1
    prepStopDuration[idsReinitiate] <- at - prepStopTime[idsReinitiate]
    prepStartTime[idsPrepStart] <- at
    prepLastRisk[idsPrepStart] <- at
    prepStopTime[idsPrepStart] <- NA

    prep.start.counter[idsPrepStart] <- ifelse(
      is.na(prep.start.counter[idsPrepStart]),
      yes = 1,
      no = prep.start.counter[idsPrepStart] + 1
    )

    # PrEP adherence class
    needPC <- which(is.na(prepClass[idsPrepStart]))
    prepClass[idsPrepStart[needPC]] <- sample(
      x = 1:3,
      size = length(needPC),
      replace = TRUE,
      prob = prep.adhr.dist
    )
  }


  ## Output --------------------------------------------------------------------

  # Random discontinuation
  dat <- set_epi(dat, "prep.rand.stop", at, length(idsStpRand))
  dat <- set_epi(dat, "prep.rand.stop.B", at, length(idsStpRand[race == 1]))
  dat <- set_epi(dat, "prep.rand.stop.H", at, length(idsStpRand[race == 2]))
  dat <- set_epi(dat, "prep.rand.stop.W", at, length(idsStpRand[race == 3]))
  dat <- set_epi(dat, "prep.rand.stop.O", at, length(idsStpRand[race == 4]))

  # Attributes
  dat <- set_attr(dat, "prepElig", prepElig)
  dat <- set_attr(dat, "prepStat", prepStat)
  dat <- set_attr(dat, "prepClass", prepClass)

  dat <- set_attr(dat, "prepStartTime", prepStartTime)
  dat <- set_attr(dat, "prepLastRisk", prepLastRisk)

  # New Attributes
  dat <- set_attr(dat, "prep.start.counter", prep.start.counter)
  dat <- set_attr(dat, "prep.stop.counter", prep.stop.counter)
  dat <- set_attr(dat, "prepStopTime", prepStopTime)
  dat <- set_attr(dat, "prepDuration", prepDuration)
  dat <- set_attr(dat, "prepStopDuration", prepStopDuration)

  # New Trackers - fast vs. slow stoppers by race
  dat <- set_epi(dat, "fastStopper", at,
                 sum(stopper.speed == 1, na.rm = TRUE))
  dat <- set_epi(dat, "slowStopper", at,
                 sum(stopper.speed == 0, na.rm = TRUE))

  dat <- set_epi(dat, "fastStopper.B", at,
                 sum(stopper.speed == 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "fastStopper.H", at,
                 sum(stopper.speed == 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "fastStopper.W", at,
                 sum(stopper.speed == 1 & race == 3, na.rm = TRUE))
  dat <- set_epi(dat, "fastStopper.O", at,
                 sum(stopper.speed == 1 & race == 4, na.rm = TRUE))

  dat <- set_epi(dat, "slowStopper.B", at,
                 sum(stopper.speed == 0 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "slowStopper.H", at,
                 sum(stopper.speed == 0 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "slowStopper.W", at,
                 sum(stopper.speed == 0 & race == 3, na.rm = TRUE))
  dat <- set_epi(dat, "slowStopper.O", at,
                 sum(stopper.speed == 0 & race == 4, na.rm = TRUE))

  # Fast and slow among PrEP users
  dat <- set_epi(dat, "prep.fastStopper", at,
                 sum(stopper.speed == 1 & prepStat == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.slowStopper", at,
                 sum(stopper.speed == 0 & prepStat == 1, na.rm = TRUE))

  dat <- set_epi(dat, "prep.fastStopper.B", at,
                 sum(stopper.speed == 1 & race == 1 & prepStat == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.fastStopper.H", at,
                 sum(stopper.speed == 1 & race == 2 & prepStat == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.fastStopper.W", at,
                 sum(stopper.speed == 1 & race == 3 & prepStat == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.fastStopper.O", at,
                 sum(stopper.speed == 1 & race == 4 & prepStat == 1, na.rm = TRUE))

  dat <- set_epi(dat, "prep.slowStopper.B", at,
                 sum(stopper.speed == 0 & race == 1 & prepStat == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.slowStopper.H", at,
                 sum(stopper.speed == 0 & race == 2 & prepStat == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.slowStopper.W", at,
                 sum(stopper.speed == 0 & race == 3 & prepStat == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.slowStopper.O", at,
                 sum(stopper.speed == 0 & race == 4 & prepStat == 1, na.rm = TRUE))

  # New tracker - PrEP initiation and reinitiation by race
  dat <- set_epi(dat, "prepStarters", at,
                 length(idsPrepStart))

  dat <- set_epi(dat, "prep.start.counter", at,
                 sum(prep.start.counter > 1, na.rm = TRUE))

  dat <- set_epi(dat, "prep.start.counter.B", at,
                 sum(prep.start.counter > 1 & race == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.start.counter.H", at,
                 sum(prep.start.counter > 1 & race == 2, na.rm = TRUE))
  dat <- set_epi(dat, "prep.start.counter.W", at,
                 sum(prep.start.counter > 1 & race == 3, na.rm = TRUE))
  dat <- set_epi(dat, "prep.start.counter.O", at,
                 sum(prep.start.counter > 1 & race == 4, na.rm = TRUE))

  # by stopper.speed
  dat <- set_epi(dat, "prep.slow.start.counter.B", at,
                 sum(prep.start.counter > 1 & race == 1 & stopper.speed == 0, na.rm = TRUE))
  dat <- set_epi(dat, "prep.slow.start.counter.H", at,
                 sum(prep.start.counter > 1 & race == 2 & stopper.speed == 0, na.rm = TRUE))
  dat <- set_epi(dat, "prep.slow.start.counter.W", at,
                 sum(prep.start.counter > 1 & race == 3 & stopper.speed == 0, na.rm = TRUE))
  dat <- set_epi(dat, "prep.slow.start.counter.O", at,
                 sum(prep.start.counter > 1 & race == 4 & stopper.speed == 0, na.rm = TRUE))

  dat <- set_epi(dat, "prep.fast.start.counter.B", at,
                 sum(prep.start.counter > 1 & race == 1 & stopper.speed == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.fast.start.counter.H", at,
                 sum(prep.start.counter > 1 & race == 2 & stopper.speed == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.fast.start.counter.W", at,
                 sum(prep.start.counter > 1 & race == 3 & stopper.speed == 1, na.rm = TRUE))
  dat <- set_epi(dat, "prep.fast.start.counter.O", at,
                 sum(prep.start.counter > 1 & race == 4 & stopper.speed == 1, na.rm = TRUE))

  # PrEP duration trackers
  dat <- set_epi(dat, "prepDuration", at, mean(prepDuration, na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.B", at, mean(prepDuration[race == 1], na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.H", at, mean(prepDuration[race == 2], na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.W", at, mean(prepDuration[race == 3], na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.O", at, mean(prepDuration[race == 4], na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.fast", at, mean(prepDuration[stopper.speed == 1], na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.slow", at, mean(prepDuration[stopper.speed == 0], na.rm = TRUE))

  dat <- set_epi(dat, "prepDuration.slow.B", at, mean(prepDuration[race == 1 & stopper.speed == 0], na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.slow.H", at, mean(prepDuration[race == 2 & stopper.speed == 0], na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.slow.W", at, mean(prepDuration[race == 3 & stopper.speed == 0], na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.slow.O", at, mean(prepDuration[race == 4 & stopper.speed == 0], na.rm = TRUE))

  dat <- set_epi(dat, "prepDuration.fast.B", at, mean(prepDuration[race == 1 & stopper.speed == 1], na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.fast.H", at, mean(prepDuration[race == 2 & stopper.speed == 1], na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.fast.W", at, mean(prepDuration[race == 3 & stopper.speed == 1], na.rm = TRUE))
  dat <- set_epi(dat, "prepDuration.fast.O", at, mean(prepDuration[race == 4 & stopper.speed == 1], na.rm = TRUE))

  # PrEP stop duration trackers
  dat <- set_epi(dat, "prepStopDuration", at, mean(prepStopDuration, na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.B", at, mean(prepStopDuration[race == 1], na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.H", at, mean(prepStopDuration[race == 2], na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.W", at, mean(prepStopDuration[race == 3], na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.O", at, mean(prepStopDuration[race == 4], na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.fast", at, mean(prepStopDuration[stopper.speed == 1], na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.slow", at, mean(prepStopDuration[stopper.speed == 0], na.rm = TRUE))

  dat <- set_epi(dat, "prepStopDuration.slow.B", at, mean(prepStopDuration[race == 1 & stopper.speed == 0], na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.slow.H", at, mean(prepStopDuration[race == 2 & stopper.speed == 0], na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.slow.W", at, mean(prepStopDuration[race == 3 & stopper.speed == 0], na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.slow.O", at, mean(prepStopDuration[race == 4 & stopper.speed == 0], na.rm = TRUE))

  dat <- set_epi(dat, "prepStopDuration.fast.B", at, mean(prepStopDuration[race == 1 & stopper.speed == 1], na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.fast.H", at, mean(prepStopDuration[race == 2 & stopper.speed == 1], na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.fast.W", at, mean(prepStopDuration[race == 3 & stopper.speed == 1], na.rm = TRUE))
  dat <- set_epi(dat, "prepStopDuration.fast.O", at, mean(prepStopDuration[race == 4 & stopper.speed == 1], na.rm = TRUE))

  ## Record attribute history for discontinuation and re-initiation-------------
  #dat <- attr_logger(dat, at)

  return(dat)
}

#' @title Attribute History Tracker Module
#'
#' @description Tracks time-varying attributes over time for PrEP users
#'
#' @inheritParams aging_msm
#'
#' @export
#'

attr_logger <- function(dat, at) {

  prep.start <- get_param(dat, "prep.start")

  # Run after PrEP is introduced
  if (at >= prep.start) {

    # Attributes
    prepClass          <- get_attr(dat, "prepClass")
    prepStat           <- get_attr(dat, "prepStat")
    prepStartTime      <- get_attr(dat, "prepStartTime")
    prepStopTime       <- get_attr(dat, "prepStopTime")
    prep.start.counter <- get_attr(dat, "prep.start.counter")
    prep.stop.counter  <- get_attr(dat, "prep.stop.counter")
    prepDuration       <- get_attr(dat, "prepDuration")
    prepStopDuration   <- get_attr(dat, "prepStopDuration")

    ever.prep <- which(!is.na(prepClass))

    dat <- EpiModel::record_attr_history(
      dat, at,
      "prepStat",
      ever.prep,
      prepStat[ever.prep]
    )

    dat <- EpiModel::record_attr_history(
      dat, at,
      "prepStartTime",
      ever.prep,
      prepStartTime[ever.prep]
    )

    dat <- EpiModel::record_attr_history(
      dat, at,
      "prepStopTime",
      ever.prep,
      prepStopTime[ever.prep]
    )

    dat <- EpiModel::record_attr_history(
      dat, at,
      "prep.start.counter",
      ever.prep,
      prep.start.counter[ever.prep]
    )

    dat <- EpiModel::record_attr_history(
      dat, at,
      "prep.stop.counter",
      ever.prep,
      prep.stop.counter[ever.prep]
    )

    dat <- EpiModel::record_attr_history(
      dat, at,
      "prepDuration",
      ever.prep,
      prepDuration[ever.prep]
    )

    dat <- EpiModel::record_attr_history(
      dat, at,
      "prepStopDuration",
      ever.prep,
      prepStopDuration[ever.prep]
    )

  }

  # Output
  return(dat)
}

#' @title Risk History Sub-Module
#'
#' @description Sub-Module function to track the behavioral risk and STI diagnosis history of
#'              HIV-negative MSM persons for purpose of PrEP indications in [`prep_msm`].
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_msm <- function(dat, at) {

  ## Input
  # Attributes
  active <- get_attr(dat, "active")
  last.neg.test <- get_attr(dat, "last.neg.test")
  dx <- get_attr(dat, "diag.status")

  gc.last.dxpos.time <- get_attr(dat, "gc.last.dxpos.time")
  ct.last.dxpos.time <- get_attr(dat, "ct.last.dxpos.time")
  syph.last.dxpos.time <- get_attr(dat, "syph.last.dxpos.time")
  prep.ind.uai.mono <- get_attr(dat, "prep.ind.uai.mono")
  prep.ind.uai.conc <- get_attr(dat, "prep.ind.uai.conc")


  n <- length(active)
  since.test <- at - last.neg.test

  # Parameters
  time.unit <- get_param(dat, "time.unit")

  ## Edgelist, adds uai summation per partnership from act list
  # nolint start
  pid <- NULL # For R CMD Check
  # nolint end
  uai <- as.vector(rowsum(
    dat[["temp"]][["al"]][, "uai"],
    dat[["temp"]][["al"]][, "pid"]
  ))
  el <- as.data.frame(cbind(dat[["temp"]][["el"]], uai))

  if (max(el[, 1:2]) > n) stop("riskhist max(el) > n")

  # Remove concordant positive edges
  el2 <- el[el[["st2"]] == 0, ]

  # Initialize attributes
  if (is.null(prep.ind.uai.mono)) {
    dat <- set_attr(dat, "prep.ind.uai.mono", rep(NA, n))
    dat <- set_attr(dat, "prep.ind.uai.nmain", rep(NA, n))
    dat <- set_attr(dat, "prep.ind.sti", rep(NA, n))
  }
  if (is.null(prep.ind.uai.conc)) {
    dat <- set_attr(dat, "prep.ind.uai.conc", rep(NA, n))
  }

  ## Degree ##
  main.deg <- get_degree(dat$run$el[[1]])
  casl.deg <- get_degree(dat$run$el[[2]])
  inst.deg <- get_degree(dat$run$el[[3]])


  ## Preconditions ##

  # Any UAI
  uai.any <- unique(c(
    el2[["p1"]][el2[["uai"]] > 0],
    el2[["p2"]][el2[["uai"]] > 0]
  ))

  # Monogamous partnerships: 1-sided
  tot.deg <- main.deg + casl.deg + inst.deg
  uai.mono1 <- intersect(which(tot.deg == 1), uai.any)

  # "Negative" partnerships
  tneg <- unique(c(
    el2[["p1"]][el2[["st1"]] == 0],
    el2[["p2"]][el2[["st1"]] == 0]
  ))
  fneg <- unique(c(
    el2[["p1"]][which(dx[el2[["p1"]]] == 0)],
    el2[["p2"]][which(dx[el2[["p1"]]] == 0)]
  ))
  all.neg <- c(tneg, fneg)

  ## Condition 1b: UAI in 1-sided "monogamous" "negative" partnership,
  ##               partner not tested in past 6 months
  uai.mono1.neg <- intersect(uai.mono1, all.neg)
  part.id1 <- union(
    el2[el2[["p1"]] %in% uai.mono1.neg, 2],
    el2[el2[["p2"]] %in% uai.mono1.neg, 1]
  )
  not.tested.6mo <- since.test[part.id1] > (182 / time.unit)
  not.tested.6mo <- which(not.tested.6mo)
  part.not.tested.6mo <- union(
    el2[el2[["p1"]] %in% part.id1[not.tested.6mo], 2],
    el2[el2[["p2"]] %in% part.id1[not.tested.6mo], 1]
  )

  dat <- set_attr(dat, "prep.ind.uai.mono", at, posit_ids = part.not.tested.6mo)

  ## Condition 2a: UAI + concurrency
  el2.uai <- el2[el2[["uai"]] > 0, ]
  vec <- c(el2.uai[, 1], el2.uai[, 2])
  uai.conc <- unique(vec[duplicated(vec)])
  dat <- set_attr(dat, "prep.ind.uai.conc", at, posit_ids = uai.conc)

  ## Condition 2b: UAI in non-main partnerships
  uai.nmain <- unique(c(
    el2[["p1"]][
      el2[["st1"]] == 0 &
        el2[["uai"]] > 0 &
        el2[["ptype"]] %in% 2:3
    ],
    el2[["p2"]][
      el2[["uai"]] > 0 &
        el2[["ptype"]] %in% 2:3
    ]
  ))
  dat <- set_attr(dat, "prep.ind.uai.nmain", at, posit_ids = uai.nmain)

  ## Condition 4, any STI diagnosis
  last.dx.time <- pmax(gc.last.dxpos.time, ct.last.dxpos.time,
                       syph.last.dxpos.time, na.rm = TRUE)
  dat <- set_attr(dat, "prep.ind.sti", last.dx.time)

  return(dat)
}


#' @title Proportionally Reallocate PrEP Adherence Class Probability
#'
#' @description Shifts probabilities from the high-adherence category to the lower
#'              three adherence categories while maintaining the proportional
#'              distribution of those lower categories.
#'
#' @param in.pcp Input vector of length four for the \code{prep.adhr.dist}
#'        parameters.
#' @param reall The pure percentage points to shift from the high adherence
#'        group to the lower three groups.
#'
#' @export
#'
reallocate_pcp <- function(in.pcp = c(0.089, 0.127, 0.784), reall = 0) {

  dist <- in.pcp[1] / sum(in.pcp[1:2])
  dist[2] <- in.pcp[2] / sum(in.pcp[1:2])

  out.pcp <- rep(NA, 3)
  out.pcp[1:2] <- in.pcp[1:2] - (dist * reall)
  out.pcp[3] <- 1 - sum(out.pcp[1:2])

  return(out.pcp)
}
