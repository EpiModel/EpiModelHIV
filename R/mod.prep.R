
#' @title PrEP Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
prep_msm <- function(dat, at) {

  # Function Selection ------------------------------------------------------

  if (at >= dat$param$riskh.start) {
    dat <- riskhist_msm(dat, at)
  } else {
    return(dat)
  }

  if (at < dat$param$prep.start) {
    return(dat)
  }


  # Attributes --------------------------------------------------------------

  # Core Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test

  # Oral PrEP Attributes
  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  # LA PrEP Attributes
  prepElig.la <- dat$attr$prepElig.la
  prepStat.la <- dat$attr$prepStat.la
  prepClass.la <- dat$attr$prepClass.la
  prepTimeLastInj <- dat$attr$prepTimeLastInj
  prepLA.dlevel <- dat$attr$prepLA.dlevel
  prepLA.dlevel.int <- dat$attr$prepLA.dlevel.int

  # Shared PrEP Attributes
  prepLastRisk <- dat$attr$prepLastRisk
  prepStartTime <- dat$attr$prepStartTime
  prepLastStiScreen <- dat$attr$prepLastStiScreen


  # Parameters --------------------------------------------------------------

  # Oral PrEP parameters
  prep.start.prob <- dat$param$prep.start.prob
  prep.prob.oral <- dat$param$prep.prob.oral
  prep.adhr.dist <- dat$param$prep.adhr.dist

  # LA PrEP parameters
  prep.adhr.dist.la <- dat$param$prep.adhr.dist.la
  prep.inj.int <- dat$param$prep.inj.int
  icept <- dat$param$prepla.dlevel.icpt
  icept.err <- dat$param$prepla.dlevel.icpt.err
  half.life <- dat$param$prepla.dlevel.halflife.int

  # Shared PrEP parameters
  prep.risk.reassess.method <- dat$param$prep.risk.reassess.method
  prep.discont.rate <- dat$param$prep.discont.rate
  prepla.discont.rate <- dat$param$prepla.discont.rate

  # Indications -------------------------------------------------------------

  ind1 <- dat$attr$prep.ind.uai.mono
  ind2 <- dat$attr$prep.ind.uai.nmain
  ind3 <- dat$attr$prep.ind.sti

  twind <- at - dat$param$prep.risk.int

  # No indications in window
  idsNoIndic <- which((ind1 < twind | is.na(ind1)) &
                      (ind2 < twind | is.na(ind2)) &
                      (ind3 < twind | is.na(ind3)))

  # Indications in window
  idsIndic <- which(ind1 >= twind | ind2 >= twind | ind3 >= twind)

  # Set eligibility to 1 if indications
  prepElig[idsIndic] <- 1
  prepElig.la[idsIndic] <- 1

  # Set eligibility to 0 if no indications
  prepElig[idsNoIndic] <- 0
  prepElig.la[idsNoIndic] <- 0


  ## Stoppage ------------------------------------------------------------------

  # Indication lapse
  # Rules = None, instant, yearly (CDC guidelines)
  if (prep.risk.reassess.method == "none") {
    idsStpInd <- NULL
  } else if (prep.risk.reassess.method == "inst") {
    idsRiskAssess <- which(active == 1 &
                           (prepStat == 1 | prepStat.la == 1))
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  } else if (prep.risk.reassess.method == "year") {
    idsRiskAssess <- which(active == 1 &
                           (prepStat == 1 | prepStat.la == 1) &
                           lnt == at &
                           (at - prepLastRisk) >= 52)
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  }

  # Random discontinuation // oral
  idsEligStpRand.oral <- which(active == 1 & prepStat == 1)
  vecStpRand.oral <- rbinom(length(idsEligStpRand.oral), 1, prep.discont.rate)
  idsStpRand.oral <- idsEligStpRand.oral[which(vecStpRand.oral == 1)]

  # Random discontinuation // injectable
  idsEligStpRand.inj <- which(active == 1 & prepStat.la == 1 & prepClass.la == 1)
  vecStpRand.inj <- rbinom(length(idsEligStpRand.inj), 1, prepla.discont.rate)
  idsStpRand.inj <- idsEligStpRand.inj[which(vecStpRand.inj == 1)]

  # Diagnosis
  idsStpDx <- which(active == 1 & (prepStat == 1 | prepStat.la == 1) & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & (prepStat == 1 | prepStat.la == 1))

  # Reset PrEP status
  idsStp <- c(idsStpInd, idsStpRand.oral, idsStpRand.inj, idsStpDx, idsStpDth)
  idsStp.oral <- intersect(idsStp, which(prepStat == 1))
  idsStp.la <- intersect(idsStp, which(prepStat.la == 1))

  # Update attributes for stoppers
  prepStat[idsStp.oral] <- 0
  # prepElig[idsStp.oral] <- 0
  prepStat.la[idsStp.la] <- 0
  # prepElig.la[idsStp.la] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  prepLastStiScreen[idsStp] <- NA


  ## Initiation ----------------------------------------------------------------

  ## Eligibility ##

  # Indications to start either formulation
  idsEligStart <- which(active == 1 &
                        status == 0 &
                        prepStat == 0 &
                        prepStat.la == 0 &
                        lnt == at)
  idsEligStart <- intersect(idsIndic, idsEligStart)
  prepElig[idsEligStart] <- prepElig.la[idsEligStart] <- 1

  if (!identical(which(prepElig == 1), idsIndic)) browser()

  vecStart <- rbinom(length(idsEligStart), 1, prep.start.prob)
  idsStart <- idsEligStart[which(vecStart == 1)]


  # No LA starts if LA PrEP not yet started
  if (dat$param$prep.la.start > at) {
    prep.prob.oral <- 1
  }
  vecOral <- rbinom(length(idsStart), 1, prep.prob.oral)
  idsStartOral <- idsStart[which(vecOral == 1)]
  idsStartInj <- idsStart[which(vecOral == 0)]

  ## Oral ##

  # Set attributes for oral starters
  if (length(idsStartOral) > 0) {
    prepStat[idsStartOral] <- 1
    prepStartTime[idsStartOral] <- at
    prepLastRisk[idsStartOral] <- at

    # PrEP adherence class
    needPC <- which(is.na(prepClass[idsStartOral]))
    prepClass[idsStartOral[needPC]] <- sample(x = 1:3, size = length(needPC),
                                              replace = TRUE, prob = prep.adhr.dist)
  }


  ## Injectable ##

  # Set attributes for LA starters
  if (length(idsStartInj) > 0) {
    prepStat.la[idsStartInj] <- 1
    prepStartTime[idsStartInj] <- at
    prepLastRisk[idsStartInj] <- at

    # LA PrEP adherence class
    needPC <- which(is.na(prepClass.la[idsStartInj]))
    prepClass.la[idsStartInj[needPC]] <- sample(x = 1:2, size = length(needPC),
                                                replace = TRUE, prob = prep.adhr.dist.la)
  }


  # LA PrEP Interval Injection ----------------------------------------------

  # Current PrEP Starters
  startLAtoday <- which(prepStat.la == 1 & prepStartTime == at)
  prepTimeLastInj[startLAtoday] <- at

  # Weeks since last injection
  last.inj <- at - prepTimeLastInj

  # Needing injection
  idsLA <- which(prepStat.la == 1)
  getInjection <- intersect(idsLA, which(last.inj >= prep.inj.int))

  # Update injection time attributes
  prepTimeLastInj[getInjection] <- at
  last.inj <- at - prepTimeLastInj


  # LA PrEP Drug Levels -----------------------------------------------------

  # Set drug level intercept for those newly started
  prepLA.dlevel.int[startLAtoday] <- pmax(0.1,
                                         rnorm(length(startLAtoday),
                                               icept, icept.err))

  # Update dlevel for all active users
  prepLA.dlevel <- prepLA.dlevel.int * 0.5^(last.inj/half.life)


  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepClass <- prepClass

  dat$attr$prepElig.la <- prepElig.la
  dat$attr$prepStat.la <- prepStat.la
  dat$attr$prepClass.la <- prepClass.la

  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepLastRisk <- prepLastRisk
  dat$attr$prepLastStiScreen <- prepLastStiScreen

  dat$attr$prepTimeLastInj <- prepTimeLastInj
  dat$attr$prepLA.dlevel <- prepLA.dlevel
  dat$attr$prepLA.dlevel.int <- prepLA.dlevel.int

  return(dat)
}


#' @title Risk History Sub-Module
#'
#' @description Sub-Module function to track the risk history of uninfected persons
#'              for purpose of PrEP targeting.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_msm <- function(dat, at) {

  ## Attributes
  n <- length(dat$attr$active)
  dx <- dat$attr$diag.status
  since.test <- at - dat$attr$last.neg.test
  rGC.tx <- dat$attr$rGC.tx
  uGC.tx <- dat$attr$uGC.tx
  rCT.tx <- dat$attr$rCT.tx
  uCT.tx <- dat$attr$uCT.tx

  ## Parameters
  time.unit <- dat$param$time.unit

  ## Edgelist, adds uai summation per partnership from act list
  pid <- NULL # For R CMD Check
  al <- as.data.frame(dat$temp$al)
  by_pid <- group_by(al, pid)
  uai <- summarise(by_pid, uai = sum(uai))[, 2]
  el <- as.data.frame(cbind(dat$temp$el, uai))

  if (max(el[, 1:2]) > n) stop("riskhist max(el) > n")

  # Remove concordant positive edges
  el2 <- el[el$st2 == 0, ]

  # Initialize attributes
  if (is.null(dat$attr$prep.ind.uai.mono)) {
    dat$attr$prep.ind.uai.mono <- rep(NA, n)
    dat$attr$prep.ind.uai.nmain <- rep(NA, n)
    dat$attr$prep.ind.sti <- rep(NA, n)
  }

  ## Degree ##
  main.deg <- get_degree(dat$el[[1]])
  casl.deg <- get_degree(dat$el[[2]])
  inst.deg <- get_degree(dat$el[[3]])


  ## Preconditions ##

  # Any UAI
  uai.any <- unique(c(el2$p1[el2$uai > 0],
                      el2$p2[el2$uai > 0]))

  # Monogamous partnerships: 1-sided
  tot.deg <- main.deg + casl.deg + inst.deg
  uai.mono1 <- intersect(which(tot.deg == 1), uai.any)

  # "Negative" partnerships
  tneg <- unique(c(el2$p1[el2$st1 == 0], el2$p2[el2$st1 == 0]))
  fneg <- unique(c(el2$p1[which(dx[el2$p1] == 0)], el2$p2[which(dx[el2$p1] == 0)]))
  all.neg <- c(tneg, fneg)

  ## Condition 1b: UAI in 1-sided "monogamous" "negative" partnership,
  ##               partner not tested in past 6 months
  uai.mono1.neg <- intersect(uai.mono1, all.neg)
  part.id1 <- c(el2[el2$p1 %in% uai.mono1.neg, 2], el2[el2$p2 %in% uai.mono1.neg, 1])
  not.tested.6mo <- since.test[part.id1] > (180/time.unit)
  part.not.tested.6mo <- uai.mono1.neg[which(not.tested.6mo == TRUE)]
  dat$attr$prep.ind.uai.mono[part.not.tested.6mo] <- at

  ## Condition 2b: UAI in non-main partnerships
  uai.nmain <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$ptype %in% 2:3],
                        el2$p2[el2$uai > 0 & el2$ptype %in% 2:3]))
  dat$attr$prep.ind.uai.nmain[uai.nmain] <- at

  ## Condition 4, any STI diagnosis
  idsDx <- which(rGC.tx == 1 | uGC.tx == 1 | rCT.tx == 1 | uCT.tx == 1)
  dat$attr$prep.ind.sti[idsDx] <- at

  return(dat)
}
