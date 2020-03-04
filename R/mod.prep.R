
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

  # PrEP Attributes
  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  prepLastRisk <- dat$attr$prepLastRisk
  prepStartTime <- dat$attr$prepStartTime
  prepLastStiScreen <- dat$attr$prepLastStiScreen


  # Parameters --------------------------------------------------------------

  prep.start.prob <- dat$param$prep.start.prob
  prep.adhr.dist <- dat$param$prep.adhr.dist
  prep.require.lnt <- dat$param$prep.require.lnt
  prep.risk.reassess.method <- dat$param$prep.risk.reassess.method
  prep.discont.rate <- dat$param$prep.discont.rate


  # Indications -------------------------------------------------------------

  ind1 <- dat$attr$prep.ind.uai.mono
  # ind2 <- dat$attr$prep.ind.uai.nmain
  ind2 <- dat$attr$prep.ind.uai.conc
  ind3 <- dat$attr$prep.ind.sti

  twind <- at - dat$param$prep.risk.int

  # No indications in window
  idsNoIndic <- which((ind1 < twind | is.na(ind1)) &
                      (ind2 < twind | is.na(ind2)) &
                      (ind3 < twind | is.na(ind3)))
  base.cond.no <- which(active == 0 | diag.status == 1)
  idsNoIndic <- union(idsNoIndic, base.cond.no)

  # Indications in window
  idsIndic <- which(ind1 >= twind | ind2 >= twind | ind3 >= twind)
  base.cond.yes <- which(active == 1 & diag.status == 0)
  idsIndic <- intersect(idsIndic, base.cond.yes)

  # Set eligibility to 1 if indications
  prepElig[idsIndic] <- 1

  # Set eligibility to 0 if no indications
  prepElig[idsNoIndic] <- 0


  ## Stoppage ------------------------------------------------------------------

  # Indication lapse
  # Rules = None, instant, yearly (CDC guidelines)
  if (prep.risk.reassess.method == "none") {
    idsStpInd <- NULL
  } else if (prep.risk.reassess.method == "inst") {
    idsRiskAssess <- which(active == 1 & prepStat == 1)
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  } else if (prep.risk.reassess.method == "year") {
    idsRiskAssess <- which(active == 1 &
                           prepStat == 1 &
                           lnt == at &
                           (at - prepLastRisk) >= 52)
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  }

  # Random discontinuation
  idsEligStpRand <- which(active == 1 & prepStat == 1)
  vecStpRand <- rbinom(length(idsEligStpRand), 1, prep.discont.rate)
  idsStpRand <- idsEligStpRand[which(vecStpRand == 1)]

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1)

  # Reset PrEP status
  idsStp <- c(idsStpInd, idsStpRand, idsStpDx, idsStpDth)

  # Update attributes for stoppers
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  prepLastStiScreen[idsStp] <- NA


  ## Initiation ----------------------------------------------------------------

  ## Eligibility ##

  # Indications to start
  if (prep.require.lnt == TRUE) {
    idsEligStart <- which(prepStat == 0 & lnt == at)
  } else {
    idsEligStart <- which(prepStat == 0)
  }

  idsEligStart <- intersect(idsIndic, idsEligStart)
  prepElig[idsEligStart] <- 1

  vecStart <- rbinom(length(idsEligStart), 1, prep.start.prob)
  idsStart <- idsEligStart[which(vecStart == 1)]

  # Set attributes for starters
  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1
    prepStartTime[idsStart] <- at
    prepLastRisk[idsStart] <- at

    # PrEP adherence class
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = 1:3, size = length(needPC),
                                          replace = TRUE, prob = prep.adhr.dist)
  }


  ## Output --------------------------------------------------------------------

  # Random discontinuation
  dat$epi$prep.rand.stop[at] <- length(idsStpRand)

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepClass <- prepClass

  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepLastRisk <- prepLastRisk
  dat$attr$prepLastStiScreen <- prepLastStiScreen

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
  if (is.null(dat$attr$prep.ind.uai.conc)) {
    dat$attr$prep.ind.uai.conc <- rep(NA, n)
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
  not.tested.6mo <- since.test[part.id1] > (180/7)
  part.not.tested.6mo <- uai.mono1.neg[which(not.tested.6mo == TRUE)]
  dat$attr$prep.ind.uai.mono[part.not.tested.6mo] <- at

  ## Condition 2a: UAI + concurrency
  el2.uai <- el2[el2$uai > 0, ]
  vec <- c(el2.uai[, 1], el2.uai[, 2])
  uai.conc <- unique(vec[duplicated(vec)])
  dat$attr$prep.ind.uai.conc[uai.conc] <- at

  ## Condition 2b: UAI in non-main partnerships
  uai.nmain <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$ptype %in% 2:3],
                        el2$p2[el2$uai > 0 & el2$ptype %in% 2:3]))
  dat$attr$prep.ind.uai.nmain[uai.nmain] <- at

  ## Condition 4, any STI diagnosis
  idsDx <- which(rGC.tx == 1 | uGC.tx == 1 | rCT.tx == 1 | uCT.tx == 1)
  dat$attr$prep.ind.sti[idsDx] <- at

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

  dist <- in.pcp[1]/sum(in.pcp[1:2])
  dist[2] <- in.pcp[2]/sum(in.pcp[1:2])

  out.pcp <- rep(NA, 3)
  out.pcp[1:2] <- in.pcp[1:2] - (dist * reall)
  out.pcp[3] <- 1 - sum(out.pcp[1:2])

  return(out.pcp)
}
