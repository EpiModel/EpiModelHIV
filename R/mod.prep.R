
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

  # Pull Data ---------------------------------------------------------------

  # Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  race <- dat$attr$race
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test

  prepAccess <- dat$attr$prepAccess
  prepIndic <- dat$attr$prepIndic
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  prepLastRisk <- dat$attr$prepLastRisk
  prepStartTime <- dat$attr$prepStartTime
  prepLastStiScreen <- dat$attr$prepLastStiScreen

  # Parameters
  prep.risk.reassess.method <- dat$param$prep.risk.reassess.method
  prep.rx.B <- dat$param$prep.rx.B
  prep.rx.W <- dat$param$prep.rx.W
  prep.adhr.dist.B <- dat$param$prep.adhr.dist.B
  prep.adhr.dist.W <- dat$param$prep.adhr.dist.W
  prep.discont.rate.B <- dat$param$prep.discont.rate.B
  prep.discont.rate.W <- dat$param$prep.discont.rate.W


  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 &
                        status == 0 &
                        prepStat == 0 &
                        lnt == at &
                        prepAccess == 1)

  # Core eligiblity
  ind1 <- dat$attr$prep.ind.uai.mono
  ind2 <- dat$attr$prep.ind.uai.nmain
  ind3 <- dat$attr$prep.ind.ai.sd
  ind4 <- dat$attr$prep.ind.sti

  twind <- at - dat$param$prep.risk.int
  idsIndic <- which(ind1 >= twind | ind2 >= twind | ind3 >= twind | ind4 >= twind)

  idsEligStart <- intersect(idsIndic, idsEligStart)

  prepIndic[idsIndic] <- 1

  ## Stoppage ------------------------------------------------------------------

  # No indications
  idsNoIndic <- which(ind1 < twind & ind2 < twind & ind3 < twind & ind4 < twind)
  prepIndic[idsNoIndic] <- 0

  # Risk reassessment rules
  if (prep.risk.reassess.method == "none") {
    idsStpInd <- NULL
  } else if (prep.risk.reassess.method == "inst") {
    idsRiskAssess <- which(active == 1 & prepStat == 1)
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  } else if (prep.risk.reassess.method == "year") {
    idsRiskAssess <- which(active == 1 & prepStat == 1 & lnt == at & (at - prepLastRisk) >= 52)
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  }

  # Random (memoryless) discontinuation
  idsEligStpRand.B <- which(active == 1 & prepStat == 1 & race == "B")
  idsEligStpRand.W <- which(active == 1 & prepStat == 1 & race == "W")
  vecStpRand.B <- rbinom(length(idsEligStpRand.B), 1, prep.discont.rate.B)
  vecStpRand.W <- rbinom(length(idsEligStpRand.W), 1, prep.discont.rate.W)
  idsStpRand.B <- idsEligStpRand.B[which(vecStpRand.B == 1)]
  idsStpRand.W <- idsEligStpRand.W[which(vecStpRand.W == 1)]
  idsStpRand <- union(idsStpRand.B, idsStpRand.W)

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1)

  # Reset PrEP status
  idsStp <- c(idsStpInd, idsStpRand, idsStpDx, idsStpDth)
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  prepLastStiScreen[idsStp] <- NA


  ## Initiation ----------------------------------------------------------------

  idsEligSt.B <- intersect(idsEligStart, which(race == "B"))
  idsEligSt.W <- intersect(idsEligStart, which(race == "W"))

  prepStat[idsEligSt.B] <- rbinom(length(idsEligSt.B), 1, prep.rx.B)
  prepStat[idsEligSt.W] <- rbinom(length(idsEligSt.W), 1, prep.rx.W)

  idsStart.B <- intersect(idsEligSt.B, which(prepStat == 1))
  idsStart.W <- intersect(idsEligSt.W, which(prepStat == 1))
  idsStart <- union(idsStart.B, idsStart.W)

  # Attributes
  if (length(idsStart) > 0) {
    prepStartTime[idsStart] <- at
    prepLastRisk[idsStart] <- at

    # PrEP adherence class
    needPC.B <- which(is.na(prepClass[idsStart.B]))
    prepClass[idsStart.B[needPC.B]] <- sample(x = 1:3, size = length(needPC.B),
                                              replace = TRUE, prob = prep.adhr.dist.B)
    needPC.W <- which(is.na(prepClass[idsStart.W]))
    prepClass[idsStart.W[needPC.W]] <- sample(x = 1:3, size = length(needPC.W),
                                              replace = TRUE, prob = prep.adhr.dist.W)
  }


  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepIndic <- prepIndic
  dat$attr$prepStat <- prepStat
  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepClass <- prepClass
  dat$attr$prepLastRisk <- prepLastRisk
  dat$attr$prepLastStiScreen <- prepLastStiScreen

  # Summary stats
  if (is.null(dat$epi$prepStpInd)) {
    dat$epi$prepStpInd <- rep(NA, dat$control$nsteps)
    dat$epi$prepStpDx <- rep(NA, dat$control$nsteps)
    dat$epi$prepStpRand <- rep(NA, dat$control$nsteps)
  }
  dat$epi$prepStpInd[at] <- length(idsStpInd)
  dat$epi$prepStpDx[at] <- length(idsStpDx)
  dat$epi$prepStpRand[at] <- length(idsStpRand)

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
  uid <- dat$attr$uid
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
    dat$attr$prep.ind.ai.sd <- rep(NA, n)
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

  ## Condition 3a: AI within known serodiscordant partnerships
  el2.cond3 <- el2[el2$st1 == 1 & el2$ptype %in% 1:2, ]

  # Disclosure
  discl.list <- dat$temp$discl.list
  disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
  delt.cdl <- uid[el2.cond3[, 1]] * 1e7 + uid[el2.cond3[, 2]]
  discl <- (delt.cdl %in% disclose.cdl)
  ai.sd <- el2.cond3$p2[discl == TRUE]
  dat$attr$prep.ind.ai.sd[ai.sd] <- at

  ## Condition 4, any STI diagnosis
  idsDx <- which(rGC.tx == 1 | uGC.tx == 1 | rCT.tx == 1 | uCT.tx == 1)
  dat$attr$prep.ind.sti[idsDx] <- at

  return(dat)
}
