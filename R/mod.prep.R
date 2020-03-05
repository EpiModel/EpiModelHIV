
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

  if (at < dat$param$prep.start) {
    return(dat)
  }

  ## Variables

  # Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test
  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  prepLastRisk <- dat$attr$prepLastRisk
  prepStartTime <- dat$attr$prepStartTime
  prepLastStiScreen <- dat$attr$prepLastStiScreen

  # Parameters

  prep.coverage <- dat$param$prep.coverage
  prep.cov.rate <- dat$param$prep.cov.rate
  prep.class.prob <- dat$param$prep.class.prob


  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 & status == 0 & prepStat == 0 & (lnt == at | dat$attr$last.neg.test.opp == at))

  # Core eligiblity
  ind1 <- dat$attr$prep.ind.uai.mono
  ind2 <- dat$attr$prep.ind.uai.nmain
  ind3 <- dat$attr$prep.ind.ai.sd
  ind4 <- dat$attr$prep.ind.sti

  twind <- at - dat$param$prep.risk.int
  idsEligStart <- intersect(which(ind1 >= twind | ind2 >= twind |
                                    ind3 >= twind | ind4 >= twind),
                            idsEligStart)

  prepElig[idsEligStart] <- 1


  ## Stoppage ------------------------------------------------------------------

  # No indications
  idsRiskAssess <- which(active == 1 & prepStat == 1 & lnt == at & (at - prepLastRisk) >= 52)
  prepLastRisk[idsRiskAssess] <- at

  idsEligStop <- intersect(which(ind1 < twind & ind2 < twind &
                                   ind3 < twind & ind4 < twind),
                           idsRiskAssess)

  prepElig[idsEligStop] <- 0

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1)

  # Reset PrEP status
  idsStp <- c(idsStpDx, idsStpDth, idsEligStop)
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  prepLastStiScreen[idsStp] <- NA


  ## Initiation ----------------------------------------------------------------

  prepCov <- sum(prepStat == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  prepCov <- ifelse(is.nan(prepCov), 0, prepCov)

  idsEligSt <- which(prepElig == 1)
  nEligSt <- length(idsEligSt)

  nStart <- max(0, min(nEligSt, round((prep.coverage - prepCov) *
                                        sum(prepElig == 1, na.rm = TRUE))))
  idsStart <- NULL
  if (nStart > 0) {
    if (prep.cov.rate >= 1) {
      idsStart <- ssample(idsEligSt, nStart)
    } else {
      idsStart <- idsEligSt[rbinom(nStart, 1, prep.cov.rate) == 1]
    }
  }

  # Attributes
  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1
    prepStartTime[idsStart] <- at
    prepLastRisk[idsStart] <- at

    # PrEP class
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
  }


  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepClass <- prepClass
  dat$attr$prepLastRisk <- prepLastRisk
  dat$attr$prepLastStiScreen <- prepLastStiScreen

  # Summary Statistics
  dat$epi$prepCov[at] <- prepCov
  dat$epi$prepStart[at] <- length(idsStart)

  return(dat)
}
