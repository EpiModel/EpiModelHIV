
#' @title PrEP Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging.mard
#'
#' @export
#'
prep.mard <- function(dat, at) {

  if (at < dat$param$prep.start) {
    return(dat)
  }

  ## Variables
  active <- dat$attr$active
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test

  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepEver <- dat$attr$prepEver
  prepClass <- dat$attr$prepClass

  prep.elig.model <- dat$param$prep.elig.model
  prep.coverage <- dat$param$prep.coverage
  prep.cov.method <- dat$param$prep.cov.method
  prep.cov.rate <- dat$param$prep.cov.rate
  prep.class.prob <- dat$param$prep.class.prob
  prep.risk.reassess <- dat$param$prep.risk.reassess

# browser()
  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 & status == 0 & prepStat == 0 & lnt == at)

  idsEligStop <- NULL
  if (prep.risk.reassess == TRUE) {
    idsEligStop <- which(active == 1 & prepStat == 1 & lnt == at)
  }

  # Core eligiblity scenarios
  if (prep.elig.model != "base") {
    mat <- dat$riskh[[prep.elig.model]]
    idsEligStart <- intersect(which(rowSums(mat) > 0), idsEligStart)
    idsEligStop <- intersect(which(rowSums(mat) == 0), idsEligStop)
  }

  prepElig[idsEligStart] <- 1
  prepElig[idsEligStop] <- 0


  ## Stoppage ------------------------------------------------------------------

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1)

  # Transition to ineligibility
  idsStpInelig <- idsEligStop

  # Reset PrEP status
  idsStp <- c(idsStpDx, idsStpDth, idsStpInelig)
  prepStat[idsStp] <- 0


  ## Initiation ----------------------------------------------------------------

  if (prep.cov.method == "curr") {
    prepCov <- sum(prepStat == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  }
  if (prep.cov.method == "ever") {
    prepCov <- sum(prepEver == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  }
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
    prepEver[idsStart] <- 1
    prepClass[idsStart] <- sample(c("l", "m", "h"), length(idsStart),
                                  TRUE, prep.class.prob)
  }


  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepEver <- prepEver
  dat$attr$prepClass <- prepClass

  # Summary Statistics
  dat$epi$prepCov[at] <- prepCov
  dat$epi$prepElig[at] <- length(idsEligStart)
  dat$epi$prepStart[at] <- length(idsStart)

  return(dat)
}
