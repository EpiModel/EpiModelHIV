
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
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test

  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  prep.elig.model <- dat$param$prep.elig.model
  prep.coverage <- dat$param$prep.coverage
  prep.cov.rate <- dat$param$prep.cov.rate
  prep.class.prob <- dat$param$prep.class.prob
  prep.risk.reassess <- dat$param$prep.risk.reassess


  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(status == 0 & prepStat == 0 & lnt == at)

  idsEligStop <- NULL
  if (prep.risk.reassess == TRUE) {
    idsEligStop <- which(prepStat == 1 & lnt == at)
  }

  twind <- at - dat$param$prep.risk.int

  # Core eligiblity scenarios
  if (prep.elig.model != "base") {
    if (substr(prep.elig.model, 1, 3) == "cdc") {
      if (prep.elig.model == "cdc1") {
        c1 <- dat$attr$uai.mono2.nt.6mo
        c2 <- dat$attr$uai.nonmonog
        c3 <- dat$attr$ai.sd.mc
      } else if (prep.elig.model == "cdc2") {
        c1 <- dat$attr$uai.mono2.nt.6mo
        c2 <- dat$attr$uai.nmain
        c3 <- dat$attr$ai.sd.mc
      } else if (prep.elig.model == "cdc3") {
        c1 <- dat$attr$uai.mono1.nt.6mo
        c2 <- dat$attr$uai.nmain
        c3 <- dat$attr$ai.sd.mc
      } else if (prep.elig.model == "cdc4") {
        c1 <- dat$attr$uai.mono1.nt.6mo
        c2 <- dat$attr$uai.nmain
        c3 <- dat$attr$uai.sd.mc
      }
      idsEligStart <- intersect(which(c1 >= twind |
                                      c2 >= twind |
                                      c3 >= twind),
                                idsEligStart)
      idsEligStop <- intersect(which(c1 < twind &
                                     c2 < twind &
                                     c3 < twind),
                                idsEligStop)
    } else {
      c1 <- dat$attr[[prep.elig.model]]
      idsEligStart <- intersect(which(c1 >= twind), idsEligStart)
      idsEligStop <- intersect(which(c1 < twind), idsEligStop)
    }
  }

  prepElig[idsEligStart] <- 1
  prepElig[idsEligStop] <- 0


  ## Stoppage ------------------------------------------------------------------

  # Diagnosis
  idsStpDx <- which(prepStat == 1 & diag.status == 1)

  # Transition to ineligibility
  idsStpInelig <- idsEligStop

  # Reset PrEP status
  idsStp <- c(idsStpDx, idsStpInelig)
  prepStat[idsStp] <- 0


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

    # PrEP class is fixed over PrEP cycles
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
  }


  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepClass <- prepClass

  # Summary Statistics
  dat$epi$prepCov[at] <- prepCov
  dat$epi$prepStart[at] <- length(idsStart)

  return(dat)
}
