
#' @title PrEP MSM Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection for MSM.
#'
#' @inheritParams aging_camplc
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
  active <- dat$attr$active
  status <- dat$attr$status
  asmm <- dat$attr$asmm
  ever.adol.prep <- dat$attr$ever.adol.prep
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test

  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepEver <- dat$attr$prepEver
  prepClass <- dat$attr$prepClass

  prep.elig.model <- dat$param$prep.elig.model
  prep.coverage.adol.naive <- dat$param$prep.coverage.adol.naive
  prep.coverage.adol.exp <- dat$param$prep.coverage.adol.exp
  prep.cov.method <- dat$param$prep.cov.method
  prep.cov.rate <- dat$param$prep.cov.rate
  prep.class.prob <- dat$param$prep.class.prob
  prep.risk.reassess <- dat$param$prep.risk.reassess


  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart.adol.naive <- which(active == 1 & status == 0 & prepStat == 0 & asmm == 0 & ever.adol.prep == 0 & lnt == at)
  idsEligStart.adol.exp <- which(active == 1 & status == 0 & prepStat == 0 & asmm == 0 & ever.adol.prep == 1 & lnt == at)

  idsEligStop <- NULL
  if (prep.risk.reassess == TRUE) {
    idsEligStop <- which(active == 1 & prepStat == 1 & asmm == 0 & lnt == at)
  }

  # Core eligiblity scenarios
  if (prep.elig.model != "base") {
    if (substr(prep.elig.model, 1, 3) == "cdc") {
      if (prep.elig.model == "cdc1") {
        mat.c1 <- dat$riskh$uai.mono2.nt.6mo
        mat.c2 <- dat$riskh$uai.nonmonog
        mat.c3 <- dat$riskh$ai.sd.mc
      } else if (prep.elig.model == "cdc2") {
        mat.c1 <- dat$riskh$uai.mono2.nt.6mo
        mat.c2 <- dat$riskh$uai.nmain
        mat.c3 <- dat$riskh$ai.sd.mc
      } else if (prep.elig.model == "cdc3") {
        mat.c1 <- dat$riskh$uai.mono1.nt.6mo
        mat.c2 <- dat$riskh$uai.nmain
        mat.c3 <- dat$riskh$ai.sd.mc
      } else if (prep.elig.model == "cdc4") {
        mat.c1 <- dat$riskh$uai.mono1.nt.6mo
        mat.c2 <- dat$riskh$uai.nmain
        mat.c3 <- dat$riskh$uai.sd.mc
      }
      idsEligStart.adol.naive <- intersect(which(rowSums(mat.c1, na.rm = TRUE) > 0 |
                                        rowSums(mat.c2, na.rm = TRUE) > 0 |
                                        rowSums(mat.c3, na.rm = TRUE) > 0),
                                idsEligStart.adol.naive)
      idsEligStart.adol.exp <- intersect(which(rowSums(mat.c1, na.rm = TRUE) > 0 |
                                                   rowSums(mat.c2, na.rm = TRUE) > 0 |
                                                   rowSums(mat.c3, na.rm = TRUE) > 0),
                                           idsEligStart.adol.exp)
      
      idsEligStop <- intersect(which(rowSums(mat.c1, na.rm = TRUE) == 0 &
                                       rowSums(mat.c2, na.rm = TRUE) == 0 &
                                       rowSums(mat.c3, na.rm = TRUE) == 0),
                               idsEligStop)
    } else {
      mat <- dat$riskh[[prep.elig.model]]
      idsEligStart.adol.naive <- intersect(which(rowSums(mat, na.rm = TRUE) > 0), idsEligStart.adol.naive)
      idsEligStart.adol.exp <- intersect(which(rowSums(mat, na.rm = TRUE) > 0), idsEligStart.adol.exp)
      idsEligStop <- intersect(which(rowSums(mat, na.rm = TRUE) == 0), idsEligStop)
    }
  }

  prepElig.adol.naive[idsEligStart.adol.naive] <- 1
  prepElig.adol.exp[idsEligStart.adol.exp] <- 1
  prepElig[idsEligStart.adol.naive] <- 1
  prepElig[idsEligStart.adol.exp] <- 1
  prepElig[idsEligStop] <- 0


  ## Stoppage ------------------------------------------------------------------

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & asmm == 0 & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1 & asmm == 0)

  # Transition to ineligibility
  idsStpInelig <- idsEligStop

  # Reset PrEP status
  idsStp <- c(idsStpDx, idsStpDth, idsStpInelig)
  prepStat[idsStp] <- 0


  ## Initiation ----------------------------------------------------------------

  if (prep.cov.method == "curr") {
    prepCov.adol.naive <- sum(prepStat == 1 & ever.adol.prep == 0, na.rm = TRUE)/sum(prepElig == 1 & ever.adol.prep == 0, na.rm = TRUE)
    prepCov.adol.exp <- sum(prepStat == 1 & ever.adol.prep == 1, na.rm = TRUE)/sum(prepElig == 1 & ever.adol.prep == 1, na.rm = TRUE)
  }
  if (prep.cov.method == "ever") {
    prepCov.adol.naive <- sum(prepEver == 1 & ever.adol.prep == 0, na.rm = TRUE)/sum(prepElig == 1 & ever.adol.prep == 0, na.rm = TRUE)
    prepCov.adol.exp <- sum(prepEver == 1 & ever.adol.prep == 1, na.rm = TRUE)/sum(prepElig == 1 & ever.adol.prep == 1, na.rm = TRUE)
  }
  prepCov.adol.naive <- ifelse(is.nan(prepCov.adol.naive), 0, prepCov.adol.naive)
  prepCov.adol.exp <- ifelse(is.nan(prepCov.adol.exp), 0, prepCov.adol.exp)
  
  idsEligSt.adol.naive <- which(prepElig.adol.naive == 1)
  idsEligSt.adol.exp <- which(prepElig.adol.exp == 1)
    nEligSt.adol.naive <- length(idsEligSt.adol.naive)
    nEligSt.adol.exp <- length(idsEligSt.adol.exp)

  nStart.adol.naive <- max(0, min(nEligSt.adol.naive, round((prep.coverage.adol.naive - prepCov.adol.naive) *
                                        sum(prepElig == 1 & ever.adol.prep == 0, na.rm = TRUE))))
  
  nStart.adol.exp <- max(0, min(nEligSt.adol.exp, round((prep.coverage.adol.exp - prepCov.adol.exp) *
                                                              sum(prepElig == 1 & ever.adol.prep == 1, na.rm = TRUE))))
  

  
  idsStart.adol.naive <- NULL
  if (nStart.adol.naive > 0) {
    if (prep.cov.rate.adol.naive >= 1) {
      idsStart.adol.naive <- ssample(idsEligSt.adol.naive, nStart.adol.naive)
    } else {
      idsStart.adol.naive <- idsEligSt.adol.naive[rbinom(nStart.adol.naive, 1, prep.cov.rate.adol.naive) == 1]
    }
  }

  idsStart.adol.exp <- NULL
  if (nStart.adol.exp > 0) {
    if (prep.cov.rate.adol.exp >= 1) {
      idsStart.adol.exp <- ssample(idsEligSt.adol.exp, nStart.adol.exp)
    } else {
      idsStart.adol.exp <- idsEligSt.adol.exp[rbinom(nStart.adol.exp, 1, prep.cov.rate.adol.exp) == 1]
    }
  }
  
  # Attributes
  if (length(idsStart.adol.naive) > 0) {
    prepStat[idsStart.adol.naive] <- 1
    prepEver[idsStart.adol.naive] <- 1

    # PrEP class is fixed over PrEP cycles
    needPC <- which(is.na(prepClass[idsStart.adol.naive]))
    prepClass[idsStart.adol.naive[needPC]] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
  }
  
  if (length(idsStart.adol.exp) > 0) {
    prepStat[idsStart.adol.exp] <- 1
    prepEver[idsStart.adol.exp] <- 1
    
    # PrEP class is fixed over PrEP cycles
    needPC <- which(is.na(prepClass[idsStart.adol.exp]))
    prepClass[idsStart.adol.exp[needPC]] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
  }
  


  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepElig.adol.naive <- prepElig.adol.naive
  dat$attr$prepElig.adol.exp <- prepElig.adol.exp
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepEver <- prepEver
  dat$attr$prepClass <- prepClass

  # Summary Statistics
  dat$epi$prepCov.adol.naive[at] <- prepCov.adol.naive
  dat$epi$prepCov.adol.exp[at] <- prepCov.adol.exp
  dat$epi$prepStart[at] <- length(idsStart.adol.naive) + length(idsStart.adol.exp)

  return(dat)
}
