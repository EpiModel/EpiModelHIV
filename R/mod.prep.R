
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
  ever.adult.prep  <- dat$attr$ever.adult.prep
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test
  
  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepStat.asmm <- dat$attr$prepStat.asmm
  prepEver <- dat$attr$prepEver
  prepClass <- dat$attr$prepClass

  prep.elig.model <- dat$param$prep.elig.model
  prep.coverage <- dat$param$prep.coverage
  prep.cov.method <- dat$param$prep.cov.method
  prep.cov.rate <- dat$param$prep.cov.rate
  prep.class.prob <- dat$param$prep.class.prob
  prep.risk.reassess <- dat$param$prep.risk.reassess
  prep.retained <- dat$attr$prep.retained

  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 & status == 0 & prepStat == 0 & asmm == 0 & lnt == at)

  idsEligStop <- NULL
  if (prep.risk.reassess == TRUE) {
    idsEligStop <- which(active == 1 & prepStat == 1 & prepStat.asmm == 0 & asmm == 0 & lnt == at)
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
      idsEligStart <- intersect(which(rowSums(mat.c1, na.rm = TRUE) > 0 |
                                        rowSums(mat.c2, na.rm = TRUE) > 0 |
                                        rowSums(mat.c3, na.rm = TRUE) > 0),
                                idsEligStart)

      
      idsEligStop <- intersect(which(rowSums(mat.c1, na.rm = TRUE) == 0 &
                                       rowSums(mat.c2, na.rm = TRUE) == 0 &
                                       rowSums(mat.c3, na.rm = TRUE) == 0),
                               idsEligStop)
    } else {
      mat <- dat$riskh[[prep.elig.model]]
      idsEligStart <- intersect(which(rowSums(mat, na.rm = TRUE) > 0), idsEligStart)
      idsEligStop <- intersect(which(rowSums(mat, na.rm = TRUE) == 0), idsEligStop)
    }
  }

  prepElig[idsEligStart] <- 1
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

  ##COV AMONG ELIGIBLE AND TOTAL WITH AND WITHOUT RETAINED
  

    prepCov.msm.elig <- (sum(prepStat == 1 & prepStat.asmm == 0 & asmm == 0, na.rm = TRUE)) / (sum(prepElig == 1 & asmm == 0, na.rm = TRUE))
    prepCov.msm.elig.w.ret <- (sum(prepStat == 1 & asmm == 0, na.rm = TRUE)) / (sum(prepElig == 1 & asmm == 0, na.rm = TRUE))
    prepCov.msm.all<- (sum(prepStat == 1 & prepStat.asmm == 0 & asmm == 0, na.rm = TRUE)) / (sum(asmm == 0, na.rm = TRUE))
    prepCov.msm.all.w.ret<- (sum(prepStat == 1 & asmm == 0, na.rm = TRUE)) / (sum(asmm == 0, na.rm = TRUE))

  

  prepCov.msm.elig <- ifelse(is.nan(prepCov.msm.elig), 0, prepCov.msm.elig)
  prepCov.msm.elig.w.ret <- ifelse(is.nan(prepCov.msm.elig.w.ret), 0, prepCov.msm.elig.w.ret)
  prepCov.msm.all <- ifelse(is.nan(prepCov.msm.all), 0, prepCov.msm.all)
  prepCov.msm.all.w.ret <- ifelse(is.nan(prepCov.msm.all.w.ret), 0, prepCov.msm.all.w.ret)

  
  idsEligSt <- which(prepElig == 1)
    nEligSt <- length(idsEligSt)


  nStart <- max(0, min(nEligSt, round((prep.coverage - prepCov.msm.elig) *
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
    ever.adult.prep[idsStart] <- 1

    # PrEP class is fixed over PrEP cycles
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
  }
  



  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepEver <- prepEver
  dat$attr$prepClass <- prepClass
  dat$attr$ever.adult.prep <- ever.adult.prep

  x <- which(dat$attr$prepStat == 1)
  y <- which(dat$attr$prepStat != 1)
  dat$attr$prep.time[x] <- dat$attr$prep.time[x] + 1
  dat$attr$prep.time[y] <- 0
  
  # Summary Statistics

  dat$epi$prepCov.msm.elig[at] <- prepCov.msm.elig 
  dat$epi$prepCov.msm.elig.w.ret[at] <- prepCov.msm.elig.w.ret 
  dat$epi$prepCov.msm.all[at] <- prepCov.msm.all 
  dat$epi$prepCov.msm.all.w.ret[at] <- prepCov.msm.all.w.ret
  
  dat$epi$prepStart.msm[at] <- length(idsStart)
  
    return(dat)
}
