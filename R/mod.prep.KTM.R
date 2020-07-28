
#' @title PrEP Module for up to 5 race groups heterosexuals and MSM.
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.  Only msm and msmf are eligible.
#'              If dat$attr$immig.loc=1 PrEP can not be started or re-initiated.
#'
#' @inheritParams aging_shamp
#'
#' @keywords module shamp msm msmf
#'
#' @export
#'
prep_KTM <- function(dat, at) {

  if (at < dat$param$prep.start) {
    return(dat)
  }

  ## Variables
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test

  prepElig <- dat$attr$prepElig
  prep.elig.time <- dat$attr$prep.elig.time
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  prep.elig.model <- dat$param$prep.elig.model
  prep.coverage <- dat$param$prep.coverage
  prep.cov.rate <- dat$param$prep.cov.rate
  prep.class.prob <- dat$param$prep.class.prob
  prep.risk.reassess <- dat$param$prep.risk.reassess
  prep.start.prob <- dat$param$prep.start.prob
  prep.stop.prob <- dat$param$prep.stop.prob
  prep.window <- dat$param$prep.window


  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(prepElig == 1 & prepStat==0 & status == 0 & prep.elig.time > at - prep.window)
  selected <- rbinom(length(idsEligStop),1,prep.start.prob)
  idsStart <-idsEligStart[selected==1] 
  
  idsEligStop <- which(prepStat == 1)
  selected <- rbinom(length(idsEligStop),1,prep.stop.prob)
  idsStop <-idsEligStop[selected==1] 


  ## Stoppage ------------------------------------------------------------------

  # Diagnosis
  idsStpDx <- which(prepStat == 1 & diag.status == 1)


  # Reset PrEP status
  idsStp <- c(idsStpDx, idsStop)
  prepStat[idsStp] <- 0
  prepElig[idsStp] <- 0
  prep.elig.time[idsStp] <- NA

  ## Initiation ----------------------------------------------------------------




  # Attributes
  if (length(idsStart) > 0) {
    prepStat[idsSart] <- 1
    prepElig[idsSart] <- 0
    prep.elig.time[idsSart] <- NA
    
    # PrEP class is fixed over PrEP cycles
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
  }

  ##Timed out
  timed.out<-which(prep.elig.time < at - prep.window)
  prepElig[timed.out] <- 0
  prep.elig.time[timed.out] <- NA

  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prep.elig.time <- prep.elig.time
  dat$attr$prepStat <- prepStat
  dat$attr$prepClass <- prepClass

  # Summary Statistics
  dat$epi$prepCov[at] <- prepCov
  dat$epi$prepStart[at] <- length(idsStart)

  return(dat)
}
