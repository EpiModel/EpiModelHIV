
#' @title PrEP ASMM Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) among ASMM to prevent HIV infection.
#'
#' @inheritParams aging_camplc
#'
#' @export
#'
prep_adol <- function(dat, at) {

  if (at < dat$param$prep.start.asmm) {
    return(dat)
  }
  

  ## Variables
  ##Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  race<- dat$attr$race
  asmm <- dat$attr$asmm


  diag.status <- dat$attr$diag.status
  debuted<-dat$attr$debuted
  debuted.time<-dat$attr$has.debuted.time
  active.time<-dat$attr$active.time
  AI.time<-dat$attr$AI.time
  lnt <- dat$attr$last.neg.test
  everAI<-dat$attr$everAI
  of.age<-dat$attr$of.age
  uaicount<-dat$attr$uaicount
 
  prepElig.asmm <- dat$attr$prepElig.asmm
  prepStat <- dat$attr$prepStat
  prepEver <- dat$attr$prepEver
  prepClass <- dat$attr$prepClass
  prepStart.time <- dat$attr$prepStart.time
  ever.adol.prep <- dat$attr$ever.adol.prep
  ever.adult.prep <- dat$attr$ever.adult.prep
  

##PrEP params.
  #Uniform.
  prep.elig.model <- dat$param$prep.elig.model.asmm
  prep.cov.method <- dat$param$prep.cov.method.asmm
  prep.risk.reassess <- dat$param$prep.risk.reassess.asmm
  prepSpell<- dat$param$prepSpell.asmm
  
  prep.coverage <- dat$param$prep.coverage.asmm
  prep.cov.rate <- dat$param$prep.cov.rate.asmm
  prep.class.prob <- dat$param$prep.class.prob.asmm
  prepDrop<- dat$param$prepDrop.asmm
  prep.delay<-dat$param$prep.delay.asmm
  prep.uaicount.thresh<-dat$param$prep.uaicount.thresh.asmm
  


  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 & status == 0 & prepStat == 0 & asmm == 1)


  idsEligStop <- NULL
  if (prep.risk.reassess == TRUE) {
    idsEligStop <- which(active == 1 & prepStat == 1 & asmm == 1)
  }

  if (prep.elig.model == "none")  prepElig.asmm[idsEligStart] <-0
  
  # Core eligiblity scenarios
  if (prep.elig.model != "none") {
    if (substr(prep.elig.model, 1, 4) == "adol") {
      if (prep.elig.model == "adol.entry") {
        c1 <- active
        c2 <- active
        c3 <- active
      } else if (prep.elig.model == "adol.debuted") {
        c1 <- active
        c2 <- debuted
        c3 <- debuted
      } else if (prep.elig.model == "adol.AI") {
        c1 <- active
        c2 <- debuted
        c3 <- everAI
      } else if (prep.elig.model == "adol.entry.older") {
        c1 <- active
        c2 <- of.age
        c3 <- of.age
      } else if (prep.elig.model == "adol.debuted.older") {
        c1 <- active
        c2 <- debuted
        c3 <- of.age
      } else if (prep.elig.model == "adol.AI.older") {
        c1 <- active
        c2 <- everAI
        c3 <- of.age
      } else if (prep.elig.model == "adol.entry.time") {
        c1 <- active
        c2 <- ifelse (at - active.time > prep.delay,1,0) 
        c3 <- ifelse (at - active.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.debuted.time") {
        c1 <- debuted
        c2 <- ifelse (at - debuted.time > prep.delay,1,0) 
        c3 <- ifelse (at - debuted.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.AI.time") {
        c1 <- everAI
        c2 <- ifelse (at - AI.time > prep.delay,1,0) 
        c3 <- ifelse (at - AI.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.entry.older.time") {
        c1 <- active
        c2 <- of.age
        c3 <- ifelse (at - active.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.debuted.older.time") {
        c1 <- debuted
        c2 <- of.age
        c3 <- ifelse (at - debuted.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.AI.older.time") {
        c1 <- everAI
        c2 <- of.age
        c3 <- ifelse (at - AI.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.riskhist") {
        c1 <- ifelse (uaicount > prep.uaicount.thresh,1,0)
        c2 <- ifelse (uaicount > prep.uaicount.thresh,1,0)
        c3 <- ifelse (uaicount > prep.uaicount.thresh,1,0)
      } else if (prep.elig.model == "adol.riskhist.older") {
        c1 <- ifelse (uaicount > prep.uaicount.thresh,1,0)
        c2 <- of.age
        c3 <- of.age
      } else if (prep.elig.model == "adol.riskhist.time") {
        c1 <- ifelse (uaicount > prep.uaicount.thresh,1,0)
        c2 <- ifelse (at - AI.time > prep.delay,1,0)
        c3 <- ifelse (at - AI.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.riskhist.older.time") {
        c1 <- ifelse (uaicount > prep.uaicount.thresh,1,0)
        c2 <- of.age
        c3 <- ifelse (at - AI.time > prep.delay,1,0)
      }
      
     #Check length of c1 and idsEligStart 
      idsEligStart <- intersect(which(c1 > 0 & c2 > 0 & c3 > 0),idsEligStart)
      idsEligStop <- intersect(which(c1 == 0 | c2 == 0 | c3 == 0),idsEligStop)
    } 
    
    prepElig.asmm[idsEligStart] <- 1
    prepElig.asmm[idsEligStop] <- 0

    
  }

  

  ## Stoppage ------------------------------------------------------------------

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & asmm == 1 & diag.status == 1)


  # Death
  idsStpDth <- which(active == 0 & prepStat == 1 & asmm == 1)

  
  #Drop out
  idsStpDrop<-integer(0)
  if (prepSpell == TRUE){
  idsStpDrop <-which(active==1 & prepStat == 1 & asmm == 1)
  
  if (length(idsStpDrop>=1)){
          idsStpDrop <-ssample(idsStpDrop, prepDrop*length(idsStpDrop), replace=FALSE)
  }

  }
  
 #Drop out for MSM continuing adolecent PrEP into adulthood
 # msm.resid <-which(active==1 & prepStat == 1 & asmm == 0 & ever.adult.prep == 0)
 #idsStpDrop.msm.resid<-integer(0)
#  if (length(msm.resid>=1)){
#    msm.resid.list<-rbinom(length(msm.resid),1,prepDrop)
#    idsStpDrop.msm.resid <-msm.resid[msm.resid.list==1]
#  }

  
  # Transition to ineligibility
  idsStpInelig <- idsEligStop
  
  #aged out
  idsEligEnd <- which(active ==1 & prepElig.asmm ==1 & asmm ==0)
  
  # Reset PrEP status
 # idsStp <- c(idsStpDx, idsStpDth, idsStpDrop, idsStpInelig, idsStpDrop.msm.resid, idsEligEnd)
  idsStp <- c(idsStpDx, idsStpDth, idsStpDrop, idsStpInelig, idsEligEnd)
  prepStat[idsStp] <- 0
  prepClass[idsStp] <-NA
  

  #Drops are added back to eligible list.
  prepElig.asmm[idsStpDrop] <- 1
  
  # Remove those that have aged out from eligibility

  prepElig.asmm[idsEligEnd] <-0
  
  ## Initiation ----------------------------------------------------------------


  if (prep.cov.method == "curr") {
    prepCov <- sum(prepStat == 1 & asmm == 1, na.rm = TRUE)/sum(prepElig.asmm == 1, na.rm = TRUE)
  }
  if (prep.cov.method == "ever") {
    prepCov <- sum(prepEver == 1 & asmm == 1, na.rm = TRUE)/sum(prepElig.asmm == 1, na.rm = TRUE)
  }
  prepCov <- ifelse(is.nan(prepCov), 0, prepCov)

  
  idsEligSt <- which(prepElig.asmm == 1 & asmm == 1)
  nEligSt <- length(idsEligSt)
  

  nStart <- max(0, min(nEligSt, round((prep.coverage - prepCov) *
                                        sum(prepElig.asmm == 1, na.rm = TRUE))))
  
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
    ever.adol.prep[idsStart] <- 1
    prepStart.time[idsStart]<-at
    
    # PrEP class is fixed over PrEP cycles
    needPC <- which(is.na(prepClass[idsStart]))
    needPC <- idsStart[needPC]
    prepClass[needPC] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
  }

  ## Output --------------------------------------------------------------------
  # Attributes
  
  dat$attr$ever.adol.prep <- ever.adol.prep 
  dat$attr$prepElig.asmm <- prepElig.asmm
  dat$attr$prepStat <- prepStat
  dat$attr$prepEver <- prepEver
  dat$attr$prepClass <- prepClass
  dat$attr$prepStart.time<-prepStart.time
  
  # Summary Statistics
  dat$epi$prepCov.asmm[at] <- prepCov
  dat$epi$prepStart.asmm[at] <- length(idsStart)

  dat$epi$prepStart[at] <- dat$epi$prepStart[at] + length(idsStart)
  
 

  return(dat)
}
