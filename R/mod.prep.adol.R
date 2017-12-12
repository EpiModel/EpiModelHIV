
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
  debuted<-dat$attr$has.debuted
  debuted.time<-dat$attr$has.debuted.time
  active.time<-dat$attr$active.time
  AI.time<-dat$attr$AI.time
  lnt <- dat$attr$last.neg.test
  everAI<-dat$attr$everAI
  of.age<-dat$attr$of.age
  uaicount<-dat$attr$uaicount
 
  prepElig <- dat$attr$prepElig
  prepElig.asmm <- dat$attr$prepElig.asmm
  prepStat <- dat$attr$prepStat
  prepEver <- dat$attr$prepEver
  prepClass <- dat$attr$prepClass
  prepStart.time <- dat$attr$prepStart.time
  ever.adol.prep <- dat$attr$ever.adol.prep
  

##PrEP params.
  #Uniform.
  prep.elig.model <- dat$param$prep.elig.model.asmm
  prep.cov.method <- dat$param$prep.cov.method.asmm
  prep.risk.reassess <- dat$param$prep.risk.reassess.asmm
  prepSpell<- dat$param$prepSpell.asmm
  
  ##Race specific.
  prep.coverage.b <- dat$param$prep.coverage.b.asmm
  prep.cov.rate.b <- dat$param$prep.cov.rate.b.asmm
  prep.class.prob.b <- dat$param$prep.class.prob.b.asmm
  prepDrop.b<- dat$param$prepDrop.b.asmm
  prep.delay.b<-dat$param$prep.delay.b.asmm
  prep.uaicount.thresh.b<-dat$param$prep.uaicount.thresh.b.asmm
  
  prep.coverage.w <- dat$param$prep.coverage.w.asmm
  prep.cov.rate.w <- dat$param$prep.cov.rate.w.asmm
  prep.class.prob.w <- dat$param$prep.class.prob.w.asmm
  prepDrop.w<- dat$param$prepDrop.w.asmm
  prep.delay.w<-dat$param$prep.delay.w.asmm
  prep.uaicount.thresh.w<-dat$param$prep.uaicount.thresh.w.asmm
  

  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart.b <- which(active == 1 & status == 0 & prepStat == 0 & asmm == 1 & race=="B")
  idsEligStart.w <- which(active == 1 & status == 0 & prepStat == 0 & asmm == 1 & race=="W")
  
  idsEligStop.b <- NULL
  if (prep.risk.reassess == TRUE) {
    idsEligStop.b <- which(active == 1 & prepStat == 1 & asmm == 1 & race=="B")
  }
  
  idsEligStop.w <- NULL
  if (prep.risk.reassess == TRUE) {
    idsEligStop.w <- which(active == 1 & prepStat == 1 & asmm == 1 & race=="W")
  }

  if (prep.elig.model == "none")  prepElig[idsEligStart.b] <-0
  if (prep.elig.model == "none")  prepElig[idsEligStart.w] <-0
  
  # Core eligiblity scenarios
  if (prep.elig.model != "none") {
    if (substr(prep.elig.model, 1, 4) == "adol") {
      if (prep.elig.model == "adol.entry") {
        c1.b <- active
        c2.b <- active
        c3.b <- active
        c1.w <- active
        c2.w <- active
        c3.w <- active
      } else if (prep.elig.model == "adol.debuted") {
        c1.b <- active
        c2.b <- debuted
        c3.b <- debuted
        c1.w <- active
        c2.w <- debuted
        c3.w <- debuted
      } else if (prep.elig.model == "adol.AI") {
        c1.b <- active
        c2.b <- debuted
        c3.b <- everAI
        c1.w <- active
        c2.w <- debuted
        c3.w <- everAI
      } else if (prep.elig.model == "adol.entry.older") {
        c1.b <- active
        c2.b <- of.age
        c3.b <- of.age
        c1.w <- active
        c2.w <- of.age
        c3.w <- of.age
      } else if (prep.elig.model == "adol.debuted.older") {
        c1.b <- active
        c2.b <- debuted
        c3.b <- of.age
        c1.w <- active
        c2.w <- debuted
        c3.w <- of.age
      } else if (prep.elig.model == "adol.AI.older") {
        c1.b <- active
        c2.b <- everAI
        c3.b <- of.age
        c1.w <- active
        c2.w <- everAI
        c3.w <- of.age
      } else if (prep.elig.model == "adol.entry.time") {
        c1.b <- active
        c2.b <- ifelse (at - active.time > prep.delay.b,1,0) 
        c3.b <- ifelse (at - active.time > prep.delay.b,1,0)
        c1.w <- active
        c2.w <- ifelse (at - active.time > prep.delay.w,1,0) 
        c3.w <- ifelse (at - active.time > prep.delay.w,1,0)
      } else if (prep.elig.model == "adol.debuted.time") {
        c1.b <- debuted
        c2.b <- ifelse (at - debuted.time > prep.delay.b,1,0) 
        c3.b <- ifelse (at - debuted.time > prep.delay.b,1,0)
        c1.w <- debuted
        c2.w <- ifelse (at - debuted.time > prep.delay.w,1,0) 
        c3.w <- ifelse (at - debuted.time > prep.delay.w,1,0) 
      } else if (prep.elig.model == "adol.AI.time") {
        c1.b <- everAI
        c2.b <- ifelse (at - AI.time > prep.delay.b,1,0) 
        c3.b <- ifelse (at - AI.time > prep.delay.b,1,0)
        c1.w <- everAI
        c2.w <- ifelse (at - AI.time > prep.delay.w,1,0) 
        c3.w <- ifelse (at - AI.time > prep.delay.w,1,0)
      } else if (prep.elig.model == "adol.entry.older.time") {
        c1.b <- active
        c2.b <- of.age
        c3.b <- ifelse (at - active.time > prep.delay.b,1,0)
        c1.w <- active
        c2.w <- of.age
        c3.w <- ifelse (at - active.time > prep.delay.w,1,0)
      } else if (prep.elig.model == "adol.debuted.older.time") {
        c1.b <- debuted
        c2.b <- of.age
        c3.b <- ifelse (at - debuted.time > prep.delay.b,1,0)
        c1.w <- debuted
        c2.w <- of.age
        c3.w <- ifelse (at - debuted.time > prep.delay.w,1,0)
      } else if (prep.elig.model == "adol.AI.older.time") {
        c1.b <- everAI
        c2.b <- of.age
        c3.b <- ifelse (at - AI.time > prep.delay.b,1,0)
        c1.w <- everAI
        c2.w <- of.age
        c3.w <- ifelse (at - AI.time > prep.delay.w,1,0)
      } else if (prep.elig.model == "adol.riskhist") {
        c1.b <- ifelse (uaicount > prep.uaicount.thresh.b,1,0)
        c2.b <- ifelse (uaicount > prep.uaicount.thresh.b,1,0)
        c3.b <- ifelse (uaicount > prep.uaicount.thresh.b,1,0)
        c1.w <- ifelse (uaicount > prep.uaicount.thresh.w,1,0)
        c2.w <- ifelse (uaicount > prep.uaicount.thresh.w,1,0)
        c3.w <- ifelse (uaicount > prep.uaicount.thresh.w,1,0)
      } else if (prep.elig.model == "adol.riskhist.older") {
        c1.b <- ifelse (uaicount > prep.uaicount.thresh.b,1,0)
        c2.b <- of.age
        c3.b <- of.age
        c1.w <- ifelse (uaicount > prep.uaicount.thresh.w,1,0)
        c2.w <- of.age
        c3.w <- of.age
      } else if (prep.elig.model == "adol.riskhist.time") {
        c1.b <- ifelse (uaicount > prep.uaicount.thresh.b,1,0)
        c2.b <- ifelse (at - AI.time > prep.delay.b,1,0)
        c3.b <- ifelse (at - AI.time > prep.delay.b,1,0)
        c1.w <- ifelse (uaicount > prep.uaicount.thresh.w,1,0)
        c2.w <- ifelse (at - AI.time > prep.delay.w,1,0)
        c3.w <- ifelse (at - AI.time > prep.delay.w,1,0)
      } else if (prep.elig.model == "adol.riskhist.older.time") {
        c1.b <- ifelse (uaicount > prep.uaicount.thresh.b,1,0)
        c2.b <- of.age
        c3.b <- ifelse (at - AI.time > prep.delay.b,1,0)
        c1.w <- ifelse (uaicount > prep.uaicount.thresh.w,1,0)
        c2.w <- of.age
        c3.w <- ifelse (at - AI.time > prep.delay.w,1,0)
      }
      
     #Check length of c1.b and idsEligStart.b 
      idsEligStart.b <- intersect(which(c1.b > 0 & c2.b > 0 & c3.b > 0),idsEligStart.b)
      idsEligStop.b <- intersect(which(c1.b == 0 | c2.b == 0 | c3.b == 0),idsEligStop.b)
      
      idsEligStart.w <- intersect(which(c1.w > 0 & c2.w > 0 & c3.w > 0),idsEligStart.w)
      idsEligStop.w <- intersect(which(c1.w == 0 | c2.w == 0 | c3.w == 0),idsEligStop.w)
    } 
    
    prepElig[idsEligStart.b] <- 1
    prepElig[idsEligStop.b] <- 0  
    
    prepElig[idsEligStart.w] <- 1
    prepElig[idsEligStop.w] <- 0 
    
    prepElig.asmm[idsEligStart.b] <- 1
    prepElig.asmm[idsEligStop.b] <- 0
    prepElig.asmm[idsEligStart.w] <- 1
    prepElig.asmm[idsEligStop.w] <- 0 
    
  }

  

  ## Stoppage ------------------------------------------------------------------

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & asmm == 1, diag.status == 1)


  # Death
  idsStpDth <- which(active == 0 & prepStat == 1 & asmm == 1)

  
  #Drop out
  idsStpDrop.b<-integer(0)
  idsStpDrop.w<-integer(0)
  if (prepSpell == TRUE){
  idsStpDrop.b <-which(active==1 & prepStat == 1 & race=="B" & asmm == 1)
  idsStpDrop.w <-which(active==1 & prepStat == 1 & race=="W" & asmm == 1)
  
  if (length(idsStpDrop.b>=1)){
          idsStpDrop.b <-ssample(idsStpDrop.b, prepDrop.b*length(idsStpDrop.b), replace=FALSE)
  }
  if (length(idsStpDrop.w>=1)){
          idsStpDrop.w <-ssample(idsStpDrop.w, prepDrop.w*length(idsStpDrop.w), replace=FALSE)
  }
  }
  

  
  # Transition to ineligibility
  idsStpInelig.b <- idsEligStop.b
  idsStpInelig.w <- idsEligStop.w

  # Reset PrEP status
  idsStp <- c(idsStpDx, idsStpDth, idsStpDrop.b, idsStpDrop.w, idsStpInelig.b, idsStpInelig.w)
  prepStat[idsStp] <- 0
  prepClass[idsStp] <-NA
  

  #Drops are added back to eligible list.
  prepElig[idsStpDrop.b] <- 1
  prepElig[idsStpDrop.w] <- 1
  prepElig.asmm[idsStpDrop.b] <- 1
  prepElig.asmm[idsStpDrop.w] <- 1

  ## Initiation ----------------------------------------------------------------
  blk<-which(race=="B")
  wht<-which(race=="W")
  

  if (prep.cov.method == "curr") {
    prepCov.b <- sum(prepStat[blk] == 1 & asmm[blk] == 1, na.rm = TRUE)/sum(prepElig[blk] == 1 & asmm[blk] == 1, na.rm = TRUE)
    prepCov.w <- sum(prepStat[wht] == 1 & asmm[wht] == 1, na.rm = TRUE)/sum(prepElig[wht] == 1 & asmm[wht] == 1, na.rm = TRUE)
  }
  if (prep.cov.method == "ever") {
    prepCov.b <- sum(prepEver[blk] == 1 & asmm[blk] == 1, na.rm = TRUE)/sum(prepElig[blk] == 1 & asmm[blk] == 1, na.rm = TRUE)
    prepCov.w <- sum(prepEver[wht] == 1 & asmm[wht] == 1, na.rm = TRUE)/sum(prepElig[wht] == 1 & asmm[wht] == 1, na.rm = TRUE)
  }
  prepCov.b <- ifelse(is.nan(prepCov.b), 0, prepCov.b)
  prepCov.w <- ifelse(is.nan(prepCov.w), 0, prepCov.w)
  
  idsEligSt.b <- which(prepElig == 1 & race == "B" & asmm == 1)
  nEligSt.b <- length(idsEligSt.b)
  
  idsEligSt.w <- which(prepElig == 1 & race== "W" & asmm == 1)
  nEligSt.w <- length(idsEligSt.w)

  nStart.b <- max(0, min(nEligSt.b, round((prep.coverage.b - prepCov.b) *
                                             sum(prepElig[blk] == 1 & asmm[blk] == 1, na.rm = TRUE))))
  
  nStart.w <- max(0, min(nEligSt.w, round((prep.coverage.w - prepCov.w) *
                                        sum(prepElig[wht] == 1 & asmm[wht] == 1, na.rm = TRUE))))
  
  idsStart.b <- NULL
  if (nStart.b > 0) {
    if (prep.cov.rate.b >= 1) {
      idsStart.b <- ssample(idsEligSt.b, nStart.b)
    } else {
      idsStart.b <- idsEligSt.b[rbinom(nStart.b, 1, prep.cov.rate.b) == 1]
    }
  }
  
  idsStart.w <- NULL
  if (nStart.w > 0) {
    if (prep.cov.rate.w >= 1) {
      idsStart.w <- ssample(idsEligSt.w, nStart.w)
    } else {
      idsStart.w <- idsEligSt.w[rbinom(nStart.w, 1, prep.cov.rate.w) == 1]
    }
  } 
  
  
  
  

  # Attributes
  if (length(idsStart.b) > 0) {
    prepStat[idsStart.b] <- 1
    prepEver[idsStart.b] <- 1
    ever.adol.prep[idsStart.b] <- 1
    prepStart.time[idsStart.b]<-at

    # PrEP class is fixed over PrEP cycles
    needPC.b <- which(is.na(prepClass[idsStart.b]))
    prepClass[needPC.b] <- sample(x = 0:3, size = length(needPC.b),
                                          replace = TRUE, prob = prep.class.prob.b)
  }

  if (length(idsStart.w) > 0) {
    prepStat[idsStart.w] <- 1
    prepEver[idsStart.w] <- 1
    ever.adol.prep[idsStart.w] <- 1
    prepStart.time[idsStart.w]<-at
    
    # PrEP class is fixed over PrEP cycles
    needPC.w <- which(is.na(prepClass[idsStart.w]))
    prepClass[needPC.w] <- sample(x = 0:3, size = length(needPC.w),
                                          replace = TRUE, prob = prep.class.prob.w)
  }

  ## Output --------------------------------------------------------------------
  # Attributes
  
  dat$attr$ever.adol.prep <- ever.adol.prep 
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepEver <- prepEver
  dat$attr$prepClass <- prepClass
  dat$attr$prepStart.time<-prepStart.time
  
  # Summary Statistics
  dat$epi$prepCov.B.asmm[at] <- prepCov.b
  dat$epi$prepCov.W.asmm[at] <- prepCov.w
  dat$epi$prepStart.B.asmm[at] <- length(idsStart.b)
  dat$epi$prepStart.W.asmm[at] <- length(idsStart.w)
  dat$epi$prepStart.asmm[at] <- length(idsStart.b) + length(idsStart.w)
  dat$epi$prepStart[at] <- dat$epi$prepStart[at] + length(idsStart.b) + length(idsStart.w) 
  
  dat$epi$prepCov.asmm[at] <- (prepCov.b*(sum(active[blk] & asmm[blk] == 1))+prepCov.w*(sum(active[wht] & asmm[wht] == 1)))/(sum(active  & asmm == 1))

 

  return(dat)
}
