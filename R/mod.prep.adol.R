
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
  age<- dat$attr$age


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
  prepStat.asmm <- dat$attr$prepStat.asmm
  prepEver <- dat$attr$prepEver
  prepClass <- dat$attr$prepClass
  prepStart.time <- dat$attr$prepStart.time
  ever.adol.prep <- dat$attr$ever.adol.prep

  
##PrEP params.
  #Uniform.
  prep.elig.model <- dat$param$prep.elig.model.asmm
  prep.delay<-dat$param$prep.delay.asmm
  prep.uaicount.thresh<-dat$param$prep.uaicount.thresh.asmm

  prep.uptake.asmm <- dat$param$prep.uptake.asmm
  prep.disc.asmm <- dat$param$prep.disc.asmm
  prep.class.prob.asmm <- dat$param$prep.class.prob.asmm

  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 & status == 0 & prepStat == 0 & asmm == 1)

  idsEligStop <- which(active == 1 & prepStat.asmm == 1)


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
    } 
    
    notElig <- which(asmm != 1)
    prepElig.asmm[idsEligStart] <- 1
    prepElig.asmm[notElig] <- 0

    
  }

  

  ## Stoppage ------------------------------------------------------------------

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat.asmm == 1 & diag.status == 1)


  # Death
  idsStpDth <- which(active == 0 & prepStat.asmm == 1)

  
  #Drop out
  idsStpDrop<-integer(0)

  idsStpDrop <-which(active==1 & prepStat.asmm == 1 & (diag.status != 1 | is.na(diag.status)==TRUE))
  
  if (length(idsStpDrop>=1)){
          drops <- rbinom(length(idsStpDrop),1,prep.disc.asmm) 
          idsStpDrop <-idsStpDrop[drops==1]
  }

 
  

  

  # Reset PrEP status
  ##Looks like those eding eligibility are reset to not prep here.
  idsStp <- c(idsStpDx, idsStpDth, idsStpDrop)
  prepStat[idsStp] <- 0
  prepStat.asmm[idsStp] <- 0
  prepClass[idsStp] <-NA
  

  ## Initiation ----------------------------------------------------------------



 prepCov <- sum(prepStat.asmm == 1 & asmm == 1, na.rm = TRUE)/sum(prepElig.asmm == 1 & asmm == 1, na.rm = TRUE)



  
  idsStart <- NULL
 if (length(idsEligStart) > 0){
      ids <- rbinom(length(idsEligStart),1,prep.uptake.asmm)
      idsStart <- idsEligStart[ids==1]
  } 
  
  
  
  

  # Attributes

  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1
    prepStat.asmm[idsStart] <- 1
    prepEver[idsStart] <- 1
    ever.adol.prep[idsStart] <- 1
    prepStart.time[idsStart]<-at
    
    # PrEP class is fixed over PrEP cycles
    needPC <- which(is.na(prepClass[idsStart]))
    needPC <- idsStart[needPC]
    prepClass[needPC] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob.asmm)
  }

  ## Output --------------------------------------------------------------------
  # Attributes
  
  dat$attr$ever.adol.prep <- ever.adol.prep 
  dat$attr$prepElig.asmm <- prepElig.asmm
  dat$attr$prepStat <- prepStat
  dat$attr$prepStat.asmm <- prepStat.asmm
  dat$attr$prepEver <- prepEver
  dat$attr$prepClass <- prepClass
  dat$attr$prepStart.time<-prepStart.time

  
  # Summary Statistics
  dat$epi$prepCov.asmm[at] <- prepCov
  dat$epi$prepStart.asmm[at] <- length(idsStart)

 

  return(dat)
}
