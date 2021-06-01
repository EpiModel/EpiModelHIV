
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

  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  prep.ps <- dat$param$prep.ps
  prep.risk <- dat$param$prep.risk
  prep.disc <- dat$param$prep.disc
  prep.class.prob <- dat$param$prep.class.prob
  prep.start.prob.ps <- dat$param$prep.start.prob.ps
  prep.start.prob.parts <- dat$param$prep.start.prob.parts
  prep.stop.prob <- dat$param$prep.stop.prob
  prep.window <- dat$param$prep.window
  
  prep.ind.ps <- dat$attr$prep.ind.ps
  prep.ind.disc <- dat$attr$prep.ind.disc
  prep.ind.parts <- dat$attr$prep.ind.parts



  ## Eligibility ---------------------------------------------------------------
  idsEligStart.ps <-NULL
  idsStart.ps <-NULL
  idsEligStart.disc <-NULL
  idsStart.disc <-NULL
  idsEligStart.parts <-NULL
  idsStart.parts <-NULL
  
  # Base eligibility
  if(prep.ps == TRUE){
  idsEligStart.ps <- which(prep.ind.ps == 1 & prepStat==0 & status == 0)
  selected <- rbinom(length(idsEligStart.ps),1,prep.start.prob.ps)
  idsStart.ps <-idsEligStart.ps[selected==1]}
  
  #Use the same Prep start prob as ps since both are a form of known discordant rel
  if(prep.disc == TRUE){
    idsEligStart.disc <- which(dat$attr$tested.negative == 1 & prep.ind.disc == 1 & prepStat==0 & status == 0)
    selected <- rbinom(length(idsEligStart.disc),1,prep.start.prob.ps)
    idsStart.disc <-idsEligStart.disc[selected==1]}
  
  if(prep.risk == "RISK"){
    idsEligStart.parts <- which(dat$attr$tested.negative == 1 & prep.ind.parts == 1 & prepStat==0 & status == 0)
    selected <- rbinom(length(idsEligStart.parts),1,prep.start.prob.parts)
    idsStart.parts <-idsEligStart.parts[selected==1]}
  
  if(prep.risk == "ALL"){
    idsEligStart.parts <- which(dat$attr$tested.negative == 1 & prepStat==0 & status == 0)
    selected <- rbinom(length(idsEligStart.parts),1,prep.start.prob.parts)
    idsStart.parts <-idsEligStart.parts[selected==1]}
  
  idsEligStop <- which(prepStat == 1)
  selected <- rbinom(length(idsEligStop),1,prep.stop.prob)
  idsStop <-idsEligStop[selected==1] 


  ## Stoppage ------------------------------------------------------------------

  # Diagnosis
  idsStpDx <- which(prepStat == 1 & diag.status == 1)


  # Reset PrEP status
  idsStp <- c(idsStpDx, idsStop)
  prepStat[idsStp] <- 0

  ## Initiation ----------------------------------------------------------------


  # Attributes
  ids.Elig <- c(idsEligStart.ps, idsEligStart.disc, idsEligStart.parts)
  dat$epi$prepElig[at] <- max(0,length(ids.Elig))
  idsStart <- c(idsStart.ps, idsStart.disc, idsStart.parts)
  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1

    
    # PrEP class is fixed over PrEP cycles
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
  }


  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepStat <- prepStat
  dat$attr$prepClass <- prepClass

  # Summary Statistics
  dat$epi$prepCov[at] <- sum(prepStat)/(sum(dat$epi$num.poi,na.rm = TRUE)-sum(dat$epi$i.num.poi,na.rm = TRUE))
  dat$epi$prepStart[at] <- length(idsStart)
  
  #Clear PREP indication attributed for the next time step
  dat$epi$prep.ind.ps[at] <-max(0,sum(dat$attr$prep.ind.ps,na.rm = TRUE))
  dat$epi$prep.ind.disc[at] <-max(0,sum(dat$attr$prep.ind.disc,na.rm=TRUE))
  dat$epi$prep.ind.parts[at] <-max(0,sum(dat$attr$prep.ind.parts,na.rm=TRUE))
  
  dat$attr$prep.ind.ps <-rep(0,length(dat$attr$prep.ind.ps))
  dat$attr$prep.ind.disc <-rep(0,length(dat$attr$prep.ind.disc))
  dat$attr$prep.ind.parts <-rep(0,length(dat$attr$prep.ind.parts))
  
  return(dat)
}
