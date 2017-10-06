

#' @title Transmission Module for those immigrating out of the country and returning
#'
#' @description Stochastically simulates disease transmission from the heat bath given the current
#'              state of the active HIV negative population.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This function takes all active HIV negative MSMF nodes and calculates a
#' transmission probability for each individual based on race specific exposure to the heatbath
#' which incorperates race specific prevalence estiamtes for determining risk.
#' After transmission events, individual-level attributes for the infected
#' persons are updated and summary statistics for incidence calculated.
#'
#' The input values \code{dat$param$msm.aq.prob.'R'} were R is a race specific probability taking values (B,BI,H,HI,W).
#' 
#'
#' @return
#' For each new infection, the disease status, infection time, and related
#' HIV attributes are updated for the infected node. Summary statistics for
#' disease incidence overall, and by race are calculated and
#' stored on \code{dat$epi}.
#'
#' @keywords module msmf msm heatbath
#' @export
#'

immigration_shamp <- function(dat, at){
  
  
  ## Infected from the heat bath
  #Attributes

  race <- dat$attr$race
  sex <- dat$attr$sex
  active <- dat$attr$active
  status <- dat$attr$status
  immig.loc <- dat$attr$immig.loc
  sex.ident <- dat$attr$sex.ident

  #Parameters

  depart.BI.f <- dat$param$immig.depart.BI.f
  depart.HI.f <- dat$param$immig.depart.HI.f
  depart.BI.m <- dat$param$immig.depart.BI.m
  depart.HI.m <- dat$param$immig.depart.HI.m
  return.BI.f <- dat$param$immig.return.BI.f
  return.HI.f <- dat$param$immig.return.HI.f
  return.BI.m <- dat$param$immig.return.BI.m
  return.HI.m <- dat$param$immig.return.HI.m
  aq.prob.BI.f <- dat$param$immig.aq.prob.BI.f
  aq.prob.HI.f <- dat$param$immig.aq.prob.HI.f
  aq.prob.BI.m <- dat$param$immig.aq.prob.BI.m
  aq.prob.HI.m <- dat$param$immig.aq.prob.HI.m
  
  infected <- NULL
  ids.BI.f.inf <- NULL 
  ids.BI.m.inf <- NULL
  ids.HI.f.inf <- NULL
  ids.HI.m.inf <- NULL
  
  departing <- NULL
  ids.BI.f.d <- NULL
  ids.BI.m.d <- NULL
  ids.HI.f.d <- NULL
  ids.HI.m.d <- NULL
  
  retruning <- NULL
  ids.BI.f.r <- NULL
  ids.BI.m.r <- NULL
  ids.HI.f.r <- NULL
  ids.HI.m.r <- NULL
  
  #Of those away who becomes infected
  ids.BI.f.inf<-which(active==1 & immig.loc==1 & sex=="F" & race =="BI" & status==0)
  ids.BI.m.inf<-which(active==1 & immig.loc==1 & sex=="M" & race =="BI" & status==0)
  ids.H.f.inf<-which(active==1 & immig.loc==1 & sex=="F" & race =="HI" & status==0)
  ids.HI.m.inf<-which(active==1 & immig.loc==1 & sex=="M" & race =="HI" & status==0)
  
  if(length(ids.BI.f.inf) > 0) {
  ids.BI.f.inf.t<-rbinom(length(ids.BI.f.inf),1,prob = aq.prob.BI.f)
  ids.BI.f.inf<-ids.BI.f.inf[ids.BI.f.inf.t==1]}

  if(length(ids.BI.m.inf) > 0) {  
  ids.BI.m.inf.t<-rbinom(length(ids.BI.m.inf),1,prob = aq.prob.BI.m)
  ids.BI.m.inf<-ids.BI.m.inf[ids.BI.m.inf.t==1]}
  
  if(length(ids.HI.f.inf) > 0) {
  ids.HI.f.inf.t<-rbinom(length(ids.HI.f.inf),1,prob = aq.prob.HI.f)
  ids.HI.f.inf<-ids.HI.f.inf[ids.HI.f.inf.t==1]}
  
  if(length(ids.HI.m.inf) > 0) {
  ids.HI.m.inf.t<-rbinom(length(ids.HI.m.inf),1,prob = aq.prob.HI.m)
  ids.HI.m.inf<-ids.HI.m.inf[ids.HI.m.inf.t==1]}
  
  infected <- c(ids.BI.f.inf, ids.BI.m.inf, ids.HI.f.inf, ids.HI.m.inf)
 


  #Who is leaving
  ids.BI.f.d<-which(active==1 & immig.loc==0 & sex=="F" & race =="BI")
  ids.BI.m.d<-which(active==1 & immig.loc==0 & sex=="M" & race =="BI")
  ids.HI.f.d<-which(active==1 & immig.loc==0 & sex=="F" & race =="HI")
  ids.HI.m.d<-which(active==1 & immig.loc==0 & sex=="M" & race =="HI")
  
  if(length(ids.BI.f.d) > 0){
  ids.BI.f.d.t<-rbinom(length(ids.BI.f.d),1,prob = depart.BI.f)
  ids.BI.f.d<-ids.BI.f.d[ids.BI.f.d.t==1]}
  
  if(length(ids.BI.m.d) > 0){
  ids.BI.m.d.t<-rbinom(length(ids.BI.m.d),1,prob = depart.BI.m)
  ids.BI.m.d<-ids.BI.m.d[ids.BI.m.d.t==1]}
  
  if(length(ids.HI.f.d) > 0){
  ids.HI.f.d.t<-rbinom(length(ids.HI.f.d),1,prob = depart.HI.f)
  ids.HI.f.d<-ids.HI.f.d[ids.HI.f.d.t==1]}
  
  if(length(ids.HI.m.d) > 0){
  ids.HI.m.d.t<-rbinom(length(ids.HI.m.d),1,prob = depart.HI.m)
  ids.HI.m.d<-ids.HI.m.d[ids.HI.m.d.t==1]}

  departing <- c(ids.BI.f.d, ids.BI.m.d, ids.HI.f.d, ids.HI.m.d)

  #Who is returning
  ids.BI.f.r<-which(active==1 & immig.loc==1 & sex=="F" & race =="BI")
  ids.BI.m.r<-which(active==1 & immig.loc==1 & sex=="M" & race =="BI")
  ids.HI.f.r<-which(active==1 & immig.loc==1 & sex=="F" & race =="HI")
  ids.HI.m.r<-which(active==1 & immig.loc==1 & sex=="M" & race =="HI")
  
  if(length(ids.BI.f.r) > 0){
  ids.BI.f.r.t<-rbinom(length(ids.BI.f.r),1,prob = return.BI.f)
  ids.BI.f.r<-ids.BI.f.r[ids.BI.f.r.t==1]}
  
  if(length(ids.BI.m.r) > 0){
  ids.BI.m.r.t<-rbinom(length(ids.BI.m.r),1,prob = return.BI.m)
  ids.BI.m.r<-ids.BI.m.r[ids.BI.m.r.t==1]}
  
  if(length(ids.HI.f.r) > 0){
  ids.HI.f.r.t<-rbinom(length(ids.HI.f.r),1,prob = return.HI.f)
  ids.HI.f.r<-ids.HI.f.r[ids.HI.f.r.t==1]}
  
  if(length(ids.HI.m.r) > 0){
  ids.HI.m.r.t<-rbinom(length(ids.HI.m.r),1,prob = return.HI.m)
  ids.HI.m.r<-ids.HI.m.r[ids.HI.m.r.t==1]}
  
  returning <- c(ids.BI.f.r, ids.BI.m.r, ids.HI.f.r, ids.HI.m.r)
  
  
  ##Move the people
  
  if (length(returning) > 0) {dat$attr$immig.loc[returning]<-0}
  if (length(departing) > 0) {dat$attr$immig.loc[departing]<-1} 

  ##Track infections
 if (length(infected) >= 1){
    dat$attr$status[infected] <- 1
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- "AR"
    dat$attr$stage.time[infected] <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected] <- 0
    dat$attr$inf.class[infected] <- "FA"
      
    
    dat$attr$infector[infected] <- "FA"
    dat$attr$inf.role[infected] <- "FA"
    dat$attr$inf.type[infected] <- "FA"
    dat$attr$inf.diag[infected] <- "FA" 
    dat$attr$inf.tx[infected] <- "FA"
    dat$attr$inf.stage[infected] <- "FA"

    dat$attr$cum.time.on.tx[infected] <- 0
    dat$attr$cum.time.off.tx[infected] <- 0
    
    
    # Summary Output
    dat$epi$incid[at] <- dat$epi$incid[at] + length(infected)
    
    dat$epi$incid.BI[at] <- sum(dat$epi$incid.BI[at] , sum(race[infected] == "BI", na.rm = TRUE), na.rm = TRUE)
    dat$epi$incid.HI[at] <- sum(dat$epi$incid.HI[at] , sum(race[infected] == "HI", na.rm = TRUE), na.rm = TRUE)

    dat$epi$incid.f[at] <- sum(dat$epi$incid.f[at] , sum(sex[infected] == "F", na.rm = TRUE), na.rm = TRUE) 
    dat$epi$incid.m[at] <- sum(dat$epi$incid.m[at] , sum(sex[infected] == "M", na.rm = TRUE), na.rm = TRUE)
    
    dat$epi$incid.msf[at] <- sum(dat$epi$incid.msf[at] , sum(sex.ident[infected] == "msf", na.rm = TRUE), na.rm = TRUE)
    dat$epi$incid.msm[at] <- sum(dat$epi$incid.msm[at] , sum(sex.ident[infected] == "msm", na.rm = TRUE), na.rm = TRUE)
    dat$epi$incid.msmf[at] <- sum(dat$epi$incid.msmf[at] , sum(sex.ident[infected] == "msmf", na.rm = TRUE), na.rm = TRUE)
    
    dat$epi$incid.FA[at] <- length(infected)
     }
    
    if (length(infected) < 1){ dat$epi$incid.FA[at] <-0}
  
  
    
  return(dat)
}





