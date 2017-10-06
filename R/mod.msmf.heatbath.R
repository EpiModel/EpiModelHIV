

#' @title Transmission Module for Heat Bath Exposure of MSMF to an MSM heatbath
#'
#' @description Stochastically simulates disease transmission from the heat bath given the current
#'              state of the active HIV negative population.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This function takes all active HIV negative MSMF nodes and calculates a
#' transmission probability for each individual based race specific exposure to the heatbath
#' which incorperates race specific prevalence estiamtes for determining risk.
#' After transmission events, individual-level attributes for the infected
#' persons are updated and summary statistics for incidence calculated.
#'
#' The input values \code{dat$param$msm.aq.prob.'R'} 
#' where R is a race specific probability taking values (B,BI,H,HI,W).
#' 
#'
#' @return
#' For each new infection, the disease status, infection time, and related
#' HIV attributes are updated for the infected node. Summary statistics for
#' disease incidence overall, and by race are calculated and
#' stored on \code{dat$epi}.
#'
#' @keywords module SHAMP msmf msm heatbath
#' @export
#'

heatbath_msmf_shamp <- function(dat, at){
  
  
  ## Infected from the heat bath
  #Attributes

  race<-dat$attr$race
  active<-dat$attr$active
  status<-dat$attr$status
  prepStat<-dat$attr$prepStat
  prepClass<-dat$attr$prepClass
  sex<-dat$attr$sex
  sex.ident<-dat$attr$sex.ident
  

  #Parameters
  prep.hr <- dat$param$prep.class.hr
  prep.class.prob <- dat$param$prep.class.prob
  
  msm.aq.prob.B <- dat$param$msm.aq.prob.B
  msm.aq.prob.BI <- dat$param$msm.aq.prob.BI
  msm.aq.prob.H <- dat$param$msm.aq.prob.H
  msm.aq.prob.HI <- dat$param$msm.aq.prob.HI
  msm.aq.prob.W <- dat$param$msm.aq.prob.W


  #Select active negative msmf by race
  ids.B<-which(active==1 & sex.ident=="msmf" & race =="B" & status == 0)
  ids.BI<-which(active==1 & sex.ident=="msmf" & race =="BI" & status == 0)
  ids.H<-which(active==1 & sex.ident=="msmf" & race =="H" & status == 0)
  ids.HI<-which(active==1 & sex.ident=="msmf" & race =="HI" & status == 0)
  ids.W<-which(active==1 & sex.ident=="msmf" & race =="W" & status == 0)
  
  #Determine who becomes infected 
 

  infected.B <- rep(0,length(ids.B))
  infected.BI <- rep(0,length(ids.BI))
  infected.H <- rep(0,length(ids.H))
  infected.HI <- rep(0,length(ids.HI))
  infected.W <- rep(0,length(ids.W))
  

  
#Calculate new infections.


    infected.B <- rbinom(length(infected.B),1,msm.aq.prob.B)
    infected.BI <- rbinom(length(infected.BI),1,msm.aq.prob.BI)
    infected.H <- rbinom(length(infected.H),1,msm.aq.prob.H)
    infected.HI <- rbinom(length(infected.HI),1,msm.aq.prob.HI)
    infected.W <- rbinom(length(infected.W),1,msm.aq.prob.W)
    
    infected.B <- ids.B[infected.B==1]
    infected.BI <- ids.BI[infected.BI==1]
    infected.H <- ids.H[infected.H==1]
    infected.HI <- ids.HI[infected.HI==1]
    infected.W <- ids.W[infected.W==1]


    infected<-c(as.numeric(infected.B),as.numeric(infected.BI),as.numeric(infected.H),as.numeric(infected.HI),as.numeric(infected.W))
    
 if (length(infected) >= 1){
    dat$attr$status[infected] <- 1
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- "AR"
    dat$attr$stage.time[infected] <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected] <- 0
    dat$attr$inf.class[infected] <- "MSM"
      
    
    dat$attr$infector[infected] <- "MSM"
    dat$attr$inf.role[infected] <- "MSM"
    dat$attr$inf.type[infected] <- "MSM"
    dat$attr$inf.diag[infected] <- "MSM" 
    dat$attr$inf.tx[infected] <- "MSM"
    dat$attr$inf.stage[infected] <- "MSM"
    
    dat$attr$cum.time.on.tx[infected] <- 0
    dat$attr$cum.time.off.tx[infected] <- 0
    
    
    # Summary Output
    dat$epi$incid[at] <- sum(dat$epi$incid[at] , length(infected), na.rm=TRUE)
    
    dat$epi$incid.B[at] <- sum(dat$epi$incid.B[at] , sum(race[infected] == "B", na.rm = TRUE), na.rm = TRUE)
    dat$epi$incid.BI[at] <- sum(dat$epi$incid.BI[at] , sum(race[infected] == "BI", na.rm = TRUE), na.rm = TRUE)
    dat$epi$incid.H[at] <- sum(dat$epi$incid.H[at] , sum(race[infected] == "H", na.rm = TRUE), na.rm = TRUE)
    dat$epi$incid.HI[at] <- sum(dat$epi$incid.HI[at] , sum(race[infected] == "HI", na.rm = TRUE), na.rm = TRUE)
    dat$epi$incid.W[at] <- sum(dat$epi$incid.W[at] , sum(race[infected] == "W", na.rm = TRUE), na.rm = TRUE)
  
    dat$epi$incid.m[at] <- sum(dat$epi$incid.m[at] , sum(sex[infected] == "M", na.rm = TRUE), na.rm = TRUE)
    
    dat$epi$incid.msf[at] <- sum(dat$epi$incid.msf[at] , sum(sex.ident[infected] == "msf", na.rm = TRUE), na.rm = TRUE)
    dat$epi$incid.msm[at] <- sum(dat$epi$incid.msm[at] , sum(sex.ident[infected] == "msm", na.rm = TRUE), na.rm = TRUE)
    dat$epi$incid.msmf[at] <- sum(dat$epi$incid.msmf[at] , sum(sex.ident[infected] == "msmf", na.rm = TRUE), na.rm = TRUE)
    
    dat$epi$incid.MSM[at] <- length(infected)
     }
    
    if (length(infected) < 1){ dat$epi$incid.H[at] <-0}
  
  
    
  return(dat)
}





