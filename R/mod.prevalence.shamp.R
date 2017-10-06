
#' @title Prevalence Calculations within Time Steps for up to 5 race groups heterosexuals and MSM.
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Summary statistic calculations are of two broad forms: prevalence and
#' incidence. This function establishes the summary statistic vectors for both
#' prevalence and incidence at time 1, and then calculates the prevalence
#' statistics for times 2 onward. Incidence statistics (e.g., number of new
#' infections or deaths) are calculated within the modules as they depend on
#' vectors that are not stored external to the module.
#'
#' @return
#' This function returns the \code{dat} object with an updated summary of current
#' attributes stored in \code{dat$epi}.
#'
#' @keywords module SHAMP msm het
#'
#' @export
#'
prevalence_shamp <- function(dat, at) {

  race <- dat$attr$race
  status <- dat$attr$status
  sex.ident<-dat$attr$sex.ident
  sex<-dat$attr$sex
  prepStat <- dat$attr$prepStat
  inf.class <- dat$attr$inf.class

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  if (at == 1) {
    dat$epi$num <- rNA

    dat$epi$num.B <- rNA
    dat$epi$num.BI <- rNA
    dat$epi$num.H <- rNA
    dat$epi$num.HI <- rNA
    dat$epi$num.W <- rNA
    dat$epi$num.f <- rNA
    dat$epi$num.m <- rNA
    dat$epi$num.msf <- rNA  
    dat$epi$num.msm <- rNA
    dat$epi$num.msmf <- rNA
    
    dat$epi$num.B.f <- rNA
    dat$epi$num.BI.f <- rNA
    dat$epi$num.H.f <- rNA
    dat$epi$num.HI.f <- rNA
    dat$epi$num.W.f <- rNA
    
    dat$epi$num.B.m <- rNA
    dat$epi$num.BI.m <- rNA
    dat$epi$num.H.m <- rNA
    dat$epi$num.HI.m <- rNA
    dat$epi$num.W.m <- rNA
    
    dat$epi$num.B.msf <- rNA
    dat$epi$num.BI.msf <- rNA
    dat$epi$num.H.msf <- rNA
    dat$epi$num.HI.msf <- rNA
    dat$epi$num.W.msf <- rNA
    
    dat$epi$num.B.msm <- rNA
    dat$epi$num.BI.msm <- rNA
    dat$epi$num.H.msm <- rNA
    dat$epi$num.HI.msm <- rNA
    dat$epi$num.W.msm <- rNA
    
    dat$epi$num.B.msmf <- rNA
    dat$epi$num.BI.msmf <- rNA
    dat$epi$num.H.msmf <- rNA
    dat$epi$num.HI.msmf <- rNA
    dat$epi$num.W.msmf <- rNA
    
    dat$epi$s.num <- rNA
    dat$epi$i.num <- rNA
    
    dat$epi$i.num.B <- rNA
    dat$epi$i.num.BI <- rNA
    dat$epi$i.num.H <- rNA
    dat$epi$i.num.HI <- rNA
    dat$epi$i.num.W <- rNA 
    dat$epi$i.num.f <- rNA 
    dat$epi$i.num.m <- rNA 
    dat$epi$i.num.msf <- rNA 
    dat$epi$i.num.msm <- rNA 
    dat$epi$i.num.msmf <- rNA 
    
    dat$epi$i.num.B.f <- rNA
    dat$epi$i.num.BI.f <- rNA
    dat$epi$i.num.H.f <- rNA
    dat$epi$i.num.HI.f <- rNA
    dat$epi$i.num.W.f <- rNA
    
    dat$epi$i.num.B.m <- rNA
    dat$epi$i.num.BI.m <- rNA
    dat$epi$i.num.H.m <- rNA
    dat$epi$i.num.HI.m <- rNA
    dat$epi$i.num.W.m <- rNA
    
    dat$epi$i.num.B.msm <- rNA
    dat$epi$i.num.BI.msm <- rNA
    dat$epi$i.num.H.msm <- rNA
    dat$epi$i.num.HI.msm <- rNA
    dat$epi$i.num.W.msm <- rNA
    
    dat$epi$i.prev <- rNA
    dat$epi$i.prev.B <- rNA
    dat$epi$i.prev.BI <- rNA
    dat$epi$i.prev.H <- rNA
    dat$epi$i.prev.HI <- rNA
    dat$epi$i.prev.W <- rNA
    dat$epi$i.prev.m <- rNA
    dat$epi$i.prev.f <- rNA
    dat$epi$i.prev.msf <- rNA
    dat$epi$i.prev.msm <- rNA
    dat$epi$i.prev.msmf <- rNA
    
    dat$epi$i.prev.B.f <- rNA
    dat$epi$i.prev.BI.f <- rNA
    dat$epi$i.prev.H.f <- rNA
    dat$epi$i.prev.HI.f <- rNA
    dat$epi$i.prev.W.f <- rNA
    
    dat$epi$i.prev.B.m <- rNA
    dat$epi$i.prev.BI.m <- rNA
    dat$epi$i.prev.H.m <- rNA
    dat$epi$i.prev.HI.m <- rNA
    dat$epi$i.prev.W.m <- rNA
    
    dat$epi$i.prev.B.msf <- rNA
    dat$epi$i.prev.BI.msf <- rNA
    dat$epi$i.prev.H.msf <- rNA
    dat$epi$i.prev.HI.msf <- rNA
    dat$epi$i.prev.W.msf <- rNA
    
    dat$epi$i.prev.B.msm <- rNA
    dat$epi$i.prev.BI.msm <- rNA
    dat$epi$i.prev.H.msm <- rNA
    dat$epi$i.prev.HI.msm <- rNA
    dat$epi$i.prev.W.msm <- rNA
    
    dat$epi$i.prev.B.msmf <- rNA
    dat$epi$i.prev.BI.msmf <- rNA
    dat$epi$i.prev.H.msmf <- rNA
    dat$epi$i.prev.HI.msmf <- rNA
    dat$epi$i.prev.W.msmf <- rNA
    
    dat$epi$nBirths <- rNA
    dat$epi$dth.gen <- rNA
    dat$epi$dth.dis <- rNA
    dat$epi$dth.age <- rNA
    dat$epi$incid <- rNA
    
    
    dat$epi$incid.B <- rNA
    dat$epi$incid.BI <- rNA
    dat$epi$incid.H <- rNA
    dat$epi$incid.HI <- rNA
    dat$epi$incid.W <- rNA
    
    dat$epi$incid.m <- rNA
    dat$epi$incid.f <- rNA
    dat$epi$incid.msf <- rNA
    dat$epi$incid.msm <- rNA
    dat$epi$incid.msmf <- rNA
    
    dat$epi$incid.FA <- rep(0, nsteps)
    dat$epi$incid.MSM <- rep(0, nsteps)
    dat$epi$incid.Lhet <- rep(0, nsteps)
    
    dat$epi$prepCurr <- rNA
    dat$epi$prepCov <- rNA
    dat$epi$prepElig <- rNA
    dat$epi$prepStart <- rNA
    dat$epi$incid.prep0 <- rNA
    dat$epi$incid.prep1 <- rNA
    dat$epi$i.num.prep0 <- rNA
    dat$epi$i.num.prep1 <- rNA

    dat$epi$cprob.always.pers <- rNA
    dat$epi$cprob.always.inst <- rNA
    
    dat$epi$i.num.MSM.inf <- rep(0, nsteps)
    dat$epi$i.num.MSMds.inf <- rep(0, nsteps)
    dat$epi$i.num.FA.inf <- rep(0, nsteps)
    dat$epi$i.num.FAds.inf <- rep(0, nsteps)
    dat$epi$i.num.Lhet.inf <- rep(0, nsteps)
    
    dat$epi$prop.MSM.inf <- rep(0, nsteps)
    dat$epi$prop.MSMds.inf <- rep(0, nsteps)
    dat$epi$prop.FA.inf <- rep(0, nsteps)
    dat$epi$prop.FAds.inf <- rep(0, nsteps)
    dat$epi$prop.Lhet.inf <- rep(0, nsteps)
  }


  dat$epi$num[at] <- length(status)
  dat$epi$num.B[at] <- sum(race == "B", na.rm = TRUE)
  dat$epi$num.BI[at] <- sum(race == "BI", na.rm = TRUE)
  dat$epi$num.H[at] <- sum(race == "H", na.rm = TRUE)
  dat$epi$num.HI[at] <- sum(race == "HI", na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == "W", na.rm = TRUE)
  dat$epi$num.f[at] <- sum(sex == "F", na.rm = TRUE)
  dat$epi$num.m[at] <- sum(sex == "M", na.rm = TRUE)
  dat$epi$num.msf[at] <- sum(sex.ident == "msf", na.rm = TRUE)
  dat$epi$num.msm[at] <- sum(sex.ident == "msm", na.rm = TRUE)
  dat$epi$num.msmf[at] <- sum(sex.ident == "msmf", na.rm = TRUE)

  
  dat$epi$num.B.f[at] <- sum(race == "B" & sex== "F", na.rm = TRUE)
  dat$epi$num.BI.f[at] <- sum(race == "BI" & sex== "F", na.rm = TRUE)
  dat$epi$num.H.f[at] <- sum(race == "H" & sex== "F", na.rm = TRUE)
  dat$epi$num.HI.f[at] <- sum(race == "HI" & sex== "F", na.rm = TRUE)
  dat$epi$num.W.f[at] <- sum(race == "W" & sex== "F", na.rm = TRUE)
  
  dat$epi$num.B.m[at] <- sum(race == "B" & sex== "M", na.rm = TRUE)
  dat$epi$num.BI.m[at] <- sum(race == "BI" & sex== "M", na.rm = TRUE)
  dat$epi$num.H.m[at] <- sum(race == "H" & sex== "M", na.rm = TRUE)
  dat$epi$num.HI.m[at] <- sum(race == "HI" & sex== "M", na.rm = TRUE)
  dat$epi$num.W.m[at] <- sum(race == "W" & sex== "M", na.rm = TRUE)
  
  dat$epi$num.B.msf[at] <- sum(race == "B" & sex.ident== "msf", na.rm = TRUE)
  dat$epi$num.BI.msf[at] <- sum(race == "BI" & sex.ident== "msf", na.rm = TRUE)
  dat$epi$num.H.msf[at] <- sum(race == "H" & sex.ident== "msf", na.rm = TRUE)
  dat$epi$num.HI.msf[at] <- sum(race == "HI" & sex.ident== "msf", na.rm = TRUE)
  dat$epi$num.W.msf[at] <- sum(race == "W" & sex.ident== "msf", na.rm = TRUE)
  
  dat$epi$num.B.msm[at] <- sum(race == "B" & sex.ident== "msm", na.rm = TRUE)
  dat$epi$num.BI.msm[at] <- sum(race == "BI" & sex.ident== "msm", na.rm = TRUE)
  dat$epi$num.H.msm[at] <- sum(race == "H" & sex.ident== "msm", na.rm = TRUE)
  dat$epi$num.HI.msm[at] <- sum(race == "HI" & sex.ident== "msm", na.rm = TRUE)
  dat$epi$num.W.msm[at] <- sum(race == "W" & sex.ident== "msm", na.rm = TRUE)
  
  dat$epi$num.B.msmf[at] <- sum(race == "B" & sex.ident== "msmf", na.rm = TRUE)
  dat$epi$num.BI.msmf[at] <- sum(race == "BI" & sex.ident== "msmf", na.rm = TRUE)
  dat$epi$num.H.msmf[at] <- sum(race == "H" & sex.ident== "msmf", na.rm = TRUE)
  dat$epi$num.HI.msmf[at] <- sum(race == "HI" & sex.ident== "msmf", na.rm = TRUE)
  dat$epi$num.W.msmf[at] <- sum(race == "W" & sex.ident== "msmf", na.rm = TRUE)

  
  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  
  dat$epi$i.num.B[at] <- sum(status == 1 & race == "B", na.rm = TRUE)
  dat$epi$i.num.BI[at] <- sum(status == 1 & race == "BI", na.rm = TRUE)
  dat$epi$i.num.H[at] <- sum(status == 1 & race == "H", na.rm = TRUE)
  dat$epi$i.num.HI[at] <- sum(status == 1 & race == "HI", na.rm = TRUE)
  dat$epi$i.num.W[at] <- sum(status == 1 & race == "W", na.rm = TRUE)
  dat$epi$i.num.f[at] <- sum(status == 1 & sex == "F", na.rm = TRUE)
  dat$epi$i.num.m[at] <- sum(status == 1 & sex == "M", na.rm = TRUE)
  dat$epi$i.num.msf[at] <- sum(status == 1 & sex.ident == "msf", na.rm = TRUE)
  dat$epi$i.num.msm[at] <- sum(status == 1 & sex.ident == "msm", na.rm = TRUE)
  dat$epi$i.num.msmf[at] <- sum(status == 1 & sex.ident == "msmf", na.rm = TRUE)
  
  dat$epi$i.num.B.f[at] <- sum(status == 1 & race == "B" & sex == "F", na.rm = TRUE)
  dat$epi$i.num.BI.f[at] <- sum(status == 1 & race == "BI" & sex == "F", na.rm = TRUE)
  dat$epi$i.num.H.f[at] <- sum(status == 1 & race == "H" & sex == "F", na.rm = TRUE)
  dat$epi$i.num.HI.f[at] <- sum(status == 1 & race == "HI" & sex == "F", na.rm = TRUE)
  dat$epi$i.num.W.f[at] <- sum(status == 1 & race == "W" & sex == "F", na.rm = TRUE)
  
  dat$epi$i.num.B.m[at] <- sum(status == 1 & race == "B" & sex == "M", na.rm = TRUE)
  dat$epi$i.num.BI.m[at] <- sum(status == 1 & race == "BI" & sex == "M", na.rm = TRUE)
  dat$epi$i.num.H.m[at] <- sum(status == 1 & race == "H" & sex == "M", na.rm = TRUE)
  dat$epi$i.num.HI.m[at] <- sum(status == 1 & race == "HI" & sex == "M", na.rm = TRUE)
  dat$epi$i.num.W.m[at] <- sum(status == 1 & race == "W" & sex == "M", na.rm = TRUE) 
  
  dat$epi$i.num.B.msf[at] <- sum(status == 1 & race == "B" & sex.ident== "msf", na.rm = TRUE)
  dat$epi$i.num.BI.msf[at] <- sum(status == 1 & race == "BI" & sex.ident== "msf", na.rm = TRUE)
  dat$epi$i.num.H.msf[at] <- sum(status == 1 & race == "H" & sex.ident== "msf", na.rm = TRUE)
  dat$epi$i.num.HI.msf[at] <- sum(status == 1 & race == "HI" & sex.ident== "msf", na.rm = TRUE)
  dat$epi$i.num.W.msf[at] <- sum(status == 1 & race == "W" & sex.ident== "msf", na.rm = TRUE)
  
  dat$epi$i.num.B.msm[at] <- sum(status == 1 & race == "B" & sex.ident== "msm", na.rm = TRUE)
  dat$epi$i.num.BI.msm[at] <- sum(status == 1 & race == "BI" & sex.ident== "msm", na.rm = TRUE)
  dat$epi$i.num.H.msm[at] <- sum(status == 1 & race == "H" & sex.ident== "msm", na.rm = TRUE)
  dat$epi$i.num.HI.msm[at] <- sum(status == 1 & race == "HI" & sex.ident== "msm", na.rm = TRUE)
  dat$epi$i.num.W.msm[at] <- sum(status == 1 & race == "W" & sex.ident== "msm", na.rm = TRUE)
  
  dat$epi$i.num.B.msmf[at] <- sum(status == 1 & race == "B" & sex.ident== "msmf", na.rm = TRUE)
  dat$epi$i.num.BI.msmf[at] <- sum(status == 1 & race == "BI" & sex.ident== "msmf", na.rm = TRUE)
  dat$epi$i.num.H.msmf[at] <- sum(status == 1 & race == "H" & sex.ident== "msmf", na.rm = TRUE)
  dat$epi$i.num.HI.msmf[at] <- sum(status == 1 & race == "HI" & sex.ident== "msmf", na.rm = TRUE)
  dat$epi$i.num.W.msmf[at] <- sum(status == 1 & race == "W" & sex.ident== "msmf", na.rm = TRUE) 
  
  dat$epi$i.num.MSM.inf[at] <- sum(inf.class=="MSM", na.rm = TRUE)
  dat$epi$i.num.MSMds.inf[at] <- sum(inf.class=="MSMds", na.rm = TRUE)
  dat$epi$i.num.FA.inf[at] <- sum(inf.class=="FA", na.rm = TRUE)
  dat$epi$i.num.FAds.inf[at] <- sum(inf.class=="FAds", na.rm = TRUE)
  dat$epi$i.num.Lhet.inf[at] <- sum(inf.class=="Lhet", na.rm = TRUE)
  
  dat$epi$prop.MSM.inf[at] <- dat$epi$i.num.MSM.inf[at]/dat$epi$i.num[at]
  dat$epi$prop.MSMds.inf[at] <- dat$epi$i.num.MSMds.inf[at]/dat$epi$i.num[at]
  dat$epi$prop.FA.inf[at] <- dat$epi$i.num.FA.inf[at]/dat$epi$i.num[at]
  dat$epi$prop.FAds.inf[at] <- dat$epi$i.num.FAds.inf[at]/dat$epi$i.num[at]
  dat$epi$prop.Lhet.inf[at] <- dat$epi$i.num.Lhet.inf[at]/dat$epi$i.num[at]
  
  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.B[at] <- dat$epi$i.num.B[at] / dat$epi$num.B[at]
  dat$epi$i.prev.BI[at] <- dat$epi$i.num.BI[at] / dat$epi$num.BI[at]
  dat$epi$i.prev.H[at] <- dat$epi$i.num.H[at] / dat$epi$num.H[at]
  dat$epi$i.prev.HI[at] <- dat$epi$i.num.HI[at] / dat$epi$num.HI[at]
  dat$epi$i.prev.W[at] <- dat$epi$i.num.W[at] / dat$epi$num.W[at]
  dat$epi$i.prev.f[at] <- dat$epi$i.num.f[at] / dat$epi$num.f[at]
  dat$epi$i.prev.m[at] <- dat$epi$i.num.m[at] / dat$epi$num.m[at]
  dat$epi$i.prev.msf[at] <- dat$epi$i.num.msf[at] / dat$epi$num.msf[at]  
  dat$epi$i.prev.msm[at] <- dat$epi$i.num.msm[at] / dat$epi$num.msm[at]
  dat$epi$i.prev.msmf[at] <- dat$epi$i.num.msmf[at] / dat$epi$num.msmf[at]
  
  dat$epi$i.prev.B.f[at] <- dat$epi$i.num.B.f[at] / dat$epi$num.B.f[at]
  dat$epi$i.prev.BI.f[at] <- dat$epi$i.num.BI.f[at] / dat$epi$num.BI.f[at]
  dat$epi$i.prev.H.f[at] <- dat$epi$i.num.H.f[at] / dat$epi$num.H.f[at]
  dat$epi$i.prev.HI.f[at] <- dat$epi$i.num.HI.f[at] / dat$epi$num.HI.f[at]
  dat$epi$i.prev.W.f[at] <- dat$epi$i.num.W.f[at] / dat$epi$num.W.f[at]
  
  dat$epi$i.prev.B.m[at] <- dat$epi$i.num.B.m[at] / dat$epi$num.B.m[at]
  dat$epi$i.prev.BI.m[at] <- dat$epi$i.num.BI.m[at] / dat$epi$num.BI.m[at]
  dat$epi$i.prev.H.m[at] <- dat$epi$i.num.H.m[at] / dat$epi$num.H.m[at]
  dat$epi$i.prev.HI.m[at] <- dat$epi$i.num.HI.m[at] / dat$epi$num.HI.m[at]
  dat$epi$i.prev.W.m[at] <- dat$epi$i.num.W.m[at] / dat$epi$num.W.m[at]
  
  dat$epi$i.prev.B.msf[at] <- dat$epi$i.num.B.msf[at] / dat$epi$num.B.msf[at]
  dat$epi$i.prev.BI.msf[at] <- dat$epi$i.num.BI.msf[at] / dat$epi$num.BI.msf[at]
  dat$epi$i.prev.H.msf[at] <- dat$epi$i.num.H.msf[at] / dat$epi$num.H.msf[at]
  dat$epi$i.prev.HI.msf[at] <- dat$epi$i.num.HI.msf[at] / dat$epi$num.HI.msf[at]
  dat$epi$i.prev.W.msf[at] <- dat$epi$i.num.W.msf[at] / dat$epi$num.W.msf[at]
  
  dat$epi$i.prev.B.msm[at] <- dat$epi$i.num.B.msm[at] / dat$epi$num.B.msm[at]
  dat$epi$i.prev.BI.msm[at] <- dat$epi$i.num.BI.msm[at] / dat$epi$num.BI.msm[at]
  dat$epi$i.prev.H.msm[at] <- dat$epi$i.num.H.msm[at] / dat$epi$num.H.msm[at]
  dat$epi$i.prev.HI.msm[at] <- dat$epi$i.num.HI.msm[at] / dat$epi$num.HI.msm[at]
  dat$epi$i.prev.W.msm[at] <- dat$epi$i.num.W.msm[at] / dat$epi$num.W.msm[at]
  
  dat$epi$i.prev.B.msmf[at] <- dat$epi$i.num.B.msmf[at] / dat$epi$num.B.msmf[at]
  dat$epi$i.prev.BI.msmf[at] <- dat$epi$i.num.BI.msmf[at] / dat$epi$num.BI.msmf[at]
  dat$epi$i.prev.H.msmf[at] <- dat$epi$i.num.H.msmf[at] / dat$epi$num.H.msmf[at]
  dat$epi$i.prev.HI.msmf[at] <- dat$epi$i.num.HI.msmf[at] / dat$epi$num.HI.msmf[at]
  dat$epi$i.prev.W.msmf[at] <- dat$epi$i.num.W.msmf[at] / dat$epi$num.W.msmf[at]
  


  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepElig[at] <- sum(dat$attr$prepElig == 1, na.rm = TRUE)
  dat$epi$i.num.prep0[at] <- sum((is.na(prepStat) | prepStat == 0) &
                                 status == 1, na.rm = TRUE)
  dat$epi$i.num.prep1[at] <- sum(prepStat == 1 & status == 1, na.rm = TRUE)
  dat$epi$i.prev.prep0[at] <- dat$epi$i.num.prep0[at] /
                              sum((is.na(prepStat) | prepStat == 0), na.rm = TRUE)
  if (at == 1) {
    dat$epi$i.prev.prep1[1] <- 0
  } else {
    dat$epi$i.prev.prep1[at] <- dat$epi$i.num.prep1[at] /
                                sum(prepStat == 1, na.rm = TRUE)
  }

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}
