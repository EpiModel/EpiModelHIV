
#' @title Prevalence Calculations within Time Steps for up to 5 race groups heterosexuals and MSM.
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging_shamp
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
  diagnosed <- dat$attr$diag.status

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  if (at == 1) {
    dat$epi$num <- rNA

    dat$epi$num.f <- rNA
    dat$epi$num.m <- rNA

    

    dat$epi$s.num <- rNA
    dat$epi$i.num <- rNA
    

    dat$epi$i.num.f <- rNA 
    dat$epi$i.num.m <- rNA 


    
    dat$epi$prev <- rNA
    dat$epi$prev.m <- rNA
    dat$epi$prev.f <- rNA


    dat$epi$nBirths <- rNA
    dat$epi$dth.gen <- rNA
    dat$epi$dth.dis <- rNA
    dat$epi$dth.age <- rNA
    dat$epi$dth.age.pos <- rNA
    dat$epi$dth.age.pos.ontx <- rNA
    dat$epi$incid <- rNA
    dat$epi$incid.m <- rNA
    dat$epi$incid.f <- rNA


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
    
    dat$epi$duplicates <-  rep(0, nsteps)
    dat$epi$duplicates.w.OT <-  rep(0, nsteps)
    
    dat$epi$nCohabs <-  rep(0, nsteps)
    dat$epi$nPers <-  rep(0, nsteps)
    dat$epi$nOT <-  rep(0, nsteps)
    
    dat$epi$diag.prevalent <- rep(0, nsteps)
    dat$epi$diag.acute <- rep(0, nsteps)
    dat$epi$n.tests <- rep(0, nsteps)
    dat$epi$n.tests.ab <- rep(0, nsteps)
    dat$epi$n.tests.rna <- rep(0, nsteps)
    
    dat$epi$partners.sought.new <- rep(0, nsteps)
    dat$epi$partners.found <- rep(0, nsteps)
    dat$epi$partners.positive <- rep(0, nsteps)
    dat$epi$partners.negative <- rep(0, nsteps)
    
  }


  dat$epi$num[at] <- length(status)

  dat$epi$num.f[at] <- sum(sex == "F", na.rm = TRUE)
  dat$epi$num.m[at] <- sum(sex == "M", na.rm = TRUE)

  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)

  dat$epi$i.num.f[at] <- sum(status == 1 & sex == "F", na.rm = TRUE)
  dat$epi$i.num.m[at] <- sum(status == 1 & sex == "M", na.rm = TRUE)


  
  dat$epi$prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]

  dat$epi$prev.f[at] <- dat$epi$i.num.f[at] / dat$epi$num.f[at]
  dat$epi$prev.m[at] <- dat$epi$i.num.m[at] / dat$epi$num.m[at]
  
  
  
  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepElig[at] <- sum(dat$attr$prepElig == 1, na.rm = TRUE)
  dat$epi$i.num.prep0[at] <- sum((is.na(prepStat) | prepStat == 0) &
                                 status == 1, na.rm = TRUE)
  dat$epi$i.num.prep1[at] <- sum(prepStat == 1 & status == 1, na.rm = TRUE)
  dat$epi$prev.prep0[at] <- dat$epi$i.num.prep0[at] /
                              sum((is.na(prepStat) | prepStat == 0), na.rm = TRUE)
  if (at == 1) {
    dat$epi$prev.prep1[1] <- 0
  } else {
    dat$epi$prev.prep1[at] <- dat$epi$i.num.prep1[at] /
                                sum(prepStat == 1, na.rm = TRUE)
  }
  
  
  dat$epi$nCohabs[at] <- sum(get_degree(dat$el[[1]]), na.rm = TRUE)/2
  dat$epi$nPers[at] <- sum(get_degree(dat$el[[2]]), na.rm = TRUE)/2
  dat$epi$nOT[at] <- sum(get_degree(dat$el[[3]]), na.rm = TRUE)/2
  
  
  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}
