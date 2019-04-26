
#' @title Prevalence Calculations within Time Steps
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging_camplc
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
#' @keywords module msm
#'
#' @export
#'
prevalence_msm <- function(dat, at) {

  active <- dat$attr$active
  race <- dat$attr$race
  status <- dat$attr$status
  prepStat <- dat$attr$prepStat
  debuted <- dat$attr$debuted
  asmm <- dat$attr$asmm
  age <- floor(dat$attr$age)
  everAI <- dat$attr$everAI
  cond.int.active <- dat$attr$cond.int.active
  
  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  if (at == 1) {
    dat$epi$num <- rNA
    dat$epi$num.B <- rNA
    dat$epi$num.W <- rNA
    dat$epi$num.asmm <- rNA
    dat$epi$num.msm <- rNA
    dat$epi$num.age18 <- rNA
    dat$epi$num.deb <- rNA
    dat$epi$num.asmm.deb <- rNA
    dat$epi$num.age18.deb <- rNA
    dat$epi$num.asmm.everAI <- rNA
    dat$epi$s.num <- rNA
    dat$epi$i.num <- rNA
    dat$epi$i.num.B <- rNA
    dat$epi$i.num.W <- rNA
    dat$epi$i.num.msm <- rNA
    dat$epi$i.num.asmm <- rNA
    dat$epi$i.num.age18 <- rNA
    dat$epi$i.prev <- rNA
    dat$epi$i.prev.B <- rNA
    dat$epi$i.prev.W <- rNA
    dat$epi$i.prev.msm <- rNA
    dat$epi$i.prev.asmm <- rNA
    dat$epi$i.prev.age18 <- rNA
    dat$epi$incid <- rNA
    dat$epi$incid.msm <- rNA
    dat$epi$incid.asmm <- rNA

    dat$epi$debuted <- rNA
    dat$epi$debuted.asmm <- rNA
    
    dat$epi$prepCurr <- rNA
    dat$epi$prepCurr.msm <- rNA
    dat$epi$prepCurr.asmm <- rNA
    dat$epi$prepCurr.ret <- rNA
    dat$epi$prepEver <- rNA
    dat$epi$prepCov.msm.elig <- rNA
    dat$epi$prepCov.msm.elig.w.ret <- rNA 
    dat$epi$prepCov.msm.all <- rNA
    dat$epi$prepCov.msm.all.w.ret <- rNA
    dat$epi$prepCov.asmm <- rNA
    dat$epi$prepCov.elig.asmm <- rNA 
    dat$epi$prepCov.all.asmm <- rNA 
    dat$epi$prepCov.deb.asmm <- rNA 
    dat$epi$prepCov.everAI.asmm <- rNA
    dat$epi$prepCov.16.AI.asmm <- rNA
    dat$epi$prepElig <- rNA
    dat$epi$prepStart <- rNA
    dat$epi$prepStart.msm <- rNA
    dat$epi$prepStart.asmm <- rNA
    dat$epi$i.num.prep0 <- rNA
    dat$epi$i.num.prep1 <- rNA

    dat$epi$cprob.always.pers <- rNA
    dat$epi$cprob.always.inst <- rNA

    dat$epi$incid.msm <- rNA
    dat$epi$incid.asmm <- rNA
    
    ##condom intervention
    dat$epi$num.cond.int.active <- rNA
    dat$epi$num.cond.int.active.asmm <- rNA
    dat$epi$num.cond.int.active.asmm.AI <- rNA
    dat$epi$num.cond.int.active.msm <- rNA
    
    dat$epi$prev.cond.int.active <- rNA
    dat$epi$prev.cond.int.active.asmm <- rNA
    dat$epi$prev.cond.int.active.asmm.AI <- rNA
    dat$epi$prev.cond.int.active.msm <- rNA
    
    dat$epi$incid.cond.int.start <- rNA
    dat$epi$incid.cond.int.stop <- rNA
    
  }


  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$num.B[at] <- sum(active == 1 & race == "B", na.rm = TRUE)
  dat$epi$num.W[at] <- sum(active == 1 & race == "W", na.rm = TRUE)
  dat$epi$num.asmm[at] <- sum(active == 1 & asmm == 1, na.rm = TRUE)
  dat$epi$num.msm[at] <- sum(active == 1 & asmm == 0, na.rm = TRUE)
  dat$epi$num.age18[at] <- sum(active == 1 & age == 18, na.rm = TRUE)
  dat$epi$num.deb[at] <- sum(active == 1 & debuted == 1, na.rm = TRUE)
  dat$epi$debuted[at] <- sum(active == 1 & debuted == 1, na.rm = TRUE)
  dat$epi$debuted.asmm[at] <- sum(active == 1 & debuted == 1 & asmm == 1, na.rm = TRUE)
  dat$epi$num.asmm.deb[at] <- sum(active == 1 & debuted == 1 & asmm == 1, na.rm = TRUE)
  dat$epi$num.age18.deb[at] <- sum(active == 1 & debuted ==1 & age == 18, na.rm = TRUE)
  dat$epi$num.asmm.everAI[at] <- sum(active == 1 & everAI == 1 & asmm == 1, na.rm = TRUE)
  dat$epi$s.num[at] <- sum(active == 1 & status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(active == 1 & status == 1, na.rm = TRUE)
  dat$epi$i.num.B[at] <- sum(active == 1 & status == 1 & race == "B", na.rm = TRUE)
  dat$epi$i.num.W[at] <- sum(active == 1 & status == 1 & race == "W", na.rm = TRUE)
  dat$epi$i.num.msm[at] <- sum(active == 1 & status == 1 & asmm == 0, na.rm = TRUE)
  dat$epi$i.num.asmm[at] <- sum(active == 1 & status == 1 & asmm == 1, na.rm = TRUE)
  dat$epi$i.num.age18[at] <- sum(active == 1 & status == 1 & age == 18, na.rm = TRUE)
  dat$epi$i.num.asmm.deb[at] <- sum(active == 1 & debuted == 1 & status == 1 & asmm == 1, na.rm = TRUE)
  dat$epi$i.num.age18.deb[at] <- sum(active == 1 & debuted ==1 & status == 1 & age == 18, na.rm = TRUE)
  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.B[at] <- dat$epi$i.num.B[at] / dat$epi$num.B[at]
  dat$epi$i.prev.W[at] <- dat$epi$i.num.W[at] / dat$epi$num.W[at]

  dat$epi$prepCurr[at] <- sum(active == 1 & prepStat == 1, na.rm = TRUE)
  dat$epi$prepCurr.msm[at] <- sum(active == 1 & prepStat == 1 & asmm == 0, na.rm = TRUE)
  dat$epi$prepCurr.asmm[at] <- sum(active == 1 & prepStat == 1 & asmm == 1, na.rm = TRUE)
  dat$epi$prepCurr.ret[at] <- sum(active == 1 & dat$attr$prepStat.asmm == 1 & asmm == 0, na.rm = TRUE)
  dat$epi$prepElig.msm[at] <- sum(active == 1 & dat$attr$prepElig == 1, na.rm = TRUE)
  dat$epi$prepElig.asmm[at] <- sum(active == 1 & dat$attr$prepElig.asmm == 1, na.rm = TRUE)
  dat$epi$prepEver[at] <- sum(active == 1 & dat$attr$prepEver == 1, na.rm = TRUE)
  dat$epi$i.num.prep0[at] <- sum(active == 1 & (is.na(prepStat) | prepStat == 0) & status == 1, na.rm = TRUE)
  dat$epi$i.num.prep1[at] <- sum(active == 1 & prepStat == 1 & status == 1, na.rm = TRUE)
  dat$epi$i.prev.prep0[at] <- dat$epi$i.num.prep0[at] /
    sum(active == 1 & (is.na(prepStat) | prepStat == 0), na.rm = TRUE)
  if (at == 1) {
    dat$epi$i.prev.prep1[1] <- 0
  } else {
    dat$epi$i.prev.prep1[at] <- dat$epi$i.num.prep1[at] / sum(active == 1 & prepStat == 1, na.rm = TRUE)
  }
  
  dat$epi$prepCov.elig.asmm[at] <- (dat$epi$prepCurr.asmm[at] / dat$epi$prepElig.asmm[at])
  dat$epi$prepCov.16.AI.asmm[at] <- (dat$epi$prepCurr.asmm[at] / sum(active == 1 & asmm == 1 & everAI ==1 & age >= 16, na.rm = TRUE ))
  dat$epi$prepCov.all.asmm[at] <- (dat$epi$prepCurr.asmm[at] / dat$epi$num.asmm[at])
  dat$epi$prepCov.deb.asmm[at] <- (dat$epi$prepCurr.asmm[at] / dat$epi$num.asmm.deb[at])
  dat$epi$prepCov.everAI.asmm[at] <- (dat$epi$prepCurr.asmm[at] / dat$epi$num.asmm.everAI[at])
  
  
  dat$epi$i.prev.msm[at] <- sum(active == 1 & status ==1 & asmm == 0, na.rm = TRUE) / dat$epi$num.msm[at]
  dat$epi$i.prev.asmm[at] <- sum(active == 1 & status ==1 & debuted == 1 & asmm == 1, na.rm = TRUE) / dat$epi$num.asmm.deb[at]
  dat$epi$i.prev.age18[at] <- sum(active == 1 & status ==1 & debuted == 1 & age == 18, na.rm = TRUE) / dat$epi$num.age18.deb[at]
  dat$epi$prepStart[at] <- dat$epi$prepStart.msm[at] + dat$epi$prepStart.asmm[at]
  
  dat$epi$num.cond.int.active[at] <- sum(active == 1 & cond.int.active == 1, na.rm = TRUE)
  dat$epi$num.cond.int.active.asmm[at] <- sum(active == 1 & cond.int.active ==1 & asmm == 1, na.rm = TRUE)
  dat$epi$num.cond.int.active.asmm.AI[at] <- sum(active == 1 & cond.int.active ==1 & everAI ==1 & asmm == 1, na.rm = TRUE)
  dat$epi$num.cond.int.active.msm[at] <- sum(active == 1 & cond.int.active ==1 & asmm == 0, na.rm = TRUE)
  
  dat$epi$prev.cond.int.active[at] <- dat$epi$num.cond.int.active[at] / sum(active == 1, na.rm = TRUE)
  dat$epi$prev.cond.int.active.asmm[at] <- dat$epi$num.cond.int.active.asmm[at] / sum(active == 1 & asmm == 1, na.rm = TRUE)
  dat$epi$prev.cond.int.active.asmm.AI[at] <-  dat$epi$num.cond.int.active.asmm.AI[at] / sum(active == 1 & everAI ==1 & asmm == 1, na.rm = TRUE)
  dat$epi$prev.cond.int.active.msm[at] <- dat$epi$num.cond.int.active.msm[at] / sum(active == 1 & asmm == 0, na.rm = TRUE)

  return(dat)
}


#' @title Prevalence Module
#'
#' @description Module function to calculate and store summary statistics for
#'              disease prevalence, demographics, and other epidemiological
#'              outcomes.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
prevalence_het <- function(dat, at) {

  status <- dat$attr$status
  male <- dat$attr$male
  age <- dat$attr$age

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  # Initialize vectors
  if (at == 1) {
    dat$epi$i.num <- rNA
    dat$epi$num <- rNA

    dat$epi$i.num.male <- rNA
    dat$epi$i.num.feml <- rNA
    dat$epi$i.prev.male <- rNA
    dat$epi$i.prev.feml <- rNA

    dat$epi$num.male <- rNA
    dat$epi$num.feml <- rNA
    dat$epi$meanAge <- rNA
    dat$epi$propMale <- rNA

    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA
  }

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$num[at] <- length(status)

  dat$epi$i.num.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.num.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE)
  dat$epi$i.prev.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE) /
    sum(male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE) /
    sum(male == 0, na.rm = TRUE)

  dat$epi$num.male[at] <- sum(male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(male == 0, na.rm = TRUE)
  dat$epi$meanAge[at] <- mean(age, na.rm = TRUE)
  dat$epi$propMale[at] <- mean(male, na.rm = TRUE)

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}
