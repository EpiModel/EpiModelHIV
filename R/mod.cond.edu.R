
#' @title Condom Intervention Module
#'
#' @description Module function selects those eligible for a Condom use intervention and updates the
#'              intervention effect timer.
#'
#' @inheritParams aging_camplc
#'
#' @details
#' ASMM may undergo a sex education intervention that increases condom use.  This module selects which ASMM 
#' are exposed to the intervention and set the duration of the intervention effect.
#'
#' @return
#' Updates attr cond.edu.timer to indicate who is exposed and the number of timesteps since exposure.
#'
#' @keywords module msm
#' @export
#'


condoms_edu_campcl <- function(dat, at) {

 if (dat$param$cond.edu == TRUE & at >= dat$param$cond.edu.start) {
  
 if (at == dat$param$cond.edu.start) {dat$attr$age.temp <- dat$attr$age}
   
   #attributes
   age <- dat$attr$age
   age.temp <- dat$attr$age.temp
   everAI <- dat$attr$everAI
   cond.int.active <- dat$attr$cond.int.active
   cond.int.active.dur <- dat$attr$cond.int.active.dur
   
   
   #parameters
   cond.edu.at.age <- dat$param$cond.edu.at.age
   cond.edu.cov <- dat$param$cond.edu.cov
   cond.edu.effect.dur <- dat$param$cond.edu.effect.dur
   
   if (dat$param$sex.exp.edu == FALSE & dat$param$sex.edu.no.sex == FALSE) {
   elig <- which(age.temp < cond.edu.at.age & age >= cond.edu.at.age)
   num.elig <- length(elig)
   if (num.elig >=1){
   selected <- rbinom(num.elig,1,cond.edu.cov)
   elig <- elig[selected==1]
   }
   }
   
   if (dat$param$sex.exp.edu == FALSE & dat$param$sex.edu.no.sex == TRUE) {
      elig <- which(age.temp < cond.edu.at.age & age >= cond.edu.at.age)
      num.elig <- length(elig)
      if (num.elig >=1){
         selected <- rbinom(num.elig,1,cond.edu.cov)
         elig <- elig[selected==1]
         AI <- everAI[elig]
         elig <- elig[AI==0]
      }
   }
   
   if (dat$param$sex.exp.edu == TRUE & dat$param$sex.edu.no.sex == FALSE) {
     elig <- which(age.temp < cond.edu.at.age & age >= cond.edu.at.age)
     num.elig <- length(elig)
     if (num.elig >=1){
     selected <- rbinom(num.elig,1,cond.edu.cov)
     elig <- elig[selected==1]
     AI <- everAI[elig]
     elig <- elig[AI==1]
   }
   }
   
   #Move the effect clock up one time step for those already active
   active <- which(cond.int.active == 1)
   cond.int.active.dur[active] <- cond.int.active.dur[active] + 1
   
   expired <- which(cond.int.active.dur > cond.edu.effect.dur)
   cond.int.active.dur[expired] <- 0
   cond.int.active[expired] <- 0
  
   
   
   cond.int.active[elig] <-1
   cond.int.active.dur[elig] <- 1

   #Set new attributes
   dat$attr$age.temp <- dat$attr$age
   dat$attr$cond.int.active <-cond.int.active
   dat$attr$cond.int.active.dur <- cond.int.active.dur
   
   #Set new outputs
   #Calculate coverage for assm, sexually experienced asmm,everyone, adult MSM 
   
   dat$epi$incid.cond.int.start[at] <- length(elig)
   dat$epi$incid.cond.int.stop[at] <- length(expired)
   

 }
  
  return(dat)
}
