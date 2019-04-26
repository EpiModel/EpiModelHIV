
#' @title Out and Debut Module
#'
#' @description Module coming out (entering the population) and debut as sexually available.
#'
#' @inheritParams aging_camplc
#'
#' @param dat Master data list object of class \code{dat} containing networks,
#'        individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @return
#' This function returns \code{dat} after updating the nodal attribute
#' \code{asmm} , \code{yamsm} ,\code{oamsm} ,\code{out} ,\code{debuted}. 
#'
#' @keywords module
#' @export
#'
out_debut_camplc <- function(dat, at) {
 

  dat$attr$asmm <- ifelse(dat$attr$age < 19,1,0)
  dat$attr$amsm <- ifelse(dat$attr$age >= 19,1,0)
  dat$attr$yamsm <- ifelse(dat$attr$age >= 19 & dat$attr$age < 26,1,0)
  dat$attr$oamsm <- ifelse(dat$attr$age >= 26,1,0)
  
  newout<-which(dat$attr$age <= dat$attr$out.age + ((dat$param$time.unit/365)/2) &  dat$attr$age >= dat$attr$out.age - ((dat$param$time.unit/365)/2) & dat$attr$out == 0)
  oldout<-which(dat$attr$age > dat$attr$out.age & dat$attr$debuted == 0)
  
  dat$attr$out[newout] <- 1
  
  dat$attr$debuted[newout] <- rbinom(newout,1,dat$param$debut.entry.prob)
  dat$attr$debuted[oldout] <- rbinom(oldout,1,dat$param$debut.prob)
  
  dat$attr$debuted <- ifelse(dat$attr$age >= 19,1,dat$attr$debuted)
  dat$attr$out <- ifelse(dat$attr$age >= 19,1,dat$attr$out)
  

  return(dat)
}
