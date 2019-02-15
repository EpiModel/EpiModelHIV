
#' @title Aging Module
#'
#' @description Module for aging over time for active nodes in the population.
#'
#' @param dat Master data list object of class \code{dat} containing networks,
#'        individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @return
#' This function returns \code{dat} after updating the nodal attribute
#' \code{age} and \code{sqrt.age}. 
#'
#' @keywords module
#' @export
#'
aging_camplc <- function(dat, at) {

  # system timer
  dat$epi$timer[at] <- proc.time()[3]

  time.unit <- dat$param$time.unit

  age <- dat$attr$age
  active <- dat$attr$active

  age[active == 1] <- age[active == 1] + time.unit / 365

  dat$attr$age <- age
  dat$attr$sqrt.age <- sqrt(age)
  dat$attr$cubert.age <- age^(1/3)
  dat$attr$of.age<-ifelse(dat$attr$age>=16 & dat$attr$active==1,1,0)
  
  dat$attr$risk.age.group<-floor(dat$attr$age)+(dat$attr$riskg/10)


  


  return(dat)
}
