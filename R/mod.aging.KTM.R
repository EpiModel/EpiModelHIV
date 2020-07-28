
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
#' \code{age} and \code{sqrt.age}. The \code{sqrt.age} vertex attribute is also
#' updated on the three networks.
#'
#' @keywords module
#' @export
#'
aging_KTM <- function(dat, at) {

  # system timer
  dat$epi$timer[at] <- proc.time()[3]

  time.unit <- dat$param$time.unit

  age <- dat$attr$age
  active <- dat$attr$active

  age[active == 1] <- age[active == 1] + time.unit / 365

  dat$attr$age <- age
  dat$attr$sqrt.age <- sqrt(age)
  dat$attr$agesq <- dat$attr$age^2

  dat$attr$age.adj<-ifelse(dat$attr$sex=="M",dat$attr$age,
                                ifelse(dat$attr$sex=="F",dat$attr$age + dat$param$age.adj,dat$attr$age))
  
  dat$attr$age_adj <- dat$attr$age.adj

  dat$attr$sqrt.age.adj<-ifelse(dat$attr$sex=="M",dat$attr$sqrt.age,
                                ifelse(dat$attr$sex=="F",sqrt(dat$attr$age-5),dat$attr$sqrt.age))
  
  
  #AGEGROUP
  dat$attr$age.group<-ifelse(dat$att$age < 20 ,1,
                          ifelse(dat$attr$age >= 20 & dat$attr$age < 25,2,
                          ifelse(dat$attr$age >= 25 & dat$attr$age < 30,3,
                          ifelse(dat$attr$age >= 30 & dat$attr$age < 35,4,
                          ifelse(dat$attr$age >= 35,5,dat$attr$age.group)))))
  
 dat$attr$age_grp <- dat$attr$age.group

 #population of interest
 
 poi<-which(dat$attr$age >= 40)
 dat$attr$poi[poi]<-1
 
    
  return(dat)
}
