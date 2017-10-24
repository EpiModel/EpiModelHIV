
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
  dat$attr$yamsm <- ifelse(dat$attr$age > 18 & dat$attr$age < 26,1,0)
  dat$attr$oamsm <- ifelse(dat$attr$age > 25,1,0)
  
  dat$attr$out <- ifelse(dat$attr$age > dat$attr$out.age,1,dat$attr$out)
  
  
  ids<-which(dat$attr$debuted == 0)
  
  dat$attr$debuted[ids] <- rbinom(ids,1,dat$init$debut.prob)
  

  dat$attr$debuted <- ifelse(dat$attr$age > 19,1,dat$attr$debuted)

  return(dat)
}
