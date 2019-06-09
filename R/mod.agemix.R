
#' @title Track age mixing
#'
#' @description Module function for tracking age on the edgelist.
#'
#' @inheritParams aging_camplc
#'
#' @details
#' The function is a diagnostic tool to determine the age mixing matrix for adult MSM and ASMM.
#'
#' @return
#' This function returns the \code{dat} object with the updated edge
#' list (\code{el.age}). Each element of \code{el.age} is a data frame with the ids of
#' all active pairs with their respective ages.
#'
#' @keywords module msm
#' @export
#'
agemix_campcl <- function(dat, at) {

if(dat$param$agemix == TRUE){
if(at < dat$param$agemix.start){dat$el.age <- list()}

if(at >= dat$param$agemix.start){
    # Attributes
    age <- dat$attr$age

    el.age <- rbind(dat$el[[1]], dat$el[[2]], dat$el[[3]], dat$el[[4]])
  

    # Add ages to the edgelist

    age1 <- age[el.age[, 1]]
    age2 <- age[el.age[, 2]]

    el.age <- cbind(el.age, age1, age2)
    colnames(el.age) <- c("p1", "p2", "age1", "age2")

  # Update edgelist in dat
  dat$el.age <- el.age
  
}
}
  
  return(dat)
}
