
#' @title Update Role Class in Main and Casual Partnerships
#'
#' @description Module function for updating act class in main and casual
#'              partnerships based on probabilities of transition.
#'
#' @inheritParams aging_msm
#'
#' @return
#' This function updates the individual-level attribute \code{role.class} on
#' \code{dat$attr}.
#'
#' @keywords module msm
#' 
#' @export
#'
roleclass_trans <- function(dat, at) {


  deg.main <- get_degree(dat$el[[1]])
  deg.pers <- get_degree(dat$el[[2]])
  deg.asmm <- get_degree(dat$el[[3]])
  
  ids.B.r <- which(dat$attr$race == "B" & dat$attr$role.class == "R" & dat$attr$age == 19 & deg.main == 0 & deg.pers == 0 & deg.asmm == 0)
  ids.B.v <- which(dat$attr$race == "B" & dat$attr$role.class == "V" & dat$attr$age == 19 & deg.main == 0 & deg.pers == 0 & deg.asmm == 0)
  
  ids.W.r <- which(dat$attr$race == "W" & dat$attr$role.class == "R" & dat$attr$age == 19 & deg.main == 0 & deg.pers == 0 & deg.asmm == 0)
  ids.W.v <- which(dat$attr$race == "W" & dat$attr$role.class == "V" & dat$attr$age == 19 & deg.main == 0 & deg.pers == 0 & deg.asmm == 0)
  
  popsize.B <- sum(dat$attr$age == 19 & dat$attr$race == "B")
  popsize.W <- sum(dat$attr$age == 19 & dat$attr$race == "W")
  
  role.class[ids.B.r] <- sample(apportion_lr(length(ids.B.r), c("I", "R"), prob=c(role.shift[1], 1-(role.shift[1])))) 
  role.class[ids.B.v] <- sample(apportion_lr(length(ids.B.v), c("I", "V"), prob=c(role.shift[2], 1-(role.shift[2])))) 
  
  role.class[ids.W.r] <- sample(apportion_lr(length(ids.W.r), c("I", "R"), prob=c(role.shift[1], 1-(role.shift[1])))) 
  role.class[ids.W.v] <- sample(apportion_lr(length(ids.W.v), c("I", "V"), prob=c(role.shift[2], 1-(role.shift[2])))) 

  
  dat$attr$role.class <- new.role.class

  return(dat)
}
