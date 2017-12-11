
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

  role.class <- dat$attr$role.class
  role.shift <- dat$param$role.shift
  
  deg.main <- get_degree(dat$el[[1]])
  deg.pers <- get_degree(dat$el[[2]])
  deg.asmm <- get_degree(dat$el[[3]])
  
  ids.B.r <- which(dat$attr$race == "B" & dat$attr$role.class == "R" & dat$attr$age >= 19 & dat$attr$age <= 19.03 & deg.main == 0 & deg.pers == 0 & deg.asmm == 0)
  ids.B.v <- which(dat$attr$race == "B" & dat$attr$role.class == "V" & dat$attr$age >= 19 & dat$attr$age <= 19.03 & deg.main == 0 & deg.pers == 0 & deg.asmm == 0)
  
  ids.W.r <- which(dat$attr$race == "W" & dat$attr$role.class == "R" & dat$attr$age >= 19 & dat$attr$age <= 19.03  & deg.main == 0 & deg.pers == 0 & deg.asmm == 0)
  ids.W.v <- which(dat$attr$race == "W" & dat$attr$role.class == "V" & dat$attr$age >= 19 & dat$attr$age <= 19.03  & deg.main == 0 & deg.pers == 0 & deg.asmm == 0)
  
  ir.list<-c(rep("I",role.shift[1]*1000),rep("R",(1-role.shift[1])*1000))
  iv.list<-c(rep("I",role.shift[2]*1000),rep("V",(1-role.shift[2])*1000))
  
  role.class[ids.B.r] <- sample(ir.list,size=length(ids.B.r)) 
  role.class[ids.B.v] <- sample(iv.list,size=length(ids.B.v)) 
  
  role.class[ids.W.r] <- sample(ir.list,size=length(ids.W.r))
  role.class[ids.W.v] <- sample(iv.list,size=length(ids.B.v)) 

  
  dat$attr$role.class <- role.class

  return(dat)
}
