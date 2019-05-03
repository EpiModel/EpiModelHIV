
#' @title Follow-up elements Module
#'
#' @description Module sets up adolcent risk history and infection age elements of dat.
#'
#' @inheritParams FU_elements_camplc
#'
#' @param dat Master data list object of class \code{dat} containing networks,
#'        individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @return
#' This function returns \code{dat} after adding elements to dat that are not tracked during the burnin
#' \code{riskhist} , \code{age.inf.vec}, \code{counts}. 
#'
#' @keywords module
#' @export
#'
FU_elements_camplc <- function(dat, at) {
 
  
  if (dat$param$FU == TRUE & dat$control$start == at) {
  # Risk history lists adol 
    
  if (dat$param$prep.risk.hist.asmm == TRUE){
    
  t1<-t2<-t3<-t4<-t5<-t6<-t7<-t8<-t9<-t10<-t11<-t12<-t13<-t14<-t15<-t16<-t17<-t18<-t19<-t20<-t21<-t22<-t23<-t24<-t25<-t26<-rep(0,length(dat$attr$uid))
  dat$riskhist<-cbind(dat$attr$uid,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26)
  dat$riskhist<-as.data.frame(dat$riskhist)
  
  }
    
   
  #Track the age at which each new infection is aquired
  dat$age.inf.vec <- list()
  
  # Risk history lists
  nc <- ceiling(dat$param$prep.risk.int)
  dat$riskh <- list()
  rh.names <- c("uai.mono2.nt.3mo", "uai.mono2.nt.6mo",
                "uai.mono1.nt.3mo", "uai.mono1.nt.6mo",
                "uai.nonmonog", "uai.nmain",
                "ai.sd.mc", "uai.sd.mc")
  for (i in 1:length(rh.names)) {
    dat$riskh[[rh.names[i]]] <- matrix(NA, ncol = nc, nrow = length(dat$attr$uid))
  }
  
  
  }

  return(dat)
}
