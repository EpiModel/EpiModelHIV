
#' @title Update PrEP adherence from asmm to MSM values at age 19
#'
#' @description Module function for updating act class in main and casual
#'              partnerships based on probabilities of transition.
#'
#' @inheritParams aging_msm
#'
#' @return
#' This function updates the individual-level attribute \code{prepClass} on
#' \code{dat$attr}.
#'
#' @keywords module msm
#' 
#' @export
#'
prep_adherence_trans <- function(dat, at) {


 
  ids.c <-NULL
  
  prepClass <- dat$attr$prepClass
  age<- dat$attr$age
  

  ids.c <- which(age >= 19 & age < 19+(1/52) & is.na(prepClass) == FALSE)
  old.class<-prepClass[ids.c]
  new.class<-old.class
  if(length(new.class) > 0){
  for (i in 1:length(new.class)){
    if(new.class[i] == 2){new.class[i] <-sample(x=0:3, size=1,prob=c(.082,.287,.41,.221))}
    if(new.class[i] == 3){new.class[i] <-4}
  }
  
  prepClass[ids.c]<-new.class
  dat$attr$prepClass <- prepClass
  }
  
  
  return(dat)
}
