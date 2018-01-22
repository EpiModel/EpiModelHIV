
# SHAMP  -----------------------------------------------------------------

#' @title Cut population size for testing
#'
#' @description Imports the data for simulations using ergm.ego and SHAMP function.  Ego and alter data are cut to a limited age range to 
#' reduce object size.
#'
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' data.
#'
#' @keywords shamp
#'
#' @export
make_small <- function(data) {
 
##Drop egos and alters younger than 20 or older than 30.
  

    ids<-which(data[[1]]$egos$age > 19 & data[[1]]$egos$age < 31 )


  data[[1]]$egos<-data[[1]]$egos[data[[1]]$egos$ego %in% ids == TRUE,]
  data[[1]]$altersCohab<-data[[1]]$altersCohab[data[[1]]$altersCohab$ego %in% ids == TRUE,]
  data[[1]]$altersPers<-data[[1]]$altersPers[data[[1]]$altersPers$ego %in% ids == TRUE,]
  data[[1]]$altersOT<-data[[1]]$altersOT[data[[1]]$altersOT$ego %in% ids == TRUE,]
  
 data[[1]]$egos$age <- sample(18:60,size=length(data[[1]]$egos$age),replace=TRUE)
 data[[1]]$altersCohab$age <- sample(18:60,size=length(data[[1]]$altersCohab$age),replace=TRUE)
 data[[1]]$altersPers$age <- sample(18:60,size=length(data[[1]]$altersPers$age),replace=TRUE)
 data[[1]]$altersOT$age <- sample(18:60,size=length(data[[1]]$altersOT$age),replace=TRUE)


  return(list(data))

}


