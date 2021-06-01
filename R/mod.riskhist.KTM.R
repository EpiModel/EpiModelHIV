
#' @title Risk History Module
#'
#' @description Module function to track the risk history of uninfected persons
#'              for purpose of intervention targeting.
#'
#' @inheritParams aging_shamp
#'
#' @keywords module SHAMP
#'
#' @export
#'
riskhist_KTM <- function(dat, at) {

  if (at < dat$param$riskh.start) {
    return(dat)
  }

  dat$attr$prep.ind.disc <- rep(0,length(dat$attr$age))
  dat$attr$prep.ind.parts <- rep(0,length(dat$attr$age))
  uid <- dat$attr$uid
  prep.risk.int <- dat$param$prep.risk.int
  tested.negative <- dat$attr$tested.negative
  

  #Partners list
  cel.temp <- dat$cel.temp
  cel.complete <- dat$cel.complete
  
  #For current negatives
  if(any(tested.negative==1)==TRUE){
    
    negs <- which(tested.negative==1)
    ids.neg <- uid[negs]

#GET THOSE WITH A POSITIVE DISCLOSED PARTNER
    disc.part <- intersect(ids.neg,dat$temp$discl.list[,2])
    x <-dat$attr$uid %in% disc.part
    dat$attr$prep.ind.disc[x==TRUE] <-1
    
    #GET PARTNERS IN PREP RISK INDICATION INTERVAL (0=concurent)
    if(length(cel.complete[,12]) > 0){
    rows <- which(cel.complete[,12] > at-prep.risk.int)
    recent <- rbind(cel.temp,cel.complete[rows,])
    
    ids.rel <- c(recent[,1], recent[,2])
    n_occur <- data.frame(table(ids.rel))
    ids.rel<-n_occur[n_occur$Freq > 1,]
    ids.rel<-as.vector(ids.rel$ids.rel)
    ids.rel<-as.numeric(ids.rel)
    ids.rel<-intersect(ids.neg,ids.rel)
    
    x <-dat$attr$uid %in% ids.rel
    dat$attr$prep.ind.parts[x==TRUE] <-1
    }
    
    
    
  }
  

  return(dat)
}
