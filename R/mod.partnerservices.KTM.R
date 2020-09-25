# SHAMP -----------------------------------------------------------------

#' @title Partner Services Module
#'
#' @description Module function for tracking partnerships over a predetermined time window.
#'
#' @inheritParams aging_shamp
#'
#' @keywords module shamp
#'
#' @export
#'
pservices_KTM <- function(dat, at) {

  
intervention_TM <- dat$param$intervention_TM
count <- 0
  
  
if(intervention_TM=="TMP" & at > 2){

    #Attributes    
index.acute <-  dat$attr$PS.index.acute
index.prev <- dat$attr$PS.index.prev
p.prev.follow.win <- dat$param$p.prev.follow.win
p.acute.follow.win <- dat$param$p.acute.follow.win
uid <- dat$attr$uid

#Parameters
partner.serv.part <- dat$attr$partner.serv.part
partner.serv.part.time <- dat$attr$partner.serv.part.time
  
#Partners list
cel.temp <- dat$cel.temp
cel.complete <- dat$cel.complete

#For prevalent index cases
if(any(index.prev==1)==TRUE){
  
  index.prev <- which(index.prev==1)
  index.prev.ids <- uid[index.prev]
  
  #GET CURRENT PARTNERS
  p1 <- cel.temp[,1]
  p2 <- cel.temp[,2]
  
  p1 <- match(p1,index.prev.ids)
  p2 <- match(p2,index.prev.ids)
  
  p1<-which(p1 >=0)
  p2<-which(p2 >=0)
  
  part.p1 <-cel.temp[p1,2]
  part.p2 <-cel.temp[p2,1]
  partners <- c(part.p1, part.p2)
  partners <- unique(partners)
  count <- length(partners)

  x <-dat$attr$uid %in% partners
  dat$attr$partner.serv.part[x==TRUE] <-1
  dat$attr$partner.serv.part.time[x==TRUE] <- 1
  
  #GET PARTNERS IN THE LAST YEAR
  rows <- which(cel.complete[,12] > at-p.prev.follow.win)
  recent <- cel.complete[rows,]
  
  p1 <- recent[,1]
  p2 <- recent[,2]
  
  p1 <- match(p1,index.prev.ids)
  p2 <- match(p2,index.prev.ids)
  
  p1<-which(p1 >=0)
  p2<-which(p2 >=0)
  
  part.p1 <-recent[p1,2]
  part.p2 <-recent[p2,1]
  partners <- c(part.p1, part.p2)
  partners <- unique(partners)
  count <- count + length(partners)
  
  x <-dat$attr$uid %in% partners
  dat$attr$partner.serv.part[x==TRUE] <-1
  dat$attr$partner.serv.part.time[x==TRUE] <- 1
  
  

}



#For acute index cases
if(any(index.acute==1)==TRUE){
  
  index.acute <- which(index.acute==1)
  index.acute.ids <- uid[index.acute]
  
  #GET CURRENT PARTNERS
  p1 <- cel.temp[,1]
  p2 <- cel.temp[,2]
  
  p1 <- match(p1,index.acute.ids)
  p2 <- match(p2,index.acute.ids)
  
  p1<-which(p1 >=0)
  p2<-which(p2 >=0)
  
  part.p1 <-cel.temp[p1,2]
  part.p2 <-cel.temp[p2,1]
  partners <- c(part.p1, part.p2)
  partners <- unique(partners)
  count <- count + length(partners)
  
  x <-dat$attr$uid %in% partners
  dat$attr$partner.serv.part[x==TRUE] <-1
  dat$attr$partner.serv.part.time[x==TRUE] <- 1
  
  #GET PARTNERS IN THE LAST 12 weeks
  rows <- which(cel.complete[,12] >= at - p.acute.follow.win)
  recent <- cel.complete[rows,]
  
  p1 <- recent[,1]
  p2 <- recent[,2]
  
  p1 <- match(p1,index.acute.ids)
  p2 <- match(p2,index.acute.ids)
  
  p1<-which(p1 >=0)
  p2<-which(p2 >=0)
  
  part.p1 <-recent[p1,2]
  part.p2 <-recent[p2,1]
  partners <- c(part.p1, part.p2)
  partners <- unique(partners)
  count <- count + length(partners)
  
  x <-dat$attr$uid %in% partners
  dat$attr$partner.serv.part[x==TRUE] <-1
  dat$attr$partner.serv.part.time[x==TRUE] <- 1
  
  
  
}


#Clear partner services index status
dat$attr$PS.index.acute <- rep(0,length(dat$attr$PS.index.acute))
dat$attr$PS.index.prev <- rep(0,length(dat$attr$PS.index.prev))


#Number of new partners
dat$epi$partners.sought.new[at] <- count


}



  return(dat)
}

