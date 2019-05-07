
#' @title Concurrency diagnostics.
#'
#' @description Module function for analysing the duration of concurrencies of particular types.
#'
#' @inheritParams aging_shamp
#'
#' @details
#' At each time step sets of concurent relationships are identified, 
#' the race of the individual in the relationship is stored and the 
#' duration of the overlap is calculated.
#'
#' @return
#' This function returns the \code{dat} object with the updated partnership list for diagnostics
#' on \code{dat$plist}.
#'
#' @keywords module SHAMP cocnurrency diagnostics
#' @export
#'
ConcDurDx_shamp <- function(dat, at){

  if(dat$param$conc_dur_dx == TRUE){

    # Variables --------------------------------------------------------------

    # Attributes

    uid <- dat$attr$uid
    race <- dat$attr$race
    sex <- dat$attr$sex
    age <- floor(dat$attr$age)
    sex.ident <- dat$attr$sex.ident
    Ecohab <- dat$attr$Ecohab
    Ecohab.timer <- dat$attr$Ecohab.timer
    cel.temp <- dat$cel.temp
    cel.complete <- dat$cel.complete
      
    ##  USE UID FOR EDGES
    ## get the race and sex 
      el1<-as.data.frame(dat$el[[1]])
      el1[,3] <- apply(el1[,c(1,2)], 1, FUN=max)
      el1[,4] <- apply(el1[,c(1,2)], 1, FUN=min)
      el1[,5]<-rep("cohab",length(el1[,1]))
      el1 <- el1[,3:5]
 
      el2<-as.data.frame(dat$el[[2]])
      el2[,3] <- apply(el2[,c(1,2)], 1, FUN=max)
      el2[,4] <- apply(el2[,c(1,2)], 1, FUN=min)
      el2[,5]<-rep("pers",length(el2[,1]))
      el2 <- el2[,3:5]
      
      cel <- rbind(el1,el2)
      cel[,4] <- sex[cel[,1]]
      cel[,5] <- sex[cel[,2]]
      cel[,6] <- race[cel[,1]]
      cel[,7] <- race[cel[,2]]
      cel[,8] <- age[cel[,1]]
      cel[,9] <- age[cel[,2]]
     
      ##Check on Ecohab dissolution 
      
      cel[,10] <- Ecohab[cel[,1]]
      cel[,11] <- Ecohab[cel[,2]]
      cel[,12] <- Ecohab.timer[cel[,1]]
      cel[,13] <- Ecohab.timer[cel[,2]]
      colnames(cel)<-c("p1","p2","type", "s1","s2","r1","r2","age1","age2","Ec1","Ec2","Ect1","Ect2")
      cel$Ec1<-ifelse(cel$type=="pers",0,cel$Ec1)
      cel$Ec2<-ifelse(cel$type=="pers",0,cel$Ec2)
      cel$Ect1<-ifelse(cel$type=="pers",0,cel$Ect1)
      cel$Ect2<-ifelse(cel$type=="pers",0,cel$Ect2)
      
      ##  USE UID FOR EDGES
      ##Covert to UID.
      cel[,1] <- uid[cel[,1]]
      cel[,2] <- uid[cel[,2]]
      
      ##Assign rel ID.    
      rels.cel <- cel[,14] <- (cel[,1] * 1000000000) + cel[,2]
      colnames(cel)<-c("p1","p2","type", "s1","s2","r1","r2","age1","age2","Ec1","Ec2","Ect1","Ect2", "ID")
      
      
      
      ##Count duplicates. 
      trelcount <- length(rels.cel)
      urelcount <- length(unique(rels.cel))
      dup_count <- trelcount - urelcount
      dat$epi$duplicates[at] <- dup_count

      ##Count duplicates including OT
      el3<-as.data.frame(dat$el[[3]])
      el3[,3] <- apply(el3[,c(1,2)], 1, FUN=max)
      el3[,4] <- apply(el3[,c(1,2)], 1, FUN=min)
      el3 <- el3[,3:4]
      el3[,1] <- uid[el3[,1]]
      el3[,2] <- uid[el3[,2]]
      
      rels.OT.cel <- (el3[,1] * 1000000000) + el3[,2]
      all.rels.cel<-c(rels.cel,rels.OT.cel)
      
      trelcount <- length(all.rels.cel)
      urelcount <- length(unique(all.rels.cel))
      dup_count <- trelcount - urelcount
      dat$epi$duplicates.w.OT[at] <- dup_count

      ##drop duplicates for counting concurrency durations (there presence interfears with the merge)
      ##Low counts will have a finite inpact of duration calculations.
      drop <- which(duplicated(rels.cel)==TRUE) 
      if(length(drop>0)){cel<-cel[-drop,]}
      

      if(at == 2){cel.temp <- cel
      cel.temp[,15] <- rep(at,length(cel.temp[,1])) 
      cel.temp[,16] <- rep(NA,length(cel.temp[,1])) 
      colnames(cel.temp)<-c("p1","p2","type", "s1","s2","r1","r2","age1","age2","Ec1","Ec2","Ect1","Ect2", "ID","start", "end")
      cel.complete <- cel.temp[0,]
      }
    
      if(at > 2){
        ##Get the rels from cel.temp that are not in cel.
        ids.old <- unique(cel.temp[,"ID"])
        
        ids.new<-NULL
        for(i in 1:length(dat$temp$new.edges[,1])){
        ids.new[i] <-  (max(dat$attr$uid[dat$temp$new.edges[i,1]],dat$attr$uid[dat$temp$new.edges[i,2]]) * 1000000000) 
                      + min(dat$attr$uid[dat$temp$new.edges[i,1]],dat$attr$uid[dat$temp$new.edges[i,2]])
        }
        
        ids.cur <- unique(cel[,"ID"])
        ids.ongoing <- intersect(ids.old,ids.cur)
        ids.ended <- setdiff(ids.old,ids.ongoing) 
        
        ended <- cel.temp[cel.temp$ID %in% ids.ended,]
        ended$end <-at
        cel.complete <- rbind(cel.complete,ended)
        
        ##look at ended Ecohab status?
        ##How many in cohabs in ended are Ecohab Vs Not
        ##compare to base.t
        
        cel.temp <- merge(x = cel, y = cel.temp, by = c("p1","p2"), all.x = TRUE)
        cel.temp <- cel.temp[,c(1:14,27,28)]
        
        
      }

      

      colnames(cel.temp)<-c("p1","p2","type", "s1","s2","r1","r2","age1","age2","Ec1","Ec2","Ect1","Ect2", "ID","start", "end")
      now <- which(is.na(cel.temp$start)==TRUE)
      cel.temp$start[now] <- at
      
      dat$cel.temp <- cel.temp
      dat$cel.complete <- cel.complete
  }
  
      
  return(dat)
}
