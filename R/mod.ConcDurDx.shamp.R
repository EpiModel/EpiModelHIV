
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
    sex.ident <- dat$attr$sex.ident
    cel.temp <- dat$cel.temp
      
    ##  USE UID FOR EDGES
    ## the race and sex now need to be in reference to UID
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
      colnames(cel)<-c("p1","p2","type", "s1","s2","r1","r2")

      if(at == 2){cel.temp <- cel
      
      cel.temp[,8] <- rep(at,length(cel.temp[,1])) 
      cel.temp[,1] <- uid[cel.temp[,1]]
      cel.temp[,2] <- uid[cel.temp[,2]]
      colnames(cel.temp)<-c("p1","p2","type", "s1","s2","r1","r2","start")
      }
    
      if(at > 2){
        
        cel[,1] <- uid[cel[,1]]
        cel[,2] <- uid[cel[,2]]
        
        cel.temp <- merge(x = cel, y = cel.temp, by = c("p1","p2"), all.x = TRUE)
        cel.temp <- cel.temp[,c(1:7,13)]
      }

      

      colnames(cel.temp)<-c("p1","p2","type", "s1","s2","r1","r2","start")
      now <- which(is.na(cel.temp$start)==TRUE)
      cel.temp$start[now] <- at
      
      dat$cel.temp <- cel.temp
      ##CHECK SIZE OF dat$cel.temp.
      
      rels <- (cel.temp[,1] * 1000000000) + cel.temp[,2]
      trelcount <- length(rels)
      urelcount <- length(unique(rels))
      duplicates <- trelcount - urelcount
      dat$epi$duplicates[at] <- duplicates
      
  }
  
      
  return(dat)
}
