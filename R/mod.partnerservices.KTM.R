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
  
 # If diag time=at-1 get there partners.  last 12 months for prev and 3 months for acute.
dat$attr$diag.time
cel.temp <- dat$cel.temp
cel.complete <- dat$cel.complete
#  follow acute partner for last three mongths and prevalent for the last 12.
  
#The partner services relationship edgelist

  psel <- dat$psel

#The current edge list
el <- rbind(dat$el[[1]],dat$el[[2]],dat$el[[3]])  

#Edges formed in this time step
new.edges <- rbind(dat$temp$new.edges, dat$el[[3]])
    #assign new edges a start time.
    start.time <- rep(at,length(new.edges[,1]))
    new.edges <- cbind(new.edges,start.time)
    
#Merge new edges with psel
psel<-rbind(psel,new.edges)
    

## drop relationships that are older than the window and not active  


  return(dat)
}

