
#' @title Sexual Acts Module for up to 5 race groups and 2 sexes.
#'
#' @description Module function for setting the number of sexual acts on the
#'              discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' The number of acts at each time step is specified as a function of the race of
#' both members in a pair and whether the relationship is male male
#' or heterosexual.  Rates are averages of the reported rates by each race.
#' For one-off partnerships, this is deterministically
#' set at 1, whereas for main and causal partnerships it is a stochastic draw
#' from a Poisson distribution. The number of total acts may further be modified
#' by the level of HIV viral suppression in an infected person.
#'
#' @return
#' This function returns the \code{dat} object with the updated discordant act
#' list (\code{dal}). Each element of \code{dal} is a data frame with the ids of the
#' discordant pair repeated the number of times they have AI.
#'
#' @keywords module SHAMP msm het
#' @export
#'
acts_shamp <- function(dat, at) {

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

    # Attributes
    status <- dat$attr$status
    race <- dat$attr$race
    sex <- dat$attr$sex
    sex.ident <- dat$attr$sex.ident
    immig.loc <- dat$attr$immig.loc

    # Parameters
    ai.scale <- dat$param$ai.scale
    vi.scale <- dat$param$vi.scale
    if (type == "main") {
      base.ai.B.rate <- dat$param$base.ai.main.B.rate
      base.ai.BI.rate <- dat$param$base.ai.main.BI.rate
      base.ai.H.rate <- dat$param$base.ai.main.H.rate
      base.ai.HI.rate <- dat$param$base.ai.main.HI.rate
      base.ai.W.rate <- dat$param$base.ai.main.W.rate
      base.vi.B.rate <- dat$param$base.vi.main.B.rate
      base.vi.BI.rate <- dat$param$base.vi.main.BI.rate
      base.vi.H.rate <- dat$param$base.vi.main.H.rate
      base.vi.HI.rate <- dat$param$base.vi.main.HI.rate
      base.vi.W.rate <- dat$param$base.vi.main.W.rate
      
      fixed.ai <- FALSE
      fixed.vi <- FALSE
      ptype <- 1
      el <- dat$el[[1]]
    }
    if (type == "pers") {
      base.ai.B.rate <- dat$param$base.ai.pers.B.rate
      base.ai.BI.rate <- dat$param$base.ai.pers.BI.rate
      base.ai.H.rate <- dat$param$base.ai.pers.H.rate
      base.ai.HI.rate <- dat$param$base.ai.pers.HI.rate
      base.ai.W.rate <- dat$param$base.ai.pers.W.rate
      base.vi.B.rate <- dat$param$base.vi.pers.B.rate
      base.vi.BI.rate <- dat$param$base.vi.pers.BI.rate
      base.vi.H.rate <- dat$param$base.vi.pers.H.rate
      base.vi.HI.rate <- dat$param$base.vi.pers.HI.rate
      base.vi.W.rate <- dat$param$base.vi.pers.W.rate
      
      fixed.ai <- FALSE
      fixed.vi <- FALSE
      ptype <- 2
      el <- dat$el[[2]]
    }
    if (type == "inst") {
      base.ai.B.rate <- 1
      base.ai.BI.rate <- 1
      base.ai.H.rate <- 1
      base.ai.HI.rate <- 1
      base.ai.W.rate <- 1
      base.vi.B.rate <- 1
      base.vi.BI.rate <- 1
      base.vi.H.rate <- 1
      base.vi.HI.rate <- 1
      base.vi.W.rate <- 1
      fixed.ai <- ifelse(ai.scale != 1, FALSE, TRUE)
      fixed.vi <- ifelse(vi.scale != 1, FALSE, TRUE)
      ptype <- 3
      el <- dat$el[[3]]
    }

    ## Processes ##

    # Construct edgelist

    st1 <- status[el[, 1]]
    st2 <- status[el[, 2]]
    disc <- abs(st1 - st2) == 1
    el[which(disc == 1 & st2 == 1), ] <- el[which(disc == 1 & st2 == 1), 2:1]
    el <- cbind(el, status[el[, 1]], status[el[, 2]])
    colnames(el) <- c("p1", "p2", "st1", "st2")

    if (nrow(el) > 0) {

      # Base AI rates
      ai.rate <- rep(0, nrow(el))
      race.p1 <- race[el[, 1]]
      race.p2 <- race[el[, 2]]
      sex.p1 <- sex[el[, 1]]
      sex.p2 <- sex[el[, 2]]
      immig.loc.p1 <- immig.loc[el[, 1]]
      immig.loc.p2 <- immig.loc[el[, 2]]
      num.het <- (sex.p1 == "F") + (sex.p2 == "F")
      num.away <- (immig.loc.p1 == 1) + (immig.loc.p2 == 1)
      
            
      ai.rate[race.p1=="B" & race.p2=="B" & num.het==0 & num.away==0] <- base.ai.B.rate
      ai.rate[race.p1=="B" & race.p2=="BI" & num.het==0 & num.away==0] <- mean(c(base.ai.B.rate,base.ai.BI.rate))
      ai.rate[race.p1=="B" & race.p2=="H" & num.het==0 & num.away==0] <- mean(c(base.ai.B.rate,base.ai.H.rate))
      ai.rate[race.p1=="B" & race.p2=="HI" & num.het==0 & num.away==0] <- mean(c(base.ai.B.rate,base.ai.HI.rate))
      ai.rate[race.p1=="B" & race.p2=="W" & num.het==0 & num.away==0] <- mean(c(base.ai.B.rate,base.ai.W.rate))
      
      ai.rate[race.p1=="BI" & race.p2=="B" & num.het==0 & num.away==0] <- mean(c(base.ai.B.rate,base.ai.BI.rate))
      ai.rate[race.p1=="H" & race.p2=="B" & num.het==0 & num.away==0] <- mean(c(base.ai.B.rate,base.ai.H.rate))
      ai.rate[race.p1=="HI" & race.p2=="B" & num.het==0 & num.away==0] <- mean(c(base.ai.B.rate,base.ai.HI.rate))
      ai.rate[race.p1=="W" & race.p2=="B" & num.het==0 & num.away==0] <- mean(c(base.ai.B.rate,base.ai.W.rate))
      
      ai.rate[race.p1=="BI" & race.p2=="BI" & num.het==0 & num.away==0] <- base.ai.BI.rate
      ai.rate[race.p1=="BI" & race.p2=="H" & num.het==0 & num.away==0] <- mean(c(base.ai.BI.rate,base.ai.H.rate))
      ai.rate[race.p1=="BI" & race.p2=="HI" & num.het==0 & num.away==0] <- mean(c(base.ai.BI.rate,base.ai.HI.rate))
      ai.rate[race.p1=="BI" & race.p2=="W" & num.het==0 & num.away==0] <- mean(c(base.ai.BI.rate,base.ai.W.rate))
      
      ai.rate[race.p1=="H" & race.p2=="BI" & num.het==0 & num.away==0] <- mean(c(base.ai.BI.rate,base.ai.H.rate))
      ai.rate[race.p1=="HI" & race.p2=="BI" & num.het==0 & num.away==0] <- mean(c(base.ai.BI.rate,base.ai.HI.rate))
      ai.rate[race.p1=="W" & race.p2=="BI" & num.het==0 & num.away==0] <- mean(c(base.ai.BI.rate,base.ai.W.rate))
      
      ai.rate[race.p1=="H" & race.p2=="H" & num.het==0 & num.away==0] <- base.ai.H.rate
      ai.rate[race.p1=="H" & race.p2=="HI" & num.het==0 & num.away==0] <- mean(c(base.ai.H.rate,base.ai.HI.rate))
      ai.rate[race.p1=="H" & race.p2=="W" & num.het==0 & num.away==0] <- mean(c(base.ai.H.rate,base.ai.W.rate))
      
      ai.rate[race.p1=="HI" & race.p2=="H" & num.het==0 & num.away==0] <- mean(c(base.ai.H.rate,base.ai.HI.rate))
      ai.rate[race.p1=="W" & race.p2=="H" & num.het==0 & num.away==0] <- mean(c(base.ai.H.rate,base.ai.W.rate))
      
      ai.rate[race.p1=="HI" & race.p2=="HI" & num.het==0 & num.away==0] <- base.ai.HI.rate
      ai.rate[race.p1=="HI" & race.p2=="W" & num.het==0 & num.away==0] <- mean(c(base.ai.HI.rate,base.ai.W.rate))
      
      ai.rate[race.p1=="W" & race.p2=="HI" & num.het==0 & num.away==0] <- mean(c(base.ai.HI.rate,base.ai.W.rate))
      
      ai.rate[race.p1=="W" & race.p2=="W" & num.het==0 & num.away==0] <- base.ai.W.rate
      
    
      ai.rate <- ai.rate * ai.scale
      
      # Base VI rates
      vi.rate <- rep(0, nrow(el))

      vi.rate[race.p1=="B" & race.p2=="B" & num.het==1 & num.away==0] <- base.vi.B.rate
      vi.rate[race.p1=="B" & race.p2=="BI" & num.het==1 & num.away==0] <- mean(c(base.vi.B.rate,base.vi.BI.rate))
      vi.rate[race.p1=="B" & race.p2=="H" & num.het==1 & num.away==0] <- mean(c(base.vi.B.rate,base.vi.H.rate))
      vi.rate[race.p1=="B" & race.p2=="HI" & num.het==1 & num.away==0] <- mean(c(base.vi.B.rate,base.vi.HI.rate))
      vi.rate[race.p1=="B" & race.p2=="W" & num.het==1 & num.away==0] <- mean(c(base.vi.B.rate,base.vi.W.rate))
      
      vi.rate[race.p1=="BI" & race.p2=="B" & num.het==1 & num.away==0] <- mean(c(base.vi.B.rate,base.vi.BI.rate))
      vi.rate[race.p1=="H" & race.p2=="B" & num.het==1 & num.away==0] <- mean(c(base.vi.B.rate,base.vi.H.rate))
      vi.rate[race.p1=="HI" & race.p2=="B" & num.het==1 & num.away==0] <- mean(c(base.vi.B.rate,base.vi.HI.rate))
      vi.rate[race.p1=="W" & race.p2=="B" & num.het==1 & num.away==0] <- mean(c(base.vi.B.rate,base.vi.W.rate))
      
      vi.rate[race.p1=="BI" & race.p2=="BI" & num.het==1 & num.away==0] <- base.vi.BI.rate
      vi.rate[race.p1=="BI" & race.p2=="H" & num.het==1 & num.away==0] <- mean(c(base.vi.BI.rate,base.vi.H.rate))
      vi.rate[race.p1=="BI" & race.p2=="HI" & num.het==1 & num.away==0] <- mean(c(base.vi.BI.rate,base.vi.HI.rate))
      vi.rate[race.p1=="BI" & race.p2=="W" & num.het==1 & num.away==0] <- mean(c(base.vi.BI.rate,base.vi.W.rate))
      
      vi.rate[race.p1=="H" & race.p2=="BI" & num.het==1 & num.away==0] <- mean(c(base.vi.BI.rate,base.vi.H.rate))
      vi.rate[race.p1=="HI" & race.p2=="BI" & num.het==1 & num.away==0] <- mean(c(base.vi.BI.rate,base.vi.HI.rate))
      vi.rate[race.p1=="W" & race.p2=="BI" & num.het==1 & num.away==0] <- mean(c(base.vi.BI.rate,base.vi.W.rate))
      
      vi.rate[race.p1=="H" & race.p2=="H" & num.het==1 & num.away==0] <- base.vi.H.rate
      vi.rate[race.p1=="H" & race.p2=="HI" & num.het==1 & num.away==0] <- mean(c(base.vi.H.rate,base.vi.HI.rate))
      vi.rate[race.p1=="H" & race.p2=="W" & num.het==1 & num.away==0] <- mean(c(base.vi.H.rate,base.vi.W.rate))
      
      vi.rate[race.p1=="HI" & race.p2=="H" & num.het==1 & num.away==0] <- mean(c(base.vi.H.rate,base.vi.HI.rate))
      vi.rate[race.p1=="W" & race.p2=="H" & num.het==1 & num.away==0] <- mean(c(base.vi.H.rate,base.vi.W.rate))
      
      vi.rate[race.p1=="HI" & race.p2=="HI" & num.het==1 & num.away==0] <- base.vi.HI.rate
      vi.rate[race.p1=="HI" & race.p2=="W" & num.het==1 & num.away==0] <- mean(c(base.vi.HI.rate,base.vi.W.rate))
      
      vi.rate[race.p1=="W" & race.p2=="HI" & num.het==1 & num.away==0] <- mean(c(base.vi.HI.rate,base.vi.W.rate))
      
      vi.rate[race.p1=="W" & race.p2=="W" & num.het==1 & num.away==0] <- base.vi.W.rate

      vi.rate <- vi.rate * vi.scale
      

      # Final act number
      ai.rate[is.na(ai.rate)] <- 0
      if (fixed.ai == FALSE) {
        ai <- rpois(length(ai.rate), ai.rate)
      } else {
        ai <- round(ai.rate)
      }
      vi.rate[is.na(vi.rate)] <- 0
      if (fixed.vi == FALSE) {
        vi <- rpois(length(vi.rate), vi.rate)
      } else {
        vi <- round(vi.rate)
      }
      
      # Full edge list
      el <- cbind(el, ptype, ai, vi)
      colnames(el)[5:7] <- c("ptype", "ai", "vi")
      
 
      if (type == "main") {
        dat$temp$el <- el
      } else {
        dat$temp$el <- rbind(dat$temp$el, el)
      }
    }

  } # loop over type end

  # Remove inactive edges and edges 
  
  acts<-dat$temp$el[,"ai"]+dat$temp$el[,"vi"]
  keep<-which(acts > 0)
  dat$temp$el <- dat$temp$el[keep, ]


  
  return(dat)
}
