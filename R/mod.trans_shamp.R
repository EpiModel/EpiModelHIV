
# SHAMP tranmission for f, msf, msm, msmf populations  -----------------------------------------------------------------

#' @title Transmission Module for AI and VI
#'
#' @description Stochastically simulates disease transmission given the current
#'              state of the discordand edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This is the final substantive function that occurs within the time loop at
#' each time step. This function takes the discordant edgelist and calculates a
#' transmission probability for each row (one sexual act) between dyads on the
#' network. After transmission events, individual-level attributes for the infected
#' persons are updated and summary statistics for incidence calculated.
#'
#' The per-act transmission probability depends on the following elements:
#' insertive versus receptive role, viral load of the infected partner, an
#' acute stage infection excess risk, condom use, and the CCR5 genetic allele.
#' Given these transmission probabilities, transmission is stochastically
#' simulating by drawing from a binomial distribution for each act conditional
#' on the per-act probability.
#'
#' @return
#' For each new infection, the disease status, infection time, and related
#' HIV attributes are updated for the infected node. Summary statistics for
#' disease incidence overall, and by race and age groups are calculated and
#' stored on \code{dat$epi}.
#'
#' @keywords module SHAMP msm het 
#'
#' @export
#'
trans_shamp <- function(dat, at){

  # Variables -----------------------------------------------------------

  # Attributes
  vl <- dat$attr$vl
  stage <- dat$attr$stage
  ccr5 <- dat$attr$ccr5
  circ <- dat$attr$circ
  diag.status <- dat$attr$diag.status
  status <- dat$attr$diag.status
  tx.status <- dat$attr$tx.status
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  sex<-dat$attr$sex
  race <- dat$attr$race
  sex.ident<-dat$attr$sex.ident
  inf.class <- dat$attr$inf.class

  # Parameters
  URAI.prob <- dat$param$URAI.prob
  UIAI.prob <- dat$param$UIAI.prob
  URVI.prob <- dat$param$URVI.prob
  UIVI.prob <- dat$param$UIVI.prob
  
  acute.rr <- dat$param$acute.rr
  condom.rr <- dat$param$condom.rr
  circ.rr <- dat$param$circ.rr
  
  
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr
  prep.hr <- dat$param$prep.class.hr

  # Data
  
  dal <- dat$temp$dal
  dal <- dal[sample(1:nrow(dal)), ]
  ncols <- dim(dal)[2]
  

  if (nrow(dal) == 0) {
    return(dat)
  }

  ## Reorder by role: ins on the left, rec on the right,
  ##                  with flippers represented twice
  disc.ip <- as.data.frame(dal[dal[, "ins"] %in% 1:2, ])
  disc.rp <- as.data.frame(dal[dal[, "ins"] %in% c(0, 2), c(2:1, 3:ncols)])
  colnames(disc.ip)[1:2] <- c("i", "r")
  colnames(disc.rp)[1:2] <- c("i", "r")
  
  sex.i <- sex[disc.ip[,"i"]]
  sex.r <- sex[disc.ip[,"r"]]
  het <- (sex.i == "F") + (sex.r == "F")
  
  disc.ip<-cbind(disc.ip,sex.i,sex.r,het)
  
  sex.i <- sex[disc.rp[,"i"]]
  sex.r <- sex[disc.rp[,"r"]]
  het <- (sex.i == "F") + (sex.r == "F")
  
  disc.rp<-cbind(disc.rp,sex.i,sex.r,het)


  # PATP: Insertive Partner Infected (Col 1) --------------------------------

  # Attributes of infected
  ip.vl <- vl[disc.ip[, 1]]
  ip.stage <- stage[disc.ip[, 1]]

  # Attributes of susceptible
  ip.ccr5 <- ccr5[disc.ip[, 2]]
  ip.prep <- prepStat[disc.ip[, 2]]
  ip.prepcl <- prepClass[disc.ip[, 2]]

  # Base TP from VL
  ip.tprob.ai <- URAI.prob * 2.45^(ip.vl - 4.5)
  ip.tprob.vi <- URVI.prob * 2.45^(ip.vl - 4.5)
  
  #set base to zero for non-event (AI/VI MM/HET)
  ai<-which(disc.ip$het!=1)
  vi<-which(disc.ip$het==1)
  
  ip.tprob.ai[vi] <- 0
  ip.tprob.vi[ai] <- 0
  
  # Transform to log odds
  ip.tlo.ai <- log(ip.tprob.ai/(1-ip.tprob.ai))
  ip.tlo.vi <- log(ip.tprob.vi/(1-ip.tprob.vi))
  
  # Condom use
  not.UAI <- which(disc.ip[, "uai"] == 0)
  not.UVI <- which(disc.ip[, "uvi"] == 0)
  
  ip.tlo.ai[not.UAI] <- ip.tlo.ai[not.UAI] + log(condom.rr)
  ip.tlo.vi[not.UVI] <- ip.tlo.vi[not.UVI] + log(condom.rr)
  
  # CCR5
  ip.tlo.ai[ip.ccr5 == "DD"] <- ip.tlo.ai[ip.ccr5 == "DD"] + -Inf
  ip.tlo.ai[ip.ccr5 == "DW"] <- ip.tlo.ai[ip.ccr5 == "DW"] + log(ccr5.heteroz.rr)
  
  ip.tlo.vi[ip.ccr5 == "DD"] <- ip.tlo.vi[ip.ccr5 == "DD"] + -Inf
  ip.tlo.vi[ip.ccr5 == "DW"] <- ip.tlo.vi[ip.ccr5 == "DW"] + log(ccr5.heteroz.rr)

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(ip.prep == 1 & ip.prepcl == i-1)
    ip.tlo.ai[temp.ids] <- ip.tlo.ai[temp.ids] + log(prep.hr[i])
    ip.tlo.vi[temp.ids] <- ip.tlo.vi[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% c(1, 2))
  ip.tlo.ai[isAcute] <- ip.tlo.ai[isAcute] + log(acute.rr)
  ip.tlo.vi[isAcute] <- ip.tlo.vi[isAcute] + log(acute.rr)

  # Retransformation to probability
  ip.tprob.ai <- plogis(ip.tlo.ai)
  ip.tprob.vi <- plogis(ip.tlo.vi)
  
  
  if(any(ip.tprob.ai < 0 | ip.tprob.ai > 1 | ip.tprob.ai == "NaN")) {warning("AI Transmission probability not 0-1, set to min0, max1, suggest reduce baseline prob.",call. = FALSE)}
  if(any(ip.tprob.vi < 0 | ip.tprob.vi > 1 | ip.tprob.vi == "NaN")) {warning("VI Transmission probability not 0-1, set to min0, max1, suggest reduce baseline prob.",call. = FALSE)}
  
  ip.tprob.ai <- ifelse(ip.tprob.ai >= 1,1,
                        ifelse(ip.tprob.ai <= 0,0,ip.tprob.ai))
  ip.tprob.ai <- ifelse(is.na(ip.tprob.ai) == TRUE,1,ip.tprob.ai)
  
  ip.tprob.vi <- ifelse(ip.tprob.vi >= 1,1,
                        ifelse(ip.tprob.vi <= 0,0,ip.tprob.vi))
  ip.tprob.vi <- ifelse(is.na(ip.tprob.vi) == TRUE,1,ip.tprob.vi)

  # PATP: Receptive Infected (Col 2) --------------------------------


  # Attributes of infected
  rp.vl <- vl[disc.rp[, 2]]
  rp.stage <- stage[disc.rp[, 2]]

  # Attributes of susceptible
  rp.circ <- circ[disc.rp[, 1]]
  rp.ccr5 <- ccr5[disc.rp[, 1]]
  rp.prep <- prepStat[disc.rp[, 1]]
  rp.prepcl <- prepClass[disc.rp[, 1]]

  # Base TP from VL
  rp.tprob.ai <- UIAI.prob * 2.45^(rp.vl - 4.5)
  rp.tprob.vi <- UIVI.prob * 2.45^(rp.vl - 4.5)
  
  #set base to zero for non-event (AI/VI MM/HET)
  ai<-which(disc.rp$het!=1)
  vi<-which(disc.rp$het==1)
  
  rp.tprob.ai[vi] <- 0
  rp.tprob.vi[ai] <- 0
  
  
  
  # Transform to log odds
  rp.tlo.ai <- log(rp.tprob.ai/(1-rp.tprob.ai))
  rp.tlo.vi <- log(rp.tprob.vi/(1-rp.tprob.vi))
  
  # Circumcision
  rp.tlo.ai[rp.circ == 1] <- rp.tlo.ai[rp.circ == 1] + log(circ.rr)
  rp.tlo.vi[rp.circ == 1] <- rp.tlo.vi[rp.circ == 1] + log(circ.rr)
  
  # Condom use
  not.UAI <- which(disc.rp[, "uai"] == 0)
  rp.tlo.ai[not.UAI] <- rp.tlo.ai[not.UAI] + log(condom.rr)
  
  not.UVI <- which(disc.rp[, "uvi"] == 0)
  rp.tlo.vi[not.UVI] <- rp.tlo.vi[not.UVI] + log(condom.rr)

  # CCR5
  rp.tlo.ai[rp.ccr5 == "DD"] <- rp.tlo.ai[rp.ccr5 == "DD"] + -Inf
  rp.tlo.ai[rp.ccr5 == "DW"] <- rp.tlo.ai[rp.ccr5 == "DW"] + log(ccr5.heteroz.rr)
  
  rp.tlo.vi[rp.ccr5 == "DD"] <- rp.tlo.vi[rp.ccr5 == "DD"] + -Inf
  rp.tlo.vi[rp.ccr5 == "DW"] <- rp.tlo.vi[rp.ccr5 == "DW"] + log(ccr5.heteroz.rr)

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(rp.prep == 1 & rp.prepcl == i-1)
    rp.tlo.ai[temp.ids] <- rp.tlo.ai[temp.ids] + log(prep.hr[i])
    rp.tlo.vi[temp.ids] <- rp.tlo.vi[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% c(1, 2))
  rp.tlo.ai[isAcute] <- rp.tlo.ai[isAcute] + log(acute.rr)
  rp.tlo.vi[isAcute] <- rp.tlo.vi[isAcute] + log(acute.rr)
  
  # Retransformation to probability
  rp.tprob.ai <- plogis(rp.tlo.ai)
  rp.tprob.vi <- plogis(rp.tlo.vi)
  
  if(any(rp.tprob.ai < 0 | rp.tprob.ai > 1 | rp.tprob.ai== "NaN")) {warning("AI Transmission probability not 0-1, set to min0, max1, suggest reduce baseline prob.",call. = FALSE)}
  if(any(rp.tprob.vi < 0 | rp.tprob.vi > 1 | rp.tprob.vi== "NaN")) {warning("VI Transmission probability not 0-1, set to min0, max1, suggest reduce baseline prob.",call. = FALSE)}
  
  rp.tprob.ai <- ifelse(rp.tprob.ai >= 1,1,
                        ifelse(rp.tprob.ai <= 0,0,rp.tprob.ai))
  rp.tprob.ai <- ifelse(is.na(rp.tprob.ai) == TRUE,1,rp.tprob.ai)
  
  rp.tprob.vi <- ifelse(rp.tprob.vi >= 1,1,
                        ifelse(rp.tprob.vi <= 0,0,rp.tprob.vi))
  rp.tprob.vi <- ifelse(is.na(rp.tprob.vi) == TRUE,1,rp.tprob.vi)

  # Transmission --------------------------------------------------------

  ## Bernoulli transmission events
  trans.ip.ai <- rbinom(length(ip.tprob.ai), 1, ip.tprob.ai)
  trans.rp.ai <- rbinom(length(rp.tprob.ai), 1, rp.tprob.ai)
  
  trans.ip.vi <- rbinom(length(ip.tprob.vi), 1, ip.tprob.vi)
  trans.rp.vi <- rbinom(length(rp.tprob.vi), 1, rp.tprob.vi)


  # Output --------------------------------------------------------------

  # Update attributes
  
## Forward transmssion from a heatbath and foreighn aquired get seperate designation.
  inf.class<-ifelse(inf.class=="MSM","MSMds",inf.class)
  inf.class<-ifelse(inf.class=="FA","FAds",inf.class)

  infected <- infector <- inf.type <- NULL
  if (sum(trans.ip.ai, trans.rp.ai,trans.ip.vi, trans.rp.vi) > 0) {

    infected <- c(disc.ip[trans.ip.ai == 1, 2],
                     disc.rp[trans.rp.ai == 1, 1],
                     disc.ip[trans.ip.vi == 1, 2],
                     disc.rp[trans.rp.vi == 1, 1])
    
    
    
    infector <- c(disc.ip[trans.ip.ai == 1, 1],
                  disc.rp[trans.rp.ai == 1, 2],
                  disc.ip[trans.ip.vi == 1, 1],
                  disc.rp[trans.rp.vi == 1, 2])
    
    inf.role <- c(rep(0, sum(trans.ip.ai)), rep(1, sum(trans.rp.ai)), rep(0, sum(trans.ip.vi)), rep(1, sum(trans.rp.vi)))
    inf.type <- c(disc.ip[trans.ip.ai == 1, "ptype"],
                  disc.ip[trans.ip.vi == 1, "ptype"],
                  disc.rp[trans.rp.ai == 1, "ptype"],
                  disc.rp[trans.rp.vi == 1, "ptype"])

    inf.stage <- stage[infector]
    inf.diag <- diag.status[infector]
    inf.tx <- tx.status[infector]

    dat$attr$status[infected] <- 1
    dat$attr$inf.class[infected] <- inf.class[infector]
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- 1
    dat$attr$stage.time[infected] <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected] <- 0

    dat$attr$infector[infected] <- infector
    dat$attr$inf.role[infected] <- inf.role
    dat$attr$inf.type[infected] <- inf.type
    dat$attr$inf.diag[infected] <- inf.diag
    dat$attr$inf.tx[infected] <- inf.tx
    dat$attr$inf.stage[infected] <- inf.stage

    dat$attr$cum.time.on.tx[infected] <- 0
    dat$attr$cum.time.off.tx[infected] <- 0
  }

  # Summary Output
  dat$epi$incid[at] <- sum(dat$epi$incid[at] , length(infected), na.rm = TRUE)
  dat$epi$incid.B[at] <- sum(dat$epi$incid.B[at] , sum(race[infected] == "B", na.rm = TRUE), na.rm = TRUE)
  dat$epi$incid.BI[at] <- sum(dat$epi$incid.BI[at] , sum(race[infected] == "BI", na.rm = TRUE), na.rm = TRUE)
  dat$epi$incid.H[at] <- sum(dat$epi$incid.H[at] , sum(race[infected] == "H", na.rm = TRUE), na.rm = TRUE)
  dat$epi$incid.HI[at] <- sum(dat$epi$incid.HI[at] , sum(race[infected] == "HI", na.rm = TRUE), na.rm = TRUE) 
  dat$epi$incid.W[at] <- sum(dat$epi$incid.W[at] , sum(race[infected] == "W", na.rm = TRUE), na.rm = TRUE) 
  
  dat$epi$incid.f[at] <- sum(dat$epi$incid.f[at] , sum(sex[infected] == "F", na.rm = TRUE), na.rm = TRUE) 
  dat$epi$incid.m[at] <- sum(dat$epi$incid.m[at] , sum(sex[infected] == "M", na.rm = TRUE), na.rm = TRUE)  
  
  dat$epi$incid.msf[at] <- sum(dat$epi$incid.msf[at] , sum(sex.ident[infected] == "msf", na.rm = TRUE), na.rm = TRUE) 
  dat$epi$incid.msm[at] <- sum(dat$epi$incid.msm[at] , sum(sex.ident[infected] == "msm", na.rm = TRUE), na.rm = TRUE)  
  dat$epi$incid.msmf[at] <- sum(dat$epi$incid.msmf[at] , sum(sex.ident[infected] == "msmf", na.rm = TRUE), na.rm = TRUE)

  dat$epi$incid.Lhet[at] <- length(infected)

  return(dat)
}

