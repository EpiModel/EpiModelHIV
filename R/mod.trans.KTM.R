
# SHAMP tranmission for f, msf, msm, msmf populations  -----------------------------------------------------------------

#' @title Transmission Module for AI and VI
#'
#' @description Stochastically simulates disease transmission given the current
#'              state of the discordand edgelist.
#'
#' @inheritParams aging_shamp
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
trans_KTM <- function(dat, at){

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
  age<-dat$attr$age
  race <- dat$attr$race
  sex.ident<-dat$attr$sex.ident
  inf.class <- dat$attr$inf.class

  # Parameters
  URVI.prob <- dat$param$URVI.prob * dat$param$VI.foi.scale
  UIVI.prob <- dat$param$UIVI.prob * dat$param$VI.foi.scale
  
  acute.rr <- dat$param$acute.rr
  condom.rr <- dat$param$condom.rr
  circ.rr <- dat$param$circ.rr
  
  
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr
  prep.hr <- dat$param$prep.class.hr

  # Data
  
  dal <- dat$temp$dal
  
  if (nrow(dal) == 0) {return(dat)}
  
  if (nrow(dal) == 1) {
    dal <- dal[sample(1:nrow(dal)), ]
    dal<-as.matrix(t(dal))
  }

  if (nrow(dal) > 1) {dal <- dal[sample(1:nrow(dal)), ]}
  
  ncols <- dim(dal)[2]
  


  ## Reorder by role: ins on the left, rec on the right,
  ##                  with flippers represented twice
  disc.ip <- as.data.frame(dal[dal[, "ins"] %in% 1:2,, drop=FALSE ])
  disc.rp <- as.data.frame(dal[dal[, "ins"] %in% c(0, 2), c(2:1, 3:ncols), drop=FALSE])
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
  ip.tprob.vi <- URVI.prob * 2.45^(ip.vl - 4.5)
  
  #set base to zero for non-event (AI/VI MM/HET)
  vi<-which(disc.ip$het==1)
  
  
  # Transform to log odds
  ip.tlo.vi <- log(ip.tprob.vi/(1-ip.tprob.vi))
  
  # Condom use
  not.UVI <- which(disc.ip[, "uvi"] == 0)
  
  ip.tlo.vi[not.UVI] <- ip.tlo.vi[not.UVI] + log(condom.rr)
  
  # CCR5
  ip.tlo.vi[ip.ccr5 == "DD"] <- ip.tlo.vi[ip.ccr5 == "DD"] + -Inf
  ip.tlo.vi[ip.ccr5 == "DW"] <- ip.tlo.vi[ip.ccr5 == "DW"] + log(ccr5.heteroz.rr)

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(ip.prep == 1 & ip.prepcl == i-1)
    ip.tlo.vi[temp.ids] <- ip.tlo.vi[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% c(1, 2))
  ip.tlo.vi[isAcute] <- ip.tlo.vi[isAcute] + log(acute.rr)

  # Retransformation to probability
  ip.tprob.vi <- plogis(ip.tlo.vi)
  
  
 if(any(ip.tprob.vi < 0 | ip.tprob.vi > 1 | ip.tprob.vi == "NaN")) {warning("VI Transmission probability not 0-1, set to min0, max1, suggest reduce baseline prob.",call. = FALSE)}
  
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
  rp.tprob.vi <- UIVI.prob * 2.45^(rp.vl - 4.5)
  
  #set base to zero for non-event (AI/VI MM/HET)
  vi<-which(disc.rp$het==1)
  
  
  
  # Transform to log odds
  rp.tlo.vi <- log(rp.tprob.vi/(1-rp.tprob.vi))
  
  # Circumcision

  rp.tlo.vi[rp.circ == 1] <- rp.tlo.vi[rp.circ == 1] + log(circ.rr)
  
  # Condom use
  not.UVI <- which(disc.rp[, "uvi"] == 0)
  rp.tlo.vi[not.UVI] <- rp.tlo.vi[not.UVI] + log(condom.rr)

  # CCR5
  rp.tlo.vi[rp.ccr5 == "DD"] <- rp.tlo.vi[rp.ccr5 == "DD"] + -Inf
  rp.tlo.vi[rp.ccr5 == "DW"] <- rp.tlo.vi[rp.ccr5 == "DW"] + log(ccr5.heteroz.rr)

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(rp.prep == 1 & rp.prepcl == i-1)
    rp.tlo.vi[temp.ids] <- rp.tlo.vi[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% c(1, 2))
  rp.tlo.vi[isAcute] <- rp.tlo.vi[isAcute] + log(acute.rr)
  
  # Retransformation to probability
  rp.tprob.vi <- plogis(rp.tlo.vi)
  
  if(any(rp.tprob.vi < 0 | rp.tprob.vi > 1 | rp.tprob.vi== "NaN")) {warning("VI Transmission probability not 0-1, set to min0, max1, suggest reduce baseline prob.",call. = FALSE)}

  
  rp.tprob.vi <- ifelse(rp.tprob.vi >= 1,1,
                        ifelse(rp.tprob.vi <= 0,0,rp.tprob.vi))
  rp.tprob.vi <- ifelse(is.na(rp.tprob.vi) == TRUE,1,rp.tprob.vi)

  # Transmission --------------------------------------------------------

  ## Bernoulli transmission events

  trans.ip.vi <- rbinom(length(ip.tprob.vi), 1, ip.tprob.vi)
  trans.rp.vi <- rbinom(length(rp.tprob.vi), 1, rp.tprob.vi)


  # Output --------------------------------------------------------------

  # Update attributes
  
## Forward transmssion from a heatbath and foreign aquired get seperate designation.
  inf.class<-ifelse(inf.class=="MSM","MSMds",inf.class)
  inf.class<-ifelse(inf.class=="FA","FAds",inf.class)

  infected <- infector <- inf.type <- NULL
  if (sum(trans.ip.vi, trans.rp.vi) > 0) {

    infected <- c(disc.ip[trans.ip.vi == 1, 2],
                     disc.rp[trans.rp.vi == 1, 1])
    
    
    
    infector <- c(disc.ip[trans.ip.vi == 1, 1],
                  disc.rp[trans.rp.vi == 1, 2])
    
    inf.role <- c(rep(0, sum(trans.ip.vi)), rep(1, sum(trans.rp.vi)))
    inf.type <- c(disc.ip[trans.ip.vi == 1, "ptype"],
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
    dat$attr$age.inf[infected] <- dat$attr$age[infected]
    
    for(k in 1:length(infected)){
    dat$attr$infected.gen[infected[k]] <- dat$attr$infected.gen[infector[k]]+1
    }
    

    
    dat$trans.el$infected <- c(dat$trans.el$infected, dat$attr$uid[infected])
    
    dat$trans.el$infected.sex <- c(dat$trans.el$infected.sex, dat$attr$sex[infected])
    dat$trans.el$infected.race <- c(dat$trans.el$infected.race, dat$attr$race[infected])
    dat$trans.el$infected.age <- c(dat$trans.el$infected.age, dat$attr$age[infected])
    dat$trans.el$infected.sex.ident <- c(dat$trans.el$infected.sex.ident, dat$attr$sex.ident[infected])
    dat$trans.el$infected.deg.tot <- c(dat$trans.el$infected.deg.tot, dat$attr$deg.tot[infected])
    
    new.gens <- dat$attr$infected.gen[infector]
    new.gens <- new.gens+1
    dat$trans.el$infected.gen <- c(dat$trans.el$infected.gen, new.gens)
    
    
    
    dat$trans.el$infector <- c(dat$trans.el$infector, dat$attr$uid[infector])
    dat$trans.el$infector.sex <- c(dat$trans.el$infector.sex, dat$attr$sex[infector])
    dat$trans.el$infector.race <- c(dat$trans.el$infector.race, dat$attr$race[infector])
    dat$trans.el$infector.age <- c(dat$trans.el$infector.age, dat$attr$age[infector])
    dat$trans.el$infector.sex.ident <- c(dat$trans.el$infector.sex.ident, dat$attr$sex.ident[infector])
    dat$trans.el$infector.gen <- c(dat$trans.el$infector.gen, dat$attr$infected.gen[infector])
    dat$trans.el$infector.deg.tot <- c(dat$trans.el$infector.deg.tot, dat$attr$deg.tot[infector])    

    dat$trans.el$time <- c(dat$trans.el$time, rep(at, length(infected)))
                           
  }

  # Summary Output

  dat$epi$incid[at] <- sum(dat$epi$incid[at] , length(infected), na.rm = TRUE)
  dat$epi$incid.f[at] <- max(0, sum(sex[infected] == "F", na.rm = TRUE), na.rm = TRUE) 
  dat$epi$incid.m[at] <- max(0 , sum(sex[infected] == "M", na.rm = TRUE), na.rm = TRUE) 
  
  dat$epi$incid.poi[at] <- max(0, sum(age[infected] < 40, na.rm = TRUE), na.rm = TRUE)
  age.poi <- age[infected] 
  poi <- which(age.poi < 40)
  age.poi <- age.poi[poi]
  dat$epi$incid.age.poi[at] <-max(0,mean(age.poi, na.rm = TRUE))
  
  

  return(dat)
}

