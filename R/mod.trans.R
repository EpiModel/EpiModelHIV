
#' @title Transmission Module
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
#' @keywords module
#' @export
#'
trans_msm <- function(dat, at){

  # Variables -----------------------------------------------------------

  # Attributes
  vl <- dat$attr$vl
  stage <- dat$attr$stage
  ccr5 <- dat$attr$ccr5
  circ <- dat$attr$circ
  diag.status <- dat$attr$diag.status
  tx.status <- dat$attr$tx.status
  race <- dat$attr$race
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  # Parameters
  URAI.prob <- dat$param$URAI.prob
  UIAI.prob <- dat$param$UIAI.prob
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
  disc.ip <- dal[dal[, "ins"] %in% 1:2, ]
  disc.rp <- dal[dal[, "ins"] %in% c(0, 2), c(2:1, 3:ncols)]
  colnames(disc.ip)[1:2] <- c("i", "r")
  colnames(disc.rp)[1:2] <- c("i", "r")


  # PATP: Insertive Man Infected (Col 1) --------------------------------

  # Attributes of infected
  ip.vl <- vl[disc.ip[, 1]]
  ip.stage <- stage[disc.ip[, 1]]

  # Attributes of susceptible
  ip.ccr5 <- ccr5[disc.ip[, 2]]
  ip.prep <- prepStat[disc.ip[, 2]]
  ip.prepcl <- prepClass[disc.ip[, 2]]

  # Base TP from VL
  ip.tprob <- URAI.prob * 2.45^(ip.vl - 4.5)
  
  # Transform to log odds
  ip.tlo <- log(ip.tprob/(1-ip.tprob))
  
  # Condom use
  not.UAI <- which(disc.ip[, "uai"] == 0)
  ip.tlo[not.UAI] <- ip.tlo[not.UAI] + log(condom.rr/(1-condom.rr))
  
  # CCR5
  ip.tlo[ip.ccr5 == "DD"] <- ip.tlo[ip.ccr5 == "DD"] + -Inf
  ip.tlo[ip.ccr5 == "DW"] <- ip.tlo[ip.ccr5 == "DW"] + 
                             log(ccr5.heteroz.rr/(1-ccr5.heteroz.rr))
  
  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(ip.prep == 1 & ip.prepcl == i-1)
    ip.tlo[temp.ids] <- ip.tlo[temp.ids] + log(prep.hr[i])
  }
  
  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% c("AR", "AF"))
  ip.tlo[isAcute] <- ip.tlo[isAcute] + log(acute.rr)
  
  # Retransformation to probability
  ip.tprob <- exp(ip.tlo)/(1+exp(ip.tlo))
  stopifnot(ip.tprob >= 0, ip.tprob <= 1)
  
  
  # PATP: Receptive Man Infected (Col 2) --------------------------------

  # Attributes of infected
  rp.vl <- vl[disc.rp[, 2]]
  rp.stage <- stage[disc.rp[, 2]]

  # Attributes of susceptible
  rp.circ <- circ[disc.rp[, 1]]
  rp.ccr5 <- ccr5[disc.rp[, 1]]
  rp.prep <- prepStat[disc.rp[, 1]]
  rp.prepcl <- prepClass[disc.rp[, 1]]

  # Base TP from VL
  rp.tprob <- UIAI.prob * 2.45^(rp.vl - 4.5)
  
  # Transform to log odds
  rp.tlo <- log(rp.tprob/(1-rp.tprob))
  
  # Circumcision
  rp.tlo[rp.circ == 1] <- rp.tlo[rp.circ == 1] + log(circ.rr/(1-circ.rr))
  
  # Condom use
  not.UAI <- which(disc.rp[, "uai"] == 0)
  rp.tlo[not.UAI] <- rp.tlo[not.UAI] + log(condom.rr/(1-condom.rr))
  
  # CCR5
  rp.tlo[rp.ccr5 == "DD"] <- rp.tlo[rp.ccr5 == "DD"] + -Inf
  rp.tlo[rp.ccr5 == "DW"] <- rp.tlo[rp.ccr5 == "DW"] +
    log(ccr5.heteroz.rr/(1-ccr5.heteroz.rr))
  
  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(rp.prep == 1 & rp.prepcl == i-1)
    rp.tlo[temp.ids] <- rp.tlo[temp.ids] + log(prep.hr[i])
  }
  
  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% c("AR", "AF"))
  rp.tlo[isAcute] <- rp.tlo[isAcute] + log(acute.rr)

  # Retransformation to probability
  rp.tprob <- exp(rp.tlo)/(1+exp(rp.tlo))
  stopifnot(rp.tprob >= 0, rp.tprob <= 1)

  # Transmission --------------------------------------------------------

  ## Bernoulli transmission events
  trans.ip <- rbinom(length(ip.tprob), 1, ip.tprob)
  trans.rp <- rbinom(length(rp.tprob), 1, rp.tprob)


  # Output --------------------------------------------------------------

  # Update attributes

  infected <- infector <- inf.type <- NULL
  if (sum(trans.ip, trans.rp) > 0) {

    infected <- c(disc.ip[trans.ip == 1, 2],
                  disc.rp[trans.rp == 1, 1])
    infector <- c(disc.ip[trans.ip == 1, 1],
                  disc.rp[trans.rp == 1, 2])
    inf.role <- c(rep(0, sum(trans.ip)), rep(1, sum(trans.rp)))
    inf.type <- c(disc.ip[trans.ip == 1, "ptype"],
                  disc.rp[trans.rp == 1, "ptype"])

    inf.stage <- stage[infector]
    inf.diag <- diag.status[infector]
    inf.tx <- tx.status[infector]

    dat$attr$status[infected] <- 1
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- "AR"
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
  dat$epi$incid[at] <- length(infected)


  return(dat)
}
