
#' @title Transmission Module
#'
#' @description Stochastically simulates disease transmission given the current
#'              state of the discordand edgelist.
#'
#' @inheritParams aging.mard
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
trans.mard <- function(dat, at){

  ## Variables

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
  pce <- dat$param$prep.class.effect

  # Data
  dal <- dat$temp$dal
  dal <- dal[sample(1:nrow(dal)), ]
  ncols <- dim(dal)[2]


  ## Processes

  ## Reorder by role: ins on the left, rec on the right,
  ##                  with flippers represented twice
  disc.ip <- dal[dal$ins %in% c("P", "B"), ]
  disc.rp <- dal[dal$ins %in% c("N", "B"), c(2:1, 3:ncols)]
  names(disc.ip)[1:2] <- c("i", "r"); names(disc.rp)[1:2] <- c("i", "r")


  ## PATP: Insertive Man Infected (Column 1)

  # Attributes of infected
  ip.vl <- vl[disc.ip[, 1]]
  ip.stage <- stage[disc.ip[, 1]]

  # Attributes of susceptible
  ip.ccr5 <- ccr5[disc.ip[, 2]]
  ip.prep <- prepStat[disc.ip[, 2]]
  ip.prepcl <- prepClass[disc.ip[, 2]]

  # Base TP from VL
  trans.ip.prob <- URAI.prob * 2.45^(ip.vl - 4.5)

  # Condom use
  trans.ip.prob[disc.ip$uai == 0] <- trans.ip.prob[disc.ip$uai == 0] * condom.rr

  # CCR5
  trans.ip.prob[ip.ccr5 == "DD"] <- trans.ip.prob[ip.ccr5 == "DD"] * 0
  trans.ip.prob[ip.ccr5 == "DW"] <- trans.ip.prob[ip.ccr5 == "DW"] * ccr5.heteroz.rr

  # PrEP
  trans.ip.prob[which(ip.prep == 1 & ip.prepcl == "l")] <-
            trans.ip.prob[which(ip.prep == 1 & ip.prepcl == "l")] * (1 - pce[1])
  trans.ip.prob[which(ip.prep == 1 & ip.prepcl == "m")] <-
            trans.ip.prob[which(ip.prep == 1 & ip.prepcl == "m")] * (1 - pce[2])
  trans.ip.prob[which(ip.prep == 1 & ip.prepcl == "h")] <-
            trans.ip.prob[which(ip.prep == 1 & ip.prepcl == "h")] * (1 - pce[3])

  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% c("AR", "AF"))
  trans.ip.prob[isAcute] <- trans.ip.prob[isAcute] * acute.rr


  ## PATP: Receptive Man Infected (Column 2)

  # Attributes of infected
  rp.vl <- vl[disc.rp[, 2]]
  rp.stage <- stage[disc.rp[, 2]]

  # Attributes of susceptible
  rp.circ <- circ[disc.rp[, 1]]
  rp.ccr5 <- ccr5[disc.rp[, 1]]
  rp.prep <- prepStat[disc.rp[, 1]]
  rp.prepcl <- prepClass[disc.rp[, 1]]

  # Base TP from VL
  trans.rp.prob <- UIAI.prob * 2.45^(rp.vl - 4.5)

  # Circumcision
  trans.rp.prob[rp.circ == 1] <- trans.rp.prob[rp.circ == 1] * circ.rr

  # Condom use
  trans.rp.prob[disc.rp$uai == 0] <- trans.rp.prob[disc.rp$uai == 0] * condom.rr

  # CCR5
  trans.rp.prob[rp.ccr5 == "DD"] <- trans.rp.prob[rp.ccr5 == "DD"] * 0
  trans.rp.prob[rp.ccr5 == "DW"] <- trans.rp.prob[rp.ccr5 == "DW"] * ccr5.heteroz.rr

  # PrEP
  trans.rp.prob[which(rp.prep == 1 & rp.prepcl == "l")] <-
            trans.rp.prob[which(rp.prep == 1 & rp.prepcl == "l")] * (1 - pce[1])
  trans.rp.prob[which(rp.prep == 1 & rp.prepcl == "m")] <-
            trans.rp.prob[which(rp.prep == 1 & rp.prepcl == "m")] * (1 - pce[2])
  trans.rp.prob[which(rp.prep == 1 & rp.prepcl == "h")] <-
            trans.rp.prob[which(rp.prep == 1 & rp.prepcl == "h")] * (1 - pce[3])

  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% c("AR", "AF"))
  trans.rp.prob[isAcute] <- trans.rp.prob[isAcute] * acute.rr


  ## Bound range of PATP
  trans.ip.prob <- pmin(trans.ip.prob, 1)
  trans.rp.prob <- pmin(trans.rp.prob, 1)

  ## Save to DAL
  disc.ip$prob <- trans.ip.prob
  disc.rp$prob <- trans.rp.prob

  ## Bernoulli transmission events
  trans.ip <- rbinom(length(trans.ip.prob), 1, trans.ip.prob)
  trans.rp <- rbinom(length(trans.rp.prob), 1, trans.rp.prob)


  ## Output

  # Update attributes

  infected <- infector <- inf.type <- NULL
  if (sum(trans.ip, trans.rp) > 0) {

    infected <- c(disc.ip[trans.ip == 1, 2],
                  disc.rp[trans.rp == 1, 1])
    infector <- c(disc.ip[trans.ip == 1, 1],
                  disc.rp[trans.rp == 1, 2])
    inf.role <- c(rep(0, sum(trans.ip)), rep(1, sum(trans.rp)))
    inf.type <- c(disc.ip[trans.ip == 1, "type"],
                  disc.rp[trans.rp == 1, "type"])

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
  dat$epi$incid.B[at] <- sum(race[infected] == "B")
  dat$epi$incid.W[at] <- sum(race[infected] == "W")
  dat$epi$incid.acte[at] <- sum(stage[infector] %in% c("AR", "AF"))
  dat$epi$incid.chrn[at] <- sum(stage[infector] == "C")
  dat$epi$incid.aids[at] <- sum(stage[infector] == "D")
  dat$epi$incid.main[at] <- sum(inf.type == "main")
  dat$epi$incid.casl[at] <- sum(inf.type == "pers")
  dat$epi$incid.inst[at] <- sum(inf.type == "inst")
  dat$epi$incid.prep0[at] <- sum(prepStat[infected] == 0)
  dat$epi$incid.prep1[at] <- sum(prepStat[infected] == 1)

  dat$epi$acts[at] <- nrow(disc.ip) + nrow(disc.rp)
  dat$epi$acts.B[at] <- sum(disc.ip$uai[race[disc.ip$r] == "B"] %in% 0:1) +
                        sum(disc.rp$uai[race[disc.ip$i] == "B"] %in% 0:1)
  dat$epi$acts.W[at] <- sum(disc.ip$uai[race[disc.ip$r] == "W"] %in% 0:1) +
                        sum(disc.rp$uai[race[disc.ip$i] == "W"] %in% 0:1)
  dat$epi$patp[at] <- mean(c(disc.ip$prob, disc.rp$prob))
  dat$epi$patp.B[at] <- mean(c(disc.ip$prob[race[disc.ip$r] == "B"],
                               disc.rp$prob[race[disc.rp$i] == "B"]))
  dat$epi$patp.W[at] <- mean(c(disc.ip$prob[race[disc.ip$r] == "W"],
                               disc.rp$prob[race[disc.rp$i] == "W"]))
  dat$epi$prob.uai[at] <- mean(c(disc.ip$uai, disc.rp$uai))
  dat$epi$prob.uai.B[at] <- mean(c(disc.ip$uai[race[disc.ip$r] == "B"],
                                   disc.rp$uai[race[disc.rp$i] == "B"]))
  dat$epi$prob.uai.W[at] <- mean(c(disc.ip$uai[race[disc.ip$r] == "W"],
                                   disc.rp$uai[race[disc.rp$i] == "W"]))

  return(dat)
}
