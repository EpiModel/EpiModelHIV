
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

  # Variables ---------------------------------------------------------------

  # Attributes
  vl <- dat$attr$vl
  stage <- dat$attr$stage
  stage.time <- dat$attr$stage.time
  ccr5 <- dat$attr$ccr5
  circ <- dat$attr$circ
  diag.status <- dat$attr$diag.status
  tx.status <- dat$attr$tx.status
  race <- dat$attr$race
  age <- dat$attr$age

  # Parameters
  betabase.URAI <- dat$param$betabase.URAI
  betabase.UIAI <- dat$param$betabase.UIAI
  betamult.acute <- dat$param$betamult.acute
  vl.acute.fall.dur <- dat$param$vl.acute.fall.dur
  betamult.condom <- dat$param$betamult.condom
  betamult.circ <- dat$param$betamult.circ
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr

  # Data
  if (dat$control$save.dal == TRUE) {
    dal <- dat$temp$dal[[at]]
  } else {
    dal <- dat$temp$dal
  }
  dal <- dal[sample(1:nrow(dal)), ]
  ncols <- dim(dal)[2]


  # Processes ---------------------------------------------------------------

  # Reorder by role: ins on the left, rec on the right,
  #                  with flippers represented twice
  disc.ip <- dal[dal$ins %in% c("P", "B"), ]
  disc.rp <- dal[dal$ins %in% c("N", "B"), c(2:1, 3:ncols)]
  names(disc.ip)[1:2] <- c("i", "r")
  names(disc.rp)[1:2] <- c("i", "r")

  # Transmission probability: insertive infected
  ip.vl <- vl[disc.ip[, 1]]
  ip.stage <- stage[disc.ip[, 1]]
  ip.stage.time <- stage.time[disc.ip[, 1]]
  ip.ccr5 <- ccr5[disc.ip[, 2]]

  tp.ip <- betabase.URAI * 2.45 ^ (ip.vl - 4.5)
  tp.ip[ip.stage == "AR"] <- tp.ip[ip.stage == "AR"] * betamult.acute
  tp.ip[ip.stage == "AF"] <- tp.ip[ip.stage == "AF"] *
    (1 + (betamult.acute - 1) * (vl.acute.fall.dur - ip.stage.time[ip.stage == "AF"]) /
                                 vl.acute.fall.dur)
  tp.ip[disc.ip$uai == 0] <- tp.ip[disc.ip$uai == 0] * betamult.condom
  tp.ip[ip.ccr5 == "DD"] <- tp.ip[ip.ccr5 == "DD"] * 0
  tp.ip[ip.ccr5 == "DW"] <- tp.ip[ip.ccr5 == "DW"] * ccr5.heteroz.rr

  # Transmission probability: receptive position
  rp.vl <- vl[disc.rp[, 2]]
  rp.stage <- stage[disc.rp[, 2]]
  rp.stage.time <- stage.time[disc.rp[, 2]]
  rp.circ <- circ[disc.rp[, 1]]
  rp.ccr5 <- ccr5[disc.rp[, 1]]

  tp.rp <- betabase.UIAI * 2.45 ^ (rp.vl - 4.5)
  tp.rp[rp.stage == "AR"] <- tp.rp[rp.stage == "AR"] * betamult.acute
  tp.rp[rp.stage == "AF"] <- tp.rp[rp.stage == "AF"] *
    (1 + (betamult.acute - 1) * (vl.acute.fall.dur - rp.stage.time[rp.stage == "AF"] ) /
                                 vl.acute.fall.dur)
  tp.rp[rp.circ == 1] <- tp.rp[rp.circ == 1] * betamult.circ
  tp.rp[disc.rp$uai == 0] <- tp.rp[disc.rp$uai == 0] * betamult.condom
  tp.rp[rp.ccr5 == "DD"] <- tp.rp[rp.ccr5 == "DD"] * 0
  tp.rp[rp.ccr5 == "DW"] <- tp.rp[rp.ccr5 == "DW"] * ccr5.heteroz.rr

  stopifnot(min(tp.ip) >= 0, max(tp.ip) <= 1,
            min(tp.rp) >= 0, max(tp.rp) <= 1)

  # Stochastic transmission
  trans.ip <- rbinom(length(tp.ip), 1, tp.ip)
  trans.rp <- rbinom(length(tp.rp), 1, tp.rp)


  # Output ------------------------------------------------------------------

  # Update attributes

  infected <- NULL
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
  dat$epi$incid.y[at] <- sum(age[infected] < 30)
  dat$epi$incid.o[at] <- sum(age[infected] >= 30)
  dat$epi$incid.B.y[at] <- sum(race[infected] == "B" & age[infected] < 30)
  dat$epi$incid.B.o[at] <- sum(race[infected] == "B" & age[infected] >= 30)
  dat$epi$incid.W.y[at] <- sum(race[infected] == "W" & age[infected] < 30)
  dat$epi$incid.W.o[at] <- sum(race[infected] == "W" & age[infected] >= 30)

  return(dat)
}
