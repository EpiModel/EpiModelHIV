
# MSM -----------------------------------------------------------------

#' @title Transmission Module
#'
#' @description Stochastically simulates disease transmission given the current
#'              state of the discordand edgelist.
#'
#' @inheritParams aging_camplc
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
#' @keywords module msm
#'
#' @export
#'
trans_msm <- function(dat, at){

  ## Variables

  # Attributes
  vl <- dat$attr$vl
  stage <- dat$attr$stage
  ccr5 <- dat$attr$ccr5
  circ <- dat$attr$circ
  diag.status <- dat$attr$diag.status
  tx.status <- dat$attr$tx.status
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


  ## Processes

  ## Reorder by role: ins on the left, rec on the right,
  ##                  with flippers represented twice
  disc.ip <- dal[dal[, "ins"] %in% 1:2, ]
  disc.rp <- dal[dal[, "ins"] %in% c(0, 2), c(2:1, 3:ncols)]
  colnames(disc.ip)[1:2] <- c("i", "r")
  colnames(disc.rp)[1:2] <- c("i", "r")


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
  trans.ip.prob[disc.ip[, "uai"] == 0] <- trans.ip.prob[disc.ip[, "uai"] == 0] * condom.rr

  # CCR5
  trans.ip.prob[ip.ccr5 == "DD"] <- trans.ip.prob[ip.ccr5 == "DD"] * 0
  trans.ip.prob[ip.ccr5 == "DW"] <- trans.ip.prob[ip.ccr5 == "DW"] * ccr5.heteroz.rr

  # PrEP
  trans.ip.prob[which(ip.prep == 1 & ip.prepcl == 0)] <-
    trans.ip.prob[which(ip.prep == 1 & ip.prepcl == 0)] * prep.hr[1]
  trans.ip.prob[which(ip.prep == 1 & ip.prepcl == 1)] <-
    trans.ip.prob[which(ip.prep == 1 & ip.prepcl == 1)] * prep.hr[2]
  trans.ip.prob[which(ip.prep == 1 & ip.prepcl == 2)] <-
    trans.ip.prob[which(ip.prep == 1 & ip.prepcl == 2)] * prep.hr[3]
  trans.ip.prob[which(ip.prep == 1 & ip.prepcl == 3)] <-
    trans.ip.prob[which(ip.prep == 1 & ip.prepcl == 3)] * prep.hr[4]

  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% c("AR", "AF"))
  trans.ip.prob[isAcute] <- trans.ip.prob[isAcute] * acute.rr

  ## TODO: multiplier for prevalent STI infection
  ## Need to do this twice, given differentials in multipliers

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
  trans.rp.prob[disc.rp[, "uai"] == 0] <- trans.rp.prob[disc.rp[, "uai"] == 0] * condom.rr

  # CCR5
  trans.rp.prob[rp.ccr5 == "DD"] <- trans.rp.prob[rp.ccr5 == "DD"] * 0
  trans.rp.prob[rp.ccr5 == "DW"] <- trans.rp.prob[rp.ccr5 == "DW"] * ccr5.heteroz.rr

  # PrEP
  trans.rp.prob[which(rp.prep == 1 & rp.prepcl == 0)] <-
    trans.rp.prob[which(rp.prep == 1 & rp.prepcl == 0)] * prep.hr[1]
  trans.rp.prob[which(rp.prep == 1 & rp.prepcl == 1)] <-
    trans.rp.prob[which(rp.prep == 1 & rp.prepcl == 1)] * prep.hr[2]
  trans.rp.prob[which(rp.prep == 1 & rp.prepcl == 2)] <-
    trans.rp.prob[which(rp.prep == 1 & rp.prepcl == 2)] * prep.hr[3]
  trans.rp.prob[which(rp.prep == 1 & rp.prepcl == 3)] <-
    trans.rp.prob[which(rp.prep == 1 & rp.prepcl == 3)] * prep.hr[4]

  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% c("AR", "AF"))
  trans.rp.prob[isAcute] <- trans.rp.prob[isAcute] * acute.rr


  ## Bound range of PATP
  trans.ip.prob <- pmin(trans.ip.prob, 1)
  trans.rp.prob <- pmin(trans.rp.prob, 1)


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
    
    dat$age.inf.vec$age <- c(dat$age.inf.vec$age,dat$attr$age[unique(infected)])
    dat$age.inf.vec$time <- c(dat$age.inf.vec$time,rep(at,length(unique(infected))))
    
  }

  # Summary Output
  dat$epi$incid[at] <- length(unique(infected))
  
  msm<-sum(dat$attr$asmm[unique(infected)] == 0)
  asmm<-sum(dat$attr$asmm[unique(infected)] == 1)
  
  dat$epi$incid.msm[at] <- msm
  dat$epi$incid.asmm[at] <- asmm
  

  return(dat)
}



# HET -----------------------------------------------------------------

#' @title Infection Module
#'
#' @description Module function to simulate transmission over an active
#'              discordant edgelist.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
#'
trans_het <- function(dat, at) {

  ## Discordant Edgelist
  del <- discord_edgelist_het(dat, at)

  nInf <- 0
  idsInf <- idsTrans <- NULL

  if (!is.null(del)) {

    ## Acts
    nedges <- length(del[[1]])

    act.rate.early <- dat$param$act.rate.early
    act.rate.late <- dat$param$act.rate.late
    act.rate.cd4 <- dat$param$act.rate.cd4

    cd4Count <- dat$attr$cd4Count[del$inf]

    isLate <- which(cd4Count < act.rate.cd4)

    rates <- rep(act.rate.early, nedges)
    rates[isLate] <- act.rate.late


    # Process
    act.rand <- dat$param$acts.rand
    if (act.rand == TRUE) {
      numActs <- rpois(nedges, rates)
    } else {
      numActs <- rates
    }

    cond.prob <- dat$param$cond.prob
    cond.prob <- rep(cond.prob, nedges)

    del$numActs <- numActs

    if (act.rand == TRUE) {
      del$protActs <- rbinom(nedges, rpois(nedges, numActs), cond.prob)
    } else {
      del$protActs <- numActs * cond.prob
    }

    del$protActs <- pmin(numActs, del$protActs)
    del$unprotActs <- numActs - del$protActs

    stopifnot(all(del$unprotActs >= 0))


    ## Transmission

    # Base transmission probability
    vlLevel <- dat$attr$vlLevel[del$inf]
    males <- dat$attr$male[del$sus]
    ages <- dat$attr$age[del$sus]
    circs <- dat$attr$circStat[del$sus]
    prop.male <- dat$epi$propMale[at - 1]
    base.tprob <- hughes_tp(vlLevel, males, ages, circs, prop.male)

    # Acute and aids stage multipliers
    acute.stage.mult <- dat$param$acute.stage.mult
    aids.stage.mult <- dat$param$aids.stage.mult

    isAcute <- which(at - dat$attr$infTime[del$inf] <
                       (dat$param$vl.acute.topeak + dat$param$vl.acute.toset))
    isAIDS <- which(dat$attr$cd4Count[del$inf] < 200)

    base.tprob[isAcute] <- base.tprob[isAcute] * acute.stage.mult
    base.tprob[isAIDS] <- base.tprob[isAIDS] * aids.stage.mult


    # Condoms
    # Probability as a mixture function of protected and unprotected acts
    cond.eff <- dat$param$cond.eff
    prob.stasis.protacts <- (1 - base.tprob*(1 - cond.eff)) ^ del$protActs
    prob.stasis.unptacts <- (1 - base.tprob) ^ del$unprotActs
    prob.stasis <- prob.stasis.protacts * prob.stasis.unptacts
    finl.tprob <- 1 - prob.stasis

    # Transmission
    del$base.tprob <- base.tprob
    del$finl.tprob <- finl.tprob

    stopifnot(length(unique(sapply(del, length))) == 1)

    # Random transmission given final trans prob
    idsTrans <- which(rbinom(nedges, 1, del$finl.tprob) == 1)

    # Subset discord edgelist to transmissions
    del <- keep.attr(del, idsTrans)


    ## Update Nodal Attr
    idsInf <- unique(del$sus)
    idsTrans <- unique(del$inf)
    nInf <- length(idsInf)

    if (nInf > 0) {
      dat$attr$status[idsInf] <- 1
      dat$attr$infTime[idsInf] <- at
      dat$attr$ageInf[idsInf] <- dat$attr$age[idsInf]
      dat$attr$dxStat[idsInf] <- 0
      dat$attr$vlLevel[idsInf] <- 0
      dat$attr$txCD4min[idsInf] <-
        pmin(rnbinom(nInf,
                     size = nbsdtosize(dat$param$tx.init.cd4.mean,
                                       dat$param$tx.init.cd4.sd),
                     mu = dat$param$tx.init.cd4.mean),
             dat$param$tx.elig.cd4)
    }

    ## Transmission data frame
    if (dat$control$save.transmat == TRUE) {
      if (nInf > 0) {
        if (at == 2) {
          dat$stats$transmat <- as.data.frame(del)
        } else {
          dat$stats$transmat <- rbind(dat$stats$transmat, as.data.frame(del))
        }
      }
    }

  }

  ## Incidence vector
  dat$epi$si.flow[at] <- nInf
  dat$epi$si.flow.male[at] <- sum(dat$attr$male[idsInf] == 1, na.rm = TRUE)
  dat$epi$si.flow.feml[at] <- sum(dat$attr$male[idsInf] == 0, na.rm = TRUE)

  return(dat)
}


discord_edgelist_het <- function(dat, at) {

  status <- dat$attr$status

  idsInft <- which(status == 1)
  nInft <- length(idsInft)

  del <- NULL

  if (nInft > 0) {

    if (is.null(dat$el)) {
      el <- get.dyads.active(dat$nw, at = at)
    } else {
      el <- dat$el
    }

    if (nrow(el) > 0) {
      el <- el[sample(1:nrow(el)), , drop = FALSE]

      disc <- which(abs(status[el[, 1]] - status[el[, 2]]) == 1)
      if (length(disc) > 0) {
        tmp.del <- el[disc, ]
        tmp.del[status[tmp.del[, 2]] == 1, ] <- tmp.del[status[tmp.del[, 2]] == 1, 2:1]

        del <- list()
        del$sus <- tmp.del[, 2]
        del$inf <- tmp.del[, 1]
      }
    }

  }

  return(del)
}


hughes_tp <- function(vls, susmales, susages, suscircs, prop.male, fmat = FALSE) {

  suscircs[is.na(suscircs)] <- 0

  sus.hsv2 <- 0.59*prop.male + 0.86*(1 - prop.male)
  sus.gud <- 0.039*prop.male + 0.053*(1 - prop.male)
  sus.tvagin <- 0.068*prop.male + 0.12*(1 - prop.male)
  sus.cerv <- 0.066*(1 - prop.male)

  interc <- -8.3067
  coef.vl <- 1.062566
  coef.male <- 0.6430989
  coef.age <- -0.0403451
  coef.hsv2 <- 0.7625081
  coef.circ <- -0.6377294
  coef.gud <- 0.9749536
  coef.vagin <- 0.9435334
  coef.cerv <- 1.288279

  tp.full <- exp(interc + coef.vl*(vls - 4) +
                   coef.male*susmales + coef.age*(susages - 35) +
                   coef.hsv2*sus.hsv2 + coef.circ*susmales*suscircs +
                   coef.gud*sus.gud + coef.vagin*sus.tvagin +
                   coef.cerv*sus.cerv)

  if (fmat == TRUE) {
    tp.full <- data.frame(tp.full, vls, susmales, susages, suscircs)
  }

  return(tp.full)
}
