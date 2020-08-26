
# MSM -----------------------------------------------------------------

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
#' acute stage infection excess risk, and condom use.
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
hivtrans_msm <- function(dat, at) {

  # Variables -----------------------------------------------------------

  # Attributes
  vl <- dat$attr$vl
  stage <- dat$attr$stage
  circ <- dat$attr$circ
  status <- dat$attr$status
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT
  race <- dat$attr$race
  tx.status <- dat$attr$tx.status

  # Parameters
  URAI.prob <- dat$param$URAI.prob
  UIAI.prob <- dat$param$UIAI.prob
  trans.scale <- dat$param$trans.scale
  acute.rr <- dat$param$acute.rr

  cond.eff <- dat$param$cond.eff
  cond.fail <- dat$param$cond.fail

  circ.rr <- dat$param$circ.rr
  prep.hr <- dat$param$prep.adhr.hr
  hiv.ugc.rr <- dat$param$hiv.ugc.rr
  hiv.uct.rr <- dat$param$hiv.uct.rr
  hiv.rgc.rr <- dat$param$hiv.rgc.rr
  hiv.rct.rr <- dat$param$hiv.rct.rr
  hiv.dual.rr <- dat$param$hiv.dual.rr


  # Data
  al <- dat$temp$al
  dal <- al[which(status[al[, 1]] == 1 & status[al[, 2]] == 0), ]
  dal <- dal[sample(1:nrow(dal)), ]
  ncols <- dim(dal)[2]

  if (nrow(dal) == 0) {
    return(dat)
  }

  ## Reorder by role: ins on the left, rec on the right
  disc.ip <- dal[dal[, "ins"] == 1, ]
  disc.rp <- dal[dal[, "ins"] == 0, c(2:1, 3:ncols)]
  colnames(disc.ip)[1:2] <- colnames(disc.rp)[1:2] <- c("ins", "rec")


  # PATP: Insertive Man Infected (Col 1) --------------------------------

  # Attributes of infected
  ip.vl <- vl[disc.ip[, 1]]
  ip.stage <- stage[disc.ip[, 1]]
  ip.txStat <- tx.status[disc.ip[, 1]]

  # Attributes of susceptible
  ip.prep <- prepStat[disc.ip[, 2]]
  ip.prepcl <- prepClass[disc.ip[, 2]]
  ip.rGC <- rGC[disc.ip[, 2]]
  ip.rCT <- rCT[disc.ip[, 2]]

  # Base TP from VL
  ip.tprob <- pmin(0.99, URAI.prob * 2.45^(ip.vl - 4.5))

  # Adjustment (based on Supervie JAIDS) for VL Suppressed, on ART
  ip.noTrans <- which(ip.vl <= log10(200) & ip.txStat == 1)
  ip.tprob[ip.noTrans] <- 2.2/1e5

  # Transform to log odds
  ip.tlo <- log(ip.tprob/(1 - ip.tprob))

  # Condom use
  not.UAI <- which(disc.ip[, "uai"] == 0)
  condom.rr <- rep(NA, nrow(disc.ip))
  races <- sort(unique(race[disc.ip[, 1]]))
  for (i in races) {
    not.UAI.race <- intersect(not.UAI, which(race[disc.ip[, 1]] == i))
    condom.rr[not.UAI.race] <- 1 - (cond.eff - cond.fail[i])
  }
  ip.tlo[not.UAI] <- ip.tlo[not.UAI] + log(condom.rr[not.UAI])

  # PrEP, by adherence class
  ip.on.prep <- which(ip.prep == 1)
  ip.tlo[ip.on.prep] <- ip.tlo[ip.on.prep] + log(prep.hr[ip.prepcl[ip.on.prep]])

  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% 1:2)
  ip.tlo[isAcute] <- ip.tlo[isAcute] + log(acute.rr)

  ## Multiplier for STI
  is.rGC <- which(ip.rGC == 1)
  is.rCT <- which(ip.rCT == 1)
  is.rect.dual <- intersect(is.rGC, is.rCT)
  is.rGC.sing <- setdiff(is.rGC, is.rect.dual)
  is.rCT.sing <- setdiff(is.rCT, is.rect.dual)
  ip.tlo[is.rGC.sing] <- ip.tlo[is.rGC.sing] + log(hiv.rgc.rr)
  ip.tlo[is.rCT.sing] <- ip.tlo[is.rCT.sing] + log(hiv.rct.rr)
  ip.tlo[is.rect.dual] <- ip.tlo[is.rect.dual] +
    max(log(hiv.rgc.rr), log(hiv.rct.rr)) +
    min(log(hiv.rgc.rr), log(hiv.rct.rr)) * hiv.dual.rr

  # Race-specific scalar for calibration
  races <- race[disc.ip[, 2]]
  ip.tlo <- ip.tlo + log(trans.scale[races])

  # Convert back to probability
  ip.tprob <- plogis(ip.tlo)
  stopifnot(ip.tprob >= 0, ip.tprob <= 1)


  # PATP: Receptive Man Infected (Col 2) --------------------------------

  # Attributes of infected
  rp.vl <- vl[disc.rp[, 2]]
  rp.stage <- stage[disc.rp[, 2]]
  rp.txStat <- tx.status[disc.rp[, 2]]

  # Attributes of susceptible
  rp.circ <- circ[disc.rp[, 1]]
  rp.prep <- prepStat[disc.rp[, 1]]
  rp.prepcl <- prepClass[disc.rp[, 1]]
  rp.uGC <- uGC[disc.rp[, 1]]
  rp.uCT <- uCT[disc.rp[, 1]]

  # Base TP from VL
  rp.tprob <- pmin(0.99, UIAI.prob * 2.45^(rp.vl - 4.5))

  # Adjustment (based on Supervie JAIDS) for VL Suppressed, on ART
  rp.noTrans <- which(rp.vl <= log10(200) & rp.txStat == 1)
  rp.tprob[rp.noTrans] <- 2.2/1e5

  # Transform to log odds
  rp.tlo <- log(rp.tprob/(1 - rp.tprob))

  # Circumcision
  rp.tlo[rp.circ == 1] <- rp.tlo[rp.circ == 1] + log(circ.rr)

  # Condom use
  not.UAI <- which(disc.rp[, "uai"] == 0)
  condom.rr <- rep(NA, nrow(disc.rp))
  races <- sort(unique(race[disc.rp[, 1]]))
  for (i in races) {
    not.UAI.race <- intersect(not.UAI, which(race[disc.rp[, 1]] == i))
    condom.rr[not.UAI.race] <- 1 - (cond.eff - cond.fail[i])
  }
  rp.tlo[not.UAI] <- rp.tlo[not.UAI] + log(condom.rr[not.UAI])

  # PrEP, by adherence class
  rp.on.prep <- which(rp.prep == 1)
  rp.tlo[rp.on.prep] <- rp.tlo[rp.on.prep] + log(prep.hr[rp.prepcl[rp.on.prep]])

  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% 1:2)
  rp.tlo[isAcute] <- rp.tlo[isAcute] + log(acute.rr)

  ## Multiplier for STI
  is.uGC <- which(rp.uGC == 1)
  is.uCT <- which(rp.uCT == 1)
  is.ureth.dual <- intersect(is.uGC, is.uCT)
  is.uGC.sing <- setdiff(is.uGC, is.ureth.dual)
  is.uCT.sing <- setdiff(is.uCT, is.ureth.dual)
  rp.tlo[is.uGC.sing] <- rp.tlo[is.uGC.sing] + log(hiv.ugc.rr)
  rp.tlo[is.uCT.sing] <- rp.tlo[is.uCT.sing] + log(hiv.uct.rr)
  rp.tlo[is.ureth.dual] <- rp.tlo[is.ureth.dual] +
    max(log(hiv.ugc.rr), log(hiv.uct.rr)) +
    min(log(hiv.ugc.rr), log(hiv.uct.rr)) * hiv.dual.rr

  # Race-specific scalar for calibration
  races <- race[disc.rp[, 1]]
  rp.tlo <- rp.tlo + log(trans.scale[races])

  # Convert back to probability
  rp.tprob <- plogis(rp.tlo)
  stopifnot(rp.tprob >= 0, rp.tprob <= 1)


  # Transmission --------------------------------------------------------

  trans.ip <- rbinom(length(ip.tprob), 1, ip.tprob)
  trans.rp <- rbinom(length(rp.tprob), 1, rp.tprob)


  # Output --------------------------------------------------------------

  infected <- NULL
  if (sum(trans.ip, trans.rp) > 0) {
    infected <- c(disc.ip[trans.ip == 1, 2],
                  disc.rp[trans.rp == 1, 1])

    # Attributes of newly infected
    dat$attr$status[infected] <- 1
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- 1
    dat$attr$stage.time[infected] <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected] <- 0
    dat$attr$cuml.time.on.tx[infected] <- 0
    dat$attr$cuml.time.off.tx[infected] <- 0

    # Attributes of transmitter
    transmitter <- as.numeric(c(disc.ip[trans.ip == 1, 1],
                                disc.rp[trans.rp == 1, 2]))
    tab.trans <- table(transmitter)
    uni.trans <- as.numeric(names(tab.trans))
    dat$attr$count.trans[uni.trans] <- dat$attr$count.trans[uni.trans] +
                                       as.numeric(tab.trans)
   }

  # Summary Output
  dat$epi$incid[at] <- length(infected)
  dat$epi$incid.B[at] <- sum(dat$attr$race[infected] == 1)
  dat$epi$incid.H[at] <- sum(dat$attr$race[infected] == 2)
  dat$epi$incid.W[at] <- sum(dat$attr$race[infected] == 3)

  if (length(infected) > 0) {
    dat$epi$incid.undx[at] <- sum(dat$attr$diag.status[transmitter] == 0)
    dat$epi$incid.dx[at] <- sum(dat$attr$diag.status[transmitter] == 1 &
                                dat$attr$cuml.time.on.tx[transmitter] == 0)
    dat$epi$incid.linked[at] <- sum(dat$attr$diag.status[transmitter] == 1 &
                                    dat$attr$cuml.time.on.tx[transmitter] > 0 &
                                    dat$attr$vl[transmitter] > log10(200))
    dat$epi$incid.vsupp[at] <- sum(dat$attr$diag.status[transmitter] == 1 &
                                   dat$attr$cuml.time.on.tx[transmitter] > 0 &
                                   dat$attr$vl[transmitter] <= log10(200))
  } else {
    dat$epi$incid.undx[at] <- 0
    dat$epi$incid.dx[at] <- 0
    dat$epi$incid.linked[at] <- 0
    dat$epi$incid.vsupp[at] <- 0
  }

  return(dat)
}


#' @export
#' @rdname hivtrans_msm
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

    el <- dat$el[[1]]

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

keep.attr <- function(attrList, keep) {
  lapply(attrList, function(x) x[keep])
}

nbsdtosize <- function(mu, sd) {
  mu ^ 2 / (sd ^ 2 - mu)
}
