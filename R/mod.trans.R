
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
hiv_trans_msm <- function(dat, at) {

  # Variables -----------------------------------------------------------

  # Attributes
  vl <- dat$attr$vl
  stage <- dat$attr$stage
  ccr5 <- dat$attr$ccr5
  circ <- dat$attr$circ
  status <- dat$attr$status
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT
  stage.syph <- dat$attr$stage.syph

  # Parameters
  URAI.prob <- dat$param$URAI.prob
  UIAI.prob <- dat$param$UIAI.prob

  acute.rr <- dat$param$acute.rr
  condom.rr <- dat$param$condom.rr
  circ.rr <- dat$param$circ.rr
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr
  prep.hr <- dat$param$prep.class.hr
  hiv.ugc.rr <- dat$param$hiv.ugc.rr
  hiv.uct.rr <- dat$param$hiv.uct.rr
  hiv.rgc.rr <- dat$param$hiv.rgc.rr
  hiv.rct.rr <- dat$param$hiv.rct.rr
  hiv.syph.rr <- dat$param$hiv.syph.rr
  hiv.rgc.rct.rr <- dat$param$hiv.rgc.rct.rr
  hiv.rgc.syph.rr <- dat$param$hiv.rgc.syph.rr
  hiv.rct.syph.rr <- dat$param$hiv.rct.syph.rr
  hiv.ugc.uct.rr <- dat$param$hiv.ugc.uct.rr
  hiv.ugc.syph.rr <- dat$param$hiv.ugc.syph.rr
  hiv.uct.syph.rr <- dat$param$hiv.uct.syph.rr
  hiv.all.ureth.rr <- dat$param$hiv.all.ureth.rr
  hiv.all.rect.rr <- dat$param$hiv.all.rect.rr

  hiv.trans.syph.rr <- dat$param$hiv.trans.syph.rr
  hiv.trans.gc.rr <- dat$param$hiv.trans.gc.rr
  hiv.trans.ct.rr <- dat$param$hiv.trans.ct.rr
  hiv.trans.gc.ct.rr <- dat$param$hiv.trans.ct.rr
  hiv.trans.gc.syph.rr <- dat$param$hiv.trans.gc.syph.rr
  hiv.trans.ct.syph.rr <- dat$param$hiv.trans.ct.syph.rr
  hiv.trans.allsti.rr <- dat$param$hiv.trans.allsti.rr


  # Data
  al <- dat$temp$al
  dal <- al[which(status[al[, 1]] == 1 & status[al[, 2]] == 0), ]
  dal <- dal[sample(1:nrow(dal)), ]
  ncols <- dim(dal)[2]

  al <- cbind(al, st1 = as.vector(dat$attr$status[al[ ,"p1"]]))
  al <- cbind(al, st2 = as.vector(dat$attr$status[al[ ,"p2"]]))

  al.negneg <- al[al[, "st1"] == 0 & al[, "st2"] == 0, , drop = FALSE]
  al.negpos <- al[(al[, "st1"] == 1 & al[, "st2"] == 0) |
                    (al[, "st1"] == 0 & al[, "st2"] == 1), , drop = FALSE]
  al.pospos <- al[al[, "st1"] == 1 & al[, "st2"] == 1, , drop = FALSE]

  num.acts.negneg <- nrow(al.negneg)
  num.acts.negpos <- nrow(al.negpos)
  num.acts.pospos <- nrow(al.pospos)

  prop.uai.negneg <- sum(al.negneg[, "uai"] == 1) / nrow(al.negneg)
  prop.uai.negpos <- sum(al.negpos[, "uai"] == 1) / nrow(al.negpos)
  prop.uai.pospos <- sum(al.pospos[, "uai"] == 1) / nrow(al.pospos)

  prop.acts.negneg <- nrow(al.negneg) / (nrow(al))
  prop.acts.negpos <- nrow(al.negpos) / (nrow(al))
  prop.acts.pospos <- nrow(al.pospos) / (nrow(al))


  if (nrow(dal) == 0) {
    return(dat)
  }

  ## Reorder by role: ins on the left, rec on the right, flippers represented twice
  disc.ip <- dal[dal[, "ins"] %in% 1:2, ]
  disc.rp <- dal[dal[, "ins"] %in% c(0, 2), c(2:1, 3:ncols)]
  colnames(disc.ip)[1:2] <- colnames(disc.rp)[1:2] <- c("ins", "rec")


  # PATP: Insertive Man Infected (Col 1) --------------------------------

  # Attributes of infected
  ip.vl <- vl[disc.ip[, 1]]
  ip.stage <- stage[disc.ip[, 1]]
  ip.stage.syph.infector <- stage.syph[disc.ip[, 1]]
  ip.uGC.infector <- uGC[disc.ip[, 1]]
  ip.uCT.infector <- uCT[disc.ip[, 1]]

  # Attributes of susceptible
  ip.ccr5 <- ccr5[disc.ip[, 2]]
  ip.prep <- prepStat[disc.ip[, 2]]
  ip.prepcl <- prepClass[disc.ip[, 2]]
  ip.rGC <- rGC[disc.ip[, 2]]
  ip.rCT <- rCT[disc.ip[, 2]]
  ip.stage.syph.infectee <- stage.syph[disc.ip[, 2]]

  # Base TP from VL
  ip.tprob <- URAI.prob * 2.45^(ip.vl - 4.5)

  # Transform to log odds
  ip.tlo <- log(ip.tprob/(1 - ip.tprob))

  # Condom use
  not.UAI <- which(disc.ip[, "uai"] == 0)
  ip.tlo[not.UAI] <- ip.tlo[not.UAI] + log(condom.rr)

  # CCR5
  ip.tlo[ip.ccr5 == "DD"] <- ip.tlo[ip.ccr5 == "DD"] + -Inf
  ip.tlo[ip.ccr5 == "DW"] <- ip.tlo[ip.ccr5 == "DW"] + log(ccr5.heteroz.rr)

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(ip.prep == 1 & ip.prepcl == i - 1)
    ip.tlo[temp.ids] <- ip.tlo[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% 1:2)
  ip.tlo[isAcute] <- ip.tlo[isAcute] + log(acute.rr)

  ## Multiplier for HIV acquisition due to rectal STI in HIV-negative partner
  is.rGC <- which(ip.rGC == 1)
  is.rCT <- which(ip.rCT == 1)
  is.syph.infectee <- which(ip.stage.syph.infectee %in%  c(1, 2, 3))

  ### Single infections
  # NG
  is.rGC.sing <- setdiff(is.rGC, is.rCT)
  is.rGC.sing <- setdiff(is.rGC.sing, is.syph.infectee)

  # CT
  is.rCT.sing <- setdiff(is.rCT, is.rGC)
  is.rCT.sing <- setdiff(is.rCT.sing, is.syph.infectee)

  # Syph
  is.syph.sing <- setdiff(is.syph.infectee, is.rGC)
  is.syph.sing <- setdiff(is.syph.sing, is.rCT)

  ### Coinfections
  # NG and CT
  is.rGC.rCT <- intersect(is.rGC, is.rCT)
  is.rGC.rCT <- setdiff(is.rGC.rCT, is.syph.infectee)

  # NG and Syph
  is.rGC.syph <- intersect(is.rGC, is.syph.infectee)
  is.rGC.syph <- setdiff(is.rGC.syph, is.rCT)

  # CT and Syph
  is.rCT.syph <- intersect(is.rCT, is.syph.infectee)
  is.rCT.syph <- setdiff(is.rCT.syph, is.rGC)

  # All three infections
  is.all <- intersect(is.rGC.rCT, is.syph.infectee)

  ## Add relative risks
  # Single infections
  ip.tlo[is.rGC.sing] <- ip.tlo[is.rGC.sing] + log(hiv.rgc.rr)
  ip.tlo[is.rCT.sing] <- ip.tlo[is.rCT.sing] + log(hiv.rct.rr)
  ip.tlo[is.syph.sing] <- ip.tlo[is.syph.sing] + log(hiv.syph.rr)

  # Two infections
  ip.tlo[is.rGC.rCT] <- ip.tlo[is.rGC.rCT] +
    max(log(hiv.rgc.rr), log(hiv.rct.rr)) +
    min(log(hiv.rgc.rr), log(hiv.rct.rr)) * hiv.rgc.rct.rr

  ip.tlo[is.rGC.syph] <- ip.tlo[is.rGC.syph] +
    max(log(hiv.rgc.rr), log(hiv.syph.rr)) +
    min(log(hiv.rgc.rr), log(hiv.syph.rr)) * hiv.rgc.syph.rr

  ip.tlo[is.rCT.syph] <- ip.tlo[is.rCT.syph] +
    max(log(hiv.rct.rr), log(hiv.syph.rr)) +
    min(log(hiv.rct.rr), log(hiv.syph.rr)) * hiv.rct.syph.rr

  # Three infections
  ip.tlo[is.all] <- ip.tlo[is.all] +
    max(log(hiv.rct.rr), log(hiv.rgc.rr), log(hiv.syph.rr)) +
    min(log(hiv.rct.rr), log(hiv.rgc.rr), log(hiv.syph.rr)) * hiv.all.rect.rr

  ## Multiplier for HIV transmission due to urethral STI in HIV-positive partner
  is.syph.infector <- which(ip.stage.syph.infector %in% c(1, 2, 3))
  is.uGC.infector <- which(ip.uGC.infector == 1)
  is.uCT.infector <- which(ip.uCT.infector == 1)

  ### Single infections
  # NG
  is.uGC.sing <- setdiff(is.uGC.infector, is.uCT.infector)
  is.uGC.sing <- setdiff(is.uGC.sing, is.syph.infector)

  # CT
  is.uCT.sing <- setdiff(is.uCT.infector, is.uGC.infector)
  is.uCT.sing <- setdiff(is.uCT.sing, is.syph.infector)

  # Syph
  is.syph.sing <- setdiff(is.syph.infector, is.uGC.infector)
  is.syph.sing <- setdiff(is.syph.sing, is.uCT.infector)

  ### Coinfections
  # NG and CT
  is.uGC.uCT <- intersect(is.uGC.infector, is.uCT.infector)
  is.uGC.uCT <- setdiff(is.uGC.uCT, is.syph.infector)

  # NG and Syph
  is.uGC.syph <- intersect(is.uGC.infector, is.syph.infector)
  is.uGC.syph <- setdiff(is.uGC.syph, is.uCT.infector)

  # CT and Syph
  is.uCT.syph <- intersect(is.uCT.infector, is.syph.infector)
  is.uCT.syph <- setdiff(is.uCT.syph, is.uGC.infector)

  # All three infections
  is.all <- intersect(is.uGC.uCT, is.syph.infector)

  ## Add relative risks
  # Single infections
  ip.tlo[is.uGC.sing] <- ip.tlo[is.uGC.sing] + log(hiv.trans.gc.rr)
  ip.tlo[is.uCT.sing] <- ip.tlo[is.uCT.sing] + log(hiv.trans.ct.rr)
  ip.tlo[is.syph.sing] <- ip.tlo[is.syph.sing] + log(hiv.trans.syph.rr)

  # Two infections
  ip.tlo[is.uGC.uCT] <- ip.tlo[is.uGC.uCT] +
    max(log(hiv.trans.gc.rr), log(hiv.trans.ct.rr)) +
    min(log(hiv.trans.gc.rr), log(hiv.trans.ct.rr)) * hiv.trans.gc.ct.rr

  ip.tlo[is.uGC.syph] <- ip.tlo[is.uGC.syph] +
    max(log(hiv.trans.gc.rr), log(hiv.trans.syph.rr)) +
    min(log(hiv.trans.gc.rr), log(hiv.trans.syph.rr)) * hiv.trans.gc.syph.rr

  ip.tlo[is.uCT.syph] <- ip.tlo[is.uCT.syph] +
    max(log(hiv.trans.ct.rr), log(hiv.trans.syph.rr)) +
    min(log(hiv.trans.ct.rr), log(hiv.trans.syph.rr)) * hiv.trans.ct.syph.rr

  # Three infections
  ip.tlo[is.all] <- ip.tlo[is.all] +
    max(log(hiv.trans.ct.rr), log(hiv.trans.gc.rr), log(hiv.trans.syph.rr)) +
    min(log(hiv.trans.ct.rr), log(hiv.trans.gc.rr), log(hiv.trans.syph.rr)) * hiv.trans.allsti.rr


  ## Re-transform back to probability
  ip.tprob <- plogis(ip.tlo)
  stopifnot(ip.tprob >= 0, ip.tprob <= 1)


  # PATP: Receptive Man Infected (Col 2) --------------------------------

  # Attributes of infected
  rp.vl <- vl[disc.rp[, 2]]
  rp.stage <- stage[disc.rp[, 2]]
  rp.stage.syph.infector <- stage.syph[disc.rp[, 2]]
  rp.rGC.infector <- rGC[disc.rp[, 2]]
  rp.rCT.infector <- rCT[disc.rp[, 2]]

  # Attributes of susceptible
  rp.circ <- circ[disc.rp[, 1]]
  rp.ccr5 <- ccr5[disc.rp[, 1]]
  rp.prep <- prepStat[disc.rp[, 1]]
  rp.prepcl <- prepClass[disc.rp[, 1]]
  rp.uGC <- uGC[disc.rp[, 1]]
  rp.uCT <- uCT[disc.rp[, 1]]
  rp.stage.syph.infectee <- stage.syph[disc.rp[, 1]]

  # Base TP from VL
  rp.tprob <- UIAI.prob * 2.45^(rp.vl - 4.5)

  # Transform to log odds
  rp.tlo <- log(rp.tprob/(1 - rp.tprob))

  # Circumcision
  rp.tlo[rp.circ == 1] <- rp.tlo[rp.circ == 1] + log(circ.rr)

  # Condom use
  not.UAI <- which(disc.rp[, "uai"] == 0)
  rp.tlo[not.UAI] <- rp.tlo[not.UAI] + log(condom.rr)

  # CCR5
  rp.tlo[rp.ccr5 == "DD"] <- rp.tlo[rp.ccr5 == "DD"] + -Inf
  rp.tlo[rp.ccr5 == "DW"] <- rp.tlo[rp.ccr5 == "DW"] + log(ccr5.heteroz.rr)

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(rp.prep == 1 & rp.prepcl == i - 1)
    rp.tlo[temp.ids] <- rp.tlo[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% 1:2)
  rp.tlo[isAcute] <- rp.tlo[isAcute] + log(acute.rr)

  ## Multiplier for HIV acquisition due to urethral STI in HIV-negative partner
  is.uGC <- which(rp.uGC == 1)
  is.uCT <- which(rp.uCT == 1)
  is.syph.infectee <- which(rp.stage.syph.infectee %in% c(1, 2, 3))

  ### Single infections
  # NG
  is.uGC.sing <- setdiff(is.uGC, is.uCT)
  is.uGC.sing <- setdiff(is.uGC.sing, is.syph.infectee)

  # CT
  is.uCT.sing <- setdiff(is.uCT, is.uGC)
  is.uCT.sing <- setdiff(is.uCT.sing, is.syph.infectee)

  # Syph
  is.syph.sing <- setdiff(is.syph.infectee, is.uGC)
  is.syph.sing <- setdiff(is.syph.sing, is.uCT)

  ### Coinfections
  # NG and CT
  is.uGC.uCT <- intersect(is.uGC, is.uCT)
  is.uGC.uCT <- setdiff(is.uGC.uCT, is.syph.infectee)

  # NG and Syph
  is.uGC.syph <- intersect(is.uGC, is.syph.infectee)
  is.uGC.syph <- setdiff(is.uGC.syph, is.uCT)

  # CT and Syph
  is.uCT.syph <- intersect(is.uCT, is.syph.infectee)
  is.uCT.syph <- setdiff(is.uCT.syph, is.uGC)

  # All three infections
  is.all <- intersect(is.uGC.uCT, is.syph.infectee)

  # Single infections
  rp.tlo[is.uGC.sing] <- rp.tlo[is.uGC.sing] + log(hiv.ugc.rr)
  rp.tlo[is.uCT.sing] <- rp.tlo[is.uCT.sing] + log(hiv.uct.rr)
  rp.tlo[is.syph.sing] <- rp.tlo[is.syph.sing] + log(hiv.syph.rr)

  # Two infections
  rp.tlo[is.uGC.uCT] <- rp.tlo[is.uGC.uCT] +
    max(log(hiv.ugc.rr), log(hiv.uct.rr)) +
    min(log(hiv.ugc.rr), log(hiv.uct.rr)) * hiv.ugc.uct.rr

  rp.tlo[is.uGC.syph] <- rp.tlo[is.uGC.syph] +
    max(log(hiv.ugc.rr), log(hiv.syph.rr)) +
    min(log(hiv.ugc.rr), log(hiv.syph.rr)) * hiv.ugc.syph.rr

  rp.tlo[is.uCT.syph] <- rp.tlo[is.uCT.syph] +
    max(log(hiv.uct.rr), log(hiv.syph.rr)) +
    min(log(hiv.uct.rr), log(hiv.syph.rr)) * hiv.uct.syph.rr

  # Three infections
  rp.tlo[is.all] <- rp.tlo[is.all] +
     max(log(hiv.uct.rr), log(hiv.ugc.rr), log(hiv.syph.rr)) +
     min(log(hiv.uct.rr), log(hiv.ugc.rr), log(hiv.syph.rr)) * hiv.all.ureth.rr

  ## Multiplier for HIV transmission due to rectal STI in HIV-positive partner
  is.syph.infector <- which(rp.stage.syph.infector %in% c(1, 2, 3))
  is.rGC.infector <- which(rp.rGC.infector == 1)
  is.rCT.infector <- which(rp.rCT.infector == 1)

  ### Single infections
  # NG
  is.rGC.sing <- setdiff(is.rGC.infector, is.rCT.infector)
  is.rGC.sing <- setdiff(is.rGC.sing, is.syph.infector)

  # CT
  is.rCT.sing <- setdiff(is.rCT.infector, is.rGC.infector)
  is.rCT.sing <- setdiff(is.rGC.sing, is.syph.infector)

  # Syph
  is.syph.sing <- setdiff(is.syph.infector, is.rGC.infector)
  is.syph.sing <- setdiff(is.syph.sing, is.rCT.infector)

  ### Coinfections
  # NG and CT
  is.rGC.rCT <- intersect(is.rGC.infector, is.rCT.infector)
  is.rGC.rCT <- setdiff(is.rGC.rCT, is.syph.infector)

  # NG and Syph
  is.rGC.syph <- intersect(is.rGC.infector, is.syph.infector)
  is.rGC.syph <- setdiff(is.rGC.syph, is.rCT.infector)

  # CT and Syph
  is.rCT.syph <- intersect(is.rCT.infector, is.syph.infector)
  is.rCT.syph <- setdiff(is.rCT.syph, is.rGC.infector)

  # All three infections
  is.all <- intersect(is.rGC.rCT, is.syph.infector)

  # Single infections
  rp.tlo[is.rGC.sing] <- rp.tlo[is.rGC.sing] + log(hiv.trans.gc.rr)
  rp.tlo[is.rCT.sing] <- rp.tlo[is.rCT.sing] + log(hiv.trans.ct.rr)
  rp.tlo[is.syph.sing] <- rp.tlo[is.syph.sing] + log(hiv.trans.syph.rr)

  # Two infections
  rp.tlo[is.rGC.rCT] <- rp.tlo[is.rGC.rCT] +
    max(log(hiv.trans.gc.rr), log(hiv.trans.ct.rr)) +
    min(log(hiv.trans.gc.rr), log(hiv.trans.ct.rr)) * hiv.trans.gc.ct.rr

  rp.tlo[is.rGC.syph] <- rp.tlo[is.rGC.syph] +
    max(log(hiv.trans.gc.rr), log(hiv.trans.syph.rr)) +
    min(log(hiv.trans.gc.rr), log(hiv.trans.syph.rr)) * hiv.trans.gc.syph.rr

  rp.tlo[is.rCT.syph] <- rp.tlo[is.rCT.syph] +
    max(log(hiv.trans.ct.rr), log(hiv.trans.syph.rr)) +
    min(log(hiv.trans.ct.rr), log(hiv.trans.syph.rr)) * hiv.trans.ct.syph.rr

  # Three infections
  rp.tlo[is.all] <- rp.tlo[is.all] +
    max(log(hiv.trans.ct.rr), log(hiv.trans.gc.rr), log(hiv.trans.syph.rr)) +
    min(log(hiv.trans.ct.rr), log(hiv.trans.gc.rr), log(hiv.trans.syph.rr)) * hiv.trans.allsti.rr

  ## Retransformation to probability
  rp.tprob <- plogis(rp.tlo)
  stopifnot(rp.tprob >= 0, rp.tprob <= 1)


  # Transmission --------------------------------------------------------

  ## Bernoulli transmission events
  trans.ip <- rbinom(length(ip.tprob), 1, ip.tprob)
  trans.rp <- rbinom(length(rp.tprob), 1, rp.tprob)


  # Output --------------------------------------------------------------


  # Update attributes

  infected <- inf.type <- NULL
  if (sum(trans.ip, trans.rp) > 0) {

    infected <- unique(c(disc.ip[trans.ip == 1, 2],
                  disc.rp[trans.rp == 1, 1]))
    inf.role <- c(rep(0, sum(trans.ip)), rep(1, sum(trans.rp)))

    inf.type <- c(disc.ip[trans.ip == 1, "ptype"],
                  disc.rp[trans.rp == 1, "ptype"])

    dat$attr$status[infected] <- 1
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- 1
    dat$attr$stage.time[infected] <- 0
    dat$attr$stage.time.ar.ndx[infected] <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected] <- 0

    dat$attr$inf.role[infected] <- inf.role[infected]
    dat$attr$inf.type[infected] <- inf.type[infected]

    dat$attr$cum.time.on.tx[infected] <- 0
    dat$attr$cum.time.off.tx[infected] <- 0
  }
  dat$attr$time.hivneg[status == 0] <- dat$attr$time.hivneg[status == 0] + 1

  trans <- rbind(disc.ip[trans.ip == 1, ], disc.rp[trans.rp == 1, ])
  dat$epi$sum_GC[at] <- length(which(((rGC[trans[, 2]] == 1 | uGC[trans[, 1]] == 1) & trans[, 6] == 1) |
                                       ((uGC[trans[, 2]] == 1 | rGC[trans[, 1]] == 1) & trans[, 6] == 0)))

  dat$epi$sum_CT[at] <- length(which(((rCT[trans[, 2]] == 1 | uCT[trans[, 1]] == 1) & trans[, 6] == 1) |
                                       ((uCT[trans[, 2]] == 1 | rCT[trans[, 1]] == 1) & trans[, 6] == 0)))

  dat$epi$sum_syph[at] <- length(which(stage.syph[trans[, 2]] %in% c(1,2,3) | stage.syph[trans[, 1]] %in% c(1,2,3)))

  dat$epi$sum_urethral[at] <- length(which(((uGC[trans[, 1]] == 1 | uCT[trans[, 1]] == 1) & trans[, 6] == 1) |
                                             ((uGC[trans[, 2]] == 1 | uCT[trans[, 2]] == 1) & trans[, 6] == 0)))

  dat$epi$sum_rectal[at] <- length(which(((rGC[trans[, 2]] == 1 | rCT[trans[, 2]] == 1) & trans[, 6] == 1) |
                                           ((rGC[trans[, 1]] == 1 | rCT[trans[, 1]] == 1) & trans[, 6] == 0)))
  #2x2 for PAF
  #               HIV+
  #             STI+  STI-
  #HIV-   STI +  1    2
  #       STI -  3    4
  dat$epi$cell1_gc[at] <- length(which(
    # P1 is infected, p1 has urethral and p2 has rectal OR
    (status[trans[, 1]] == 1 & rGC[trans[, 2]] == 1 & uGC[trans[, 1]] == 1) |
      # P2 is infected, p1 has urethral and p2 has rectal
    (status[trans[, 2]] == 1 & rGC[trans[, 2]] == 1 & uGC[trans[, 1]] == 1)))

  dat$epi$cell2_gc[at] <- length(which(
    # P1 is infected, p1 does not have urethral GC, p2 has rectal GC
    (status[trans[, 1]] == 1 & uGC[trans[, 1]] == 0 & rGC[trans[, 2]] == 1) |
    # P2 is infected, p1 has urethral GC, p2 does not have rectal GC
    (status[trans[, 2]] == 1 & uGC[trans[, 1]] == 1 & rGC[trans[, 2]] == 0)))

  dat$epi$cell3_gc[at] <- length(which(
    # P1 is infected, p1 has urethral GC, p2 does not have rectal GC
    (status[trans[, 1]] == 1 & uGC[trans[, 1]] == 1 & rGC[trans[, 2]] == 0) |
    # P2 is infected, p1 does not have urethral GC, p2 does have rectal GC
    (status[trans[, 2]] == 1 & uGC[trans[, 1]] == 0 & rGC[trans[, 2]] == 1)))

  dat$epi$cell4_gc[at] <- length(which(
    # P1 is infected, p1 does not have urethral GC, p2 does not have rectal GC
    (status[trans[, 1]] == 1 & uGC[trans[, 1]] == 0 & rGC[trans[, 2]] == 0) |
    # P2 is infected, p1 does not have urethral GC, p2 does not have rectal GC
    (status[trans[, 2]] == 1 & uGC[trans[, 1]] == 0 & rGC[trans[, 2]] == 0)))

  dat$epi$cell1_ct[at] <- length(which(
    # P1 is infected, p1 has urethral and p2 has rectal OR
    (status[trans[, 1]] == 1 & rCT[trans[, 2]] == 1 & uCT[trans[, 1]] == 1) |
      # P2 is infected, p1 has urethral and p2 has rectal
      (status[trans[, 2]] == 1 & rCT[trans[, 2]] == 1 & uCT[trans[, 1]] == 1)))

  dat$epi$cell2_ct[at] <- length(which(
    # P1 is infected, p1 does not have urethral CT, p2 has rectal CT
    (status[trans[, 1]] == 1 & uCT[trans[, 1]] == 0 & rCT[trans[, 2]] == 1) |
    # P2 is infected, p1 has urethral CT, p2 does not have rectal CT
    (status[trans[, 2]] == 1 & uCT[trans[, 1]] == 1 & rCT[trans[, 2]] == 0)))

  dat$epi$cell3_ct[at] <- length(which(
    # P1 is infected, p1 has urethral CT, p2 does not have rectal CT
    (status[trans[, 1]] == 1 & uCT[trans[, 1]] == 1 & rCT[trans[, 2]] == 0) |
    # P2 is infected, p1 does not have urethral CT, p2 does have rectal CT
    (status[trans[, 2]] == 1 & uCT[trans[, 1]] == 0 & rCT[trans[, 2]] == 1)))

  dat$epi$cell4_ct[at] <- length(which(
    # P1 is infected, p1 does not have urethral CT, p2 does not have rectal CT
    (status[trans[, 1]] == 1 & uCT[trans[, 1]] == 0 & rCT[trans[, 2]] == 0) |
    # P2 is infected, p1 does not have urethral CT, p2 does not have rectal CT
    (status[trans[, 2]] == 1 & uCT[trans[, 1]] == 0 & rCT[trans[, 2]] == 0)))

  dat$epi$cell1_syph[at] <- length(which(
    #P1 is infected, P1 has syphilis, P2 has syphilis
    status[trans[, 1]] == 1 & stage.syph[trans[, 2]] %in% c(1,2,3) & stage.syph[trans[, 1]] %in% c(1,2,3) |
    #P2 is infected, P1 has syphilis, P2 has syphilis
    status[trans[, 2]] == 1 & stage.syph[trans[, 2]] %in% c(1,2,3) & stage.syph[trans[, 1]] %in% c(1,2,3)))

  dat$epi$cell2_syph[at] <- length(which(
    #P1 is infected, P1 does not have syphilis, P2 has syphilis
    status[trans[, 1]] == 1 & stage.syph[trans[, 2]] %in% c(1,2,3) & !(stage.syph[trans[, 1]] %in% c(1,2,3)) |
    #P2 is infected, P1 has syphilis, P2 does not have syphilis
    status[trans[, 2]] == 1 & !(stage.syph[trans[, 2]] %in% c(1,2,3)) & stage.syph[trans[, 1]] %in% c(1,2,3)))

  dat$epi$cell3_syph[at] <- length(which(
    #P1 is infected, P1 has syphilis, P2 does not have syphilis OR
    status[trans[, 1]] == 1 & stage.syph[trans[, 1]] %in% c(1,2,3) & !(stage.syph[trans[, 2]] %in% c(1,2,3)) |
    #P2 is infected, P1 does not have syphilis, P2 has syphilis
    status[trans[, 2]] == 1 & !(stage.syph[trans[, 1]] %in% c(1,2,3)) & stage.syph[trans[, 2]] %in% c(1,2,3)))

  dat$epi$cell4_syph[at] <- length(which(
    #P1 is infected, P1 does not have syphilis, P2 does not have syphilis OR
    status[trans[, 1]] == 1 & !(stage.syph[trans[, 1]] %in% c(1,2,3)) & !(stage.syph[trans[, 2]] %in% c(1,2,3)) |
    #P2 is infected, P1 does not have syphilis, P2 does not have syphilis
    status[trans[, 2]] == 1 & !(stage.syph[trans[, 1]] %in% c(1,2,3)) & !(stage.syph[trans[, 2]] %in% c(1,2,3))))

  dat$epi$cell1_sti[at] <- length(which(
    #P1 is infected, P1 has urethral STI, P2 has rectal STI OR
    (status[trans[, 1]] == 1 & ((uGC[trans[, 1]] == 1 | uCT[trans[ , 1]] == 1 | stage.syph[trans[, 1]] %in% c(1,2,3))) &
                                  (rGC[trans[, 2]] == 1 | rCT[trans[ , 2]] == 1 | stage.syph[trans[, 2]] %in% c(1,2,3))) |
    #P2 is infected, P1 has urethral STI, P2 has rectal STI OR
    (status[trans[, 2]] == 1 & ((uGC[trans[, 1]] == 1 | uCT[trans[ , 1]] == 1 | stage.syph[trans[, 1]] %in% c(1,2,3))) &
                                  (rGC[trans[, 2]] == 1 | rCT[trans[ , 2]] == 1 | stage.syph[trans[, 2]] %in% c(1,2,3)))))

  dat$epi$cell2_sti[at] <- length(which(
    #P1 is infected, P1 does not have urethral STI, P2 has rectal STI OR
    (status[trans[, 1]] == 1 & ((uGC[trans[, 1]] == 0 & uCT[trans[, 1]] == 0 & !(stage.syph[trans[, 1]] %in% c(1,2,3)))) &
                                  (rGC[trans[, 2]] == 1 | rCT[trans[, 2]] == 1 | stage.syph[trans[, 2]] %in% c(1,2,3))) |
    #P2 is infected, P1 has urethral STI, P2 does not have rectal STI OR
    (status[trans[, 2]] == 1 & ((uGC[trans[, 1]] == 1 | uCT[trans[, 1]] == 1 | stage.syph[trans[, 1]] %in% c(1,2,3))) &
                                    (rGC[trans[, 2]] == 0 & rCT[trans[, 2]] == 0 & !(stage.syph[trans[, 2]] %in% c(1,2,3))))))

  dat$epi$cell3_sti[at] <- length(which(
    #P1 is infected, P1 has urethral STI, P2 does not have rectal STI OR
    (status[trans[, 1]] == 1 & ((uGC[trans[, 1]] == 1 | uCT[trans[, 1]] == 1 | stage.syph[trans[, 1]] %in% c(1,2,3))) &
                                  (rGC[trans[, 2]] == 0 & rCT[trans[, 2]] == 0 & !(stage.syph[trans[, 2]] %in% c(1,2,3)))) |
    #P2 is infected, P1 does not have urethral STI, P2 has rectal STI OR
    (status[trans[, 2]] == 1 & ((uGC[trans[, 1]] == 0 & uCT[trans[, 1]] == 0 & !(stage.syph[trans[, 1]] %in% c(1,2,3)))) &
                                    (rGC[trans[, 2]] == 1 | rCT[trans[, 2]] == 1 | stage.syph[trans[, 2]] %in% c(1,2,3)))))

  dat$epi$cell4_sti[at] <- length(which(
    #P1 is infected, P1 does not have urethral STI, P2 does not have rectal STI OR
    (status[trans[, 1]] == 1 & ((uGC[trans[, 1]] == 0 & uCT[trans[, 1]] == 0 & !(stage.syph[trans[, 1]] %in% c(1,2,3)))) &
                                  (rGC[trans[, 2]] == 0 & rCT[trans[, 2]] == 0 & !(stage.syph[trans[, 2]] %in% c(1,2,3)))) |
    #P2 is infected, P1 does not have urethral STI, P2 does not have rectal STI OR
    (status[trans[, 2]] == 1 & ((uGC[trans[, 1]] == 0 & uCT[trans[, 1]] == 0 & !(stage.syph[trans[, 1]] %in% c(1,2,3)))) &
                                  (rGC[trans[, 2]] == 0 & rCT[trans[, 2]] == 0 & !(stage.syph[trans[, 2]] %in% c(1,2,3))))))

  # Summary Output
  dat$epi$incid[at] <- length(infected)

  dat$epi$trans.main[at] <- sum(inf.type == 1)
  dat$epi$trans.pers[at] <- sum(inf.type == 2)
  dat$epi$trans.inst[at] <- sum(inf.type == 3)

  dat$epi$num.acts.negneg[at] <- num.acts.negneg
  dat$epi$num.acts.negpos[at] <- num.acts.negpos
  dat$epi$num.acts.pospos[at] <- num.acts.pospos

  dat$epi$prop.uai.negneg[at] <- prop.uai.negneg
  dat$epi$prop.uai.negpos[at] <- prop.uai.negpos
  dat$epi$prop.uai.pospos[at] <- prop.uai.pospos

  dat$epi$prop.acts.negneg[at] <- prop.acts.negneg
  dat$epi$prop.acts.negpos[at] <- prop.acts.negpos
  dat$epi$prop.acts.pospos[at] <- prop.acts.pospos

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
