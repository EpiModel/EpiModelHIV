
#' @title STI Transmission Module
#'
#' @description Stochastically simulates GC/CT transmission given the current
#'              state of the edgelist.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
sti_trans_msm <- function(dat, at) {

  # Parameters ----------------------------------------------------------

  # Acquisition probability | exposure
  rgc.tprob <- dat$param$rgc.tprob
  ugc.tprob <- dat$param$ugc.tprob
  rct.tprob <- dat$param$rct.tprob
  uct.tprob <- dat$param$uct.tprob
  syph.tprob <- dat$param$syph.tprob

  # Relative risk by syphilis stage
  syph.incub.rr <- dat$param$syph.incub.rr
  syph.earlat.rr <- dat$param$syph.earlat.rr
  syph.late.rr <- dat$param$syph.late.rr

  # Probability of symptoms | infection
  rgc.sympt.prob <- dat$param$rgc.sympt.prob
  ugc.sympt.prob <- dat$param$ugc.sympt.prob
  rct.sympt.prob <- dat$param$rct.sympt.prob
  uct.sympt.prob <- dat$param$uct.sympt.prob
  syph.incub.sympt.prob <- dat$param$syph.incub.sympt.prob

  # Relative risk of infection | condom use
  sti.cond.rr <- dat$param$sti.cond.rr


  # Attributes ----------------------------------------------------------

  # Current infection state
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT
  syphilis <- dat$attr$syphilis
  stage.syph <- dat$attr$stage.syph
  syph.sympt <- dat$attr$syph.sympt
  stage.time.syph <- dat$attr$stage.time.syph

  # Infection time
  rGC.infTime <- dat$attr$rGC.infTime
  uGC.infTime <- dat$attr$uGC.infTime
  rCT.infTime <- dat$attr$rCT.infTime
  uCT.infTime <- dat$attr$uCT.infTime
  syph.infTime <- dat$attr$syph.infTime
  last.rGC.infTime <- dat$attr$last.rGC.infTime
  last.uGC.infTime <- dat$attr$last.uGC.infTime
  last.rCT.infTime <- dat$attr$last.rCT.infTime
  last.uCT.infTime <- dat$attr$last.uCT.infTime
  last.syph.infTime <- dat$attr$last.syph.infTime

  # GC/CT Infection symptoms (non-varying)
  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  rCT.sympt <- dat$attr$rCT.sympt
  uCT.sympt <- dat$attr$uCT.sympt

  # Diagnosis status
  diag.status.gc <- dat$attr$diag.status.gc
  diag.status.ct <- dat$attr$diag.status.ct
  diag.status.syph <- dat$attr$diag.status.syph

  # Pull act list
  al <- dat$temp$al

  ## ins variable coding
  # ins = 0 : p2 is insertive
  # ins = 1 : p1 is insertive
  # ins = 2 : both p1 and p2 are insertive

  # Rectal GC -----------------------------------------------------------

  # Requires: uGC in insertive man, and no rGC in receptive man
  p1Inf_rgc <- which(uGC[al[, "p1"]] == 1 & uGC.infTime[al[, "p1"]] < at &
                     rGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(1, 2))
  p2Inf_rgc <- which(uGC[al[, "p2"]] == 1 & uGC.infTime[al[, "p2"]] < at &
                     rGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(0, 2))
  allActs_rgc <- c(p1Inf_rgc, p2Inf_rgc)

  # UAI modifier
  uai_rgc <- al[, "uai"][allActs_rgc]
  tprob_rgc <- rep(rgc.tprob, length(allActs_rgc))
  tprob_rgc[uai_rgc == 0] <- tprob_rgc[uai_rgc == 0] * sti.cond.rr

  # Stochastic transmission
  trans_rgc <- rbinom(length(allActs_rgc), 1, tprob_rgc)

  # Determine the infected partner
  idsInf_rgc <- NULL
  if (sum(trans_rgc) > 0) {
    transAL_rgc <- al[allActs_rgc[trans_rgc == 1], , drop = FALSE]
    idsInf_rgc <- unique(ifelse(uGC[transAL_rgc[, "p1"]] == 1,
                                transAL_rgc[, "p2"], transAL_rgc[, "p1"]))
  }

  # Update attributes
  rGC[idsInf_rgc] <- 1
  rGC.infTime[idsInf_rgc] <- last.rGC.infTime[idsInf_rgc] <- at
  rGC.sympt[idsInf_rgc] <- rbinom(length(idsInf_rgc), 1, rgc.sympt.prob)
  diag.status.gc[idsInf_rgc] <- 0


  # Urethral GC ---------------------------------------------------------

  # Requires: rGC in receptive man, and no uGC in insertive man
  p1Inf_ugc <- which(rGC[al[, "p1"]] == 1 & rGC.infTime[al[, "p1"]] < at &
                     uGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf_ugc <- which(rGC[al[, "p2"]] == 1 & rGC.infTime[al[, "p2"]] < at &
                     uGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs_ugc <- c(p1Inf_ugc, p2Inf_ugc)

  # UAI modifier
  uai_ugc <- al[, "uai"][allActs_ugc]
  tprob_ugc <- rep(ugc.tprob, length(allActs_ugc))
  tprob_ugc[uai_ugc == 0] <- tprob_ugc[uai_ugc == 0] * sti.cond.rr

  # Stochastic transmission
  trans_ugc <- rbinom(length(allActs_ugc), 1, tprob_ugc)

  # Determine the newly infected partner
  idsInf_ugc <- NULL
  if (sum(trans_ugc) > 0) {
    transAL_ugc <- al[allActs_ugc[trans_ugc == 1],  , drop = FALSE]
    idsInf_ugc <- unique(ifelse(uGC[transAL_ugc[, "p1"]] == 1,
                                transAL_ugc[, "p2"], transAL_ugc[, "p1"]))
  }

  # Update attributes
  uGC[idsInf_ugc] <- 1
  uGC.infTime[idsInf_ugc] <- last.uGC.infTime[idsInf_ugc] <- at
  uGC.sympt[idsInf_ugc] <- rbinom(length(idsInf_ugc), 1, ugc.sympt.prob)
  diag.status.gc[idsInf_ugc] <- 0


  # Rectal CT -----------------------------------------------------------

  # Requires: uCT in insertive man, and no rCT in receptive man
  p1Inf_rct <- which(uCT[al[, "p1"]] == 1 & uCT.infTime[al[, "p1"]] < at &
                     rCT[al[, "p2"]] == 0 & al[, "ins"] %in% c(1, 2))
  p2Inf_rct <- which(uCT[al[, "p2"]] == 1 & uCT.infTime[al[, "p2"]] < at &
                     rCT[al[, "p1"]] == 0 & al[, "ins"] %in% c(0, 2))
  allActs_rct <- c(p1Inf_rct, p2Inf_rct)

  # UAI modifier
  uai_rct <- al[, "uai"][allActs_rct]
  tprob_rct <- rep(rct.tprob, length(allActs_rct))
  tprob_rct[uai_rct == 0] <- tprob_rct[uai_rct == 0] * sti.cond.rr

  # Stochastic transmission
  trans_rct <- rbinom(length(allActs_rct), 1, tprob_rct)

  # Determine the newly infected partner
  idsInf_rct <- NULL
  if (sum(trans_rct) > 0) {
    transAL_rct <- al[allActs_rct[trans_rct == 1],  , drop = FALSE]
    idsInf_rct <- unique(ifelse(uCT[transAL_rct[, "p1"]] == 1,
                                transAL_rct[, "p2"], transAL_rct[, "p1"]))
  }

  # Update attributes
  rCT[idsInf_rct] <- 1
  rCT.infTime[idsInf_rct] <- last.rCT.infTime[idsInf_rct] <- at
  rCT.sympt[idsInf_rct] <- rbinom(length(idsInf_rct), 1, rct.sympt.prob)
  diag.status.ct[idsInf_rct] <- 0


  # Urethral CT ---------------------------------------------------------

  # Requires: rCT in receptive man, and no uCT in insertive man
  p1Inf_uct <- which(rCT[al[, "p1"]] == 1 & rCT.infTime[al[, "p1"]] < at &
                     uCT[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf_uct <- which(rCT[al[, "p2"]] == 1 & rCT.infTime[al[, "p2"]] < at &
                     uCT[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs_uct <- c(p1Inf_uct, p2Inf_uct)

  # UAI modifier
  uai_uct <- al[, "uai"][allActs_uct]
  tprob_uct <- rep(uct.tprob, length(allActs_uct))
  tprob_uct[uai_uct == 0] <- tprob_uct[uai_uct == 0] * sti.cond.rr

  # Transmission
  trans_uct <- rbinom(length(allActs_uct), 1, tprob_uct)

  # Determine the newly infected partner
  idsInf_uct <- NULL
  if (sum(trans_uct) > 0) {
    transAL_uct <- al[allActs_uct[trans_uct == 1],  , drop = FALSE]
    idsInf_uct <- unique(ifelse(uCT[transAL_uct[, "p1"]] == 1,
                                transAL_uct[, "p2"], transAL_uct[, "p1"]))
  }

  # Update attributes
  uCT[idsInf_uct] <- 1
  uCT.infTime[idsInf_uct] <- last.uCT.infTime[idsInf_uct] <- at
  uCT.sympt[idsInf_uct] <- rbinom(length(idsInf_uct), 1, uct.sympt.prob)
  diag.status.ct[idsInf_uct] <- 0


  # Syphilis ---------------------------------------------------------

  # Find the syphilis discordant pairs from the act list
  p1Inf_syph <- al[which(syphilis[al[, "p1"]] == 1 &
                           syph.infTime[al[, "p1"]] < at &
                           syphilis[al[, "p2"]] == 0), ,
                   drop = FALSE]
  p2Inf_syph <- al[which(syphilis[al[, "p2"]] == 1 &
                           syph.infTime[al[, "p2"]] < at &
                           syphilis[al[, "p1"]] == 0), ,
                   drop = FALSE]

  # Invert ordering of p2Inf_syph matrix: this puts the infected person always
  # in column one, then changes the categorization of who is R/I because the "ins"
  # variable is expressed in terms of "p1" in the original ordering. This allows
  # using the HIV-related, position-specific transmission code as intended.
  idsP1Ins <- which(p2Inf_syph[, "ins"] == 1)
  idsP1Rec <- which(p2Inf_syph[, "ins"] == 0)
  p2Inf_syph[, 1:2] <- p2Inf_syph[, 2:1]
  p2Inf_syph[idsP1Ins, "ins"] <- 0
  p2Inf_syph[idsP1Rec, "ins"] <- 1
  stopifnot(!any(dat$attr$role.class[p2Inf_syph[, 1]] == "R" &
                   p2Inf_syph[, "ins"] != 0))

  allActs_syph <- rbind(p1Inf_syph, p2Inf_syph)
  ncols <- ncol(allActs_syph)

  # Reorder by role: ins on the left, rec on the right, flippers represented twice
  disc.syph.ip <- allActs_syph[allActs_syph[, "ins"] %in% 1:2, , drop = FALSE]
  disc.syph.rp <- allActs_syph[allActs_syph[, "ins"] %in% c(0, 2), c(2:1, 3:ncols), drop = FALSE]
  colnames(disc.syph.ip)[1:2] <- colnames(disc.syph.rp)[1:2] <- c("ins", "rec")

  ## Insertive Man Infected with Syphilis (Col 1)
  if (nrow(disc.syph.ip) == 0) {
    trans.syph.ip <- NULL
  } else {
    # Syphilis stage of infected partner
    ip.stage.syph <- stage.syph[disc.syph.ip[, 1]]
    stopifnot(all(!is.na(ip.stage.syph)))

    # Base TP from VL
    ip.syph.tprob <- rep(syph.tprob, length(ip.stage.syph))

    # Transform to log odds
    ip.syph.tlo <- log(ip.syph.tprob/(1 - ip.syph.tprob))

    # Condom use multiplier
    not.syph.ip.UAI <- which(disc.syph.ip[, "uai"] == 0)
    ip.syph.tlo[not.syph.ip.UAI] <- ip.syph.tlo[not.syph.ip.UAI] + log(sti.cond.rr)

    # Incubating stage multiplier
    isincub <- which(ip.stage.syph %in% 1)
    ip.syph.tlo[isincub] <- ip.syph.tlo[isincub] + log(syph.incub.rr)

    # Early latent-stage multiplier
    isearlat <- which(ip.stage.syph %in% 4)
    ip.syph.tlo[isearlat] <- ip.syph.tlo[isearlat] + log(syph.earlat.rr)

    # Retransformation to probability
    ip.syph.tprob <- plogis(ip.syph.tlo)

    # Late stage multiplier (not log odds b/c log 0 = undefined)
    islate <- which(ip.stage.syph %in% c(5, 6))
    ip.syph.tprob[islate] <- ip.syph.tprob[islate] * syph.late.rr

    # Check for valid probabilities
    stopifnot(ip.syph.tprob >= 0, ip.syph.tprob <= 1)

    ## Bernoulli Transmission Events
    trans.syph.ip <- rbinom(length(ip.syph.tprob), 1, ip.syph.tprob)
  }

  ## Receptive Man Infected with Syphilis (Col 2)
  if (nrow(disc.syph.rp) == 0) {
    trans.syph.rp <- NULL
  } else {
    # Syphilis stage of infected partner
    rp.stage.syph <- stage.syph[disc.syph.rp[, 2]]
    stopifnot(all(!is.na(rp.stage.syph)))

    # Base TP from VL
    rp.syph.tprob <- rep(syph.tprob, length(rp.stage.syph))

    # Transform to log odds
    rp.syph.tlo <- log(rp.syph.tprob/(1 - rp.syph.tprob))

    # Condom use multiplier
    not.syph.rp.UAI <- which(disc.syph.rp[, "uai"] == 0)
    rp.syph.tlo[not.syph.rp.UAI] <- rp.syph.tlo[not.syph.rp.UAI] + log(sti.cond.rr)

    # Incubating stage multiplier
    isincub <- which(rp.stage.syph %in% 1)
    rp.syph.tlo[isincub] <- rp.syph.tlo[isincub] + log(syph.incub.rr)

    # Early latent stage multipliers
    isearlat <- which(rp.stage.syph %in% 4)
    rp.syph.tlo[isearlat] <- rp.syph.tlo[isearlat] + log(syph.earlat.rr)

    # Retransformation to probability
    rp.syph.tprob <- plogis(rp.syph.tlo)

    # Late stage multiplier (not log odds b/c log 0 = undefined)
    islate <- which(rp.stage.syph %in% c(5, 6))
    rp.syph.tprob[islate] <- rp.syph.tprob[islate] * (syph.late.rr)

    # Check for valid probabilities
    stopifnot(rp.syph.tprob >= 0, rp.syph.tprob <= 1)

    # Bernoulli transmission events
    trans.syph.rp <- rbinom(length(rp.syph.tprob), 1, rp.syph.tprob)
  }

  # Update attributes for newly infected
  idsInf_syph <- NULL
  if (sum(trans.syph.ip, trans.syph.rp, na.rm = TRUE) > 0) {
    idsInf_syph <- unique(c(disc.syph.ip[trans.syph.ip == 1, 2],
                       disc.syph.rp[trans.syph.rp == 1, 1]))
    syphilis[idsInf_syph] <- 1
    syph.infTime[idsInf_syph] <- last.syph.infTime[idsInf_syph] <- at
    stage.syph[idsInf_syph] <- 1
    stage.time.syph[idsInf_syph] <- 0
    diag.status.syph[idsInf_syph] <- 0
    syph.sympt[idsInf_syph] <- rbinom(length(idsInf_syph), 1, syph.incub.sympt.prob)
  }

  # Output --------------------------------------------------------------

  # Gonorrhea attributes
  dat$attr$rGC <- rGC
  dat$attr$uGC <- uGC
  dat$attr$rGC.infTime <- rGC.infTime
  dat$attr$uGC.infTime <- uGC.infTime
  dat$attr$last.rGC.infTime <- last.rGC.infTime
  dat$attr$last.uGC.infTime <- last.uGC.infTime
  dat$attr$rGC.sympt <- rGC.sympt
  dat$attr$uGC.sympt <- uGC.sympt
  dat$attr$diag.status.gc <- diag.status.gc

  # Chlamydia attributes
  dat$attr$rCT <- rCT
  dat$attr$uCT <- uCT
  dat$attr$rCT.infTime <- rCT.infTime
  dat$attr$uCT.infTime <- uCT.infTime
  dat$attr$last.rCT.infTime <- last.rCT.infTime
  dat$attr$last.uCT.infTime <- last.uCT.infTime
  dat$attr$rCT.sympt <- rCT.sympt
  dat$attr$uCT.sympt <- uCT.sympt
  dat$attr$diag.status.ct <- diag.status.ct

  # Syphilis attributes
  dat$attr$syphilis <- syphilis
  dat$attr$syph.infTime <- syph.infTime
  dat$attr$last.syph.infTime <- last.syph.infTime
  dat$attr$stage.syph <- stage.syph
  dat$attr$syph.sympt <- syph.sympt
  dat$attr$stage.time.syph <- stage.time.syph
  dat$attr$diag.status.syph <- diag.status.syph

  # Summary incidence statistics
  dat$epi$incid.rgc[at] <- length(idsInf_rgc)
  dat$epi$incid.ugc[at] <- length(idsInf_ugc)
  dat$epi$incid.gc[at] <- length(unique(c(idsInf_rgc,idsInf_ugc)))
  dat$epi$incid.rct[at] <- length(idsInf_rct)
  dat$epi$incid.uct[at] <- length(idsInf_uct)
  dat$epi$incid.ct[at] <- length(unique(c(idsInf_rct,idsInf_uct)))
  dat$epi$incid.syph[at] <- length(idsInf_syph)
  dat$epi$incid.sti[at] <- dat$epi$incid.gc[at] + dat$epi$incid.ct[at] + dat$epi$incid.syph[at]
  dat$epi$incid.gcct.prep[at] <- length(intersect(unique(c(idsInf_rgc, idsInf_ugc,
                                                           idsInf_rct, idsInf_uct)),
                                                  which(dat$attr$prepStat == 1)))
  dat$epi$incid.syph.prep[at] <- length(intersect(idsInf_syph,
                                                  which(dat$attr$prepStat == 1)))

  dat$epi$incid.gcct[at] <- length(unique(c(idsInf_rgc,idsInf_ugc))) + length(unique(c(idsInf_rct,idsInf_uct)))


  # Risk group-specific
  dat$epi$incid.gc.tttraj1[at] <- length(which(dat$attr$tt.traj.gc.hivneg[unique(c(idsInf_rgc,idsInf_ugc))] == 1)) + length(which(dat$attr$tt.traj.gc.hivpos[unique(c(idsInf_rgc,idsInf_ugc))] == 1))
  dat$epi$incid.ct.tttraj1[at] <- length(which(dat$attr$tt.traj.ct.hivneg[unique(c(idsInf_rct,idsInf_uct))] == 1)) + length(which(dat$attr$tt.traj.ct.hivpos[unique(c(idsInf_rct,idsInf_uct))] == 1))
  dat$epi$incid.syph.tttraj1[at] <- length(which(dat$attr$tt.traj.syph.hivneg[unique(c(idsInf_syph))] == 1)) + length(which(dat$attr$tt.traj.syph.hivpos[unique(c(idsInf_syph))] == 1))
  dat$epi$incid.sti.tttraj1[at] <- dat$epi$incid.gc.tttraj1[at] + dat$epi$incid.ct.tttraj1[at] + dat$epi$incid.syph.tttraj1[at]
  dat$epi$incid.gcct.tttraj1[at] <- dat$epi$incid.gc.tttraj1[at] + dat$epi$incid.ct.tttraj1[at]


  dat$epi$incid.gc.tttraj2[at] <- length(which(dat$attr$tt.traj.gc.hivneg[unique(c(idsInf_rgc,idsInf_ugc))] == 2)) + length(which(dat$attr$tt.traj.gc.hivpos[unique(c(idsInf_rgc,idsInf_ugc))] == 2))
  dat$epi$incid.ct.tttraj2[at] <- length(which(dat$attr$tt.traj.ct.hivneg[unique(c(idsInf_rct,idsInf_uct))] == 2)) + length(which(dat$attr$tt.traj.ct.hivpos[unique(c(idsInf_rct,idsInf_uct))] == 2))
  dat$epi$incid.syph.tttraj2[at] <- length(which(dat$attr$tt.traj.syph.hivneg[unique(c(idsInf_syph))] == 2)) + length(which(dat$attr$tt.traj.syph.hivpos[unique(c(idsInf_syph))] == 2))
  dat$epi$incid.sti.tttraj2[at] <- dat$epi$incid.gc.tttraj2[at] + dat$epi$incid.ct.tttraj2[at] + dat$epi$incid.syph.tttraj2[at]
  dat$epi$incid.gcct.tttraj2[at] <- dat$epi$incid.gc.tttraj2[at] + dat$epi$incid.ct.tttraj2[at]

  # Stop check for STI attributes
  stopifnot(all(!is.na(dat$attr$rGC.infTime[dat$attr$rGC == 1])),
            all(!is.na(dat$attr$rGC.sympt[dat$attr$rGC == 1])),
            all(!is.na(dat$attr$uGC.infTime[dat$attr$uGC == 1])),
            all(!is.na(dat$attr$uGC.sympt[dat$attr$uGC == 1])),
            all(!is.na(dat$attr$rCT.infTime[dat$attr$rCT == 1])),
            all(!is.na(dat$attr$rCT.sympt[dat$attr$rCT == 1])),
            all(!is.na(dat$attr$uCT.infTime[dat$attr$uCT == 1])),
            all(!is.na(dat$attr$uCT.sympt[dat$attr$uCT == 1])),
            all(!is.na(dat$attr$syph.infTime[dat$attr$syphilis == 1])))

  return(dat)
}


#' @title STI Recovery Module
#'
#' @description Stochastically simulates GC/CT recovery.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
sti_recov_msm <- function(dat, at) {

  # Parameters ----------------------------------------------------------

  rgc.asympt.int <- dat$param$rgc.asympt.int
  ugc.asympt.int <- dat$param$ugc.asympt.int
  gc.tx.int <- dat$param$gc.tx.int
  gc.ntx.int <- dat$param$gc.ntx.int

  rct.asympt.int <- dat$param$rct.asympt.int
  uct.asympt.int <- dat$param$uct.asympt.int
  ct.tx.int <- dat$param$ct.tx.int
  ct.ntx.int <- dat$param$ct.ntx.int

  syph.early.tx.int <- dat$param$syph.early.tx.int
  syph.late.tx.int <- dat$param$syph.late.tx.int

  # Attributes ----------------------------------------------------------

  # Infection status
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT
  syphilis <- dat$attr$syphilis
  stage.syph <- dat$attr$stage.syph

  # Infection time
  rGC.infTime <- dat$attr$rGC.infTime
  uGC.infTime <- dat$attr$uGC.infTime
  rCT.infTime <- dat$attr$rCT.infTime
  uCT.infTime <- dat$attr$uCT.infTime
  syph.infTime <- dat$attr$syph.infTime

  # Symptoms
  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  rCT.sympt <- dat$attr$rCT.sympt
  uCT.sympt <- dat$attr$uCT.sympt

  # Tx
  uGC.tx <- dat$attr$uGC.tx
  uGC.tx.prep <- dat$attr$uGC.tx.prep
  uGC.tx.ept <- dat$attr$uGC.tx.ept
  rGC.tx <- dat$attr$rGC.tx
  rGC.tx.prep <- dat$attr$rGC.tx.prep
  rGC.tx.ept <- dat$attr$rGC.tx.ept
  uCT.tx <- dat$attr$uCT.tx
  uCT.tx.prep <- dat$attr$uCT.tx.prep
  uCT.tx.ept <- dat$attr$uCT.tx.ept
  rCT.tx <- dat$attr$rCT.tx
  rCT.tx.prep <- dat$attr$rCT.tx.prep
  rCT.tx.ept <- dat$attr$rCT.tx.ept
  syph.incub.tx <- dat$attr$syph.incub.tx
  syph.prim.tx <- dat$attr$syph.prim.tx
  syph.seco.tx <- dat$attr$syph.seco.tx
  syph.earlat.tx <- dat$attr$syph.earlat.tx
  syph.latelat.tx <- dat$attr$syph.latelat.tx
  syph.tert.tx <- dat$attr$syph.tert.tx
  syph.tx.prep <- dat$attr$syph.tx.prep

  # GC Recovery ---------------------------------------------------------

  ## Recovery for asymptomatic untreated (natural clearance)
  idsRGC_asympt_ntx <- which(rGC == 1 &
                             rGC.infTime < at &
                             rGC.sympt == 0 &
                             (is.na(rGC.tx) | rGC.tx == 0) &
                             (is.na(rGC.tx.prep) | rGC.tx.prep == 0) &
                             (is.na(rGC.tx.ept) | rGC.tx.ept == 0))
  idsUGC_asympt_ntx <- which(uGC == 1 &
                             uGC.infTime < at &
                             uGC.sympt == 0 &
                             (is.na(uGC.tx) | uGC.tx == 0) &
                             (is.na(uGC.tx.prep) | uGC.tx.prep == 0) &
                             (is.na(uGC.tx.ept) | uGC.tx.ept == 0))

  recovRGC_asympt_ntx <- idsRGC_asympt_ntx[which(rbinom(length(idsRGC_asympt_ntx), 1, 1/rgc.asympt.int) == 1)]
  recovUGC_asympt_ntx <- idsUGC_asympt_ntx[which(rbinom(length(idsUGC_asympt_ntx), 1, 1/ugc.asympt.int) == 1)]

  ## Recovery for symptomatic untreated (natural clearance)
  idsRGC_sympt_ntx <- which(rGC == 1 &
                            rGC.infTime < at &
                            rGC.sympt == 1 &
                            (is.na(rGC.tx) | rGC.tx == 0) &
                            (is.na(rGC.tx.prep) | rGC.tx.prep == 0) &
                            (is.na(rGC.tx.ept) | rGC.tx.ept == 0))
  idsUGC_sympt_ntx <- which(uGC == 1 &
                            uGC.infTime < at &
                            uGC.sympt == 1 &
                            (is.na(uGC.tx) | uGC.tx == 0) &
                            (is.na(uGC.tx.prep) | uGC.tx.prep == 0) &
                            (is.na(uGC.tx.ept) | uGC.tx.ept == 0))

  # If NA, recovery rate for symptomatic untreated = rate for asymptomatic untreated
  if (!is.na(gc.ntx.int)) {
    recovRGC_sympt_ntx <- idsRGC_sympt_ntx[which(rbinom(length(idsRGC_sympt_ntx), 1, 1/gc.ntx.int) == 1)]
    recovUGC_sympt_ntx <- idsUGC_sympt_ntx[which(rbinom(length(idsUGC_sympt_ntx), 1, 1/gc.ntx.int) == 1)]
  } else {
    recovRGC_sympt_ntx <- idsRGC_sympt_ntx[which(rbinom(length(idsRGC_sympt_ntx), 1, 1/rgc.asympt.int) == 1)]
    recovUGC_sympt_ntx <- idsUGC_sympt_ntx[which(rbinom(length(idsUGC_sympt_ntx), 1, 1/ugc.asympt.int) == 1)]
  }

  ## Recovery for treated (both asymptomatic and symptomatic)
  idsRGC_tx <- which(rGC == 1 &
                     rGC.infTime < at &
                     (rGC.tx == 1 | rGC.tx.prep == 1 | rGC.tx.ept == 1))
  idsUGC_tx <- which(uGC == 1 &
                     uGC.infTime < at &
                     (uGC.tx == 1 | uGC.tx.prep == 1 | uGC.tx.ept == 1))

  recovRGC_tx <- idsRGC_tx[which(rbinom(length(idsRGC_tx), 1, 1/gc.tx.int) == 1)]
  recovUGC_tx <- idsUGC_tx[which(rbinom(length(idsUGC_tx), 1, 1/gc.tx.int) == 1)]

  recovRGC <- c(recovRGC_asympt_ntx, recovRGC_sympt_ntx, recovRGC_tx)
  recovUGC <- c(recovUGC_asympt_ntx, recovUGC_sympt_ntx, recovUGC_tx)

  rgc.infect.dur <- at - rGC.infTime[recovRGC]
  ugc.infect.dur <- at - uGC.infTime[recovUGC]
  gc.infect.dur <- mean(c(rgc.infect.dur, ugc.infect.dur), na.rm = TRUE)


  # CT Recovery ---------------------------------------------------------

  ## Recovery for asymptomatic untreated (natural clearance)
  idsRCT_asympt_ntx <- which(rCT == 1 &
                             rCT.infTime < at &
                             rCT.sympt == 0 &
                             (is.na(rCT.tx) | rCT.tx == 0) &
                             (is.na(rCT.tx.prep) | rCT.tx.prep == 0) &
                             (is.na(rCT.tx.ept) | rCT.tx.ept == 0))
  idsUCT_asympt_ntx <- which(uCT == 1 &
                             uCT.infTime < at &
                             uCT.sympt == 0 &
                             (is.na(uCT.tx) | uCT.tx == 0) &
                             (is.na(uCT.tx.prep) | uCT.tx.prep == 0) &
                             (is.na(uCT.tx.ept) | uCT.tx.ept == 0))

  recovRCT_asympt_ntx <- idsRCT_asympt_ntx[which(rbinom(length(idsRCT_asympt_ntx), 1, 1/rct.asympt.int) == 1)]
  recovUCT_asympt_ntx <- idsUCT_asympt_ntx[which(rbinom(length(idsUCT_asympt_ntx), 1, 1/uct.asympt.int) == 1)]

  ## Recovery for symptomatic untreated (natural clearance)
  idsRCT_sympt_ntx <- which(rCT == 1 &
                            rCT.infTime < at &
                            rCT.sympt == 1 &
                            (is.na(rCT.tx) | rCT.tx == 0) &
                            (is.na(rCT.tx.prep) | rCT.tx.prep == 0) &
                            (is.na(rCT.tx.ept) | rCT.tx.ept == 0))
  idsUCT_sympt_ntx <- which(uCT == 1 &
                            uCT.infTime < at &
                            uCT.sympt == 1 &
                            (is.na(uCT.tx) | uCT.tx == 0) &
                            (is.na(uCT.tx.prep) | uCT.tx.prep == 0) &
                            (is.na(rCT.tx.ept) | rCT.tx.ept == 0))

  # If NA, recovery rate for symptomatic untreated = rate for asymptomatic untreated
  if (!is.na(ct.ntx.int)) {
    recovRCT_sympt_ntx <- idsRCT_sympt_ntx[which(rbinom(length(idsRCT_sympt_ntx), 1, 1/ct.ntx.int) == 1)]
    recovUCT_sympt_ntx <- idsUCT_sympt_ntx[which(rbinom(length(idsUCT_sympt_ntx), 1, 1/ct.ntx.int) == 1)]
  } else {
    recovRCT_sympt_ntx <- idsRCT_sympt_ntx[which(rbinom(length(idsRCT_sympt_ntx), 1, 1/rct.asympt.int) == 1)]
    recovUCT_sympt_ntx <- idsUCT_sympt_ntx[which(rbinom(length(idsUCT_sympt_ntx), 1, 1/uct.asympt.int) == 1)]
  }

  ## Recovery for treated (both asymptomatic and symptomatic)
  idsRCT_tx <- which(rCT == 1 &
                     rCT.infTime < at &
                     (rCT.tx == 1 | rCT.tx.prep == 1 | rCT.tx.ept == 1))
  idsUCT_tx <- which(uCT == 1 &
                     uCT.infTime < at &
                     (uCT.tx == 1 | uCT.tx.prep == 1 | uCT.tx.ept == 1))

  recovRCT_tx <- idsRCT_tx[which(rbinom(length(idsRCT_tx), 1, 1/ct.tx.int) == 1)]
  recovUCT_tx <- idsUCT_tx[which(rbinom(length(idsUCT_tx), 1, 1/ct.tx.int) == 1)]

  recovRCT <- c(recovRCT_asympt_ntx, recovRCT_sympt_ntx, recovRCT_tx)
  recovUCT <- c(recovUCT_asympt_ntx, recovUCT_sympt_ntx, recovUCT_tx)

  rct.infect.dur <- at - rCT.infTime[recovRCT]
  uct.infect.dur <- at - uCT.infTime[recovUCT]
  ct.infect.dur <- mean(c(rct.infect.dur, uct.infect.dur), na.rm = TRUE)

  gcct.infect.dur <- mean(c(rgc.infect.dur, ugc.infect.dur, rct.infect.dur, uct.infect.dur), na.rm = TRUE)


  # Syphilis Recovery -------------------------------------------------

  ## Recovery for treated
  idssyph_early_tx <- which(syphilis == 1 &
                            stage.syph %in% c(1:4) &
                            syph.infTime < at & (syph.incub.tx == 1 |
                            syph.prim.tx == 1 | syph.seco.tx == 1 |
                            syph.earlat.tx == 1 | syph.tx.prep == 1))
  idssyph_late_tx <- which(syphilis == 1 &
                           stage.syph %in% c(5:6) &
                           syph.infTime < at &
                           (syph.latelat.tx == 1 | syph.tert.tx == 1 |
                              syph.tx.prep == 1))

  ## Move stage-specific treated to recovered
  recovsyph_early_tx <- idssyph_early_tx[which(rbinom(length(idssyph_early_tx), 1, 1/syph.early.tx.int) == 1)]
  recovsyph_late_tx <- idssyph_late_tx[which(rbinom(length(idssyph_late_tx), 1, 1/syph.late.tx.int) == 1)]

  ## Aggregate recovery by stage
  recovsyph <- c(recovsyph_early_tx, recovsyph_late_tx)
  syph.infect.dur <- mean(at - syph.infTime[recovsyph], na.rm = TRUE)


    # Output -----------------------------------------------------------

  ## All recovered
  recovGCCT <- c(recovUCT, recovRCT, recovRGC, recovUGC)

  # Reset EPT attributes
  dat$attr$eptindexEligdate[recovGCCT] <- NA
  dat$attr$eptpartEligReceive[recovGCCT] <- NA
  dat$attr$eptpartEligTx_GC[recovGCCT] <- NA
  dat$attr$eptpartEligTx_CT[recovGCCT] <- NA
  dat$attr$eptpartEligTxdate[recovGCCT] <- NA
  dat$attr$eptpartTx[recovGCCT] <- NA


  # Syphilis
  dat$attr$syphilis[recovsyph] <- 0
  dat$attr$stage.syph[recovsyph] <- NA
  dat$attr$stage.time.syph[recovsyph] <- NA
  dat$attr$syph.sympt[recovsyph] <- NA
  dat$attr$syph.infTime[recovsyph] <- NA
  dat$attr$diag.status.syph[recovsyph] <- NA
  dat$attr$syph.incub.tx[recovsyph] <- NA
  dat$attr$syph.prim.tx[recovsyph] <- NA
  dat$attr$syph.seco.tx[recovsyph] <- NA
  dat$attr$syph.earlat.tx[recovsyph] <- NA
  dat$attr$syph.latelat.tx[recovsyph] <- NA
  dat$attr$syph.tert.tx[recovsyph] <- NA
  dat$attr$syph.tx.prep[recovsyph] <- NA

  # Gonorrhea
  dat$attr$rGC[recovRGC] <- 0
  dat$attr$rGC.sympt[recovRGC] <- NA
  dat$attr$rGC.infTime[recovRGC] <- NA
  dat$attr$rGC.tx[recovRGC] <- NA
  dat$attr$rGC.tx.prep[recovRGC] <- NA
  dat$attr$rGC.tx.ept[recovRGC] <- NA
  dat$attr$diag.status.gc[recovRGC] <- NA
  dat$attr$uGC[recovUGC] <- 0
  dat$attr$uGC.sympt[recovUGC] <- NA
  dat$attr$uGC.infTime[recovUGC] <- NA
  dat$attr$uGC.tx[recovUGC] <- NA
  dat$attr$uGC.tx.prep[recovUGC] <- NA
  dat$attr$uGC.tx.ept[recovUGC] <- NA
  dat$attr$diag.status.gc[recovUGC] <- NA

  # Chlamydia
  dat$attr$rCT[recovRCT] <- 0
  dat$attr$rCT.sympt[recovRCT] <- NA
  dat$attr$rCT.infTime[recovRCT] <- NA
  dat$attr$rCT.tx[recovRCT] <- NA
  dat$attr$rCT.tx.prep[recovRCT] <- NA
  dat$attr$rCT.tx.ept[recovRCT] <- NA
  dat$attr$diag.status.ct[recovRCT] <- NA
  dat$attr$uCT[recovUCT] <- 0
  dat$attr$uCT.sympt[recovUCT] <- NA
  dat$attr$uCT.infTime[recovUCT] <- NA
  dat$attr$uCT.tx[recovUCT] <- NA
  dat$attr$uCT.tx.prep[recovUCT] <- NA
  dat$attr$uCT.tx.ept[recovUCT] <- NA
  dat$attr$diag.status.ct[recovUCT] <- NA

  # Summary stats
  dat$epi$recov.rgc[at] <- length(unique(recovRGC))
  dat$epi$recov.ugc[at] <- length(unique(recovUGC))
  dat$epi$recov.rct[at] <- length(unique(recovRCT))
  dat$epi$recov.uct[at] <- length(unique(recovUCT))
  dat$epi$recov.syphilis[at] <- length(unique(recovsyph))
  dat$epi$recov.earlysyph[at] <- length(unique(recovsyph_early_tx))

  dat$epi$gc.infect.dur[at] <- gc.infect.dur
  dat$epi$ct.infect.dur[at] <- ct.infect.dur
  dat$epi$gcct.infect.dur[at] <- gcct.infect.dur
  dat$epi$syph.infect.dur[at] <- syph.infect.dur

  return(dat)
}


#' @title STI Treatment Module
#'
#' @description Stochastically simulates GC/CT and syphilis diagnosis and
#'              treatment.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
sti_tx_msm <- function(dat, at) {

  # Parameters ------------------------------------------------------------
  gc.sympt.prob.tx <- dat$param$gc.sympt.prob.tx
  ct.sympt.prob.tx <- dat$param$ct.sympt.prob.tx

  gc.asympt.prob.tx <- dat$param$gc.asympt.prob.tx
  ct.asympt.prob.tx <- dat$param$ct.asympt.prob.tx

  syph.incub.sympt.prob.tx <- dat$param$syph.incub.sympt.prob.tx
  syph.incub.asympt.prob.tx <- dat$param$syph.incub.asympt.prob.tx
  syph.prim.sympt.prob.tx <- dat$param$syph.prim.sympt.prob.tx
  syph.prim.asympt.prob.tx <- dat$param$syph.prim.asympt.prob.tx
  syph.seco.sympt.prob.tx <- dat$param$syph.seco.sympt.prob.tx
  syph.seco.asympt.prob.tx <- dat$param$syph.seco.asympt.prob.tx
  syph.earlat.sympt.prob.tx <- dat$param$syph.earlat.sympt.prob.tx
  syph.earlat.asympt.prob.tx <- dat$param$syph.earlat.asympt.prob.tx
  syph.latelat.sympt.prob.tx <- dat$param$syph.latelat.sympt.prob.tx
  syph.latelat.asympt.prob.tx <- dat$param$syph.latelat.asympt.prob.tx
  syph.tert.sympt.prob.tx <- dat$param$syph.tert.sympt.prob.tx
  syph.tert.asympt.prob.tx <- dat$param$syph.tert.asympt.prob.tx

  prep.sti.screen.int <- dat$param$prep.sti.screen.int
  prep.sti.prob.tx <- dat$param$prep.sti.prob.tx
  ept.gc.success <- dat$param$ept.gc.success
  ept.ct.success <- dat$param$ept.ct.success
  ept.coverage <- dat$param$ept.coverage
  ept.cov.rate <- dat$param$ept.cov.rate

  sti.correlation.time <- dat$param$sti.correlation.time

  # Attributes ------------------------------------------------------------

  # Infection status
  role.class <- dat$attr$role.class
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT
  syphilis <- dat$attr$syphilis
  stage.syph <- dat$attr$stage.syph

  # Infection time
  rGC.infTime <- dat$attr$rGC.infTime
  uGC.infTime <- dat$attr$uGC.infTime
  rCT.infTime <- dat$attr$rCT.infTime
  uCT.infTime <- dat$attr$uCT.infTime
  syph.infTime <- dat$attr$syph.infTime

  # Symptoms
  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  rCT.sympt <- dat$attr$rCT.sympt
  uCT.sympt <- dat$attr$uCT.sympt
  syph.sympt <- dat$attr$syph.sympt

  # Tx
  uGC.tx <- dat$attr$uGC.tx
  uGC.tx.prep <- dat$attr$uGC.tx.prep
  uGC.tx.ept <- dat$attr$uGC.tx.ept
  rGC.tx <- dat$attr$rGC.tx
  rGC.tx.prep <- dat$attr$rGC.tx.prep
  rGC.tx.ept <- dat$attr$rGC.tx.ept
  uCT.tx <- dat$attr$uCT.tx
  uCT.tx.prep <- dat$attr$uCT.tx.prep
  uCT.tx.ept <- dat$attr$uCT.tx.ept
  rCT.tx <- dat$attr$rCT.tx
  rCT.tx.prep <- dat$attr$rCT.tx.prep
  rCT.tx.ept <- dat$attr$rCT.tx.ept
  syph.incub.tx <- dat$attr$syph.incub.tx
  syph.prim.tx <- dat$attr$syph.prim.tx
  syph.seco.tx <- dat$attr$syph.seco.tx
  syph.earlat.tx <- dat$attr$syph.earlat.tx
  syph.latelat.tx <- dat$attr$syph.latelat.tx
  syph.tert.tx <- dat$attr$syph.tert.tx
  syph.tx.prep <- dat$attr$syph.tx.prep
  prepStartTime <- dat$attr$prepStartTime
  prepLastStiScreen <- dat$attr$prepLastStiScreen

  # tt.traj
  tt.traj.gc.hivpos <- dat$attr$tt.traj.gc.hivpos
  tt.traj.gc.hivneg <- dat$attr$tt.traj.gc.hivneg
  tt.traj.ct.hivpos <- dat$attr$tt.traj.ct.hivpos
  tt.traj.ct.hivneg <- dat$attr$tt.traj.ct.hivneg
  tt.traj.syph.hivpos <- dat$attr$tt.traj.syph.hivpos
  tt.traj.syph.hivneg <- dat$attr$tt.traj.syph.hivneg

  # Diagnosis/ Attributes from testing
  diag.status.syph <- dat$attr$diag.status.syph
  diag.status.gc <- dat$attr$diag.status.gc
  diag.status.ct <- dat$attr$diag.status.ct
  tsinceltst.syph <- dat$attr$time.since.last.test.syph
  tsinceltst.rgc <- dat$attr$time.since.last.test.rgc
  tsinceltst.ugc <- dat$attr$time.since.last.test.ugc
  tsinceltst.rct <- dat$attr$time.since.last.test.rct
  tsinceltst.uct <- dat$attr$time.since.last.test.uct

  # EPT
  eptpartEligTx_GC <- dat$attr$eptpartEligTx_GC
  eptpartEligTx_CT <- dat$attr$eptpartEligTx_CT
  eptpartEligTxdate <- dat$attr$eptpartEligTxdate

  # Syphilis --------------------------------------------------------------

  ## Symptomatic syphilis treatment
  # Select those in incubating stage who are eligible to be treated
  idssyph_tx_sympt_incub <- which(syphilis == 1 &
                                   syph.infTime < at &
                                   stage.syph == 1 &
                                   syph.sympt == 1 &
                                   is.na(syph.incub.tx))

  # Select those who will be treated based on eligibility to be treated
  txsyph_sympt_incub <- idssyph_tx_sympt_incub[which(rbinom(length(idssyph_tx_sympt_incub),
                                                            1, syph.incub.sympt.prob.tx) == 1)]

  # Select those in primary stage who are eligible to be treated
  idssyph_tx_sympt_prim <- which(syphilis == 1 &
                                 syph.infTime < at &
                                 stage.syph == 2 &
                                 syph.sympt == 1 &
                                 is.na(syph.prim.tx))

  # Select those who will be treated based on eligibility to be treated
  txsyph_sympt_prim <- idssyph_tx_sympt_prim[which(rbinom(length(idssyph_tx_sympt_prim),
                                                          1, syph.prim.sympt.prob.tx) == 1)]

  # Select those in secondary stage who are eligible to be treated
  idssyph_tx_sympt_seco <- which(syphilis == 1 &
                                 syph.infTime < at &
                                 stage.syph == 3 &
                                 syph.sympt == 1 &
                                 is.na(syph.seco.tx))

  # Select those who will be treated based on eligibility to be treated
  txsyph_sympt_seco <- idssyph_tx_sympt_seco[which(rbinom(length(idssyph_tx_sympt_seco),
                                                          1, syph.seco.sympt.prob.tx) == 1)]

  # Select those in early latent stage who are eligible to be treated
  idssyph_tx_sympt_earlat <- which(syphilis == 1 &
                                   syph.infTime < at &
                                   stage.syph == 4 &
                                   syph.sympt == 1 &
                                   is.na(syph.earlat.tx))

  # Select those who will be treated based on eligibility to be treated
  txsyph_sympt_earlat <- idssyph_tx_sympt_earlat[which(rbinom(length(idssyph_tx_sympt_earlat),
                                                              1, syph.earlat.sympt.prob.tx) == 1)]

  # Select those in late latent stage who are eligible to be treated
  idssyph_tx_sympt_latelat <- which(syphilis == 1 &
                                    syph.infTime < at &
                                    (stage.syph == 5) &
                                    (syph.sympt == 1) &
                                    is.na(syph.latelat.tx))

  # Select those who will be treated based on eligibility to be treated
  txsyph_sympt_latelat <- idssyph_tx_sympt_latelat[which(rbinom(length(idssyph_tx_sympt_latelat),
                                                                1, syph.latelat.sympt.prob.tx) == 1)]

  # Select those in tertiary stage who are eligible to be treated
  idssyph_tx_sympt_tert <- which(syphilis == 1 &
                                 syph.infTime < at &
                                 stage.syph == 6 &
                                 syph.sympt == 1 &
                                 is.na(syph.tert.tx))

  # Select those who will be treated based on eligibility to be treated
  txsyph_sympt_tert <- idssyph_tx_sympt_tert[which(rbinom(length(idssyph_tx_sympt_tert),
                                                          1, syph.tert.sympt.prob.tx) == 1)]

  # Aggregate all those eligible to be treated
  idssyph_tx_sympt <- c(idssyph_tx_sympt_incub, idssyph_tx_sympt_prim, idssyph_tx_sympt_seco, idssyph_tx_sympt_earlat,
                        idssyph_tx_sympt_latelat, idssyph_tx_sympt_tert)

  # Aggregate all those selected to be treated
  txsyph_sympt <- c(txsyph_sympt_incub, txsyph_sympt_prim, txsyph_sympt_seco, txsyph_sympt_earlat,
                    txsyph_sympt_latelat, txsyph_sympt_tert)


  ## Asymptomatic syphilis treatment
  # Select those in primary stage who are eligible to be treated
  idssyph_tx_asympt_incub <- which(syph.infTime < at &
                                    stage.syph == 1 &
                                    syph.sympt == 0 &
                                    diag.status.syph == 1 &
                                    (is.na(syph.incub.tx) | syph.incub.tx == 0))

  # Select those to be treated
  txsyph_asympt_incub <- idssyph_tx_asympt_incub[which(rbinom(length(idssyph_tx_asympt_incub),
                                                            1, syph.incub.asympt.prob.tx) == 1)]
  # Select those in primary stage who are eligible to be treated
  idssyph_tx_asympt_prim <- which(syph.infTime < at &
                                  stage.syph == 2 &
                                  syph.sympt == 0 &
                                  diag.status.syph == 1 &
                                  (is.na(syph.prim.tx) | syph.prim.tx == 0))

  # Select those to be treated
  txsyph_asympt_prim <- idssyph_tx_asympt_prim[which(rbinom(length(idssyph_tx_asympt_prim),
                                                            1, syph.prim.asympt.prob.tx) == 1)]

  # Select those in secondary stage who are eligible to be treated
  idssyph_tx_asympt_seco <- which(syph.infTime < at &
                                  stage.syph == 3 &
                                  syph.sympt == 0 &
                                  diag.status.syph == 1 &
                                  (is.na(syph.seco.tx) | syph.seco.tx == 0))

  # Select those to be treated
  txsyph_asympt_seco <- idssyph_tx_asympt_seco[which(rbinom(length(idssyph_tx_asympt_seco),
                                                            1, syph.seco.asympt.prob.tx) == 1)]


  # Select those in early latent stage who are eligible to be treated
  idssyph_tx_asympt_earlat <- which(syph.infTime < at &
                                    stage.syph == 4 &
                                    syph.sympt == 0 &
                                    diag.status.syph == 1 &
                                    (is.na(syph.earlat.tx) | syph.earlat.tx == 0))

  # Select those to be treated
  txsyph_asympt_earlat <- idssyph_tx_asympt_earlat[which(rbinom(length(idssyph_tx_asympt_earlat),
                                                                1, syph.earlat.asympt.prob.tx) == 1)]

  # Select those in late latent stage who are eligible to be treated
  idssyph_tx_asympt_latelat <- which(syph.infTime < at &
                                     (stage.syph == 5) &
                                     (syph.sympt == 0) &
                                     diag.status.syph == 1 &
                                     (is.na(syph.latelat.tx) | syph.latelat.tx == 0))

  # Select those to be treated
  txsyph_asympt_latelat <- idssyph_tx_asympt_latelat[which(rbinom(length(idssyph_tx_asympt_latelat),
                                                                  1, syph.latelat.asympt.prob.tx) == 1)]
  # Select those in tertiary stage who are eligible to be treated
  idssyph_tx_asympt_tert <- which(syph.infTime < at &
                                  stage.syph == 6 &
                                  syph.sympt == 0 &
                                  diag.status.syph == 1 &
                                  (is.na(syph.tert.tx) | syph.tert.tx == 0))

  # Select those to be treated
  txsyph_asympt_tert <- idssyph_tx_asympt_tert[which(rbinom(length(idssyph_tx_asympt_tert),
                                                            1, syph.tert.asympt.prob.tx) == 1)]

  # Aggregate all those eligible to be treated
  idssyph_tx_asympt <- c(idssyph_tx_asympt_incub, idssyph_tx_asympt_prim, idssyph_tx_asympt_seco, idssyph_tx_asympt_earlat,
                         idssyph_tx_asympt_latelat, idssyph_tx_asympt_tert)

  # Aggregate all those selected to be treated
  txsyph_asympt <- c(txsyph_asympt_incub, txsyph_asympt_prim, txsyph_asympt_seco, txsyph_asympt_earlat,
                     txsyph_asympt_latelat, txsyph_asympt_tert)

  # All treated syphilis
  txsyph <- union(txsyph_sympt, txsyph_asympt)
  idssyph_tx <- union(idssyph_tx_sympt, idssyph_tx_asympt)

  # By stage
  idssyph_incub_tx <- union(idssyph_tx_asympt_incub, idssyph_tx_sympt_incub)
  idssyph_prim_tx <- union(idssyph_tx_asympt_prim, idssyph_tx_sympt_prim)
  idssyph_seco_tx <- union(idssyph_tx_asympt_seco, idssyph_tx_sympt_seco)
  idssyph_earlat_tx <- union(idssyph_tx_asympt_earlat, idssyph_tx_sympt_earlat)
  idssyph_latelat_tx <- union(idssyph_tx_asympt_latelat, idssyph_tx_sympt_latelat)
  idssyph_tert_tx <- union(idssyph_tx_asympt_tert, idssyph_tx_sympt_tert)

  syph_incub_tx <- union(txsyph_asympt_incub, txsyph_sympt_incub)
  syph_prim_tx <- union(txsyph_asympt_prim, txsyph_sympt_prim)
  syph_seco_tx <- union(txsyph_asympt_seco, txsyph_sympt_seco)
  syph_earlat_tx <- union(txsyph_asympt_earlat, txsyph_sympt_earlat)
  syph_latelat_tx <- union(txsyph_asympt_latelat, txsyph_sympt_latelat)
  syph_tert_tx <- union(txsyph_asympt_tert, txsyph_sympt_tert)


  # Gonorrhea ---------------------------------------------------------------

  # symptomatic gc treatment
  idsRGC_tx_sympt <- which(rGC == 1 &
                           rGC.infTime < at &
                           rGC.sympt == 1 &
                           is.na(rGC.tx))

  idsUGC_tx_sympt <- which(uGC == 1 &
                           uGC.infTime < at &
                           uGC.sympt == 1 &
                           is.na(uGC.tx))

  idsGC_tx_sympt <- c(idsRGC_tx_sympt, idsUGC_tx_sympt)

  txGC_sympt <- idsGC_tx_sympt[which(rbinom(length(idsGC_tx_sympt), 1, gc.sympt.prob.tx) == 1)]
  txRGC_sympt <- intersect(idsRGC_tx_sympt, txGC_sympt)
  txUGC_sympt <- intersect(idsUGC_tx_sympt, txGC_sympt)

  # asymptomatic gc treatment
  idsRGC_tx_asympt <- which(rGC == 1 &
                            rGC.infTime < at &
                            rGC.sympt == 0 &
                            diag.status.gc == 1 &
                            is.na(rGC.tx))

  idsUGC_tx_asympt <- which(uGC == 1 &
                            uGC.infTime < at &
                            uGC.sympt == 0 &
                            diag.status.gc == 1 &
                            is.na(uGC.tx))

  idsGC_tx_asympt <- c(idsRGC_tx_asympt, idsUGC_tx_asympt)
  txGC_asympt <- idsGC_tx_asympt[which(rbinom(length(idsGC_tx_asympt), 1, gc.asympt.prob.tx) == 1)]
  txRGC_asympt <- intersect(idsRGC_tx_asympt, txGC_asympt)
  txUGC_asympt <- intersect(idsUGC_tx_asympt, txGC_asympt)

  # all treated GC
  txRGC <- union(txRGC_sympt, txRGC_asympt)
  txUGC <- union(txUGC_sympt, txUGC_asympt)

  idsRGC_tx <- union(idsRGC_tx_sympt, idsRGC_tx_asympt)
  idsUGC_tx <- union(idsUGC_tx_sympt, idsUGC_tx_asympt)


  # Chlamydia ---------------------------------------------------------------

  # symptomatic ct treatment
  idsRCT_tx_sympt <- which(rCT == 1 &
                           rCT.infTime < at &
                           rCT.sympt == 1 &
                           is.na(rCT.tx))

  idsUCT_tx_sympt <- which(uCT == 1 &
                           uCT.infTime < at &
                           uCT.sympt == 1 &
                           is.na(uCT.tx))
  idsCT_tx_sympt <- c(idsRCT_tx_sympt, idsUCT_tx_sympt)

  txCT_sympt <- idsCT_tx_sympt[which(rbinom(length(idsCT_tx_sympt), 1, ct.sympt.prob.tx) == 1)]
  txRCT_sympt <- intersect(idsRCT_tx_sympt, txCT_sympt)
  txUCT_sympt <- intersect(idsUCT_tx_sympt, txCT_sympt)

  # asymptomatic ct treatment
  idsRCT_tx_asympt <- which(rCT == 1 &
                            rCT.infTime < at &
                            rCT.sympt == 0 &
                            diag.status.ct == 1 &
                            is.na(rCT.tx))

  idsUCT_tx_asympt <- which(uCT == 1 &
                            uCT.infTime < at &
                            uCT.sympt == 0 &
                            diag.status.ct == 1 &
                            is.na(uCT.tx))

  idsCT_tx_asympt <- c(idsRCT_tx_asympt, idsUCT_tx_asympt)
  txCT_asympt <- idsCT_tx_asympt[which(rbinom(length(idsCT_tx_asympt), 1, ct.asympt.prob.tx) == 1)]
  txRCT_asympt <- intersect(idsRCT_tx_asympt, txCT_asympt)
  txUCT_asympt <- intersect(idsUCT_tx_asympt, txCT_asympt)

  # all treated CT
  txRCT <- union(txRCT_sympt, txRCT_asympt)
  txUCT <- union(txUCT_sympt, txUCT_asympt)

  idsRCT_tx <- union(idsRCT_tx_sympt, idsRCT_tx_asympt)
  idsUCT_tx <- union(idsUCT_tx_sympt, idsUCT_tx_asympt)


  # PrEP-related treatment for all STIs --------------------------------------

  # Interval-based treatment for MSM on PrEP
  idsSTI_screen <- which(prepStartTime == at | (at - prepLastStiScreen >= prep.sti.screen.int))

  prepLastStiScreen[idsSTI_screen] <- at

  idsRGC_prep_tx <- intersect(idsSTI_screen, which(rGC == 1 & rGC.infTime < at & is.na(rGC.tx.prep)))
  idsUGC_prep_tx <- intersect(idsSTI_screen, which(uGC == 1 & uGC.infTime < at & is.na(uGC.tx.prep)))
  idsRCT_prep_tx <- intersect(idsSTI_screen, which(rCT == 1 & rCT.infTime < at & is.na(rCT.tx.prep)))
  idsUCT_prep_tx <- intersect(idsSTI_screen, which(uCT == 1 & uCT.infTime < at & is.na(uCT.tx.prep)))
  idssyph_prep_tx <- intersect(idsSTI_screen, which(syphilis == 1 & syph.infTime < at & is.na(syph.tx.prep)))

  txRGC_prep <- idsRGC_prep_tx[which(rbinom(length(idsRGC_prep_tx), 1, prep.sti.prob.tx) == 1)]
  txUGC_prep <- idsUGC_prep_tx[which(rbinom(length(idsUGC_prep_tx), 1, prep.sti.prob.tx) == 1)]
  txRCT_prep <- idsRCT_prep_tx[which(rbinom(length(idsRCT_prep_tx), 1, prep.sti.prob.tx) == 1)]
  txUCT_prep <- idsUCT_prep_tx[which(rbinom(length(idsUCT_prep_tx), 1, prep.sti.prob.tx) == 1)]
  txsyph_prep <- idssyph_prep_tx[which(rbinom(length(idssyph_prep_tx), 1, prep.sti.prob.tx) == 1)]

  # Summarize all treated for each STI - EPT-treated not eligible to provide EPT to their partners
  txRGC_all <- c(txRGC, txRGC_prep)
  txUGC_all <- c(txUGC, txUGC_prep)
  txRCT_all <- c(txRCT, txRCT_prep)
  txUCT_all <- c(txUCT, txUCT_prep)
  txsyph_all <- c(txsyph, txsyph_prep)

  # Subset all treated for GC/CT to treated with partners (for EPT)
  ept_txRGC_all <- txRGC_all[dat$attr$recentpartners[txRGC_all] > 0]
  ept_txUGC_all <- txUGC_all[dat$attr$recentpartners[txUGC_all] > 0]
  ept_txRCT_all <- txRCT_all[dat$attr$recentpartners[txRCT_all] > 0]
  ept_txUCT_all <- txUCT_all[dat$attr$recentpartners[txUCT_all] > 0]
  ept_tx_all <- unique(c(ept_txRGC_all, ept_txUGC_all, ept_txRCT_all, ept_txUCT_all))

  # Update EPT index status and eligibility for GC/CT treated with partners
  dat$attr$eptindexElig[ept_tx_all] <- 1
  dat$attr$eptindexStat[ept_tx_all] <- 0
  dat$attr$eptindexEligdate[ept_tx_all] <- at

  # EPT Treatment for Non-index (no test is done) ------------------------------

  # Have prevalent infection, are eligible for tx through EPT, are untreated,
  # are not previously assigned for EPT tx, and were provided/uptake EPT last
  # time step
  idsRGC_tx_ept <- which(rGC == 1 &
                             rGC.infTime < at &
                             eptpartEligTx_GC == 1 &
                             eptpartEligTxdate == (at - 1) &
                             (is.na(rGC.tx) | rGC.tx == 0) &
                             (is.na(rGC.tx.prep) | rGC.tx.prep == 0) &
                             is.na(rGC.tx.ept))

  idsUGC_tx_ept <- which(uGC == 1 &
                             uGC.infTime < at &
                             eptpartEligTx_GC == 1 &
                             eptpartEligTxdate == (at - 1) &
                             (is.na(uGC.tx) | uGC.tx == 0) &
                             (is.na(uGC.tx.prep) | uGC.tx.prep == 0) &
                             is.na(uGC.tx.ept))

  idsGC_tx_ept <- c(idsRGC_tx_ept, idsUGC_tx_ept)
  txGC_ept <- idsGC_tx_ept[which(rbinom(length(idsGC_tx_ept), 1, ept.gc.success) == 1)]
  txRGC_ept <- intersect(idsRGC_tx_ept, txGC_ept)
  txUGC_ept <- intersect(idsUGC_tx_ept, txGC_ept)

  idsRCT_tx_ept <- which(rCT == 1 &
                             rCT.infTime < at &
                             eptpartEligTx_CT == 1 &
                             eptpartEligTxdate == (at - 1) &
                             (is.na(rCT.tx)  | rCT.tx == 0) &
                             (is.na(rCT.tx.prep) | rCT.tx.prep == 0) &
                             is.na(rCT.tx.ept))

  idsUCT_tx_ept <- which(uCT == 1 &
                             uCT.infTime < at &
                             eptpartEligTx_CT == 1 &
                             eptpartEligTxdate == (at - 1) &
                             (is.na(uCT.tx)  | uCT.tx == 0) &
                             (is.na(uCT.tx.prep) | uCT.tx.prep == 0) &
                             is.na(uCT.tx.ept))
  idsCT_tx_ept <- c(idsRCT_tx_ept, idsUCT_tx_ept)

  txCT_ept <- idsCT_tx_ept[which(rbinom(length(idsCT_tx_ept), 1, ept.ct.success) == 1)]
  txRCT_ept <- intersect(idsRCT_tx_ept, txCT_ept)
  txUCT_ept <- intersect(idsUCT_tx_ept, txCT_ept)

  # All EPT-treated index ids
  allidsept <- unique(c(idsCT_tx_ept, idsGC_tx_ept))

  # Summarize all successfully treated for each STI, now including EPT
  alltxRGC <- c(txRGC, txRGC_prep, txRGC_ept)
  alltxUGC <- c(txUGC, txUGC_prep, txUGC_ept)
  alltxRCT <- c(txRCT, txRCT_prep, txRCT_ept)
  alltxUCT <- c(txUCT, txUCT_prep, txUCT_ept)
  alltxEPT <- c(txRGC_ept, txUGC_ept, txRCT_ept, txUCT_ept)

  # EPT Initiation for Index Partner -------------------------------------------

  # Eligibility only lasts one time step - so coverage is 0 for current eligibles
  eptCov <- 0
  idsEligSt <- which(dat$attr$eptindexElig == 1 & dat$attr$eptindexEligdate == at)
  nEligSt <- length(idsEligSt)

  nStart <- max(0, min(nEligSt, round((ept.coverage - eptCov) * nEligSt)))
  ept_idsStart <- NULL
  if (nStart > 0) {
    if (ept.cov.rate >= 1) {
      ept_idsStart <- ssample(idsEligSt, nStart)
    } else {
      ept_idsStart <- idsEligSt[rbinom(nStart, 1, ept.cov.rate) == 1]
    }
  }
  eptCov <- (length(ept_idsStart)) / nEligSt

  # Update EPT index status for those selected to receive EPT for their partners
  dat$attr$eptindexStat[ept_idsStart] <- 1

  # Correlated testing for other STIs if symptomatic for one -------------------
  if (dat$param$sti.stitx.correlation == "true") {

    # All treated for other STIs, minus those getting treated for STI through EPT
    tst.rgc <- setdiff(c(txsyph_sympt, txUGC_sympt, txRCT_sympt, txUCT_sympt), txRGC_ept)
    tst.ugc <- setdiff(c(txsyph_sympt, txRGC_sympt, txRCT_sympt, txUCT_sympt), txUGC_ept)
    tst.rct <- setdiff(c(txsyph_sympt, txUGC_sympt, txRGC_sympt, txUCT_sympt), txRCT_ept)
    tst.uct <- setdiff(c(txsyph_sympt, txUGC_sympt, txRGC_sympt, txRCT_sympt), txUCT_ept)
    tst.syph <- c(txRGC_sympt, txUGC_sympt, txRCT_sympt, txUCT_sympt)

    # Remove those just treated for STI (either sympt/asympt) from testing
    tst.rgc <- setdiff(tst.rgc, txRGC)
    tst.ugc <- setdiff(tst.ugc, txUGC)
    tst.rct <- setdiff(tst.rct, txRCT)
    tst.uct <- setdiff(tst.uct, txUCT)
    tst.syph <- setdiff(tst.syph, txsyph)

    # Subset to those not tested for particular STI recently (in last 3 weeks)
    tst.rgc <- tst.rgc[which(tsinceltst.rgc[tst.rgc] > sti.correlation.time &
                               (is.na(diag.status.gc[tst.rgc]) | diag.status.gc[tst.rgc]) &
                         role.class[tst.rgc] %in% c("R", "V"))]
    tst.ugc <- tst.ugc[which(tsinceltst.ugc[tst.ugc] > sti.correlation.time &
                               (is.na(diag.status.gc[tst.ugc]) | diag.status.gc[tst.ugc]) &
                         role.class[tst.ugc] %in% c("I", "V"))]
    tst.rct <- tst.rct[which(tsinceltst.rct[tst.rct] > sti.correlation.time &
                               (is.na(diag.status.ct[tst.rct]) | diag.status.ct[tst.rct]) &
                         role.class[tst.rct] %in% c("R", "V"))]
    tst.uct <- tst.uct[which(tsinceltst.uct[tst.uct] > sti.correlation.time &
                               (is.na(diag.status.ct[tst.uct]) | diag.status.ct[tst.uct]) &
                         role.class[tst.uct] %in% c("I", "V"))]
    tst.syph <- tst.syph[which(tsinceltst.syph[tst.syph] > sti.correlation.time &
                               is.na(diag.status.syph[tst.syph]) | diag.status.syph[tst.syph])]

    tst.syph.pos <- tst.syph[which(syphilis[tst.syph] == 1 &
                                     stage.syph[tst.syph] %in% c(2, 3, 4, 5, 6))]
    tst.syph.neg <- setdiff(tst.syph, tst.syph.pos)
    tst.earlysyph.pos <- tst.syph[which(syphilis[tst.syph] == 1 &
                                          stage.syph[tst.syph] %in% c(2, 3, 4))]
    tst.latesyph.pos <- tst.syph[which(syphilis[tst.syph] == 1 &
                                         stage.syph[tst.syph] %in% c(5, 6))]

    tst.rgc.pos <- tst.rgc[which(rGC[tst.rgc] == 1)]
    tst.ugc.pos <- tst.ugc[which(uGC[tst.ugc] == 1)]
    tst.rgc.neg <- setdiff(tst.rgc, tst.ugc.pos)
    tst.ugc.neg <- setdiff(tst.ugc, tst.ugc.pos)
    tst.gc.pos <- c(tst.rgc.pos, tst.ugc.pos)

    tst.rct.pos <- tst.rct[which(rCT[tst.rct] == 1)]
    tst.uct.pos <- tst.uct[which(uCT[tst.uct] == 1)]
    tst.rct.neg <- setdiff(tst.rct, tst.uct.pos)
    tst.uct.neg <- setdiff(tst.uct, tst.uct.pos)
    tst.ct.pos <- c(tst.rct.pos, tst.uct.pos)

    # Syphilis Attributes
    dat$attr$last.neg.test.syph[tst.syph.neg] <- at
    dat$attr$last.neg.test.syph[tst.syph.pos] <- NA
    dat$attr$diag.status.syph[tst.syph.pos] <- 1
    dat$attr$last.diag.time.syph[tst.syph.pos] <- at
    dat$attr$time.since.last.test.syph[tst.syph] <- 0

    # GC Attributes
    dat$attr$last.neg.test.rgc[tst.rgc.neg] <- at
    dat$attr$last.neg.test.ugc[tst.ugc.neg] <- at
    dat$attr$last.neg.test.rgc[tst.rgc.pos] <- NA
    dat$attr$last.neg.test.ugc[tst.ugc.pos] <- NA
    dat$attr$diag.status.gc[tst.gc.pos] <- 1
    dat$attr$last.diag.time.gc[tst.gc.pos] <- at
    dat$attr$time.since.last.test.rgc[tst.rgc] <- 0
    dat$attr$time.since.last.test.ugc[tst.ugc] <- 0

    # CT Attributes
    dat$attr$last.neg.test.rct[tst.rct.neg] <- at
    dat$attr$last.neg.test.uct[tst.uct.neg] <- at
    dat$attr$last.neg.test.rct[tst.rct.pos] <- NA
    dat$attr$last.neg.test.uct[tst.uct.pos] <- NA
    dat$attr$diag.status.ct[tst.ct.pos] <- 1
    dat$attr$last.diag.time.ct[tst.ct.pos] <- at
    dat$attr$time.since.last.test.rct[tst.rct] <- 0
    dat$attr$time.since.last.test.uct[tst.uct] <- 0

    dat$epi$rGC_symptstidxtime[at] <- length(tst.rgc)
    dat$epi$uGC_symptstidxtime[at] <- length(tst.ugc)
    dat$epi$rCT_symptstidxtime[at] <- length(tst.rct)
    dat$epi$uCT_symptstidxtime[at] <- length(tst.uct)
    dat$epi$syph_symptstidxtime[at] <- length(tst.syph)

    dat$epi$rGC_pos_symptstidxtime[at] <- length(tst.rgc.pos)
    dat$epi$uGC_pos_symptstidxtime[at] <- length(tst.ugc.pos)
    dat$epi$rCT_pos_symptstidxtime[at] <- length(tst.rct.pos)
    dat$epi$uCT_pos_symptstidxtime[at] <- length(tst.uct.pos)
    dat$epi$syph_pos_symptstidxtime[at] <- length(tst.syph.pos)
    dat$epi$syph_earlypos_symptstidxtime[at] <- length(tst.earlysyph.pos)
    dat$epi$syph_latepos_symptstidxtime[at] <- length(tst.latesyph.pos)

    # Update total tests to now include symptomatic STI-correlated testing
    dat$epi$rGCasympttests[at] <- dat$epi$rGCasympttests[at] + dat$epi$rGC_symptstidxtime[at]
    dat$epi$uGCasympttests[at] <- dat$epi$uGCasympttests[at] + dat$epi$uGC_symptstidxtime[at]
    dat$epi$GCasympttests[at] <- dat$epi$rGCasympttests[at] + dat$epi$uGCasympttests[at] +
      dat$epi$rGC_symptstidxtime[at] + dat$epi$uGC_symptstidxtime[at]

    dat$epi$rGCasympttests.pos[at] <- dat$epi$rGCasympttests.pos[at] + dat$epi$rGC_pos_symptstidxtime[at]
    dat$epi$uGCasympttests.pos[at] <- dat$epi$uGCasympttests.pos[at] + dat$epi$uGC_pos_symptstidxtime[at]
    dat$epi$GCasympttests.pos[at] <- dat$epi$rGCasympttests.pos[at] + dat$epi$uGCasympttests.pos[at] +
      dat$epi$rGC_pos_symptstidxtime[at] + dat$epi$uGC_pos_symptstidxtime[at]

    dat$epi$rCTasympttests[at] <- dat$epi$rCTasympttests[at] + dat$epi$rCT_symptstidxtime[at]
    dat$epi$uCTasympttests[at] <- dat$epi$uCTasympttests[at] + dat$epi$uCT_symptstidxtime[at]
    dat$epi$CTasympttests[at] <- dat$epi$rCTasympttests[at] + dat$epi$uCTasympttests[at] +
      dat$epi$rCT_symptstidxtime[at] + dat$epi$uCT_symptstidxtime[at]

    dat$epi$rCTasympttests.pos[at] <-  dat$epi$rCTasympttests.pos[at] +
      dat$epi$rCT_pos_symptstidxtime[at]
    dat$epi$uCTasympttests.pos[at] <- dat$epi$uCTasympttests.pos[at] +
      dat$epi$uCT_pos_symptstidxtime[at]
    dat$epi$CTasympttests.pos[at] <- dat$epi$rCTasympttests.pos[at] + dat$epi$uCTasympttests.pos[at] +
      dat$epi$rCT_pos_symptstidxtime[at] + dat$epi$uCT_pos_symptstidxtime[at]

    dat$epi$syphasympttests[at] <- dat$epi$syphasympttests[at] +
      dat$epi$syph_symptstidxtime[at]
    dat$epi$syphasympttests.pos[at] <- dat$epi$syphasympttests.pos[at] +
      dat$epi$syph_pos_symptstidxtime[at]
    dat$epi$syphearlyasympttests.pos[at] <- dat$epi$syphearlyasympttests.pos[at] +
      dat$epi$syph_earlypos_symptstidxtime[at]
    dat$epi$syphlateasympttests.pos[at] <- dat$epi$syphlateasympttests.pos[at] +
      dat$epi$syph_latepos_symptstidxtime[at]

    dat$epi$stiasympttests[at] <- dat$epi$rGCasympttests[at] + dat$epi$uGCasympttests[at] +
      dat$epi$rCTasympttests[at] + dat$epi$uCTasympttests[at] + dat$epi$syphasympttests[at] +
      dat$epi$rGC_symptstidxtime[at] + dat$epi$uGC_symptstidxtime[at] +
      dat$epi$rCT_symptstidxtime[at] + dat$epi$uCT_symptstidxtime[at] +
      dat$epi$syph_symptstidxtime[at]
    dat$epi$stiasympttests[at] <- dat$epi$rGCasympttests.pos[at] + dat$epi$uGCasympttests.pos[at] +
      dat$epi$rCTasympttests.pos[at] + dat$epi$uCTasympttests.pos[at] + dat$epi$syphasympttests.pos[at] +
      dat$epi$rGC_pos_symptstidxtime[at] + dat$epi$uGC_pos_symptstidxtime[at] +
      dat$epi$rCT_pos_symptstidxtime[at] + dat$epi$uCT_pos_symptstidxtime[at] +
      dat$epi$syph_pos_symptstidxtime[at]

  }


  # Output ---------------------------------------------------------------------
  # PrEP
  dat$attr$prepLastStiScreen <- prepLastStiScreen

  # Syphilis
  dat$attr$syph.incub.tx[idssyph_incub_tx] <- 0
  dat$attr$syph.prim.tx[idssyph_prim_tx] <- 0
  dat$attr$syph.seco.tx[idssyph_seco_tx] <- 0
  dat$attr$syph.earlat.tx[idssyph_earlat_tx] <- 0
  dat$attr$syph.latelat.tx[idssyph_latelat_tx] <- 0
  dat$attr$syph.tert.tx[idssyph_tert_tx] <- 0
  dat$attr$syph.incub.tx[syph_incub_tx] <- 1
  dat$attr$syph.prim.tx[syph_prim_tx] <- 1
  dat$attr$syph.seco.tx[syph_seco_tx] <- 1
  dat$attr$syph.earlat.tx[syph_earlat_tx] <- 1
  dat$attr$syph.latelat.tx[syph_latelat_tx] <- 1
  dat$attr$syph.tert.tx[syph_tert_tx] <- 1
  dat$attr$last.tx.time.syph[txsyph_all] <- at
  dat$attr$syph.tx.prep[idssyph_prep_tx] <- 0
  dat$attr$syph.tx.prep[txsyph_prep] <- 1
  dat$attr$last.tx.time.syph.prep[txsyph_prep] <- at
  dat$attr$last.diag.time.syph[txsyph_sympt] <- at
  dat$attr$last.neg.test.syph[txsyph_sympt] <- NA
  dat$attr$time.since.last.test.syph[txsyph_sympt] <- 0
  dat$attr$diag.status.syph[idssyph_tx_sympt] <- 0
  dat$attr$diag.status.syph[txsyph_sympt] <- 1

  # Gonorrhea
  dat$attr$rGC.tx[idsRGC_tx] <- 0
  dat$attr$rGC.tx[txRGC] <- 1
  dat$attr$uGC.tx[idsUGC_tx] <- 0
  dat$attr$uGC.tx[txUGC] <- 1
  dat$attr$rGC.tx.prep[idsRGC_prep_tx] <- 0
  dat$attr$rGC.tx.prep[txRGC_prep] <- 1
  dat$attr$rGC.tx.ept[idsRGC_tx_ept] <- 0
  dat$attr$rGC.tx.ept[txRGC_ept] <- 1
  dat$attr$last.tx.time.rgc[txRGC_all] <- at
  dat$attr$last.tx.time.ugc[txUGC_all] <- at
  dat$attr$last.tx.time.rgc.prep[txRGC_prep] <- at
  dat$attr$last.tx.time.ugc.prep[txUGC_prep] <- at
  dat$attr$uGC.tx.prep[idsUGC_prep_tx] <- 0
  dat$attr$uGC.tx.prep[txUGC_prep] <- 1
  dat$attr$uGC.tx.ept[idsUGC_tx_ept] <- 0
  dat$attr$uGC.tx.ept[txUGC_ept] <- 1
  dat$attr$rGC.tx[which((dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1 | dat$attr$uGC.tx.ept == 1) & dat$attr$rGC == 1)] <- 1
  dat$attr$uGC.tx[which((dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1 | dat$attr$rGC.tx.ept == 1) & dat$attr$uGC == 1)] <- 1
  dat$attr$last.diag.time.gc[txGC_sympt] <- at
  dat$attr$last.neg.test.rgc[txRGC_sympt] <- NA
  dat$attr$last.neg.test.ugc[txUGC_sympt] <- NA
  dat$attr$time.since.last.test.rgc[txRGC_sympt] <- 0
  dat$attr$time.since.last.test.ugc[txUGC_sympt] <- 0
  dat$attr$diag.status.gc[idsGC_tx_sympt] <- 0
  dat$attr$diag.status.gc[txGC_sympt] <- 1

  # Chlamydia
  dat$attr$rCT.tx[idsRCT_tx] <- 0
  dat$attr$rCT.tx[txRCT] <- 1
  dat$attr$uCT.tx[idsUCT_tx] <- 0
  dat$attr$uCT.tx[txUCT] <- 1
  dat$attr$rCT.tx.prep[idsRCT_prep_tx] <- 0
  dat$attr$rCT.tx.prep[txRCT_prep] <- 1
  dat$attr$rCT.tx.ept[idsRCT_tx_ept] <- 0
  dat$attr$rCT.tx.ept[txRCT_ept] <- 1
  dat$attr$last.tx.time.rct[txRCT_all] <- at
  dat$attr$last.tx.time.uct[txUCT_all] <- at
  dat$attr$last.tx.time.rct.prep[txRCT_prep] <- at
  dat$attr$last.tx.time.uct.prep[txUCT_prep] <- at
  dat$attr$uCT.tx.prep[idsUCT_prep_tx] <- 0
  dat$attr$uCT.tx.prep[txUCT_prep] <- 1
  dat$attr$uCT.tx.ept[idsUCT_tx_ept] <- 0
  dat$attr$uCT.tx.ept[txUCT_ept] <- 1
  dat$attr$rCT.tx[which((dat$attr$uCT.tx == 1 | dat$attr$uCT.tx.prep == 1 | dat$attr$uCT.tx.ept == 1) & dat$attr$rCT == 1)] <- 1
  dat$attr$uCT.tx[which((dat$attr$rCT.tx == 1 | dat$attr$rCT.tx.prep == 1 | dat$attr$rCT.tx.ept == 1) & dat$attr$uCT == 1)] <- 1
  dat$attr$last.diag.time.ct[txCT_sympt] <- at
  dat$attr$last.neg.test.rct[txRCT_sympt] <- NA
  dat$attr$last.neg.test.uct[txUCT_sympt] <- NA
  dat$attr$time.since.last.test.rct[txRCT_sympt] <- 0
  dat$attr$time.since.last.test.uct[txUCT_sympt] <- 0
  dat$attr$diag.status.ct[idsCT_tx_sympt] <- 0
  dat$attr$diag.status.ct[txCT_sympt] <- 1

  # Proportion of infections treated in past year
  dat$epi$tx.gc.prop[at] <- (length(which(dat$attr$last.tx.time.rgc <= 52 & dat$attr$last.rGC.infTime <= 52)) + length(which(dat$attr$last.tx.time.ugc <= 52 & dat$attr$last.uGC.infTime <= 52))) / (length(which(dat$attr$last.rGC.infTime <= 52)) + length(which(dat$attr$last.uGC.infTime <= 52)))

  dat$epi$tx.ct.prop[at] <- (length(which(dat$attr$last.tx.time.rct <= 52 & dat$attr$last.rCT.infTime <= 52)) + length(which(dat$attr$last.tx.time.uct <= 52 & dat$attr$last.uCT.infTime <= 52))) / (length(which(dat$attr$last.rCT.infTime <= 52)) + length(which(dat$attr$last.uCT.infTime <= 52)))

  dat$epi$tx.gcct.prop[at] <- (length(which(dat$attr$last.tx.time.rgc <= 52 & dat$attr$last.rGC.infTime <= 52)) + length(which(dat$attr$last.tx.time.ugc <= 52 & dat$attr$last.uGC.infTime <= 52)) + length(which(dat$attr$last.tx.time.rct <= 52 & dat$attr$last.rCT.infTime <= 52)) + length(which(dat$attr$last.tx.time.uct <= 52 & dat$attr$last.uCT.infTime <= 52))) /
    (length(which(dat$attr$last.rGC.infTime <= 52)) + length(which(dat$attr$last.uGC.infTime <= 52)) + length(which(dat$attr$last.rCT.infTime <= 52)) + length(which(dat$attr$last.uCT.infTime <= 52)))

  dat$epi$tx.syph.prop[at] <- length(which(dat$attr$last.tx.time.syph <= 52 & dat$attr$last.syph.infTime <= 52)) / length(which(dat$attr$last.syph.infTime <= 52))


  # Non-index EPT-treated
  dat$attr$eptpartEligTx_GC[txGC_ept] <- NA
  dat$attr$eptpartEligTx_CT[txCT_ept] <- NA
  dat$attr$eptpartEligTxdate[alltxEPT] <- NA
  dat$attr$eptpartTx[allidsept] <- 0
  dat$attr$eptpartTx[alltxEPT] <- 1

  # summary statistics
  if (is.null(dat$epi$num.asympt.tx)) {
    dat$epi$num.asympt.tx <- rep(NA, length(dat$control$nsteps))
    dat$epi$num.asympt.cases <- rep(NA, length(dat$control$nsteps))
    dat$epi$num.asympt.tx.prep <- rep(NA, length(dat$control$nsteps))
    dat$epi$num.asympt.cases.prep <- rep(NA, length(dat$control$nsteps))
    dat$epi$num.rect.tx <- rep(NA, length(dat$epi$control$nsteps))
    dat$epi$num.rect.cases <- rep(NA, length(dat$epi$control$nsteps))
    dat$epi$num.rect.tx.prep <- rep(NA, length(dat$epi$control$nsteps))
    dat$epi$num.rect.cases.prep <- rep(NA, length(dat$epi$control$nsteps))
  }

  # Number of tests for symptomatic
  dat$epi$rGCsympttests[at] <- length(unique(txRGC_sympt))
  dat$epi$uGCsympttests[at] <- length(unique(txUGC_sympt))
  dat$epi$GCsympttests[at] <- length(unique(txRGC_sympt)) + length(unique(txUGC_sympt))

  dat$epi$rCTsympttests[at] <- length(unique(txRCT_sympt))
  dat$epi$uCTsympttests[at] <- length(unique(txUCT_sympt))
  dat$epi$CTsympttests[at] <- length(unique(txRCT_sympt)) + length(unique(txUCT_sympt))

  dat$epi$syphsympttests[at] <- length(unique(txsyph_sympt))

  dat$epi$stisympttests[at] <- length(unique(txRGC_sympt)) + length(unique(txUGC_sympt)) +
    length(unique(txRCT_sympt)) + length(unique(txUCT_sympt)) + length(unique(txsyph_sympt))

  asympt.tx <- c(intersect(txRGC_all, which(dat$attr$rGC.sympt == 0)),
                 intersect(txUGC_all, which(dat$attr$uGC.sympt == 0)),
                 intersect(txRCT_all, which(dat$attr$rCT.sympt == 0)),
                 intersect(txUCT_all, which(dat$attr$uCT.sympt == 0)),
                 intersect(txsyph_all, which(dat$attr$syph.sympt == 0)))
  dat$epi$num.asympt.tx[at] <- length(unique(asympt.tx))

  asympt.cases <- c(idsRGC_tx_asympt, intersect(idsRGC_prep_tx, which(dat$attr$rGC.sympt == 0)),
                    idsUGC_tx_asympt, intersect(idsUGC_prep_tx, which(dat$attr$uGC.sympt == 0)),
                    idsRCT_tx_asympt, intersect(idsRCT_prep_tx, which(dat$attr$rCT.sympt == 0)),
                    idsUCT_tx_asympt, intersect(idsUCT_prep_tx, which(dat$attr$uCT.sympt == 0)),
                    idssyph_tx_asympt, intersect(idssyph_prep_tx, which(dat$attr$syph.sympt == 0)))
  dat$epi$num.asympt.cases[at] <- length(unique(asympt.cases))

  asympt.tx.prep <- c(intersect(txRGC_prep, which(dat$attr$rGC.sympt == 0)),
                      intersect(txUGC_prep, which(dat$attr$uGC.sympt == 0)),
                      intersect(txRCT_prep, which(dat$attr$rCT.sympt == 0)),
                      intersect(txUCT_prep, which(dat$attr$uCT.sympt == 0)),
                      intersect(txsyph_prep, which(dat$attr$syph.sympt == 0)))
  dat$epi$num.asympt.tx.prep[at] <- length(unique(asympt.tx.prep))

  asympt.cases.prep <- c(intersect(idsRGC_prep_tx, which(dat$attr$rGC.sympt == 0)),
                         intersect(idsUGC_prep_tx, which(dat$attr$uGC.sympt == 0)),
                         intersect(idsRCT_prep_tx, which(dat$attr$rCT.sympt == 0)),
                         intersect(idsUCT_prep_tx, which(dat$attr$uCT.sympt == 0)),
                         intersect(idssyph_prep_tx, which(dat$attr$syph.sympt == 0)))
  dat$epi$num.asympt.cases.prep[at] <- length(unique(asympt.cases.prep))

  rect.tx <- c(txRGC_all, txRCT_all)
  dat$epi$num.rect.tx[at] <- length(unique(rect.tx))
  rect.cases <- c(idsRGC_tx, idsRGC_prep_tx, idsRCT_tx, idsRCT_prep_tx)
  dat$epi$num.rect.cases[at] <- length(unique(rect.cases))

  rect.tx.prep <- c(txRGC_prep, txRCT_prep)
  dat$epi$num.rect.tx.prep[at] <- length(unique(rect.tx.prep))
  rect.cases.prep <- c(idsRGC_prep_tx, idsRCT_prep_tx)
  dat$epi$num.rect.cases.prep[at] <- length(unique(rect.cases.prep))

  # Track total number treated
  dat$epi$txGC[at] <- length(unique(txRGC)) + length(unique(txUGC))
  dat$epi$txGC_asympt[at] <- length(unique(txGC_asympt))
  dat$epi$txCT[at] <- length(unique(txRCT)) + length(unique(txUCT))
  dat$epi$txCT_asympt[at] <- length(unique(txCT_asympt))
  dat$epi$txsyph[at] <- length(unique(c(txsyph)))
  dat$epi$txsyph_asympt[at] <- length(unique(txsyph_asympt))
  dat$epi$txearlysyph[at] <- length(unique(c(txsyph_sympt_prim, txsyph_sympt_seco,
                                             txsyph_asympt_prim, txsyph_asympt_seco,
                                             txsyph_asympt_earlat)))
  dat$epi$txlatesyph[at] <- length(unique(c(txsyph_asympt_latelat, txsyph_asympt_tert,
                                            txsyph_sympt_tert)))
  dat$epi$txSTI_asympt[at] <- dat$epi$txGC_asympt[at] + dat$epi$txCT_asympt[at] + dat$epi$txsyph_asympt[at]
  dat$epi$txSTI[at] <- dat$epi$txGC[at] + dat$epi$txCT[at] + dat$epi$txsyph[at]

  # Risk group-specific treatment and test counters
  dat$epi$txGC.tttraj1[at] <- length(which(tt.traj.gc.hivneg[unique(c(txRGC, txUGC))] == 1 | tt.traj.gc.hivpos[unique(c(txRGC, txUGC))] == 1))
  dat$epi$txGC_asympt.tttraj1[at] <- length(which(tt.traj.gc.hivneg[unique(txGC_asympt)] == 1 | tt.traj.gc.hivpos[unique(txGC_asympt)] == 1))
  dat$epi$txGC.tttraj2[at] <- length(which(tt.traj.gc.hivneg[unique(c(txRGC, txUGC))] == 2 | tt.traj.gc.hivpos[unique(c(txRGC, txUGC))] == 2))
  dat$epi$txGC_asympt.tttraj2[at] <- length(which(tt.traj.gc.hivneg[unique(txGC_asympt)] == 2 | tt.traj.gc.hivpos[unique(txGC_asympt)] == 2))

  dat$epi$txCT.tttraj1[at] <- length(which(tt.traj.ct.hivpos[unique(c(txRCT, txUCT))] == 1 | tt.traj.ct.hivneg[unique(c(txRCT, txUCT))] == 1))
  dat$epi$txCT_asympt.tttraj1[at] <- length(which(tt.traj.ct.hivpos[unique(txCT_asympt)] == 1 | tt.traj.ct.hivneg[unique(txCT_asympt)] == 1))
  dat$epi$txCT.tttraj2[at] <- length(which(tt.traj.ct.hivpos[unique(c(txRCT, txUCT))] == 2 | tt.traj.ct.hivneg[unique(c(txRCT, txUCT))] == 2))
  dat$epi$txCT_asympt.tttraj2[at] <- length(which(tt.traj.ct.hivpos[unique(txCT_asympt)] == 2 | tt.traj.ct.hivneg[unique(txCT_asympt)] == 2))

  dat$epi$txsyph.tttraj1[at] <- length(which(tt.traj.syph.hivneg[unique(txsyph)] == 1 | tt.traj.syph.hivpos[unique(txsyph)] == 1))
  dat$epi$txsyph_asympt.tttraj1[at] <- length(which(tt.traj.syph.hivneg[unique(txsyph_asympt)] == 1 | tt.traj.syph.hivpos[unique(txsyph_asympt)] == 1))
  dat$epi$txsyph.tttraj2[at] <- length(which(tt.traj.syph.hivneg[unique(txsyph)] == 2 | tt.traj.syph.hivpos[unique(txsyph)] == 2))
  dat$epi$txsyph_asympt.tttraj2[at] <- length(which(tt.traj.syph.hivneg[unique(txsyph_asympt)] == 2 | tt.traj.syph.hivpos[unique(txsyph_asympt)] == 2))

  dat$epi$txearlysyph.tttraj1[at] <- length(which(tt.traj.syph.hivneg[unique(c(txsyph_sympt_prim, txsyph_sympt_seco,
                                             txsyph_asympt_prim, txsyph_asympt_seco,
                                             txsyph_asympt_earlat))] == 1 | tt.traj.syph.hivpos[unique(c(txsyph_sympt_prim, txsyph_sympt_seco,
                                                                                                       txsyph_asympt_prim, txsyph_asympt_seco,
                                                                                                       txsyph_asympt_earlat))] == 1))
  dat$epi$txlatesyph.tttraj1[at] <- length(which(tt.traj.syph.hivneg[unique(c(txsyph_asympt_latelat, txsyph_asympt_tert,
                                                                            txsyph_sympt_tert))] == 1 |
                                                   tt.traj.syph.hivpos[unique(c(txsyph_asympt_latelat, txsyph_asympt_tert,
                                                                              txsyph_sympt_tert))] == 1))

  dat$epi$txearlysyph.tttraj2[at] <- length(which(tt.traj.syph.hivneg[unique(c(txsyph_sympt_prim, txsyph_sympt_seco,
                                                                             txsyph_asympt_prim, txsyph_asympt_seco,
                                                                             txsyph_asympt_earlat))] == 2 |
                                                    tt.traj.syph.hivpos[unique(c(txsyph_sympt_prim, txsyph_sympt_seco, txsyph_asympt_prim,
                                                                               txsyph_asympt_seco, txsyph_asympt_earlat))] == 2))
  dat$epi$txlatesyph.tttraj2[at] <- length(which(tt.traj.syph.hivneg[unique(c(txsyph_asympt_latelat, txsyph_asympt_tert,
                                                                            txsyph_sympt_tert))] == 2 |
                                                   tt.traj.syph.hivpos[unique(c(txsyph_asympt_latelat, txsyph_asympt_tert,
                                                                              txsyph_sympt_tert))] == 2))

  dat$epi$txSTI_asympt.tttraj1[at] <- dat$epi$txGC_asympt.tttraj1[at] + dat$epi$txCT_asympt.tttraj1[at] + dat$epi$txsyph_asympt.tttraj1[at]

  dat$epi$txSTI_asympt.tttraj2[at] <- dat$epi$txGC_asympt.tttraj2[at] + dat$epi$txCT_asympt.tttraj2[at] + dat$epi$txsyph_asympt.tttraj2[at]

  dat$epi$txSTI.tttraj1[at] <- dat$epi$txGC.tttraj1[at] + dat$epi$txCT.tttraj1[at] + dat$epi$txsyph.tttraj1[at]

  dat$epi$txSTI.tttraj2[at] <- dat$epi$txGC.tttraj2[at] + dat$epi$txCT.tttraj2[at] + dat$epi$txsyph.tttraj2[at]

  dat$epi$rGCsympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 1)) +
                                        length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 1))
  dat$epi$uGCsympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 1)) +
                                        length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 1))
  dat$epi$GCsympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 1)) +
                                      length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 1)) +
                                      length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 1)) +
                                      length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 1))

  dat$epi$rGCsympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 2)) +
                                        length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 2))
  dat$epi$uGCsympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 2)) +
                                        length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 2))
  dat$epi$GCsympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 2)) +
                                      length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 2)) +
                                      length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 2)) +
                                      length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 2))

  dat$epi$rCTsympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 1)) +
                                        length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 1))
  dat$epi$uCTsympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 1)) +
                                        length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 1))
  dat$epi$CTsympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 1)) +
                                      length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 1)) +
                                      length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 1)) +
                                      length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 1))

  dat$epi$rCTsympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 2)) +
                                        length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 2))
  dat$epi$uCTsympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 2)) +
                                        length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 2))
  dat$epi$CTsympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 2)) +
                                      length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 2)) +
                                      length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 2)) +
                                      length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 2))

  dat$epi$syphsympttests.tttraj1[at] <- length(which(tt.traj.syph.hivpos[unique(txsyph_sympt)] == 1)) +
                                        length(which(tt.traj.syph.hivneg[unique(txsyph_sympt)] == 1))

  dat$epi$syphsympttests.tttraj2[at] <- length(which(tt.traj.syph.hivpos[unique(txsyph_sympt)] == 2)) +
                                        length(which(tt.traj.syph.hivneg[unique(txsyph_sympt)] == 2))

  dat$epi$stisympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 1)) +
                                        length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 1)) +
                                        length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 1)) +
                                        length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 1)) +
                                        length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 1)) +
                                        length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 1)) +
                                        length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 1)) +
                                        length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 1)) +
                                        length(which(tt.traj.syph.hivpos[unique(txsyph_sympt)] == 1)) +
                                        length(which(tt.traj.syph.hivneg[unique(txsyph_sympt)] == 1))

  dat$epi$stisympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 2)) +
                                        length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 2)) +
                                        length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 2)) +
                                        length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 2)) +
                                        length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 2)) +
                                        length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 2)) +
                                        length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 2)) +
                                        length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 2)) +
                                        length(which(tt.traj.syph.hivpos[unique(txsyph_sympt)] == 2)) +
                                        length(which(tt.traj.syph.hivneg[unique(txsyph_sympt)] == 2))

  # EPT
  # Proportion of treated GC/CT index who have current partners - e.g. eligibility for EPT
  dat$epi$propindexeptElig[at] <- ifelse(length(unique(c(txRGC_all, txUGC_all, txRCT_all, txUCT_all))) > 0,
                                         length(unique(ept_tx_all)) /
                                           length(unique(c(txRGC_all, txUGC_all, txRCT_all, txUCT_all))),
                                         0)

  # Proportion of eligible index who will receive EPT - varies by sim scenario
  dat$epi$eptCov[at] <- eptCov

  # Number of non-index treated due to EPT - mix of provision, uptake, tx success
  dat$epi$eptTx[at] <- length(unique(alltxEPT))

  # Proportion of all non-index eligible to be treated who had treatment success
  dat$epi$eptprop_tx[at] <- ifelse(length(unique(allidsept)) > 0,
                                   length(unique(alltxEPT)) / length(unique(allidsept)),
                                   0)

  return(dat)
}
