
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
sti_trans <- function(dat, at) {

  # Parameters ----------------------------------------------------------

  # Acquisition probabilities given contact with infected man
  rgc.tprob <- dat$param$rgc.tprob
  ugc.tprob <- dat$param$ugc.tprob
  rct.tprob <- dat$param$rct.tprob
  uct.tprob <- dat$param$uct.tprob

  # Probability of symptoms given infection
  rgc.sympt.prob <- dat$param$rgc.sympt.prob
  ugc.sympt.prob <- dat$param$ugc.sympt.prob
  rct.sympt.prob <- dat$param$rct.sympt.prob
  uct.sympt.prob <- dat$param$uct.sympt.prob

  # Relative risk of infection given condom use during act
  sti.cond.eff <- dat$param$sti.cond.eff
  sti.cond.fail.B <- dat$param$sti.cond.fail.B
  sti.cond.fail.W <- dat$param$sti.cond.fail.W


  # Attributes ----------------------------------------------------------

  race <- dat$attr$race

  # Current infection state
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT

  # n Times infected
  rGC.timesInf <- dat$attr$rGC.timesInf
  uGC.timesInf <- dat$attr$uGC.timesInf
  rCT.timesInf <- dat$attr$rCT.timesInf
  uCT.timesInf <- dat$attr$uCT.timesInf

  # Set disease status to 0 for new births
  newBirths <- which(dat$attr$arrival.time == at)
  rGC[newBirths] <- rGC.timesInf[newBirths] <- 0
  uGC[newBirths] <- uGC.timesInf[newBirths] <- 0
  rCT[newBirths] <- rCT.timesInf[newBirths] <- 0
  uCT[newBirths] <- uCT.timesInf[newBirths] <- 0

  # Infection time
  rGC.infTime <- dat$attr$rGC.infTime
  uGC.infTime <- dat$attr$uGC.infTime
  rCT.infTime <- dat$attr$rCT.infTime
  uCT.infTime <- dat$attr$uCT.infTime

  # Infection symptoms (non-varying)
  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  rCT.sympt <- dat$attr$rCT.sympt
  uCT.sympt <- dat$attr$uCT.sympt


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
  uai_rgc <- al[allActs_rgc, "uai"]
  tprob_rgc <- rep(rgc.tprob, length(allActs_rgc))

  # Transform to log odds
  tlo_rgc <- log(tprob_rgc/(1 - tprob_rgc))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_rgc, "p1"]], race[al[p2Inf_rgc, "p2"]])
  condom.rr <- rep(NA, length(races))
  condom.rr[races == "B"] <- 1 - (sti.cond.eff - sti.cond.fail.B)
  condom.rr[races == "W"] <- 1 - (sti.cond.eff - sti.cond.fail.W)

  tlo_rgc[uai_rgc == 0] <- tlo_rgc[uai_rgc == 0] + log(condom.rr[uai_rgc == 0])

  # Back-transform to probability
  tprob_rgc <- plogis(tlo_rgc)

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
  rGC.infTime[idsInf_rgc] <- at
  rGC.sympt[idsInf_rgc] <- rbinom(length(idsInf_rgc), 1, rgc.sympt.prob)
  rGC.timesInf[idsInf_rgc] <- rGC.timesInf[idsInf_rgc] + 1


  # Urethral GC ---------------------------------------------------------

  # Requires: rGC in receptive man, and no uGC in insertive man
  p1Inf_ugc <- which(rGC[al[, "p1"]] == 1 & rGC.infTime[al[, "p1"]] < at &
                     uGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf_ugc <- which(rGC[al[, "p2"]] == 1 & rGC.infTime[al[, "p2"]] < at &
                     uGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs_ugc <- c(p1Inf_ugc, p2Inf_ugc)

  # UAI modifier
  uai_ugc <- al[allActs_ugc, "uai"]
  tprob_ugc <- rep(ugc.tprob, length(allActs_ugc))

  # Transform to log odds
  tlo_ugc <- log(tprob_ugc/(1 - tprob_ugc))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_ugc, "p2"]], race[al[p2Inf_ugc, "p1"]])
  condom.rr <- rep(NA, length(races))
  condom.rr[races == "B"] <- 1 - (sti.cond.eff - sti.cond.fail.B)
  condom.rr[races == "W"] <- 1 - (sti.cond.eff - sti.cond.fail.W)

  tlo_ugc[uai_ugc == 0] <- tlo_ugc[uai_ugc == 0] + log(condom.rr[uai_ugc == 0])

  # Back-transform to probability
  tprob_ugc <- plogis(tlo_ugc)

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
  uGC.infTime[idsInf_ugc] <- at
  uGC.sympt[idsInf_ugc] <- rbinom(length(idsInf_ugc), 1, ugc.sympt.prob)
  uGC.timesInf[idsInf_ugc] <- uGC.timesInf[idsInf_ugc] + 1


  # Rectal CT -----------------------------------------------------------

  # Requires: uCT in insertive man, and no rCT in receptive man
  p1Inf_rct <- which(uCT[al[, "p1"]] == 1 & uCT.infTime[al[, "p1"]] < at &
                     rCT[al[, "p2"]] == 0 & al[, "ins"] %in% c(1, 2))
  p2Inf_rct <- which(uCT[al[, "p2"]] == 1 & uCT.infTime[al[, "p2"]] < at &
                     rCT[al[, "p1"]] == 0 & al[, "ins"] %in% c(0, 2))
  allActs_rct <- c(p1Inf_rct, p2Inf_rct)

  # UAI modifier
  uai_rct <- al[allActs_rct, "uai"]
  tprob_rct <- rep(rct.tprob, length(allActs_rct))

  # Transform to log odds
  tlo_rct <- log(tprob_rct/(1 - tprob_rct))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_rct, "p1"]], race[al[p2Inf_rct, "p2"]])
  condom.rr <- rep(NA, length(races))
  condom.rr[races == "B"] <- 1 - (sti.cond.eff - sti.cond.fail.B)
  condom.rr[races == "W"] <- 1 - (sti.cond.eff - sti.cond.fail.W)

  tlo_rct[uai_rct == 0] <- tlo_rct[uai_rct == 0] + log(condom.rr[uai_rct == 0])

  # Back-transform to probability
  tprob_rct <- plogis(tlo_rct)

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
  rCT.infTime[idsInf_rct] <- at
  rCT.sympt[idsInf_rct] <- rbinom(length(idsInf_rct), 1, rct.sympt.prob)
  rCT.timesInf[idsInf_rct] <- rCT.timesInf[idsInf_rct] + 1


  # Urethral CT ---------------------------------------------------------

  # Requires: rCT in receptive man, and no uCT in insertive man
  p1Inf_uct <- which(rCT[al[, "p1"]] == 1 & rCT.infTime[al[, "p1"]] < at &
                     uCT[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf_uct <- which(rCT[al[, "p2"]] == 1 & rCT.infTime[al[, "p2"]] < at &
                     uCT[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs_uct <- c(p1Inf_uct, p2Inf_uct)

  # UAI modifier
  uai_uct <- al[allActs_uct, "uai"]
  tprob_uct <- rep(uct.tprob, length(allActs_uct))

  # Transform to log odds
  tlo_uct <- log(tprob_uct/(1 - tprob_uct))

  # Modify log odds by race-specific condom effectiveness
  races <- c(race[al[p1Inf_uct, "p2"]], race[al[p2Inf_uct, "p1"]])
  condom.rr <- rep(NA, length(races))
  condom.rr[races == "B"] <- 1 - (sti.cond.eff - sti.cond.fail.B)
  condom.rr[races == "W"] <- 1 - (sti.cond.eff - sti.cond.fail.W)

  tlo_uct[uai_uct == 0] <- tlo_uct[uai_uct == 0] + log(condom.rr[uai_uct == 0])

  # Back-transform to probability
  tprob_uct <- plogis(tlo_uct)

  # Stochastic transmission
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
  uCT.infTime[idsInf_uct] <- at
  uCT.sympt[idsInf_uct] <- rbinom(length(idsInf_uct), 1, uct.sympt.prob)
  uCT.timesInf[idsInf_uct] <- uCT.timesInf[idsInf_uct] + 1


  # Output --------------------------------------------------------------

  # attributes
  dat$attr$rGC <- rGC
  dat$attr$uGC <- uGC
  dat$attr$rCT <- rCT
  dat$attr$uCT <- uCT

  dat$attr$rGC.infTime <- rGC.infTime
  dat$attr$uGC.infTime <- uGC.infTime
  dat$attr$rCT.infTime <- rCT.infTime
  dat$attr$uCT.infTime <- uCT.infTime

  dat$attr$rGC.timesInf <- rGC.timesInf
  dat$attr$uGC.timesInf <- uGC.timesInf
  dat$attr$rCT.timesInf <- rCT.timesInf
  dat$attr$uCT.timesInf <- uCT.timesInf

  dat$attr$rGC.sympt <- rGC.sympt
  dat$attr$uGC.sympt <- uGC.sympt
  dat$attr$rCT.sympt <- rCT.sympt
  dat$attr$uCT.sympt <- uCT.sympt


  # Summary stats
  dat$epi$incid.gc[at] <- length(idsInf_rgc) + length(idsInf_ugc)
  dat$epi$incid.gc.B[at] <- length(intersect(union(idsInf_rgc, idsInf_ugc), which(race == "B")))
  dat$epi$incid.gc.W[at] <- length(intersect(union(idsInf_rgc, idsInf_ugc), which(race == "W")))

  dat$epi$incid.ct[at] <- length(idsInf_rct) + length(idsInf_uct)
  dat$epi$incid.ct.B[at] <- length(intersect(union(idsInf_rct, idsInf_uct), which(race == "B")))
  dat$epi$incid.ct.W[at] <- length(intersect(union(idsInf_rct, idsInf_uct), which(race == "W")))


  # Check all infected have all STI attributes
  stopifnot(all(!is.na(rGC.infTime[rGC == 1])),
            all(!is.na(rGC.sympt[rGC == 1])),
            all(!is.na(uGC.infTime[uGC == 1])),
            all(!is.na(uGC.sympt[uGC == 1])),
            all(!is.na(rCT.infTime[rCT == 1])),
            all(!is.na(rCT.sympt[rCT == 1])),
            all(!is.na(uCT.infTime[uCT == 1])),
            all(!is.na(uCT.sympt[uCT == 1])))

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
sti_recov <- function(dat, at) {

  # Parameters ----------------------------------------------------------

  rgc.ntx.int <- dat$param$rgc.ntx.int
  ugc.ntx.int <- dat$param$ugc.ntx.int
  gc.tx.int <- dat$param$gc.tx.int

  rct.ntx.int <- dat$param$rct.ntx.int
  uct.ntx.int <- dat$param$uct.ntx.int
  ct.tx.int <- dat$param$ct.tx.int


  # GC Recovery ---------------------------------------------------------

  # Untreated (asymptomatic and symptomatic)
  idsRGC_ntx <- which(dat$attr$rGC == 1 &
                      dat$attr$rGC.infTime < at &
                      (is.na(dat$attr$rGC.tx) | dat$attr$rGC.tx == 0) &
                      (is.na(dat$attr$rGC.tx.prep) | dat$attr$rGC.tx.prep == 0))
  idsUGC_ntx <- which(dat$attr$uGC == 1 &
                             dat$attr$uGC.infTime < at &
                             (is.na(dat$attr$uGC.tx) | dat$attr$uGC.tx == 0) &
                             (is.na(dat$attr$uGC.tx.prep) | dat$attr$uGC.tx.prep == 0))

  recovRGC_ntx <- idsRGC_ntx[which(rbinom(length(idsRGC_ntx), 1,
                                          1/rgc.ntx.int) == 1)]
  recovUGC_ntx <- idsUGC_ntx[which(rbinom(length(idsUGC_ntx), 1,
                                          1/ugc.ntx.int) == 1)]


  # Treated (asymptomatic and symptomatic)
  idsRGC_tx <- which(dat$attr$rGC == 1 &
                     dat$attr$rGC.infTime < at &
                     (dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1))
  idsUGC_tx <- which(dat$attr$uGC == 1 &
                     dat$attr$uGC.infTime < at &
                     (dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1))

  recovRGC_tx <- idsRGC_tx[which(rbinom(length(idsRGC_tx), 1,
                                        1/gc.tx.int) == 1)]
  recovUGC_tx <- idsUGC_tx[which(rbinom(length(idsUGC_tx), 1,
                                        1/gc.tx.int) == 1)]

  recovRGC <- c(recovRGC_ntx, recovRGC_tx)
  recovUGC <- c(recovUGC_ntx, recovUGC_tx)

  dat$attr$rGC[recovRGC] <- 0
  dat$attr$rGC.sympt[recovRGC] <- NA
  dat$attr$rGC.infTime[recovRGC] <- NA
  dat$attr$rGC.tx[recovRGC] <- NA
  dat$attr$rGC.tx.prep[recovRGC] <- NA

  dat$attr$uGC[recovUGC] <- 0
  dat$attr$uGC.sympt[recovUGC] <- NA
  dat$attr$uGC.infTime[recovUGC] <- NA
  dat$attr$uGC.tx[recovUGC] <- NA
  dat$attr$uGC.tx.prep[recovUGC] <- NA



  # CT Recovery ---------------------------------------------------------

  # Untreated (asymptomatic and symptomatic)
  idsRCT_ntx <- which(dat$attr$rCT == 1 &
                      dat$attr$rCT.infTime < at &
                      (is.na(dat$attr$rCT.tx) | dat$attr$rCT.tx == 0) &
                      (is.na(dat$attr$rCT.tx.prep) | dat$attr$rCT.tx.prep == 0))
  idsUCT_ntx <- which(dat$attr$uCT == 1 &
                      dat$attr$uCT.infTime < at &
                      (is.na(dat$attr$uCT.tx) | dat$attr$uCT.tx == 0) &
                      (is.na(dat$attr$uCT.tx.prep) | dat$attr$uCT.tx.prep == 0))

  recovRCT_ntx <- idsRCT_ntx[which(rbinom(length(idsRCT_ntx),
                                          1, 1/rct.ntx.int) == 1)]
  recovUCT_ntx <- idsUCT_ntx[which(rbinom(length(idsUCT_ntx),
                                          1, 1/uct.ntx.int) == 1)]


  # Treated (asymptomatic and symptomatic)
  idsRCT_tx <- which(dat$attr$rCT == 1 &
                     dat$attr$rCT.infTime < at &
                     (dat$attr$rCT.tx == 1 | dat$attr$rCT.tx.prep == 1))
  idsUCT_tx <- which(dat$attr$uCT == 1 &
                     dat$attr$uCT.infTime < at &
                     (dat$attr$uCT.tx == 1 | dat$attr$uCT.tx.prep == 1))

  recovRCT_tx <- idsRCT_tx[which(rbinom(length(idsRCT_tx),
                                        1, 1/ct.tx.int) == 1)]
  recovUCT_tx <- idsUCT_tx[which(rbinom(length(idsUCT_tx),
                                        1, 1/ct.tx.int) == 1)]


  recovRCT <- c(recovRCT_ntx, recovRCT_tx)
  recovUCT <- c(recovUCT_ntx, recovUCT_tx)

  dat$attr$rCT[recovRCT] <- 0
  dat$attr$rCT.sympt[recovRCT] <- NA
  dat$attr$rCT.infTime[recovRCT] <- NA
  dat$attr$rCT.tx[recovRCT] <- NA
  dat$attr$rCT.tx.prep[recovRCT] <- NA

  dat$attr$uCT[recovUCT] <- 0
  dat$attr$uCT.sympt[recovUCT] <- NA
  dat$attr$uCT.infTime[recovUCT] <- NA
  dat$attr$uCT.tx[recovUCT] <- NA
  dat$attr$uCT.tx.prep[recovUCT] <- NA

  return(dat)
}


#' @title STI Treatment Module
#'
#' @description Stochastically simulates GC/CT diagnosis and treatment.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
sti_tx <- function(dat, at) {

  # Parameters
  gc.sympt.prob.tx.B <- dat$param$gc.sympt.prob.tx.B
  gc.sympt.prob.tx.W <- dat$param$gc.sympt.prob.tx.W
  ct.sympt.prob.tx.W <- dat$param$ct.sympt.prob.tx.B
  ct.sympt.prob.tx.B <- dat$param$ct.sympt.prob.tx.W

  gc.asympt.prob.tx.B <- dat$param$gc.asympt.prob.tx.B
  gc.asympt.prob.tx.W <- dat$param$gc.asympt.prob.tx.W
  ct.asympt.prob.tx.W <- dat$param$ct.asympt.prob.tx.B
  ct.asympt.prob.tx.B <- dat$param$ct.asympt.prob.tx.W

  prep.sti.screen.int <- dat$param$prep.sti.screen.int
  prep.sti.prob.tx <- dat$param$prep.sti.prob.tx

  prep.cont.stand.tx <- dat$param$prep.continue.stand.tx
  if (prep.cont.stand.tx == TRUE) {
    prep.stand.tx.grp <- 0:1
  } else {
    prep.stand.tx.grp <- 0
  }

  race <- dat$attr$race

  ## Symptomatic GC Treatment ##
  idsRGC_tx_sympt <- which(dat$attr$rGC == 1 &
                           dat$attr$rGC.infTime < at &
                           dat$attr$rGC.sympt == 1 &
                           is.na(dat$attr$rGC.tx) &
                           dat$attr$prepStat %in% prep.stand.tx.grp)
  idsUGC_tx_sympt <- which(dat$attr$uGC == 1 &
                           dat$attr$uGC.infTime < at &
                           dat$attr$uGC.sympt == 1 &
                           is.na(dat$attr$uGC.tx) &
                           dat$attr$prepStat %in% prep.stand.tx.grp)

  # Subset by race
  idsRGC_tx_sympt_B <- intersect(idsRGC_tx_sympt, which(race == "B"))
  idsRGC_tx_sympt_W <- intersect(idsRGC_tx_sympt, which(race == "W"))
  idsUGC_tx_sympt_B <- intersect(idsUGC_tx_sympt, which(race == "B"))
  idsUGC_tx_sympt_W <- intersect(idsUGC_tx_sympt, which(race == "W"))

  # Collect over site
  idsGC_tx_sympt_B <- union(idsRGC_tx_sympt_B, idsUGC_tx_sympt_B)
  idsGC_tx_sympt_W <- union(idsRGC_tx_sympt_W, idsUGC_tx_sympt_W)

  # Treatment by race
  txGC_sympt_B <- idsGC_tx_sympt_B[which(rbinom(length(idsGC_tx_sympt_B), 1,
                                                gc.sympt.prob.tx.B) == 1)]
  txGC_sympt_W <- idsGC_tx_sympt_W[which(rbinom(length(idsGC_tx_sympt_W), 1,
                                                gc.sympt.prob.tx.W) == 1)]
  txGC_sympt <- union(txGC_sympt_B, txGC_sympt_W)

  # Subset by site
  txRGC_sympt <- intersect(idsRGC_tx_sympt, txGC_sympt)
  txUGC_sympt <- intersect(idsUGC_tx_sympt, txGC_sympt)

  ## Asymptomatic GC Treatment ##
  idsRGC_tx_asympt <- which(dat$attr$rGC == 1 &
                            dat$attr$rGC.infTime < at &
                            dat$attr$rGC.sympt == 0 &
                            is.na(dat$attr$rGC.tx) &
                            dat$attr$prepStat %in% prep.stand.tx.grp)
  idsUGC_tx_asympt <- which(dat$attr$uGC == 1 &
                            dat$attr$uGC.infTime < at &
                            dat$attr$uGC.sympt == 0 &
                            is.na(dat$attr$uGC.tx) &
                            dat$attr$prepStat %in% prep.stand.tx.grp)

  # Subset by race
  idsRGC_tx_asympt_B <- intersect(idsRGC_tx_asympt, which(race == "B"))
  idsRGC_tx_asympt_W <- intersect(idsRGC_tx_asympt, which(race == "W"))
  idsUGC_tx_asympt_B <- intersect(idsUGC_tx_asympt, which(race == "B"))
  idsUGC_tx_asympt_W <- intersect(idsUGC_tx_asympt, which(race == "W"))

  # Collect over site
  idsGC_tx_asympt_B <- union(idsRGC_tx_asympt_B, idsUGC_tx_asympt_B)
  idsGC_tx_asympt_W <- union(idsRGC_tx_asympt_W, idsUGC_tx_asympt_W)

  # Treatment by race
  txGC_asympt_B <- idsGC_tx_asympt_B[which(rbinom(length(idsGC_tx_asympt_B), 1,
                                                  gc.asympt.prob.tx.B) == 1)]
  txGC_asympt_W <- idsGC_tx_asympt_W[which(rbinom(length(idsGC_tx_asympt_W), 1,
                                                  gc.asympt.prob.tx.W) == 1)]
  txGC_asympt <- union(txGC_asympt_B, txGC_asympt_W)

  # Subset by site
  txRGC_asympt <- intersect(idsRGC_tx_asympt, txGC_asympt)
  txUGC_asympt <- intersect(idsUGC_tx_asympt, txGC_asympt)

  ## All Treated GC ##

  # IDs of men sucessfully treated
  txRGC <- union(txRGC_sympt, txRGC_asympt)
  txUGC <- union(txUGC_sympt, txUGC_asympt)

  # IDs of men eligible for treatment
  idsRGC_tx <- union(idsRGC_tx_sympt, idsRGC_tx_asympt)
  idsUGC_tx <- union(idsUGC_tx_sympt, idsUGC_tx_asympt)


  ## Symptomatic CT Treatment ##
  idsRCT_tx_sympt <- which(dat$attr$rCT == 1 &
                           dat$attr$rCT.infTime < at &
                           dat$attr$rCT.sympt == 1 &
                           is.na(dat$attr$rCT.tx) &
                           dat$attr$prepStat %in% prep.stand.tx.grp)
  idsUCT_tx_sympt <- which(dat$attr$uCT == 1 &
                           dat$attr$uCT.infTime < at &
                           dat$attr$uCT.sympt == 1 &
                           is.na(dat$attr$uCT.tx) &
                           dat$attr$prepStat %in% prep.stand.tx.grp)

  # Subset by race
  idsRCT_tx_sympt_B <- intersect(idsRCT_tx_sympt, which(race == "B"))
  idsRCT_tx_sympt_W <- intersect(idsRCT_tx_sympt, which(race == "W"))
  idsUCT_tx_sympt_B <- intersect(idsUCT_tx_sympt, which(race == "B"))
  idsUCT_tx_sympt_W <- intersect(idsUCT_tx_sympt, which(race == "W"))

  # Collect over site
  idsCT_tx_sympt_B <- union(idsRCT_tx_sympt_B, idsUCT_tx_sympt_B)
  idsCT_tx_sympt_W <- union(idsRCT_tx_sympt_W, idsUCT_tx_sympt_W)

  # Treatment by race
  txCT_sympt_B <- idsCT_tx_sympt_B[which(rbinom(length(idsCT_tx_sympt_B), 1,
                                                ct.sympt.prob.tx.B) == 1)]
  txCT_sympt_W <- idsCT_tx_sympt_W[which(rbinom(length(idsCT_tx_sympt_W), 1,
                                                ct.sympt.prob.tx.W) == 1)]
  txCT_sympt <- union(txCT_sympt_B, txCT_sympt_W)

  # Subset by site
  txRCT_sympt <- intersect(idsRCT_tx_sympt, txCT_sympt)
  txUCT_sympt <- intersect(idsUCT_tx_sympt, txCT_sympt)


  ## Asymptomatic CT Treatment ##
  idsRCT_tx_asympt <- which(dat$attr$rCT == 1 &
                            dat$attr$rCT.infTime < at &
                            dat$attr$rCT.sympt == 0 &
                            is.na(dat$attr$rCT.tx) &
                            dat$attr$prepStat == 0)
  idsUCT_tx_asympt <- which(dat$attr$uCT == 1 &
                            dat$attr$uCT.infTime < at &
                            dat$attr$uCT.sympt == 0 &
                            is.na(dat$attr$uCT.tx) &
                            dat$attr$prepStat == 0)

  # Subset by race
  idsRCT_tx_asympt_B <- intersect(idsRCT_tx_asympt, which(race == "B"))
  idsRCT_tx_asympt_W <- intersect(idsRCT_tx_asympt, which(race == "W"))
  idsUCT_tx_asympt_B <- intersect(idsUCT_tx_asympt, which(race == "B"))
  idsUCT_tx_asympt_W <- intersect(idsUCT_tx_asympt, which(race == "W"))

  # Collect over site
  idsCT_tx_asympt_B <- union(idsRCT_tx_asympt_B, idsUCT_tx_asympt_B)
  idsCT_tx_asympt_W <- union(idsRCT_tx_asympt_W, idsUCT_tx_asympt_W)

  # Treatment by race
  txCT_asympt_B <- idsCT_tx_asympt_B[which(rbinom(length(idsCT_tx_asympt_B), 1,
                                                ct.asympt.prob.tx.B) == 1)]
  txCT_asympt_W <- idsCT_tx_asympt_W[which(rbinom(length(idsCT_tx_asympt_W), 1,
                                                ct.asympt.prob.tx.W) == 1)]
  txCT_asympt <- union(txCT_asympt_B, txCT_asympt_W)

  # Subset by site
  txRCT_asympt <- intersect(idsRCT_tx_asympt, txCT_asympt)
  txUCT_asympt <- intersect(idsUCT_tx_asympt, txCT_asympt)

  ## All Treated CT ##
  txRCT <- union(txRCT_sympt, txRCT_asympt)
  txUCT <- union(txUCT_sympt, txUCT_asympt)

  idsRCT_tx <- union(idsRCT_tx_sympt, idsRCT_tx_asympt)
  idsUCT_tx <- union(idsUCT_tx_sympt, idsUCT_tx_asympt)


  ## Interval-based treatment for MSM on PrEP ##
  idsSTI_screen <- which(dat$attr$prepStartTime == at |
                           (at - dat$attr$prepLastStiScreen >= prep.sti.screen.int))

  dat$attr$prepLastStiScreen[idsSTI_screen] <- at


  idsRGC_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$rGC == 1 &
                                      dat$attr$rGC.infTime < at &
                                      is.na(dat$attr$rGC.tx.prep)))
  idsUGC_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$uGC == 1 &
                                      dat$attr$uGC.infTime < at &
                                      is.na(dat$attr$uGC.tx.prep)))
  idsRCT_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$rCT == 1 &
                                      dat$attr$rCT.infTime < at &
                                      is.na(dat$attr$rCT.tx.prep)))
  idsUCT_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$uCT == 1 &
                                      dat$attr$uCT.infTime < at &
                                      is.na(dat$attr$uCT.tx.prep)))

  txRGC_prep <- idsRGC_prep_tx[which(rbinom(length(idsRGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txUGC_prep <- idsUGC_prep_tx[which(rbinom(length(idsUGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txRCT_prep <- idsRCT_prep_tx[which(rbinom(length(idsRCT_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txUCT_prep <- idsUCT_prep_tx[which(rbinom(length(idsUCT_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]


  ## Update Attributes ##
  dat$attr$rGC.tx[idsRGC_tx] <- 0
  dat$attr$rGC.tx[txRGC] <- 1

  dat$attr$uGC.tx[idsUGC_tx] <- 0
  dat$attr$uGC.tx[txUGC] <- 1

  dat$attr$rCT.tx[idsRCT_tx] <- 0
  dat$attr$rCT.tx[txRCT] <- 1

  dat$attr$uCT.tx[idsUCT_tx] <- 0
  dat$attr$uCT.tx[txUCT] <- 1

  dat$attr$rGC.tx.prep[idsRGC_prep_tx] <- 0
  dat$attr$rGC.tx.prep[txRGC_prep] <- 1

  dat$attr$uGC.tx.prep[idsUGC_prep_tx] <- 0
  dat$attr$uGC.tx.prep[txUGC_prep] <- 1

  dat$attr$rCT.tx.prep[idsRCT_prep_tx] <- 0
  dat$attr$rCT.tx.prep[txRCT_prep] <- 1

  dat$attr$uCT.tx.prep[idsUCT_prep_tx] <- 0
  dat$attr$uCT.tx.prep[txUCT_prep] <- 1


  ## Add tx at other anatomical site ##
  dat$attr$rGC.tx[which((dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1) &
                          dat$attr$rGC == 1)] <- 1
  dat$attr$uGC.tx[which((dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1) &
                          dat$attr$uGC == 1)] <- 1

  dat$attr$rCT.tx[which((dat$attr$uCT.tx == 1 | dat$attr$uCT.tx.prep == 1) &
                          dat$attr$rCT == 1)] <- 1
  dat$attr$uCT.tx[which((dat$attr$rCT.tx == 1 | dat$attr$rCT.tx.prep == 1) &
                          dat$attr$uCT == 1)] <- 1

  return(dat)
}
