
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
stitrans_msm <- function(dat, at) {

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
  sti.cond.fail <- dat$param$sti.cond.fail


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
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.rr[ids.race] <- 1 - (sti.cond.eff - sti.cond.fail[i])
  }

  tlo_rgc[uai_rgc == 0] <- tlo_rgc[uai_rgc == 0] + log(condom.rr[uai_rgc == 0])

  # Back-transform to probability
  tprob_rgc <- plogis(tlo_rgc)

  # Stochastic transmission
  trans_rgc <- rbinom(length(allActs_rgc), 1, tprob_rgc)

  # Determine the infected partner
  idsInf_rgc <- NULL
  if (sum(trans_rgc) > 0) {
    transAL_rgc <- al[allActs_rgc[trans_rgc == 1], , drop = FALSE]
    idsInf_rgc <- c(intersect(al[p1Inf_rgc, "p2"], transAL_rgc[, "p2"]),
                    intersect(al[p2Inf_rgc, "p1"], transAL_rgc[, "p1"]))
    stopifnot(all(rGC[idsInf_rgc] == 0))
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
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.rr[ids.race] <- 1 - (sti.cond.eff - sti.cond.fail[i])
  }

  tlo_ugc[uai_ugc == 0] <- tlo_ugc[uai_ugc == 0] + log(condom.rr[uai_ugc == 0])

  # Back-transform to probability
  tprob_ugc <- plogis(tlo_ugc)

  # Stochastic transmission
  trans_ugc <- rbinom(length(allActs_ugc), 1, tprob_ugc)

  # Determine the newly infected partner
  idsInf_ugc <- NULL
  if (sum(trans_ugc) > 0) {
    transAL_ugc <- al[allActs_ugc[trans_ugc == 1],  , drop = FALSE]
    idsInf_ugc <- c(intersect(al[p1Inf_ugc, "p2"], transAL_ugc[, "p2"]),
                    intersect(al[p2Inf_ugc, "p1"], transAL_ugc[, "p1"]))
    stopifnot(all(uGC[idsInf_ugc] == 0))
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
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.rr[ids.race] <- 1 - (sti.cond.eff - sti.cond.fail[i])
  }

  tlo_rct[uai_rct == 0] <- tlo_rct[uai_rct == 0] + log(condom.rr[uai_rct == 0])

  # Back-transform to probability
  tprob_rct <- plogis(tlo_rct)

  # Stochastic transmission
  trans_rct <- rbinom(length(allActs_rct), 1, tprob_rct)

  # Determine the newly infected partner
  idsInf_rct <- NULL
  if (sum(trans_rct) > 0) {
    transAL_rct <- al[allActs_rct[trans_rct == 1],  , drop = FALSE]
    idsInf_rct <- c(intersect(al[p1Inf_rct, "p2"], transAL_rct[, "p2"]),
                    intersect(al[p2Inf_rct, "p1"], transAL_rct[, "p1"]))
    stopifnot(all(rCT[idsInf_rct] == 0))
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
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.rr[ids.race] <- 1 - (sti.cond.eff - sti.cond.fail[i])
  }

  tlo_uct[uai_uct == 0] <- tlo_uct[uai_uct == 0] + log(condom.rr[uai_uct == 0])

  # Back-transform to probability
  tprob_uct <- plogis(tlo_uct)

  # Stochastic transmission
  trans_uct <- rbinom(length(allActs_uct), 1, tprob_uct)

  # Determine the newly infected partner
  idsInf_uct <- NULL
  if (sum(trans_uct) > 0) {
    transAL_uct <- al[allActs_uct[trans_uct == 1],  , drop = FALSE]
    idsInf_uct <- c(intersect(al[p1Inf_uct, "p2"], transAL_uct[, "p2"]),
                    intersect(al[p2Inf_uct, "p1"], transAL_uct[, "p1"]))
    stopifnot(all(uCT[idsInf_uct] == 0))
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
  dat$epi$incid.rgc[at] <- length(idsInf_rgc)
  dat$epi$incid.ugc[at] <- length(idsInf_ugc)
  dat$epi$incid.gc.B[at] <- length(intersect(union(idsInf_rgc, idsInf_ugc), which(race == 1)))
  dat$epi$incid.gc.H[at] <- length(intersect(union(idsInf_rgc, idsInf_ugc), which(race == 2)))
  dat$epi$incid.gc.W[at] <- length(intersect(union(idsInf_rgc, idsInf_ugc), which(race == 3)))

  dat$epi$incid.ct[at] <- length(idsInf_rct) + length(idsInf_uct)
  dat$epi$incid.rct[at] <- length(idsInf_rct)
  dat$epi$incid.uct[at] <- length(idsInf_uct)
  dat$epi$incid.ct.B[at] <- length(intersect(union(idsInf_rct, idsInf_uct), which(race == 1)))
  dat$epi$incid.ct.B[at] <- length(intersect(union(idsInf_rct, idsInf_uct), which(race == 2)))
  dat$epi$incid.ct.W[at] <- length(intersect(union(idsInf_rct, idsInf_uct), which(race == 3)))

  dat$epi$incid.rct.B[at] <- length(intersect(idsInf_rct, which(race == 1)))
  dat$epi$incid.rct.H[at] <- length(intersect(idsInf_rct, which(race == 2)))
  dat$epi$incid.rct.W[at] <- length(intersect(idsInf_rct, which(race == 3)))
  dat$epi$incid.uct.B[at] <- length(intersect(idsInf_uct, which(race == 1)))
  dat$epi$incid.uct.H[at] <- length(intersect(idsInf_uct, which(race == 2)))
  dat$epi$incid.uct.W[at] <- length(intersect(idsInf_uct, which(race == 3)))

  dat$epi$incid.rgc.B[at] <- length(intersect(idsInf_rgc, which(race == 1)))
  dat$epi$incid.rgc.H[at] <- length(intersect(idsInf_rgc, which(race == 2)))
  dat$epi$incid.rgc.W[at] <- length(intersect(idsInf_rgc, which(race == 3)))
  dat$epi$incid.ugc.B[at] <- length(intersect(idsInf_ugc, which(race == 1)))
  dat$epi$incid.ugc.H[at] <- length(intersect(idsInf_ugc, which(race == 2)))
  dat$epi$incid.ugc.W[at] <- length(intersect(idsInf_ugc, which(race == 3)))

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
