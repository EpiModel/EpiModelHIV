
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
  rsyph.tprob <- dat$param$rsyph.tprob
  usyph.tprob <- dat$param$usyph.tprob
  
  #Multiplier for syphilis infection
  syph.earlat.rr <- dat$param$syph.earlat.rr
  syph.late.rr <- dat$param$syph.late.rr
  syph.rhiv.rr <- dat$param$syph.rhiv.rr
  syph.uhiv.rr <- dat$param$syph.uhiv.rr
  
  # Probability of symptoms given infection
  rgc.sympt.prob <- dat$param$rgc.sympt.prob
  ugc.sympt.prob <- dat$param$ugc.sympt.prob
  rct.sympt.prob <- dat$param$rct.sympt.prob
  uct.sympt.prob <- dat$param$uct.sympt.prob
  
  # Relative risk of infection given condom use during act
  sti.cond.rr <- dat$param$sti.cond.rr

  # Cessation
  gc.prob.cease <- dat$param$gc.prob.cease
  ct.prob.cease <- dat$param$ct.prob.cease
  syph.prob.cease <- dat$param$syph.prob.cease

  # Attributes ----------------------------------------------------------

  # Current infection state
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT
  syphilis <- dat$attr$syphilis
  status <- dat$attr$status
  stage.syph <- dat$attr$stage.syph

  # n Times infected
  rGC.timesInf <- dat$attr$rGC.timesInf
  uGC.timesInf <- dat$attr$uGC.timesInf
  rCT.timesInf <- dat$attr$rCT.timesInf
  uCT.timesInf <- dat$attr$uCT.timesInf
  syph.timesInf <- dat$attr$syph.timesInf

  # Set disease status to 0 for new births
  newBirths <- which(dat$attr$arrival.time == at)
  rGC[newBirths] <- rGC.timesInf[newBirths] <- 0
  uGC[newBirths] <- uGC.timesInf[newBirths] <- 0
  rCT[newBirths] <- rCT.timesInf[newBirths] <- 0
  uCT[newBirths] <- uCT.timesInf[newBirths] <- 0
  syphilis[newBirths] <- syph.timesInf[newBirths] <- 0

  # Infection time
  rGC.infTime <- dat$attr$rGC.infTime
  uGC.infTime <- dat$attr$uGC.infTime
  rCT.infTime <- dat$attr$rCT.infTime
  uCT.infTime <- dat$attr$uCT.infTime
  syph.infTime <- dat$attr$syph.infTime
  rGC.lastinfTime <- dat$attr$rGC.lastinfTime
  uGC.lastinfTime <- dat$attr$uGC.lastinfTime
  rCT.lastinfTime <- dat$attr$rCT.lastinfTime
  uCT.lastinfTime <- dat$attr$uCT.lastinfTime
  syph.lastinfTime <- dat$attr$syph.lastinfTime

  # Infection symptoms (non-varying)
  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  rCT.sympt <- dat$attr$rCT.sympt
  uCT.sympt <- dat$attr$uCT.sympt
  stage.prim.sympt <- dat$attr$stage.prim.sympt
  stage.seco.sympt <- dat$attr$stage.seco.sympt
  stage.earlat.sympt <- dat$attr$stage.earlat.sympt
  stage.latelat.sympt <- dat$attr$stage.latelat.sympt
  stage.latelatelat.sympt <- dat$attr$stage.latelatelat.sympt  
  stage.tert.sympt <- dat$attr$stage.tert.sympt
  
  # Men who cease sexual activity during symptomatic infection
  GC.cease <- dat$attr$GC.cease
  CT.cease <- dat$attr$CT.cease
  syph.cease <- dat$attr$syph.cease
  
  # Diagnosis status
  diag.status.gc <- dat$attr$diag.status.gc
  diag.status.ct <- dat$attr$diag.status.ct

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
  rGC.infTime[idsInf_rgc] <- at
  rGC.lastinfTime[idsInf_rgc] <- at
  rGC.sympt[idsInf_rgc] <- rbinom(length(idsInf_rgc), 1, rgc.sympt.prob)
  rGC.timesInf[idsInf_rgc] <- rGC.timesInf[idsInf_rgc] + 1
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
  uGC.infTime[idsInf_ugc] <- at
  uGC.lastinfTime[idsInf_ugc] <- at
  uGC.sympt[idsInf_ugc] <- rbinom(length(idsInf_ugc), 1, ugc.sympt.prob)
  uGC.timesInf[idsInf_ugc] <- uGC.timesInf[idsInf_ugc] + 1
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
  rCT.infTime[idsInf_rct] <- at
  rCT.lastinfTime[idsInf_rct] <- at
  rCT.sympt[idsInf_rct] <- rbinom(length(idsInf_rct), 1, rct.sympt.prob)
  rCT.timesInf[idsInf_rct] <- rCT.timesInf[idsInf_rct] + 1
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
  uCT.infTime[idsInf_uct] <- at
  uCT.lastinfTime[idsInf_uct] <- at
  uCT.sympt[idsInf_uct] <- rbinom(length(idsInf_uct), 1, uct.sympt.prob)
  uCT.timesInf[idsInf_uct] <- uCT.timesInf[idsInf_uct] + 1
  diag.status.ct[idsInf_uct] <- 0

  
  # Syphilis ---------------------------------------------------------
  p1Inf_syph <- al[which(syphilis[al[, "p1"]] == 1 & syph.infTime[al[, "p1"]] < at &
                          syphilis[al[, "p2"]] == 0), , drop = FALSE]
  p2Inf_syph <- al[which(syphilis[al[, "p2"]] == 1 & syph.infTime[al[, "p2"]] < at &
                         syphilis[al[, "p1"]] == 0), , drop = FALSE]
  allActs_syph <- rbind(p1Inf_syph, p2Inf_syph)
  ncols <- dim(allActs_syph)[2]

  # Reorder by role: ins on the left, rec on the right, flippers represented twice
  disc.syph.ip <- allActs_syph[allActs_syph[, "ins"] %in% 1:2, , drop = FALSE]
  disc.syph.rp <- allActs_syph[allActs_syph[, "ins"] %in% c(0, 2), c(2:1, 3:ncols), drop = FALSE]

  ###### Insertive Man Infected with Syphilis (Col 1)
  if (is.null(dim(disc.syph.ip)[1])) {
      trans.syph.ip <- NULL
  }  else {
  colnames(disc.syph.ip)[1:2] <- c("ins","rec")

  # Attributes of infected
  ip.stage.syph <- stage.syph[disc.syph.ip[, 1]]

  # Base TP from VL
  ip.syph.tprob <- rep(rsyph.tprob, length(ip.stage.syph))

  # Transform to log odds
  ip.syph.tlo <- log(ip.syph.tprob/(1 - ip.syph.tprob))

  # Condom use multiplier
  not.syph.ip.UAI <- which(disc.syph.ip[, "uai"] == 0)
  ip.syph.tlo[not.syph.ip.UAI] <- ip.syph.tlo[not.syph.ip.UAI] + log(sti.cond.rr)

  # Early latent-stage multipliers
  isearlat <- which(ip.stage.syph %in% 4)
  ip.syph.tlo[isearlat] <- ip.syph.tlo[isearlat] + log(syph.earlat.rr)

  # Multiplier for syphilis transmission if receptive partner is HIV+ (infectee)
  is.HIV.infectee <- which(status[disc.syph.ip[, 2]] == 1)
  
  # Multiplier for syphilis transmission if insertive partner is HIV+ (infector)
  #is.HIV.infector <- which(status[disc.syph.ip[, 1]] == 1)
  
  ip.syph.tlo[is.HIV.infectee] <- ip.syph.tlo[is.HIV.infectee] + log(syph.rhiv.rr)
  #ip.syph.tlo[is.HIV.infector] <- ip.syph.tlo[is.HIV.infector] + log(syph.hiv.rr)

  # Retransformation to probability
  ip.syph.tprob <- plogis(ip.syph.tlo)

  # Late stage multiplier (not log odds b/c log 0 = undefined)
  islate <- which(ip.stage.syph %in% c(5, 6, 7))
  ip.syph.tprob[islate] <- ip.syph.tprob[islate] * syph.late.rr


  # Check for valid probabilities
  stopifnot(ip.syph.tprob >= 0, ip.syph.tprob <= 1)

  ## Bernoulli Transmission Events
  trans.syph.ip <- rbinom(length(ip.syph.tprob), 1, ip.syph.tprob)
  
  }
  
  ##### Receptive Man Infected with Syphilis (Col 2)
  if (is.null(dim(disc.syph.rp)[1])) {
          trans.syph.rp <- NULL
      }  else {
  colnames(disc.syph.rp)[1:2] <- c("ins", "rec")
  # Attributes of infected
  rp.stage.syph <- stage.syph[disc.syph.rp[, 2]]
   
  # Base TP from VL
  rp.syph.tprob <- rep(usyph.tprob, length(rp.stage.syph))
   
  # Transform to log odds
  rp.syph.tlo <- log(rp.syph.tprob/(1 - rp.syph.tprob))
   
  # Condom use multiplier
  not.syph.rp.UAI <- which(disc.syph.rp[, "uai"] == 0)
  rp.syph.tlo[not.syph.rp.UAI] <- rp.syph.tlo[not.syph.rp.UAI] + log(sti.cond.rr)
  
  # Early latent stage multipliers
  isearlat <- which(rp.stage.syph %in% 4)
  rp.syph.tlo[isearlat] <- rp.syph.tlo[isearlat] + log(syph.earlat.rr)

  # Multiplier for syphilis transmission if insertive partner is HIV+ (infectee)
  is.HIV.infectee <- which(status[disc.syph.rp[, 1]] == 1)
  
  # Multiplier for syphilis transmission if receptive partner is HIV+ (infector)
  is.HIV.infector <- which(status[disc.syph.rp[, 2]] == 1)
  
  rp.syph.tlo[is.HIV.infectee] <- rp.syph.tlo[is.HIV.infectee] + log(syph.uhiv.rr)
  #rp.syph.tlo[is.HIV.infector] <- rp.syph.tlo[is.HIV.infector] + log(syph.hiv.rr)
  
  # Retransformation to probability
  rp.syph.tprob <- plogis(rp.syph.tlo)
  
  # Late stage multiplier (not log odds b/c log 0 = undefined)
  islate <- which(rp.stage.syph %in% c(5, 6, 7))
  rp.syph.tprob[islate] <- rp.syph.tprob[islate] * (syph.late.rr)
  
  # Check for valid probabilities
  stopifnot(rp.syph.tprob >= 0, rp.syph.tprob <= 1)
  
  # Bernoulli transmission events
  trans.syph.rp <- rbinom(length(rp.syph.tprob), 1, rp.syph.tprob)
  
  }
  
  
  # Update attributes
  infected.syph <- inf.type.syph <- inf.role.syph <- NULL
  if (sum(trans.syph.ip, trans.syph.rp, na.rm = TRUE) > 0) {
      
      infected.syph <- c(disc.syph.ip[trans.syph.ip == 1, 2],
                    disc.syph.rp[trans.syph.rp == 1, 1])
      inf.role.syph <- c(rep(0, sum(trans.syph.ip)), rep(1, sum(trans.syph.rp)))
      inf.type.syph <- c(disc.syph.ip[trans.syph.ip == 1, "ptype"],
                    disc.syph.rp[trans.syph.rp == 1, "ptype"])
      
      dat$attr$syphilis[infected.syph] <- 1
      dat$attr$syph.infTime[infected.syph] <- at
      dat$attr$syph.lastinfTime[infected.syph] <- at
      dat$attr$stage.syph[infected.syph] <- 1
      dat$attr$stage.time.syph[infected.syph] <- 0
      dat$attr$diag.status.syph[infected.syph] <- 0
      syph.timesInf[infected.syph] <- syph.timesInf[infected.syph] + 1
      dat$attr$inf.role.syph[infected.syph] <- inf.role.syph
      dat$attr$inf.type.syph[infected.syph] <- inf.type.syph

  }
  
  
  # Set activity cessation attribute for newly infected -----------------

  # Symptomatic GC
  GC.sympt <- which(is.na(GC.cease) & (rGC.sympt == 1 | uGC.sympt == 1))
  idsGC.cease <- GC.sympt[which(rbinom(length(GC.sympt),
                                       1, gc.prob.cease) == 1)]
  GC.cease[GC.sympt] <- 0
  GC.cease[idsGC.cease] <- 1

  # Symptomatic CT
  CT.sympt <- which(is.na(CT.cease) & (rCT.sympt == 1 | uCT.sympt == 1))
  idsCT.cease <- CT.sympt[which(rbinom(length(CT.sympt),
                                       1, ct.prob.cease) == 1)]
  CT.cease[CT.sympt] <- 0
  CT.cease[idsCT.cease] <- 1

  
  # Symptomatic syphilis
  syphilis.sympt <- which(is.na(syph.cease) & (dat$attr$stage.prim.sympt == 1 |
                                                dat$attr$stage.seco.sympt == 1 |
                                                dat$attr$stage.earlat.sympt == 1 |
                                                dat$attr$stage.latelat.sympt == 1 |
                                                dat$attr$stage.latelatelat.sympt == 1 |
                                                dat$attr$stage.tert.sympt == 1))
  idssyph.cease <- syphilis.sympt[which(rbinom(length(syphilis.sympt),
                                        1, syph.prob.cease) == 1)]
  syph.cease[syphilis.sympt] <- 0
  syph.cease[idssyph.cease] <- 1

  # Output --------------------------------------------------------------

  # attributes
  dat$attr$rGC <- rGC
  dat$attr$uGC <- uGC
  dat$attr$rCT <- rCT
  dat$attr$uCT <- uCT

  dat$attr$rGC.infTime <- dat$attr$rGC.lastinfTime <- rGC.infTime
  dat$attr$uGC.infTime <- dat$attr$uGC.infTime <- uGC.infTime
  dat$attr$rCT.infTime <- dat$attr$rCT.infTime <- rCT.infTime
  dat$attr$uCT.infTime <- dat$attr$uCT.infTime <- uCT.infTime

  dat$attr$rGC.timesInf <- rGC.timesInf
  dat$attr$uGC.timesInf <- uGC.timesInf
  dat$attr$rCT.timesInf <- rCT.timesInf
  dat$attr$uCT.timesInf <- uCT.timesInf
  
  dat$attr$rGC.sympt <- rGC.sympt
  dat$attr$uGC.sympt <- uGC.sympt
  dat$attr$rCT.sympt <- rCT.sympt
  dat$attr$uCT.sympt <- uCT.sympt

  dat$attr$GC.cease <- GC.cease
  dat$attr$CT.cease <- CT.cease
  #dat$attr$syph.cease <- syph.cease
  
  dat$attr$diag.status.gc <- diag.status.gc
  dat$attr$diag.status.ct <- diag.status.ct
  
  # Summary stats
  dat$epi$incid.rgc[at] <- length(idsInf_rgc)
  dat$epi$incid.ugc[at] <- length(idsInf_ugc)
  dat$epi$incid.gc[at] <- length(idsInf_rgc) + length(idsInf_ugc)
  dat$epi$incid.rct[at] <- length(idsInf_rct)
  dat$epi$incid.uct[at] <- length(idsInf_uct)
  dat$epi$incid.ct[at] <- length(idsInf_rct) + length(idsInf_uct)
  dat$epi$incid.syph[at] <- length(infected.syph)

  dat$epi$incid.gcct.prep[at] <- length(intersect(unique(c(idsInf_rgc, idsInf_ugc,
                                                    idsInf_rct, idsInf_uct)),
                                           which(dat$attr$prepStat == 1)))
  dat$epi$incid.syph.prep[at] <- length(intersect(unique(infected.syph), which(dat$attr$prepStat == 1)))

  # Check all infected have all STI attributes
  stopifnot(all(!is.na(rGC.infTime[rGC == 1])),
            all(!is.na(rGC.sympt[rGC == 1])),
            all(!is.na(uGC.infTime[uGC == 1])),
            all(!is.na(uGC.sympt[uGC == 1])),
            all(!is.na(rCT.infTime[rCT == 1])),
            all(!is.na(rCT.sympt[rCT == 1])),
            all(!is.na(uCT.infTime[uCT == 1])),
            all(!is.na(uCT.sympt[uCT == 1])),
            all(!is.na(syph.infTime[syphilis == 1]))
            )

  if (is.null(dat$epi$times.rgc)) {
    dat$epi$times.rgc <- rep(NA, length(dat$epi$num))
    dat$epi$times.ugc <- rep(NA, length(dat$epi$num))
    dat$epi$times.rct <- rep(NA, length(dat$epi$num))
    dat$epi$times.uct <- rep(NA, length(dat$epi$num))
    dat$epi$times.syph <- rep(NA, length(dat$epi$num))
  }
  dat$epi$times.rgc[at] <- mean(rGC.timesInf, na.rm = TRUE)
  dat$epi$times.ugc[at] <- mean(uGC.timesInf, na.rm = TRUE)
  dat$epi$times.rct[at] <- mean(rCT.timesInf, na.rm = TRUE)
  dat$epi$times.uct[at] <- mean(uCT.timesInf, na.rm = TRUE)
  dat$epi$times.syph[at] <- mean(syph.timesInf, na.rm = TRUE)
  
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

  # GC Recovery ---------------------------------------------------------

  # Asymptomatic untreated
  idsRGC_asympt_ntx <- which(dat$attr$rGC == 1 &
                             dat$attr$rGC.infTime < at &
                             dat$attr$rGC.sympt == 0 &
                             (is.na(dat$attr$rGC.tx) | dat$attr$rGC.tx == 0) &
                             (is.na(dat$attr$rGC.tx.prep) | dat$attr$rGC.tx.prep == 0))
  idsUGC_asympt_ntx <- which(dat$attr$uGC == 1 &
                             dat$attr$uGC.infTime < at &
                             dat$attr$uGC.sympt == 0 &
                             (is.na(dat$attr$uGC.tx) | dat$attr$uGC.tx == 0) &
                             (is.na(dat$attr$uGC.tx.prep) | dat$attr$uGC.tx.prep == 0))

  recovRGC_asympt_ntx <- idsRGC_asympt_ntx[which(rbinom(length(idsRGC_asympt_ntx), 1,
                                                        1/rgc.asympt.int) == 1)]
  recovUGC_asympt_ntx <- idsUGC_asympt_ntx[which(rbinom(length(idsUGC_asympt_ntx), 1,
                                                        1/ugc.asympt.int) == 1)]

  # Symptomatic untreated
  idsRGC_sympt_ntx <- which(dat$attr$rGC == 1 &
                            dat$attr$rGC.infTime < at &
                            dat$attr$rGC.sympt == 1 &
                            (is.na(dat$attr$rGC.tx) | dat$attr$rGC.tx == 0) &
                            (is.na(dat$attr$rGC.tx.prep) | dat$attr$rGC.tx.prep == 0))
  idsUGC_sympt_ntx <- which(dat$attr$uGC == 1 &
                            dat$attr$uGC.infTime < at &
                            dat$attr$uGC.sympt == 1 &
                            (is.na(dat$attr$uGC.tx) | dat$attr$uGC.tx == 0) &
                            (is.na(dat$attr$uGC.tx.prep) | dat$attr$uGC.tx.prep == 0))

  # If parameter is null, uses recovery rate of asymptomatic untreated
  if (!is.na(gc.ntx.int)) {
    recovRGC_sympt_ntx <- idsRGC_sympt_ntx[which(rbinom(length(idsRGC_sympt_ntx), 1,
                                                        1/gc.ntx.int) == 1)]
    recovUGC_sympt_ntx <- idsUGC_sympt_ntx[which(rbinom(length(idsUGC_sympt_ntx), 1,
                                                        1/gc.ntx.int) == 1)]
  } else {
    recovRGC_sympt_ntx <- idsRGC_sympt_ntx[which(rbinom(length(idsRGC_sympt_ntx), 1,
                                                        1/rgc.asympt.int) == 1)]
    recovUGC_sympt_ntx <- idsUGC_sympt_ntx[which(rbinom(length(idsUGC_sympt_ntx), 1,
                                                        1/ugc.asympt.int) == 1)]
  }

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

  recovRGC <- c(recovRGC_asympt_ntx, recovRGC_sympt_ntx, recovRGC_tx)
  recovUGC <- c(recovUGC_asympt_ntx, recovUGC_sympt_ntx, recovUGC_tx)

  dat$attr$rGC[recovRGC] <- 0
  dat$attr$rGC.sympt[recovRGC] <- NA
  dat$attr$rGC.infTime[recovRGC] <- NA
  dat$attr$rGC.tx[recovRGC] <- NA
  dat$attr$rGC.tx.prep[recovRGC] <- NA
  dat$attr$diag.status.gc[recovRGC] <- NA

  dat$attr$uGC[recovUGC] <- 0
  dat$attr$uGC.sympt[recovUGC] <- NA
  dat$attr$uGC.infTime[recovUGC] <- NA
  dat$attr$uGC.tx[recovUGC] <- NA
  dat$attr$uGC.tx.prep[recovUGC] <- NA
  dat$attr$diag.status.gc[recovUGC] <- NA

  dat$attr$GC.cease[c(recovRGC, recovUGC)] <- NA



  # CT Recovery ---------------------------------------------------------

  # Asymptomatic untreated
  idsRCT_asympt_ntx <- which(dat$attr$rCT == 1 &
                             dat$attr$rCT.infTime < at &
                             dat$attr$rCT.sympt == 0 &
                             (is.na(dat$attr$rCT.tx) | dat$attr$rCT.tx == 0) &
                             (is.na(dat$attr$rCT.tx.prep) | dat$attr$rCT.tx.prep == 0))
  idsUCT_asympt_ntx <- which(dat$attr$uCT == 1 &
                             dat$attr$uCT.infTime < at &
                             dat$attr$uCT.sympt == 0 &
                             (is.na(dat$attr$uCT.tx) | dat$attr$uCT.tx == 0) &
                             (is.na(dat$attr$uCT.tx.prep) | dat$attr$uCT.tx.prep == 0))

  recovRCT_asympt_ntx <- idsRCT_asympt_ntx[which(rbinom(length(idsRCT_asympt_ntx),
                                                        1, 1/rct.asympt.int) == 1)]
  recovUCT_asympt_ntx <- idsUCT_asympt_ntx[which(rbinom(length(idsUCT_asympt_ntx),
                                                        1, 1/uct.asympt.int) == 1)]

  # Symptomatic untreated
  idsRCT_sympt_ntx <- which(dat$attr$rCT == 1 &
                            dat$attr$rCT.infTime < at &
                            dat$attr$rCT.sympt == 1 &
                            (is.na(dat$attr$rCT.tx) | dat$attr$rCT.tx == 0) &
                            (is.na(dat$attr$rCT.tx.prep) | dat$attr$rCT.tx.prep == 0))
  idsUCT_sympt_ntx <- which(dat$attr$uCT == 1 &
                            dat$attr$uCT.infTime < at &
                            dat$attr$uCT.sympt == 1 &
                            (is.na(dat$attr$uCT.tx) | dat$attr$uCT.tx == 0) &
                            (is.na(dat$attr$uCT.tx.prep) | dat$attr$uCT.tx.prep == 0))

  if (!is.na(ct.ntx.int)) {
    recovRCT_sympt_ntx <- idsRCT_sympt_ntx[which(rbinom(length(idsRCT_sympt_ntx),
                                                        1, 1/ct.ntx.int) == 1)]
    recovUCT_sympt_ntx <- idsUCT_sympt_ntx[which(rbinom(length(idsUCT_sympt_ntx),
                                                        1, 1/ct.ntx.int) == 1)]
  } else {
    recovRCT_sympt_ntx <- idsRCT_sympt_ntx[which(rbinom(length(idsRCT_sympt_ntx),
                                                        1, 1/rct.asympt.int) == 1)]
    recovUCT_sympt_ntx <- idsUCT_sympt_ntx[which(rbinom(length(idsUCT_sympt_ntx),
                                                        1, 1/uct.asympt.int) == 1)]
  }

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


  recovRCT <- c(recovRCT_asympt_ntx, recovRCT_sympt_ntx, recovRCT_tx)
  recovUCT <- c(recovUCT_asympt_ntx, recovUCT_sympt_ntx, recovUCT_tx)

  dat$attr$rCT[recovRCT] <- 0
  dat$attr$rCT.sympt[recovRCT] <- NA
  dat$attr$rCT.infTime[recovRCT] <- NA
  dat$attr$rCT.tx[recovRCT] <- NA
  dat$attr$rCT.tx.prep[recovRCT] <- NA
  dat$attr$diag.status.ct[recovRCT] <- NA
  
  dat$attr$uCT[recovUCT] <- 0
  dat$attr$uCT.sympt[recovUCT] <- NA
  dat$attr$uCT.infTime[recovUCT] <- NA
  dat$attr$uCT.tx[recovUCT] <- NA
  dat$attr$uCT.tx.prep[recovUCT] <- NA
  dat$attr$diag.status.ct[recovUCT] <- NA
  
  dat$attr$CT.cease[c(recovRCT, recovUCT)] <- NA

  
  
  # Syphilis Recovery -------------------------------------------------

  # Spontaneous untreated recovery during early stages (symptomatic or asymptomatic)
  
  # Treated (asymptomatic and symptomatic)
  idssyph_prim_tx <- which(dat$attr$syphilis == 1 &
                            dat$attr$stage.syph == 2 &
                            dat$attr$syph.infTime < at &
                            (dat$attr$syph.tx == 1 | dat$attr$syph.tx.prep == 1))
  idssyph_seco_tx <- which(dat$attr$syphilis == 1 &
                               dat$attr$stage.syph == 3 &
                               dat$attr$syph.infTime < at &
                               (dat$attr$syph.tx == 1 | dat$attr$syph.tx.prep == 1))
  idssyph_earlat_tx <- which(dat$attr$syphilis == 1 &
                               dat$attr$stage.syph == 4 &
                               dat$attr$syph.infTime < at &
                               (dat$attr$syph.tx == 1 | dat$attr$syph.tx.prep == 1))
  idssyph_latelat_tx <- which(dat$attr$syphilis == 1 &
                               dat$attr$stage.syph == 5 &
                               dat$attr$syph.infTime < at &
                               (dat$attr$syph.tx == 1 | dat$attr$syph.tx.prep == 1))
  idssyph_latelatelat_tx <- which(dat$attr$syphilis == 1 &
                               dat$attr$stage.syph == 6 &
                               dat$attr$syph.infTime < at &
                               (dat$attr$syph.tx == 1 | dat$attr$syph.tx.prep == 1))
  idssyph_tert_tx <- which(dat$attr$syphilis == 1 &
                               dat$attr$stage.syph == 7 &
                               dat$attr$syph.infTime < at &
                               (dat$attr$syph.tx == 1 | dat$attr$syph.tx.prep == 1))
                      
  idssyph_early_tx <- c(idssyph_prim_tx, idssyph_seco_tx, idssyph_earlat_tx)
  idssyph_late_tx <- c(idssyph_latelat_tx, idssyph_latelatelat_tx, idssyph_tert_tx)
  
  # Move stage-specific treated to recovered
  recovsyph_prim_tx <- idssyph_prim_tx[which(rbinom(length(idssyph_prim_tx), 1,
                                                      1/syph.early.tx.int) == 1)]
  recovsyph_seco_tx <- idssyph_seco_tx[which(rbinom(length(idssyph_seco_tx), 1,
                                                    1/syph.early.tx.int) == 1)]
  recovsyph_earlat_tx <- idssyph_earlat_tx[which(rbinom(length(idssyph_earlat_tx), 1,
                                                    1/syph.early.tx.int) == 1)]
  recovsyph_latelat_tx <- idssyph_latelat_tx[which(rbinom(length(idssyph_latelat_tx), 1,
                                                    1/syph.late.tx.int) == 1)]
  recovsyph_latelatelat_tx <- idssyph_latelatelat_tx[which(rbinom(length(idssyph_latelatelat_tx), 1,
                                                    1/syph.late.tx.int) == 1)]
  recovsyph_tert_tx <- idssyph_tert_tx[which(rbinom(length(idssyph_tert_tx), 1,
                                                    1/syph.late.tx.int) == 1)]
  
  
  # Recovery by early, late, and all
  recovsyph_early_tx <- c(recovsyph_prim_tx, recovsyph_seco_tx, recovsyph_earlat_tx)
  recovsyph_late_tx <- c(recovsyph_latelat_tx, recovsyph_latelatelat_tx, recovsyph_tert_tx)
  recovsyph <- c(recovsyph_early_tx, recovsyph_late_tx)
  
  # Update attributes
  dat$attr$syphilis[recovsyph] <- 0
  dat$attr$stage.syph[recovsyph_early_tx] <- NA
  dat$attr$stage.time.syph[recovsyph] <- NA
  dat$attr$stage.prim.sympt[recovsyph] <- NA
  dat$attr$stage.seco.sympt[recovsyph] <- NA
  dat$attr$stage.earlat.sympt[recovsyph] <- NA
  dat$attr$stage.latelat.sympt[recovsyph] <- NA
  dat$attr$stage.latelatelat.sympt[recovsyph] <- NA
  dat$attr$stage.tert.sympt[recovsyph] <- NA
  dat$attr$syph.infTime[recovsyph] <- NA
  dat$attr$diag.status.syph[recovsyph] <- NA
  dat$attr$syph.tx[recovsyph] <- NA
  dat$attr$syph.tx.prep[recovsyph] <- NA
  dat$attr$syph.cease[recovsyph] <- NA
  
  
  
  # Summary stats -----------------------------------------------------
  
  dat$epi$recov.rgc[at] <- length(unique(recovRGC))
  dat$epi$recov.ugc[at] <- length(unique(recovUGC))
  dat$epi$recov.rct[at] <- length(unique(recovRCT))
  dat$epi$recov.uct[at] <- length(unique(recovUCT))
  dat$epi$recov.syphilis[at] <- length(unique(recovsyph))
  dat$epi$recov.prim.syph[at] <- length(unique(recovsyph_prim_tx))
  dat$epi$recov.seco.syph[at] <- length(unique(recovsyph_seco_tx))
  dat$epi$recov.earlat.syph[at] <- length(unique(recovsyph_earlat_tx))
  dat$epi$recov.latelat.syph[at] <- length(unique(recovsyph_latelat_tx))
  dat$epi$recov.latelatelat.syph[at] <- length(unique(recovsyph_latelatelat_tx))
  dat$epi$recov.tert.syph[at] <- length(unique(recovsyph_tert_tx))
  dat$epi$recov.earlysyph[at] <- length(unique(recovsyph_early_tx))
  dat$epi$recov.latesyph[at] <- length(unique(recovsyph_late_tx))
  
  return(dat)
}


#' @title STI Treatment Module
#'
#' @description Stochastically simulates GC/CT and syphilis diagnosis and treatment.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
sti_tx <- function(dat, at) {

  # Parameters
  hivdx.syph.sympt.tx.rr <- dat$param$hivdx.syph.sympt.tx.rr
  gc.sympt.prob.tx <- dat$param$gc.sympt.prob.tx
  ct.sympt.prob.tx <- dat$param$ct.sympt.prob.tx

  gc.asympt.prob.tx <- dat$param$gc.asympt.prob.tx
  ct.asympt.prob.tx <- dat$param$ct.asympt.prob.tx
  
  syph.prim.sympt.prob.tx <- dat$param$syph.prim.sympt.prob.tx
  syph.prim.asympt.prob.tx <- dat$param$syph.prim.asympt.prob.tx
  syph.seco.sympt.prob.tx <- dat$param$syph.seco.sympt.prob.tx
  syph.seco.asympt.prob.tx <- dat$param$syph.seco.asympt.prob.tx
  syph.earlat.prob.tx <- dat$param$syph.earlat.prob.tx
  syph.latelat.prob.tx <- dat$param$syph.latelat.prob.tx
  syph.tert.sympt.prob.tx <- dat$param$syph.tert.sympt.prob.tx
  syph.tert.asympt.prob.tx <- dat$param$syph.tert.asympt.prob.tx
  
  prep.sti.screen.int <- dat$param$prep.sti.screen.int
  prep.sti.prob.tx <- dat$param$prep.sti.prob.tx

  prep.cont.stand.tx <- dat$param$prep.continue.stand.tx
  if (prep.cont.stand.tx == TRUE) {
    prep.stand.tx.grp <- 0:1
  } else {
    prep.stand.tx.grp <- 0
  }
  
  # Attributes from testing
  diag.status.syph <- dat$attr$diag.status.syph
  diag.status.gc <- dat$attr$diag.status.gc
  diag.status.ct <- dat$attr$diag.status.ct
  diag.status <- dat$attr$diag.status
  
  
  # symptomatic syphilis treatment
  idssyph_tx_sympt_prim <- which(dat$attr$syphilis == 1 & 
                                dat$attr$syph.infTime < at &
                                dat$attr$stage.syph == 2 &
                                dat$attr$stage.prim.sympt == 1 &
                                is.na(dat$attr$syph.tx) &
                                dat$attr$prepStat %in% prep.stand.tx.grp)
  idssyph_tx_sympt_prim_hivdx <- which(diag.status[idssyph_tx_sympt_prim] == 1)
  idssyph_tx_sympt_prim_hivnotdx <- setdiff(idssyph_tx_sympt_prim, idssyph_tx_sympt_prim_hivdx)
  txsyph_sympt_prim <- c(idssyph_tx_sympt_prim_hivnotdx[which(rbinom(length(idssyph_tx_sympt_prim_hivnotdx), 1, ((1 / hivdx.syph.sympt.tx.rr) * syph.prim.sympt.prob.tx)) == 1)],
                         idssyph_tx_sympt_prim_hivdx[which(rbinom(length(idssyph_tx_sympt_prim_hivdx), 1, syph.prim.sympt.prob.tx) == 1)])
  
  idssyph_tx_sympt_seco <- which(dat$attr$syphilis == 1 & 
                                     dat$attr$syph.infTime < at &
                                     dat$attr$stage.syph == 3 &
                                     dat$attr$stage.seco.sympt == 1 &
                                     is.na(dat$attr$syph.tx) &
                                     dat$attr$prepStat %in% prep.stand.tx.grp)
  
  idssyph_tx_sympt_seco_hivdx <- which(diag.status[idssyph_tx_sympt_seco] == 1)
  idssyph_tx_sympt_seco_hivnotdx <- setdiff(idssyph_tx_sympt_seco, idssyph_tx_sympt_seco_hivdx)
  txsyph_sympt_seco <- c(idssyph_tx_sympt_seco_hivnotdx[which(rbinom(length(idssyph_tx_sympt_seco_hivnotdx), 1, ((1 / hivdx.syph.sympt.tx.rr) * syph.seco.sympt.prob.tx)) == 1)],
                         idssyph_tx_sympt_seco_hivdx[which(rbinom(length(idssyph_tx_sympt_seco_hivdx), 1,  syph.seco.sympt.prob.tx) == 1)])

  idssyph_tx_sympt_tert <- which(dat$attr$syphilis == 1 & 
                                     dat$attr$syph.infTime < at &
                                     dat$attr$stage.syph == 7 &
                                     dat$attr$stage.tert.sympt == 1 &
                                     is.na(dat$attr$syph.tx) &
                                     dat$attr$prepStat %in% prep.stand.tx.grp)
  idssyph_tx_sympt_tert_hivdx <- which(diag.status[idssyph_tx_sympt_tert] == 1)
  idssyph_tx_sympt_tert_hivnotdx <- setdiff(idssyph_tx_sympt_tert, idssyph_tx_sympt_tert_hivdx)
  txsyph_sympt_tert <- c(idssyph_tx_sympt_tert_hivnotdx[which(rbinom(length(idssyph_tx_sympt_tert_hivnotdx), 1, syph.tert.sympt.prob.tx) == 1)],
                         idssyph_tx_sympt_tert_hivdx[which(rbinom(length(idssyph_tx_sympt_tert_hivdx), 1, (hivdx.syph.sympt.tx.rr * syph.tert.sympt.prob.tx)) == 1)])
 
  idssyph_tx_sympt <- c(idssyph_tx_sympt_prim_hivdx, idssyph_tx_sympt_prim_hivnotdx, idssyph_tx_sympt_seco_hivdx, idssyph_tx_sympt_seco_hivnotdx,
                        idssyph_tx_sympt_tert_hivdx, idssyph_tx_sympt_tert_hivnotdx)
  txsyph_sympt <- c(txsyph_sympt_prim, txsyph_sympt_seco, txsyph_sympt_tert)

  
  # Asymptomatic syphilis treatment
  idssyph_tx_asympt_prim <- which(dat$attr$syphilis == 1 & 
                                     dat$attr$syph.infTime < at &
                                     dat$attr$stage.syph == 2 &
                                     dat$attr$stage.prim.sympt == 0 &
                                     is.na(dat$attr$syph.tx) &
                                     dat$attr$prepStat %in% prep.stand.tx.grp)
  idssyph_tx_asympt_prim_dx <- which(diag.status.syph[idssyph_tx_asympt_prim] == 1)
  txsyph_asympt_prim <- idssyph_tx_asympt_prim_dx[which(rbinom(length(idssyph_tx_asympt_prim_dx), 1, syph.prim.asympt.prob.tx) == 1)]
  
  idssyph_tx_asympt_seco <- which(dat$attr$syphilis == 1 & 
                                     dat$attr$syph.infTime < at &
                                     dat$attr$stage.syph == 3 &
                                     dat$attr$stage.seco.sympt == 0 &
                                     is.na(dat$attr$syph.tx)  &
                                     dat$attr$prepStat %in% prep.stand.tx.grp)
  idssyph_tx_asympt_seco_dx <- which(diag.status.syph[idssyph_tx_asympt_seco] == 1)
  txsyph_asympt_seco <- idssyph_tx_asympt_seco_dx[which(rbinom(length(idssyph_tx_asympt_seco_dx), 1, syph.seco.asympt.prob.tx) == 1)]
  
  
  idssyph_tx_asympt_earlat <- which(dat$attr$syphilis == 1 & 
                                      dat$attr$syph.infTime < at &
                                      dat$attr$stage.syph == 4 &
                                      dat$attr$stage.earlat.sympt == 0 &
                                      is.na(dat$attr$syph.tx) &
                                      dat$attr$prepStat %in% prep.stand.tx.grp)
  idssyph_tx_asympt_earlat_dx <- which(diag.status.syph[idssyph_tx_asympt_earlat] == 1)
  txsyph_asympt_earlat <- idssyph_tx_asympt_earlat_dx[which(rbinom(length(idssyph_tx_asympt_earlat_dx), 1, syph.earlat.prob.tx) == 1)]
  
  idssyph_tx_asympt_latelat <- which(dat$attr$syphilis == 1 & 
                                        dat$attr$syph.infTime < at &
                                        (dat$attr$stage.syph == 5 | dat$attr$stage.syph == 6) &
                                        dat$attr$stage.latelat.sympt == 0 &
                                        is.na(dat$attr$syph.tx) &
                                        dat$attr$prepStat %in% prep.stand.tx.grp)
  idssyph_tx_asympt_latelat_dx <- which(diag.status.syph[idssyph_tx_asympt_latelat] == 1)
  txsyph_asympt_latelat <- idssyph_tx_asympt_latelat_dx[which(rbinom(length(idssyph_tx_asympt_latelat_dx), 1, syph.latelat.prob.tx) == 1)]
  
  idssyph_tx_asympt_tert <- which(dat$attr$syphilis == 1 & 
                                        dat$attr$syph.infTime < at &
                                        dat$attr$stage.syph == 7 &
                                        dat$attr$stage.tert.sympt == 0 &
                                        is.na(dat$attr$syph.tx) &
                                        dat$attr$prepStat %in% prep.stand.tx.grp)
  idssyph_tx_asympt_tert_dx <- which(diag.status.syph[idssyph_tx_asympt_tert] == 1)
  txsyph_asympt_tert <- idssyph_tx_asympt_tert_dx[which(rbinom(length(idssyph_tx_asympt_tert_dx), 1, syph.tert.asympt.prob.tx) == 1)]
  
  idssyph_tx_asympt <- c(idssyph_tx_asympt_prim, idssyph_tx_asympt_seco, idssyph_tx_asympt_earlat, idssyph_tx_asympt_latelat, idssyph_tx_asympt_tert)
  
  txsyph_asympt <- c(txsyph_asympt_prim, txsyph_asympt_seco, txsyph_asympt_earlat, txsyph_asympt_latelat, txsyph_asympt_tert)
  
  
  # All treated syphilis
  txsyph <- union(txsyph_sympt, txsyph_asympt)
  idssyph_tx <- union(idssyph_tx_sympt, idssyph_tx_asympt)
  
  # symptomatic gc treatment
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
  idsGC_tx_sympt <- c(idsRGC_tx_sympt, idsUGC_tx_sympt)

  txGC_sympt <- idsGC_tx_sympt[which(rbinom(length(idsGC_tx_sympt), 1,
                                            gc.sympt.prob.tx) == 1)]
  txRGC_sympt <- intersect(idsRGC_tx_sympt, txGC_sympt)
  txUGC_sympt <- intersect(idsUGC_tx_sympt, txGC_sympt)

  # asymptomatic gc treatment
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
  idsGC_tx_asympt <- c(idsRGC_tx_asympt, idsUGC_tx_asympt)
  idsGC_tx_asympt_dx <- which(diag.status.gc[idsGC_tx_asympt] == 1)
  
  txGC_asympt <- idsGC_tx_asympt_dx[which(rbinom(length(idsGC_tx_asympt_dx), 1,
                                              gc.asympt.prob.tx) == 1)]
  txRGC_asympt <- intersect(idsRGC_tx_asympt, txGC_asympt)
  txUGC_asympt <- intersect(idsUGC_tx_asympt, txGC_asympt)

  # all treated GC
  txRGC <- union(txRGC_sympt, txRGC_asympt)
  txUGC <- union(txUGC_sympt, txUGC_asympt)

  idsRGC_tx <- union(idsRGC_tx_sympt, idsRGC_tx_asympt)
  idsUGC_tx <- union(idsUGC_tx_sympt, idsUGC_tx_asympt)


  # symptomatic ct treatment
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
  idsCT_tx_sympt <- c(idsRCT_tx_sympt, idsUCT_tx_sympt)

  txCT_sympt <- idsCT_tx_sympt[which(rbinom(length(idsCT_tx_sympt), 1,
                                            ct.sympt.prob.tx) == 1)]
  txRCT_sympt <- intersect(idsRCT_tx_sympt, txCT_sympt)
  txUCT_sympt <- intersect(idsUCT_tx_sympt, txCT_sympt)

  # asymptomatic ct treatment - add diagnosis
  idsRCT_tx_asympt <- which(dat$attr$rCT == 1 &
                            dat$attr$rCT.infTime < at &
                            dat$attr$rCT.sympt == 0 &
                            dat$attr$diag.status.ct == 1 &
                            dat$attr$prepStat == 0)
  idsUCT_tx_asympt <- which(dat$attr$uCT == 1 &
                            dat$attr$uCT.infTime < at &
                            dat$attr$uCT.sympt == 0 &
                            is.na(dat$attr$uCT.tx) &
                            dat$attr$prepStat == 0)
  idsCT_tx_asympt <- c(idsRCT_tx_asympt, idsUCT_tx_asympt)
  idsCT_tx_asympt_dx <- which(diag.status.ct[idsCT_tx_asympt] == 1)
  
  txCT_asympt <- idsCT_tx_asympt_dx[which(rbinom(length(idsCT_tx_asympt_dx), 1,
                                              ct.asympt.prob.tx) == 1)]
  txRCT_asympt <- intersect(idsRCT_tx_asympt, txCT_asympt)
  txUCT_asympt <- intersect(idsUCT_tx_asympt, txCT_asympt)

  # all treated CT
  txRCT <- union(txRCT_sympt, txRCT_asympt)
  txUCT <- union(txUCT_sympt, txUCT_asympt)

  idsRCT_tx <- union(idsRCT_tx_sympt, idsRCT_tx_asympt)
  idsUCT_tx <- union(idsUCT_tx_sympt, idsUCT_tx_asympt)
  
  # Interval-based treatment for MSM on PrEP
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
  
  idssyph_prep_tx <- intersect(idsSTI_screen,
                                 which(dat$attr$syphilis == 1 &
                                        dat$attr$rGC.infTime < at &
                                         is.na(dat$attr$syph.tx.prep)))
 
  txRGC_prep <- idsRGC_prep_tx[which(rbinom(length(idsRGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txUGC_prep <- idsUGC_prep_tx[which(rbinom(length(idsUGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txRCT_prep <- idsRCT_prep_tx[which(rbinom(length(idsRCT_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txUCT_prep <- idsUCT_prep_tx[which(rbinom(length(idsUCT_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txsyph_prep <- idssyph_prep_tx[which(rbinom(length(idssyph_prep_tx), 1,
                                              prep.sti.prob.tx) == 1)]

  
  # update attributes
  
  dat$attr$syph.tx[idssyph_tx] <- 0
  dat$attr$syph.tx[txsyph] <- 1
  dat$attr$last.syph.tx <- at
  
  dat$attr$rGC.tx[idsRGC_tx] <- 0
  dat$attr$rGC.tx[txRGC] <- 1
  dat$attr$last.rGC.tx[txRGC] <- at
  
  dat$attr$uGC.tx[idsUGC_tx] <- 0
  dat$attr$uGC.tx[txUGC] <- 1
  dat$attr$last.uGC.tx[txUGC] <- at

  dat$attr$rCT.tx[idsRCT_tx] <- 0
  dat$attr$rCT.tx[txRCT] <- 1
  dat$attr$last.rCT.tx[txRCT] <- at

  dat$attr$uCT.tx[idsUCT_tx] <- 0
  dat$attr$uCT.tx[txUCT] <- 1
  dat$attr$last.uCT.tx[txUCT] <- at

  dat$attr$rGC.tx.prep[idsRGC_prep_tx] <- 0
  dat$attr$rGC.tx.prep[txRGC_prep] <- 1
  dat$attr$last.rGC.tx.prep[txRGC_prep] <- at

  dat$attr$uGC.tx.prep[idsUGC_prep_tx] <- 0
  dat$attr$uGC.tx.prep[txUGC_prep] <- 1
  dat$attr$last.uGC.tx.prep[txUGC_prep] <- at

  dat$attr$rCT.tx.prep[idsRCT_prep_tx] <- 0
  dat$attr$rCT.tx.prep[txRCT_prep] <- 1
  dat$attr$last.rCT.tx.prep[txRCT_prep] <- at
  
  dat$attr$uCT.tx.prep[idsUCT_prep_tx] <- 0
  dat$attr$uCT.tx.prep[txUCT_prep] <- 1
  dat$attr$last.uCT.tx.prep[txUCT_prep] <- at
  
  dat$attr$syph.tx.prep[idssyph_prep_tx] <- 0
  dat$attr$syph.tx.prep[txsyph_prep] <- 1
  dat$attr$last.syph.tx.prep[txsyph_prep] <- at
  
  # add tx at other site
  dat$attr$rGC.tx[which((dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1) & dat$attr$rGC == 1)] <- 1
  dat$attr$uGC.tx[which((dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1) & dat$attr$uGC == 1)] <- 1

  dat$attr$rCT.tx[which((dat$attr$uCT.tx == 1 | dat$attr$uCT.tx.prep == 1) & dat$attr$rCT == 1)] <- 1
  dat$attr$uCT.tx[which((dat$attr$rCT.tx == 1 | dat$attr$rCT.tx.prep == 1) & dat$attr$uCT == 1)] <- 1

  txRGC_all <- union(txRGC, txRGC_prep)
  txUGC_all <- union(txUGC, txUGC_prep)
  txRCT_all <- union(txRCT, txRCT_prep)
  txUCT_all <- union(txUCT, txUCT_prep)
  txsyph_all <- union(txsyph, txsyph_prep)
  
  # Update last treated time
  dat$attr$last.tx.time.syph[txsyph_all] <- at
  dat$attr$last.tx.time.rgc[txRGC_all] <- at
  dat$attr$last.tx.time.ugc[txUGC_all] <- at
  dat$attr$last.tx.time.rct[txRCT_all] <- at
  dat$attr$last.tx.time.uct[txUCT_all] <- at
  
  # Adding EPT eligibility here - if treated at this time step, can provide EPT but may not cause uptake
  dat$attr$eptElig[txRGC_all] <- 1
  dat$attr$eptStat[txRGC_all] <- 0
  dat$attr$eptElig[txUGC_all] <- 1
  dat$attr$eptStat[txUGC_all] <- 0
  dat$attr$eptElig[txRCT_all] <- 1
  dat$attr$eptStat[txRCT_all] <- 0
  dat$attr$eptElig[txUCT_all] <- 1
  dat$attr$eptStat[txUCT_all] <- 0
  dat$attr$eptEligdate[dat$attr$eptElig == 1] <- at

  # summary stats
  if (is.null(dat$epi$num.asympt.tx)) {
    dat$epi$rGCsympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$uGCsympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$rCTsympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$rCTsympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$syphsympttests <- rep(0, length(dat$control$nsteps))
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
  dat$epi$rGCsympttests[at] <- length(txRGC_sympt)
  dat$epi$uGCsympttests[at] <- length(txUGC_sympt)
  dat$epi$GCsympttests[at] <- length(c(txRGC_sympt, txUGC_sympt))
  
  dat$epi$rCTsympttests[at] <- length(txRCT_sympt)
  dat$epi$uCTsympttests[at] <- length(txUCT_sympt)
  dat$epi$CTsympttests[at] <- length(c(txRCT_sympt, txUCT_sympt))
  
  dat$epi$syphsympttests[at] <- length(txsyph_sympt)
  
  # Before start of STI testing, count # of tests in asymptomatic people as 0
  if (at < dat$param$stitest.start) {
      
      dat$epi$rGCasympttests[at] <- 0
      dat$epi$uGCasympttests[at] <- 0
      dat$epi$GCasympttests[at] <- 0
      dat$epi$rCTasympttests[at] <- 0
      dat$epi$rCTasympttests[at] <- 0
      dat$epi$syphasympttests[at] <- 0
  
  }
  asympt.tx <- c(intersect(txRGC_all, which(dat$attr$rGC.sympt == 0)),
                 intersect(txUGC_all, which(dat$attr$uGC.sympt == 0)),
                 intersect(txRCT_all, which(dat$attr$rCT.sympt == 0)),
                 intersect(txUCT_all, which(dat$attr$uCT.sympt == 0)),
                 intersect(txsyph_all, which(dat$attr$stage.prim.sympt == 0 |
                                              dat$attr$stage.seco.sympt == 0 | dat$attr$stage.earlat.sympt == 0 |
                                              dat$attr$stage.latelat.sympt == 0 | dat$attr$stage.latelatelat.sympt | 
                                              dat$attr$stage.tert.sympt == 0)))
                           
  dat$epi$num.asympt.tx[at] <- length(unique(asympt.tx))
  asympt.cases <- c(idsRGC_tx_asympt, intersect(idsRGC_prep_tx, which(dat$attr$rGC.sympt == 0)),
                    idsUGC_tx_asympt, intersect(idsUGC_prep_tx, which(dat$attr$uGC.sympt == 0)),
                    idsRCT_tx_asympt, intersect(idsRCT_prep_tx, which(dat$attr$rCT.sympt == 0)),
                    idsUCT_tx_asympt, intersect(idsUCT_prep_tx, which(dat$attr$uCT.sympt == 0)),
                    idssyph_tx_asympt, intersect(idssyph_prep_tx, which(dat$attr$stage.prim.sympt == 0 |
                                                                             dat$attr$stage.seco.sympt == 0 | dat$attr$stage.earlat.sympt == 0 |
                                                                             dat$attr$stage.latelat.sympt == 0 | dat$attr$stage.latelatelat.sympt |
                                                                             dat$attr$stage.tert.sympt == 0)))
  dat$epi$num.asympt.cases[at] <- length(unique(asympt.cases))


  asympt.tx.prep <- c(intersect(txRGC_prep, which(dat$attr$rGC.sympt == 0)),
                      intersect(txUGC_prep, which(dat$attr$uGC.sympt == 0)),
                      intersect(txRCT_prep, which(dat$attr$rCT.sympt == 0)),
                      intersect(txUCT_prep, which(dat$attr$uCT.sympt == 0)),
                      intersect(txsyph_prep, which(dat$attr$stage.prim.sympt == 0 |
                                                        dat$attr$stage.seco.sympt == 0 | dat$attr$stage.earlat.sympt == 0 |
                                                        dat$attr$stage.latelat.sympt == 0 | dat$attr$stage.latelatelat.sympt |
                                                        dat$attr$stage.tert.sympt == 0)))
  dat$epi$num.asympt.tx.prep[at] <- length(unique(asympt.tx.prep))
  
  asympt.cases.prep <- c(intersect(idsRGC_prep_tx, which(dat$attr$rGC.sympt == 0)),
                         intersect(idsUGC_prep_tx, which(dat$attr$uGC.sympt == 0)),
                         intersect(idsRCT_prep_tx, which(dat$attr$rCT.sympt == 0)),
                         intersect(idsUCT_prep_tx, which(dat$attr$uCT.sympt == 0)),
                         intersect(idssyph_prep_tx, which(dat$attr$stage.prim.sympt == 0 |
                                                               dat$attr$stage.seco.sympt == 0 | dat$attr$stage.earlat.sympt == 0 |
                                                               dat$attr$stage.latelat.sympt == 0 | dat$attr$stage.latelatelat.sympt | 
                                                               dat$attr$stage.tert.sympt == 0)))
  dat$epi$num.asympt.cases.prep[at] <- length(unique(asympt.cases.prep))

  rect.tx <- c(txRGC_all, txRCT_all)
  dat$epi$num.rect.tx[at] <- length(unique(rect.tx))
  rect.cases <- c(idsRGC_tx, idsRGC_prep_tx, idsRCT_tx, idsRCT_prep_tx)
  dat$epi$num.rect.cases[at] <- length(unique(rect.cases))

  rect.tx.prep <- c(txRGC_prep, txRCT_prep)
  dat$epi$num.rect.tx.prep[at] <- length(unique(rect.tx.prep))
  rect.cases.prep <- c(idsRGC_prep_tx, idsRCT_prep_tx)
  dat$epi$num.rect.cases.prep[at] <- length(unique(rect.cases.prep))

  return(dat)
}
