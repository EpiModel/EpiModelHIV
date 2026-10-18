# MUST BE BEFORE THE OTHER STIS
#   -> because of windowed counter all is calc in stitest
#
# transmission: sti_trans
#   - indep? : yes, EZ
#   - also remove checks at the end?
#
# recov: sti_recov
#   - indep? : completely
#
# screening: sti_test
#   - is it indep from other 2 STIs?: yes except for windowed counters
#
# stitx_msm
#   - is indep

#' @title Syphilis Module
#'
#' @description Stochastically simulates syphilis. Manages transmission,
#'   progression, testing and treatment.
#'
#' @inheritParams aging_msm
#'
#' @export
#'
syphilis_msm <- function(dat, at) {

  dat <- syphilis_transmission(dat, at)
  dat <- syphilis_progression(dat, at)
  dat <- syphilis_testing(dat, at)
  dat <- syphilis_treatment(dat, at)

  return(dat)
}

#' Simulate Syphilis transmission
syphilis_transmission <- function(dat, at) {
  # Setup ----------------------------------------------------------------------

  ## Attributes
  race <- get_attr(dat, "race")

  # Current infection state
  syph <- get_attr(dat, "syph")
  syph.stage <- get_attr(dat, "syph.stage")

  # n Times infected
  syph.timesInf <- get_attr(dat, "syph.timesInf")

  # Infection time
  syph.infTime <- get_attr(dat, "syph.infTime")
  syph.stageTime <- get_attr(dat, "syph.stageTime")

  # Infection symptoms (non-varying)
  syph.sympt <- get_attr(dat, "syph.sympt")

  ## Parameters
  # Acquisition probabilities given contact with infected man
  syph.ips.prob <- get_param(dat, "syph.ips.prob")
  syph.e.prob <- get_param(dat, "syph.e.prob")

  # Relative risk of infection given condom use during act
  sti.cond.eff.or    <- get_param(dat, "sti.cond.eff.or")
  sti.cond.fail.or   <- get_param(dat, "sti.cond.fail.or")

  # Process --------------------------------------------------------------------

  # Pull act list
  al <- dat[["temp"]][["al"]]

  ## ins variable coding
  # ins = 0 : p2 is insertive
  # ins = 1 : p1 is insertive
  # ins = 2 : both p1 and p2 are insertive

  # Requires: discordant partnership
  allActs_syph <- which(abs(syph[al[, 1]] - syph[al[, 2]]) == 1)

  # Stage-specific transition probability
  dal_syph <- al[allActs_syph, , drop = FALSE]
  dal_syph_rev <- which(syph[dal_syph[, "p2"]] == 1)
  dal_syph[dal_syph_rev, 1:2] <- dal_syph[dal_syph_rev, 2:1]
  dal_syph[dal_syph_rev, "ins"] <- ifelse(dal_syph[dal_syph_rev, "ins"] == 1, 0, 1)

  tprob_syph <- ifelse(syph.stage[dal_syph[, "p1"]] < 4, syph.ips.prob, syph.e.prob)
  tprob_syph <- ifelse(syph.stage[dal_syph[, "p1"]] > 4, 0, tprob_syph)

  # UAI modifier
  uai_syph <- dal_syph[, "uai"]

  # Transform to log odds
  tlo_syph <- log(tprob_syph / (1 - tprob_syph))

  # Modify log odds by condom effectiveness - specific to race of insertive partner
  races <- ifelse(dal_syph[, "ins"] == 1, race[dal_syph[, "p1"]], race[dal_syph[, "p2"]])
  condom.or <- rep(NA, length(races))
  for (i in sort(unique(races))) {
    ids.race <- which(races == i)
    condom.or[ids.race] <- 1 - (sti.cond.eff.or - sti.cond.fail.or[i])
  }

  tlo_syph[uai_syph == 0] <- tlo_syph[uai_syph == 0] + log(condom.or[uai_syph == 0])

  # Back-transform to probability
  tprob_syph <- plogis(tlo_syph)

  # Stochastic transmission
  trans_syph <- runif(length(allActs_syph)) < tprob_syph

  # Determine the infected partner
  idsInf_syph <- numeric()
  if (sum(trans_syph) > 0) {
    idsInf_syph <- unique(dal_syph[trans_syph, "p2"])
    stopifnot(all(syph[idsInf_syph] == 0))
  }

  # Update attributes
  syph[idsInf_syph] <- 1
  syph.stage[idsInf_syph] <- 1
  syph.infTime[idsInf_syph] <- at
  syph.stageTime[idsInf_syph] <- at
  syph.timesInf[idsInf_syph] <- syph.timesInf[idsInf_syph] + 1
  syph.sympt[idsInf_syph] <- 0

  # Output ---------------------------------------------------------------------

  # attributes
  dat <- set_attr(dat, "syph", syph)
  dat <- set_attr(dat, "syph.stage", syph.stage)
  dat <- set_attr(dat, "syph.infTime", syph.infTime)
  dat <- set_attr(dat, "syph.stageTime", syph.stageTime)
  dat <- set_attr(dat, "syph.timesInf", syph.timesInf)
  dat <- set_attr(dat, "syph.sympt", syph.sympt)

  # Summary stats
  dat <- set_epi(dat, "incid.syph", at, length(idsInf_syph))

  return(dat)
}

#' syphilis disease progression and recovery
syphilis_progression <- function(dat, at) {
  # Setup ----------------------------------------------------------------------

  ## Attributes
  syph <- get_attr(dat, "syph")
  syph.stage <- get_attr(dat, "syph.stage")
  syph.stageTime <- get_attr(dat, "syph.stageTime")
  syph.sympt <- get_attr(dat, "syph.sympt")
  syph.tx <- get_attr(dat, "syph.tx")

  ## Parameters
  syph.tx.rate <- 1 / get_param(dat, "syph.tx.int")
  syph.ip.rate <- get_param(dat, "syph.ip.rate")
  syph.ps.rate <- get_param(dat, "syph.ps.rate")
  syph.se.rate <- get_param(dat, "syph.se.rate")
  syph.el.rate <- get_param(dat, "syph.el.rate")
  syph.lt.rate <- get_param(dat, "syph.lt.rate")
  syph.pri.sym.prob <- get_param(dat, "syph.pri.sym.prob")
  syph.sec.sym.prob <- get_param(dat, "syph.sec.sym.prob")

  # Process --------------------------------------------------------------------
  ## Syphilis Progression

  ## Incubation to primary stage ##
  nPrim <- 0
  idsEligPrim <- which(syph == 1 & syph.stage == 1 & syph.stageTime < at)
  nEligPrim <- length(idsEligPrim)

  if (nEligPrim > 0) {
    vecPrim <- which(rbinom(nEligPrim, 1, syph.ip.rate) == 1)
    if (length(vecPrim) > 0) {
      idsPrim <- idsEligPrim[vecPrim]
      nPrim <- length(idsPrim)
      syph.stage[idsPrim] <- 2
      syph.stageTime[idsPrim] <- at
      syph.sympt[idsPrim] <- sample(
        0:1,
        size = nPrim,
        prob = c(1 - syph.pri.sym.prob, syph.pri.sym.prob),
        replace = TRUE
      )
    }
  }

  ## Primary to secondary stage ##
  nSec <- 0
  idsEligSec <- which(syph == 1 & syph.stage == 2 & syph.stageTime < at)
  nEligSec <- length(idsEligSec)

  if (nEligSec > 0) {
    vecSec <- which(rbinom(nEligSec, 1, syph.ps.rate) == 1)
    if (length(vecSec) > 0) {
      idsSec <- idsEligSec[vecSec]
      nSec <- length(idsSec)
      syph.stage[idsSec] <- 3
      syph.stageTime[idsSec] <- at
      syph.sympt[idsSec] <- sample(
        0:1,
        size = nSec,
        prob = c(1 - syph.sec.sym.prob, syph.sec.sym.prob),
        replace = TRUE
      )
    }
  }

  ## Secondary to early latent stage ##
  nEarL <- 0
  idsEligEarL <- which(syph == 1 & syph.stage == 3 & syph.stageTime < at)
  nEligEarL <- length(idsEligEarL)

  if (nEligEarL > 0) {
    vecEarL <- which(rbinom(nEligEarL, 1, syph.se.rate) == 1)
    if (length(vecEarL) > 0) {
      idsEarL <- idsEligEarL[vecEarL]
      nEarL <- length(idsEarL)
      syph.stage[idsEarL] <- 4
      syph.stageTime[idsEarL] <- at
      syph.sympt[idsEarL] <- 0
    }
  }

  ## Early latent to late latent stage ##
  nLaL <- 0
  idsEligLaL <- which(syph == 1 & syph.stage == 4 & syph.stageTime < at)
  nEligLaL <- length(idsEligLaL)

  if (nEligLaL > 0) {
    vecLal <- which(rbinom(nEligLaL, 1, syph.el.rate) == 1)
    if (length(vecLal) > 0) {
      idsLaL <- idsEligLaL[vecLal]
      nLaL <- length(idsLaL)
      syph.stage[idsLaL] <- 5
      syph.stageTime[idsLaL] <- at
      syph.sympt[idsLaL] <- 0
    }
  }

  ## Late latent to tertiary stage ##
  nTer <- 0
  idsEligTer <- which(syph == 1 & syph.stage == 5 & syph.stageTime < at)
  nEligTer <- length(idsEligTer)

  if (nEligTer > 0) {
    vecTer <- which(rbinom(nEligTer, 1, syph.lt.rate) == 1)
    if (length(vecTer) > 0) {
      idsTer <- idsEligTer[vecTer]
      nTer <- length(idsTer)
      syph.stage[idsTer] <- 6
      syph.stageTime[idsTer] <- at
      syph.sympt[idsTer] <- 1
    }
  }

  # Syphilis Recovery ----------------------------------------------------
  idsSyph_recovElig <- which(syph == 1 & syph.tx == 1)
  idsSyph_recov <- idsSyph_recovElig[runif(length(idsSyph_recovElig)) < syph.tx.rate]

  # Output ---------------------------------------------------------------------

  # Update attr
  dat <- set_attr(dat, "syph.stage", syph.stage)
  dat <- set_attr(dat, "syph.stageTime", syph.stageTime)
  dat <- set_attr(dat, "syph.sympt", syph.sympt)

  dat <- set_attr(dat, "syph", 0, idsSyph_recov)
  dat <- set_attr(dat, "syph.stage", NA, idsSyph_recov)
  dat <- set_attr(dat, "syph.sympt", NA, idsSyph_recov)
  dat <- set_attr(dat, "syph.infTime", NA, idsSyph_recov)
  dat <- set_attr(dat, "syph.stageTime", NA, idsSyph_recov)
  dat <- set_attr(dat, "syph.tx", NA, idsSyph_recov)
  dat <- set_attr(dat, "diag.status.syph", NA, idsSyph_recov)
  dat <- set_attr(dat, "diag.status.syph.prep", NA, idsSyph_recov)

  # Epi trackers
  dat <- set_epi(dat, "syph.ip.flow", at, nPrim)
  dat <- set_epi(dat, "syph.ps.flow", at, nSec)
  dat <- set_epi(dat, "syph.se.flow", at, nEarL)
  dat <- set_epi(dat, "syph.el.flow", at, nLaL)
  dat <- set_epi(dat, "syph.lt.flow", at, nTer)
  dat <- set_epi(dat, "syph.sympt.num", at, sum(syph.sympt == 1, na.rm = TRUE))

  return(dat)
}

syphilis_testing <- function(dat, at) {
  # Setup ----------------------------------------------------------------------

  ## Attributes
  syph <- get_attr(dat, "syph")
  diag.status.syph <- get_attr(dat, "diag.status.syph")
  diag.status.syph.prep <- get_attr(dat, "diag.status.syph.prep")
  syph.last.screen.time <- get_attr(dat, "syph.last.screen.time")
  syph.last.dxpos.time <- get_attr(dat, "syph.last.dxpos.time")

  diag.status <- get_attr(dat, "diag.status")

  prepStat <- get_attr(dat, "prepStat")
  prepStartTime <- get_attr(dat, "prepStartTime")

  ## Parameters
  sti.screening.dist <- get_param(dat, "sti.screening.dist")
  if (!(sti.screening.dist %in% c("deterministic", "stochastic")))
    stop("sti.screening.dist must be `deterministic` or `stochastic`")

  syph.screen.hivneg.rate <- get_param(dat, "syph.screen.hivneg.rate")
  syph.screen.hivpos.rate <- get_param(dat, "syph.screen.hivpos.rate")

  sti.screen.prep.start <- get_param(dat, "sti.screen.prep.start")
  sti.screen.prep.rate <- get_param(dat, "sti.screen.prep.rate")

  # Process --------------------------------------------------------------------

  # Screening for Non-HIV+ Dx, Off-PrEP
  syph.hivneg.nprep.elig <- which(
    (diag.status.syph == 0 | is.na(diag.status.syph)) &
      (diag.status == 0 | is.na(diag.status)) &
      prepStat == 0
  )
  syph.hivneg.nprep <- sti_screen_subfx(
    elig = syph.hivneg.nprep.elig,
    test.dist = sti.screening.dist,
    screen.rate = syph.screen.hivneg.rate,
    last.screen.time = syph.last.screen.time,
    r.screen.prob = NA,
    u.status = syph,
    r.status = NA,
    at = at
  )
  diag.status.syph[syph.hivneg.nprep$dx.pos.ids] <- 1
  syph.last.screen.time[syph.hivneg.nprep$screen.ids] <- at
  syph.last.dxpos.time[syph.hivneg.nprep$dx.pos.ids] <- at

  ## Screening for Non-HIV+ Dx, On-PrEP
  if (at >= sti.screen.prep.start) {
    syph.hivneg.prep.elig <- which(
      (diag.status.syph == 0 | is.na(diag.status.syph)) &
        (diag.status == 0 | is.na(diag.status)) &
        prepStat == 1
    )
    syph.hivneg.prep <- sti_screen_subfx(
      elig = syph.hivneg.prep.elig,
      test.dist = sti.screening.dist,
      screen.rate = sti.screen.prep.rate,
      last.screen.time = syph.last.screen.time,
      r.screen.prob = NA,
      u.status = syph,
      r.status = NA,
      prep = TRUE,
      prep.start.time = prepStartTime,
      at = at
    )
    diag.status.syph[syph.hivneg.prep$dx.pos.ids] <- 1
    diag.status.syph.prep[syph.hivneg.prep$dx.pos.ids] <- 1
    syph.last.screen.time[syph.hivneg.prep$screen.ids] <- at
    syph.last.dxpos.time[syph.hivneg.prep$dx.pos.ids] <- at
  }

  ## Screening for HIV+ Dx
  syph.hivpos.elig <- which(
    (diag.status.syph == 0 | is.na(diag.status.syph)) &
      diag.status == 1 &
      prepStat == 0
  )
  syph.hivpos <- sti_screen_subfx(
    elig = syph.hivpos.elig,
    test.dist = sti.screening.dist,
    screen.rate = syph.screen.hivpos.rate,
    last.screen.time = syph.last.screen.time,
    r.screen.prob = NA,
    u.status = syph,
    r.status = NA,
    at = at
  )
  diag.status.syph[syph.hivpos$dx.pos.ids] <- 1
  syph.last.screen.time[syph.hivpos$screen.ids] <- at
  syph.last.dxpos.time[syph.hivpos$dx.pos.ids] <- at

  ## Calculate sti.dx.n.window
  dat <- update_windowed_counter(
    dat,
    posit_ids = unique(c(
      syph.hivneg.nprep$dx.pos.ids,
      syph.hivneg.prep$dx.pos.ids,
      syph.hivpos$dx.pos.ids
    )),
    window_size = get_param(dat, "syph.dx.window"),
    counter_name = "syph.dx.count.window"
  )

  ## Update attributes
  dat <- set_attr(dat, "diag.status.syph", diag.status.syph)
  dat <- set_attr(dat, "diag.status.syph.prep", diag.status.syph.prep)
  dat <- set_attr(dat, "syph.last.screen.time", syph.last.screen.time)
  dat <- set_attr(dat, "syph.last.dxpos.time", syph.last.dxpos.time)

  # Output ---------------------------------------------------------------------
  ## Screening episodes
  syph.screens.hivneg.nprep <- length(syph.hivneg.nprep$screen.ids)
  syph.dx.hivneg.nprep <- length(syph.hivneg.nprep$dx.pos.ids)
  syph.screens.hivneg.prep <- length(syph.hivneg.prep$screen.ids)
  syph.dx.hivneg.prep <- length(syph.hivneg.prep$dx.pos.ids)
  syph.screens.hivpos <- length(syph.hivpos$screen.ids)
  syph.dx.hivpos <- length(syph.hivpos$dx.pos.ids)

  dat <- set_epi(dat, "syph.screens.hivneg.nprep", at, syph.screens.hivneg.nprep)
  dat <- set_epi(dat, "syph.dx.hivneg.nprep", at, syph.dx.hivneg.nprep)
  dat <- set_epi(dat, "syph.screens.hivneg.prep", at, syph.screens.hivneg.prep)
  dat <- set_epi(dat, "syph.dx.hivneg.prep", at, syph.dx.hivneg.prep)
  dat <- set_epi(dat, "syph.screens.hivpos", at, syph.screens.hivpos)
  dat <- set_epi(dat, "syph.dx.hivpos", at, syph.dx.hivpos)

  # Any past year screening (proportions)
  hivneg.denom <- which(diag.status == 0 & prepStat == 0)
  hivpos.denom <- which(diag.status == 1 & prepStat == 0)

  ids.syph.screened.py <- which(at - syph.last.screen.time <= 364 / get_param(dat, "time.unit"))

  syph.screened.py.hivneg <- length(intersect(ids.syph.screened.py, hivneg.denom)) / length(hivneg.denom)
  syph.screened.py.hivpos <- length(intersect(ids.syph.screened.py, hivpos.denom)) / length(hivpos.denom)

  dat <- set_epi(dat, "syph.screened.py.hivneg", at, syph.screened.py.hivneg)
  dat <- set_epi(dat, "syph.screened.py.hivpos", at, syph.screened.py.hivpos)

  return(dat)
}

# Stochastically simulates syphilis treatment.
syphilis_treatment <- function(dat, at) {
  # Setup ----------------------------------------------------------------------

  ## Attributes
  syph <- get_attr(dat, "syph")
  syph.stageTime <- get_attr(dat, "syph.stageTime")
  syph.sympt <- get_attr(dat, "syph.sympt")
  syph.tx <- get_attr(dat, "syph.tx")

  diag.status.syph <- get_attr(dat, "diag.status.syph")
  diag.status.syph.prep <- get_attr(dat, "diag.status.syph.prep")
  syph.last.screen.time <- get_attr(dat, "syph.last.screen.time")
  syph.last.dxpos.time <- get_attr(dat, "syph.last.dxpos.time")

  prepStat <- get_attr(dat, "prepStat")

  ## Parameters
  syph.sympt.tx.prob  <- get_param(dat, "syph.sympt.tx.prob")[1]
  syph.asympt.tx.prob <- get_param(dat, "syph.asympt.tx.prob")[1]
  sti.prep.tx.prob <- get_param(dat, "sti.prep.tx.prob")

  # Process --------------------------------------------------------------------

  ## Symptomatic treatment eligible
  idsSyph_txElig_sympt <- which(
    syph == 1 &
      syph.stageTime == at &
      syph.sympt == 1 &
      is.na(syph.tx) # node has not previously started symptoms or tested positive during current syphilis infection
  )
  idsSyph_tx_sympt <- idsSyph_txElig_sympt[runif(length(idsSyph_txElig_sympt)) < syph.sympt.tx.prob]

  ## Asymptomatic treatment
  idsSyph_txElig_asympt <- which(
    diag.status.syph == 1 &
      syph.last.dxpos.time == at &
      prepStat == 0
  )
  idsSyph_tx_asympt <- idsSyph_txElig_asympt[runif(length(idsSyph_txElig_asympt)) < syph.asympt.tx.prob]

  ## PrEP treatment
  idsSyph_txElig_prep <- which(
    diag.status.syph.prep == 1 &
      syph.last.dxpos.time == at &
      prepStat == 1
  )
  idsSyph_tx_prep <- idsSyph_txElig_prep[
    runif(length(idsSyph_txElig_prep)) < sti.prep.tx.prob]

  ## Combining treatment
  idsSyph_txElig <- unique(c(idsSyph_txElig_sympt, idsSyph_txElig_asympt, idsSyph_txElig_prep))
  idsSyph_tx <- unique(c(idsSyph_tx_sympt, idsSyph_tx_asympt, idsSyph_tx_prep))

  ## Update Attrs
  syph.tx[idsSyph_txElig] <- 0
  syph.tx[idsSyph_tx] <- 1

  # Set last screen/dxpos times for symptomatic Dx/Tx
  syph.last.screen.time[idsSyph_tx_sympt] <- at
  syph.last.dxpos.time[idsSyph_tx_sympt] <- at

  # Output ---------------------------------------------------------------------
  dat <- set_attr(dat, "syph.tx", syph.tx)
  dat <- set_attr(dat, "syph.last.screen.time", syph.last.screen.time)
  dat <- set_attr(dat, "syph.last.dxpos.time", syph.last.dxpos.time)

  return(dat)
}
