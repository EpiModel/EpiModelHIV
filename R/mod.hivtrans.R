# TODO: rename attributes and params with: `hiv.xxx` same as gono, chla, etc

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
#' network. After transmission events, individual-level attributes for the
#' infected persons are updated and summary statistics for incidence calculated.
#'
#' The per-act transmission probability depends on the following elements:
#' insertive versus receptive role, viral load of the infected partner, an
#' acute stage infection excess risk, and condom use. Given these transmission probabilities,
#' transmission is stochastically simulating by drawing from a binomial distribution for
#' each act conditional on the per-act probability.
#'
#' @return
#' For each new infection, the disease status, infection time, and related
#' HIV attributes are updated for the infected node. Summary statistics for
#' disease incidence overall, and by race and age groups are calculated.
#'
#' @export
#'
hivtrans_msm <- function(dat, at) {

  ## Input
  # Attributes
  vl         <- get_attr(dat, "vl")
  stage      <- get_attr(dat, "stage")
  circ       <- get_attr(dat, "circ")
  status     <- get_attr(dat, "status")
  prep       <- get_attr(dat, "prep.any")
  prep.class <- get_attr(dat, "prep.class")
  gono.rect  <- get_attr(dat, "gono.rect")
  gono.uret  <- get_attr(dat, "gono.uret")
  chla.rect  <- get_attr(dat, "chla.rect")
  chla.uret  <- get_attr(dat, "chla.uret")
  syph.stage <- get_attr(dat, "syph.stage")
  race       <- get_attr(dat, "race")
  tx.status  <- get_attr(dat, "tx.status")

  # Parameters
  hiv.urai.prob    <- get_param(dat, "hiv.urai.prob")
  hiv.uiai.prob    <- get_param(dat, "hiv.uiai.prob")
  hiv.trans.scale  <- get_param(dat, "hiv.trans.scale")
  hiv.acute.or     <- get_param(dat, "hiv.acute.or")
  hiv.cond.eff.or  <- get_param(dat, "hiv.cond.eff.or")
  cond.fail.or     <- get_param(dat, "cond.fail.or")
  hiv.circ.or      <- get_param(dat, "hiv.circ.or")
  prep.adhr.or     <- get_param(dat, "prep.adhr.or")
  hiv.gono.uret.or <- get_param(dat, "hiv.gono.uret.or")
  hiv.chla.uret.or <- get_param(dat, "hiv.chla.uret.or")
  hiv.gono.rect.or <- get_param(dat, "hiv.gono.rect.or")
  hiv.chla.rect.or <- get_param(dat, "hiv.chla.rect.or")
  hiv.dual.rr      <- get_param(dat, "hiv.dual.rr")
  hiv.syph.or      <- get_param(dat, "hiv.syph.or")

  ## Process
  # Data
  al <- dat[["temp"]][["al"]]

  isDiscordant <- status[al[, 1]] == 1 & status[al[, 2]] == 0
  if (!any(isDiscordant)) {
    dat <- set_epi(dat, "incid", at, 0)
    dat <- set_epi(dat, "incid.B", at, 0)
    dat <- set_epi(dat, "incid.H", at, 0)
    dat <- set_epi(dat, "incid.W", at, 0)
    return(dat)
  }

  dal <- al[isDiscordant, , drop = FALSE]

  ## Reorder by role: ins on the left, rec on the right
  disc.ip <- dal[dal[, "ins"] == 1, , drop = FALSE]
  disc.rp <- dal[dal[, "ins"] == 0, c(2:1, 3:ncol(dal)), drop = FALSE]
  colnames(disc.ip)[1:2] <- colnames(disc.rp)[1:2] <- c("p.ins", "p.rec")


  # PATP: Insertive Man Infected (Col 1) --------------------------------

  # Attributes of infected
  ip.vl <- vl[disc.ip[, 1]]
  ip.stage <- stage[disc.ip[, 1]]
  ip.txStat <- tx.status[disc.ip[, 1]]

  # Attributes of susceptible
  ip.prep <- prep[disc.ip[, 2]]
  ip.prepcl <- prep.class[disc.ip[, 2]]
  ip.gono.rect <- gono.rect[disc.ip[, 2]]
  ip.chla.rect <- chla.rect[disc.ip[, 2]]

  # Base TP from VL
  ip.tprob <- pmin(0.99, hiv.urai.prob * 2.45^(ip.vl - 4.5))

  # Adjustment (based on Supervie JAIDS) for VL Suppressed, on ART
  ip.noTrans <- ip.vl <= log10(200) & ip.txStat == 1
  ip.tprob[ip.noTrans] <- 2.2 / 1e5

  # Transform to log odds
  ip.tlo <- log(ip.tprob / (1 - ip.tprob))

  # Condom use
  not.UAI <- which(disc.ip[, "uai"] == 0)
  condom.or <- rep(NA, nrow(disc.ip))
  races <- sort(unique(race[disc.ip[, 1]]))
  for (i in races) {
    not.UAI.race <- intersect(not.UAI, which(race[disc.ip[, 1]] == i))
    condom.or[not.UAI.race] <- 1 - (hiv.cond.eff.or - cond.fail.or[i])
  }
  ip.tlo[not.UAI] <- ip.tlo[not.UAI] + log(condom.or[not.UAI])

  # PrEP, by adherence class
  ip.on.prep <- which(ip.prep == 1)
  ip.tlo[ip.on.prep] <- ip.tlo[ip.on.prep] + log(prep.adhr.or[ip.prepcl[ip.on.prep]])

  # Acute-stage multipliers
  isAcute <- ip.stage %in% 1:2
  ip.tlo[isAcute] <- ip.tlo[isAcute] + log(hiv.acute.or)

  ## Multiplier for STI
  is.gono.rect <- which(ip.gono.rect == 1)
  is.chla.rect <- which(ip.chla.rect == 1)
  is.rect.dual <- intersect(is.gono.rect, is.chla.rect)
  is.gono.rect.sing <- setdiff(is.gono.rect, is.rect.dual)
  is.chla.rect.sing <- setdiff(is.chla.rect, is.rect.dual)
  ip.tlo[is.gono.rect.sing] <- ip.tlo[is.gono.rect.sing] + log(hiv.gono.rect.or)
  ip.tlo[is.chla.rect.sing] <- ip.tlo[is.chla.rect.sing] + log(hiv.chla.rect.or)
  ip.tlo[is.rect.dual] <- ip.tlo[is.rect.dual] +
    max(log(hiv.gono.rect.or), log(hiv.chla.rect.or)) +
    min(log(hiv.gono.rect.or), log(hiv.chla.rect.or)) * hiv.dual.rr


  # Syph multiplier
  is_syph_active <- which(syph.stage[disc.ip[, 2]] %in% c(2, 3))
  ip.tlo[is_syph_active] <- ip.tlo[is_syph_active] + log(hiv.syph.or)

  # Race-specific scalar for calibration
  races <- race[disc.ip[, 2]]
  ip.tlo <- ip.tlo + log(hiv.trans.scale[races])

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
  rp.prep <- prep[disc.rp[, 1]]
  rp.prepcl <- prep.class[disc.rp[, 1]]
  rp.gono.uret <- gono.uret[disc.rp[, 1]]
  rp.chla.uret <- chla.uret[disc.rp[, 1]]

  # Base TP from VL
  rp.tprob <- pmin(0.99, hiv.uiai.prob * 2.45^(rp.vl - 4.5))

  # Adjustment (based on Supervie JAIDS) for VL Suppressed, on ART
  rp.noTrans <- which(rp.vl <= log10(200) & rp.txStat == 1)
  rp.tprob[rp.noTrans] <- 2.2 / 1e5

  # Transform to log odds
  rp.tlo <- log(rp.tprob / (1 - rp.tprob))

  # Circumcision
  rp.tlo[rp.circ == 1] <- rp.tlo[rp.circ == 1] + log(hiv.circ.or)

  # Condom use
  not.UAI <- which(disc.rp[, "uai"] == 0)
  condom.or <- rep(NA, nrow(disc.rp))
  races <- sort(unique(race[disc.rp[, 1]]))
  for (i in races) {
    not.UAI.race <- intersect(not.UAI, which(race[disc.rp[, 1]] == i))
    condom.or[not.UAI.race] <- 1 - (hiv.cond.eff.or - cond.fail.or[i])
  }
  rp.tlo[not.UAI] <- rp.tlo[not.UAI] + log(condom.or[not.UAI])

  # PrEP, by adherence class
  rp.on.prep <- which(rp.prep == 1)
  rp.tlo[rp.on.prep] <- rp.tlo[rp.on.prep] + log(prep.adhr.or[rp.prepcl[rp.on.prep]])

  # Acute-stage multipliers
  isAcute <- rp.stage %in% 1:2
  rp.tlo[isAcute] <- rp.tlo[isAcute] + log(hiv.acute.or)

  ## Multiplier for STI
  is.gono.uret <- which(rp.gono.uret == 1)
  is.chla.uret <- which(rp.chla.uret == 1)
  is.ureth.dual <- intersect(is.gono.uret, is.chla.uret)
  is.gono.uret.sing <- setdiff(is.gono.uret, is.ureth.dual)
  is.chla.uret.sing <- setdiff(is.chla.uret, is.ureth.dual)
  rp.tlo[is.gono.uret.sing] <- rp.tlo[is.gono.uret.sing] + log(hiv.gono.uret.or)
  rp.tlo[is.chla.uret.sing] <- rp.tlo[is.chla.uret.sing] + log(hiv.chla.uret.or)
  rp.tlo[is.ureth.dual] <- rp.tlo[is.ureth.dual] +
    max(log(hiv.gono.uret.or), log(hiv.chla.uret.or)) +
    min(log(hiv.gono.uret.or), log(hiv.chla.uret.or)) * hiv.dual.rr

  # Race-specific scalar for calibration
  races <- race[disc.rp[, 1]]
  rp.tlo <- rp.tlo + log(hiv.trans.scale[races])

  # Convert back to probability
  rp.tprob <- plogis(rp.tlo)
  stopifnot(rp.tprob >= 0, rp.tprob <= 1)


  # Transmission --------------------------------------------------------

  trans.ip <- runif(length(ip.tprob)) < ip.tprob
  trans.rp <- runif(length(rp.tprob)) < rp.tprob


  # Output --------------------------------------------------------------

  infected <- numeric()
  if (sum(trans.ip, trans.rp) > 0) {
    infected <- c(
      disc.ip[trans.ip == 1, 2],
      disc.rp[trans.rp == 1, 1]
    )

    # Attributes of newly infected
    infected <- unique(infected)
    dat <- set_attr(dat, "status", 1, posit_ids = infected)
    dat <- set_attr(dat, "inf.time", at, posit_ids = infected)
    dat <- set_attr(dat, "vl", 0, posit_ids = infected)
    dat <- set_attr(dat, "stage", 1, posit_ids = infected)
    dat <- set_attr(dat, "stage.time", 0, posit_ids = infected)
    dat <- set_attr(dat, "diag.status", 0, posit_ids = infected)
    dat <- set_attr(dat, "tx.status", 0, posit_ids = infected)
    dat <- set_attr(dat, "cuml.time.on.tx", 0, posit_ids = infected)
    dat <- set_attr(dat, "cuml.time.off.tx", 0, posit_ids = infected)

  }

  # Set Epi Trackers
  dat <- set_epi(dat, "incid", at,  length(infected))
  dat <- set_epi(dat, "incid.B", at,  sum(race[infected] == 1))
  dat <- set_epi(dat, "incid.H", at,  sum(race[infected] == 2))
  dat <- set_epi(dat, "incid.W", at,  sum(race[infected] == 3))

  return(dat)
}
