
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
#' network. After transmission events, individual-level attributes for the
#' infected persons are updated and summary statistics for incidence calculated.
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

  ## Input
  # Attributes
  vl        <- get_attr(dat, "vl")
  stage     <- get_attr(dat, "stage")
  circ      <- get_attr(dat, "circ")
  status    <- get_attr(dat, "status")
  prepStat  <- get_attr(dat, "prepStat")
  prepClass <- get_attr(dat, "prepClass")
  rGC       <- get_attr(dat, "rGC")
  uGC       <- get_attr(dat, "uGC")
  rCT       <- get_attr(dat, "rCT")
  uCT       <- get_attr(dat, "uCT")
  race      <- get_attr(dat, "race")
  tx.status <- get_attr(dat, "tx.status")

  # Parameters
  URAI.prob   <- get_param(dat, "URAI.prob")
  UIAI.prob   <- get_param(dat, "UIAI.prob")
  trans.scale <- get_param(dat, "trans.scale")
  acute.rr    <- get_param(dat, "acute.rr")

  cond.eff    <- get_param(dat, "cond.eff")
  cond.fail   <- get_param(dat, "cond.fail")

  circ.rr     <- get_param(dat, "circ.rr")
  prep.hr     <- get_param(dat, "prep.adhr.hr")
  hiv.ugc.rr  <- get_param(dat, "hiv.ugc.rr")
  hiv.uct.rr  <- get_param(dat, "hiv.uct.rr")
  hiv.rgc.rr  <- get_param(dat, "hiv.rgc.rr")
  hiv.rct.rr  <- get_param(dat, "hiv.rct.rr")
  hiv.dual.rr <- get_param(dat, "hiv.dual.rr")

  ## Process
  # Data
  al <- dat$temp$al

  isDiscordant <- status[al[, 1]] == 1 & status[al[, 2]] == 0
  if (!any(isDiscordant)) {
    return(dat)
  }

  dal <- al[isDiscordant, ]

  ## Reorder by role: ins on the left, rec on the right
  disc.ip <- dal[dal[, "ins"] == 1, ]
  disc.rp <- dal[dal[, "ins"] == 0, c(2:1, 3:ncol(dal))]
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
  ip.noTrans <- ip.vl <= log10(200) & ip.txStat == 1
  ip.tprob[ip.noTrans] <- 2.2 / 1e5

  # Transform to log odds
  ip.tlo <- log(ip.tprob / (1 - ip.tprob))

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
  isAcute <- ip.stage %in% 1:2
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
  rp.tprob[rp.noTrans] <- 2.2 / 1e5

  # Transform to log odds
  rp.tlo <- log(rp.tprob / (1 - rp.tprob))

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
  isAcute <- rp.stage %in% 1:2
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
