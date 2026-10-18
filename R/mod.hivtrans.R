
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

  # EDP Attributes
  edp.prepStat       <- get_attr(dat, "edp.prepStat")
  edp.prepClass      <- get_attr(dat, "edp.prepClass")
  sex.edp            <- get_attr(dat, "sex.edp")
  timeLastSex        <- get_attr(dat, "timeLastSex")
  prepCombo          <- get_attr(dat, "prepCombo")
  lastPrepCombo      <- get_attr(dat, "lastPrepCombo")
  pillCount          <- get_attr(dat, "pillCount")

  # Parameters
  hiv.urai.prob      <- get_param(dat, "hiv.urai.prob")
  hiv.uiai.prob      <- get_param(dat, "hiv.uiai.prob")
  hiv.trans.scale    <- get_param(dat, "hiv.trans.scale")
  hiv.acute.or       <- get_param(dat, "hiv.acute.or")

  hiv.cond.eff.or    <- get_param(dat, "hiv.cond.eff.or")
  hiv.cond.fail.or   <- get_param(dat, "hiv.cond.fail.or")

  hiv.circ.or        <- get_param(dat, "hiv.circ.or")
  prep.adhr.or       <- get_param(dat, "prep.adhr.or")
  hiv.ugc.or         <- get_param(dat, "hiv.ugc.or")
  hiv.uct.or         <- get_param(dat, "hiv.uct.or")
  hiv.rgc.or         <- get_param(dat, "hiv.rgc.or")
  hiv.rct.or         <- get_param(dat, "hiv.rct.or")
  hiv.dual.rr        <- get_param(dat, "hiv.dual.rr")

  # EDP Parameters
  prep.adhr.edp.dist  <- get_param(dat, "prep.adhr.edp.dist")
  prep.adhr.edp.rr    <- get_param(dat, "prep.adhr.edp.rr")
  edp.class.static    <- get_param(dat, "edp.class.static")

  ## Process

  # Update attributes

  # Pull acts list
  al <- dat[["temp"]][["al"]]

  # Update sex attribute for EDP users who have sex at time step
  sex.edp <- ifelse(!is.na(sex.edp), NA, NA)
  edp.user <- which(edp.prepStat == 1)
  sex.edp[edp.user] <- ifelse(edp.user %in% al[, 1] | edp.user %in% al[, 2], 1, 0)

  # Create variable for DO PrEP users who have sex at time step
  sex.do <- rep(NA, length(dat$prepStat))
  do.user <- which(prepStat == 1)
  sex.do[do.user] <- ifelse(do.user %in% al[, 1] | do.user %in% al[, 2], 1, 0)

  # Create variable for DO user who has sex AND is high adherence
  sex.do.hi <- rep(NA, length(dat$prepStat))
  sex.do.hi <- ifelse(sex.do[do.user] == 1 & prepClass[do.user] == 3, 1, 0)

  ## Apply adherence class to EDP users who have had sex
  if (edp.class.static) {
    needPC.edp <- which(is.na(edp.prepClass) & sex.edp == 1) #static adherence class
    edp.prepClass[needPC.edp] <-
      sample(
        x = 1:4,
        size = length(needPC.edp),
        replace = TRUE,
        prob = prep.adhr.edp.dist
      )
  } else {

    # Create proportional PrEP adherence classes based on prep.adhr.edp.dist
    edp.adhr1 <- c(prep.adhr.edp.dist[1], rep(prep.adhr.edp.dist[2] / 4, 4),
                   rep(prep.adhr.edp.dist[3] / 2, 2), prep.adhr.edp.dist[4])
    edp.adhr2 <- c(edp.adhr1[7] / (edp.adhr1[7] + edp.adhr1[8]),
                   edp.adhr1[8] / (edp.adhr1[7] + edp.adhr1[8]))
    edp.adhr3 <- c(edp.adhr1[5] / (edp.adhr1[5] + edp.adhr1[6] + edp.adhr1[7] + edp.adhr1[8]),
                   edp.adhr1[6] / (edp.adhr1[5] + edp.adhr1[6] + edp.adhr1[7] + edp.adhr1[8]),
                   edp.adhr1[7] / (edp.adhr1[5] + edp.adhr1[6] + edp.adhr1[7] + edp.adhr1[8]),
                   edp.adhr1[8] / (edp.adhr1[5] + edp.adhr1[6] + edp.adhr1[7] + edp.adhr1[8]))
    edp.adhr4 <- c(edp.adhr1[1] / (edp.adhr1[1] + edp.adhr1[3]),
                   edp.adhr1[3] / (edp.adhr1[1] + edp.adhr1[3]))
    edp.adhr5 <- c(edp.adhr1[2] / (edp.adhr1[2] + edp.adhr1[4]),
                   edp.adhr1[4] / (edp.adhr1[2] + edp.adhr1[4]))
    edp.adhr6 <- c(edp.adhr1[1] / (edp.adhr1[1] + edp.adhr1[2] + edp.adhr1[3] + edp.adhr1[4]),
                   edp.adhr1[2] / (edp.adhr1[1] + edp.adhr1[2] + edp.adhr1[3] + edp.adhr1[4]),
                   edp.adhr1[3] / (edp.adhr1[1] + edp.adhr1[2] + edp.adhr1[3] + edp.adhr1[4]),
                   edp.adhr1[4] / (edp.adhr1[1] + edp.adhr1[2] + edp.adhr1[3] + edp.adhr1[4]))
    edp.adhr7 <- c(edp.adhr1[5] / (edp.adhr1[5] + edp.adhr1[6]),
                   edp.adhr1[6] / (edp.adhr1[5] + edp.adhr1[6]))

    timeSinceLastSex <- ifelse(is.na(timeLastSex) & sex.edp == 1, 0, at - timeLastSex)

    needPC.na <- which(sex.edp == 1 & is.na(lastPrepCombo))
    prepCombo[needPC.na] <- sample(
      x = 1:8,
      size = length(needPC.na),
      replace = TRUE,
      prob = edp.adhr1
      )

    needPC.82 <- which(sex.edp == 1 & lastPrepCombo %in% c(8, 2))
    prepCombo[needPC.82] <-
      ifelse(timeSinceLastSex[needPC.82] == 1,
             sample(
               x = 7:8,
               size = length(needPC.82),
               replace = TRUE,
               prob = edp.adhr2
               ),
               ifelse(timeSinceLastSex[needPC.82] == 2,
                      sample(
                        x = 5:8,
                        size = length(needPC.82),
                        replace = TRUE,
                        prob = edp.adhr3
                        ),
                        sample(
                          x = 1:8,
                          size = length(needPC.82),
                          replace = TRUE,
                          prob = edp.adhr1
                          )
                      )
             )

    needPC.74 <- which(sex.edp == 1 & lastPrepCombo %in% c(7, 4))
    prepCombo[needPC.74] <-
      ifelse(timeSinceLastSex[needPC.74] == 1,
             sample(
               x = 5:6,
               size = length(needPC.74),
               replace = TRUE,
               prob = edp.adhr7
             ),
             ifelse(timeSinceLastSex[needPC.74] == 2,
                    sample(
                      x = 1:4,
                      size = length(needPC.74),
                      replace = TRUE,
                      prob = edp.adhr6
                      ),
                    sample(
                      x = 1:8,
                      size = length(needPC.74),
                      replace = TRUE,
                      prob = edp.adhr1
                      )
                    )
            )

    needPC.63 <- which(sex.edp == 1 & lastPrepCombo %in% c(6, 3))
    prepCombo[needPC.63] <-
      ifelse(timeSinceLastSex[needPC.63] == 1,
             sample(
               x = c(2, 4),
               size = length(needPC.63),
               replace = TRUE,
               prob = edp.adhr5
               ),
             ifelse(timeSinceLastSex[needPC.63] == 2,
                    sample(
                      x = 5:8,
                      size = length(needPC.63),
                      replace = TRUE,
                      prob = edp.adhr3
                      ),
                    sample(
                      x = 1:8,
                      size = length(needPC.63),
                      replace = TRUE,
                      prob = edp.adhr1
                      )
                    )
             )

    needPC.51 <- which(sex.edp == 1 & lastPrepCombo %in% c(5, 1))
    prepCombo[needPC.51] <-
      ifelse(timeSinceLastSex[needPC.51] == 1,
             sample(
               x = c(1, 3),
               size = length(needPC.51),
               replace = TRUE,
               prob = edp.adhr4
             ),
             ifelse(timeSinceLastSex[needPC.51] == 2,
                    sample(
                      x = 1:4,
                      size = length(needPC.51),
                      replace = TRUE,
                      prob = edp.adhr6
                    ),
                    sample(
                      x = 1:8,
                      size = length(needPC.51),
                      replace = TRUE,
                      prob = edp.adhr1
                    )
             )
      )

    ## Update edp.prepClass based on prepCombo
    needPC.edp <- which(sex.edp == 1)
    edp.prepClass[needPC.edp] <-
      ifelse(prepCombo[needPC.edp] == 1, 1,
             ifelse(prepCombo[needPC.edp] %in% 2:5, 2,
                    ifelse(prepCombo[needPC.edp] %in% 6:7, 3,
                           4)
                    )
             )

    ## Update edp.prepClass to NA for EDP users who do not have sex at time step
    PC.na <- which(sex.edp == 0 | is.na(sex.edp))
    edp.prepClass[PC.na] <- NA

    ## Update lastPrepCombo attribute
    lastPrepCombo[needPC.na] <- prepCombo[needPC.na]
    lastPrepCombo[needPC.82] <- prepCombo[needPC.82]
    lastPrepCombo[needPC.74] <- prepCombo[needPC.74]
    lastPrepCombo[needPC.63] <- prepCombo[needPC.63]
    lastPrepCombo[needPC.51] <- prepCombo[needPC.51]

    ## Update pill counts
    pillCount[edp.prepStat == 0] <- 0
    pillCount[edp.prepStat == 1] <- 0
    count8765 <- which(lastPrepCombo %in% 5:8 & sex.edp == 1 & (timeSinceLastSex > 6 | is.na(timeSinceLastSex)))

    timeLastSex      <- ifelse(sex.edp == 1, at, timeLastSex)
    timeSinceLastSex <- at - timeLastSex

    count8    <- which(lastPrepCombo == 8 & timeSinceLastSex %in% 0:2)
    count7    <- which(lastPrepCombo == 7 & timeSinceLastSex %in% 0:1)
    count6    <- which(lastPrepCombo == 6 & timeSinceLastSex %in% c(0, 2))
    count5    <- which(lastPrepCombo == 5 & timeSinceLastSex == 0)
    count4    <- which(lastPrepCombo == 4 & timeSinceLastSex == 1)
    count3    <- which(lastPrepCombo == 3 & timeSinceLastSex == 2)
    count2    <- which(lastPrepCombo == 2 & timeSinceLastSex %in% 1:2)

    counterIDs <- c(count8, count7, count6, count5, count4, count3, count2)

    pillCount[counterIDs] <- 1
    pillCount[count8765]  <- 2
  }

  # Subset acts list to discordant pairs
  isDiscordant <- status[al[, 1]] == 1 & status[al[, 2]] == 0
  if (!any(isDiscordant)) {
    return(dat)
  }

  dal <- al[isDiscordant, ]

  # Create variable for EDP users who have discordant sex act at time step
  sex.edp.disc <- rep(NA, length(dat$prepStat))
  sex.edp.disc[edp.user] <- ifelse(edp.user %in% dal[, 1] | edp.user %in% dal[, 2], 1, 0)

  # Create variable for discordant sex acts that correspond with excellent adherence
  sex.edp.disc.4 <- rep(NA, length(dat$prepStat))
  sex.edp.disc.4 <- ifelse(sex.edp.disc[edp.user] == 1 & edp.prepClass[edp.user] == 4, 1, 0)

  # Create variable for DO PrEP users who have discordant sex at time step
  sex.do.disc <- rep(NA, length(dat$prepStat))
  sex.do.disc[do.user] <- ifelse(do.user %in% dal[, 1] | do.user %in% dal[, 2], 1, 0)

  # Create variable for DO user who has discordant sex AND is high adherence
  sex.do.disc.hi <- rep(NA, length(dat$prepStat))
  sex.do.disc.hi <- ifelse(sex.do.disc[do.user] == 1 & prepClass[do.user] == 3, 1, 0)

  ## Reorder by role: ins on the left, rec on the right
  disc.ip <- dal[dal[, "ins"] == 1, , drop = FALSE]
  disc.rp <- dal[dal[, "ins"] == 0, c(2:1, 3:ncol(dal)), drop = FALSE]
  colnames(disc.ip)[1:2] <- colnames(disc.rp)[1:2] <- c("p.ins", "p.rec")


  # PATP: Insertive Man Infected (Col 1) --------------------------------

  # Attributes of infected
  ip.vl     <- vl[disc.ip[, 1]]
  ip.stage  <- stage[disc.ip[, 1]]
  ip.txStat <- tx.status[disc.ip[, 1]]

  # Attributes of susceptible
  ip.prep      <- prepStat[disc.ip[, 2]]
  ip.prepclass <- prepClass[disc.ip[, 2]]
  ip.rGC       <- rGC[disc.ip[, 2]]
  ip.rCT       <- rCT[disc.ip[, 2]]

  # EDP attributes of susceptible
  ip.prep.edp      <- edp.prepStat[disc.ip[, 2]]
  ip.prepclass.edp <- edp.prepClass[disc.ip[, 2]]

  # Base TP from VL
  ip.tprob <- pmin(0.99, hiv.urai.prob * 2.45^(ip.vl - 4.5))

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
    condom.rr[not.UAI.race] <- 1 - (hiv.cond.eff.or - hiv.cond.fail.or[i])
  }
  ip.tlo[not.UAI] <- ip.tlo[not.UAI] + log(condom.rr[not.UAI])

  # Daily PrEP, by adherence class
  ip.on.prep <- which(ip.prep == 1)
  ip.tlo[ip.on.prep] <- ip.tlo[ip.on.prep] + log(prep.adhr.or[ip.prepclass[ip.on.prep]])

  # EDP, by adherence class
  ip.on.edp <- which(ip.prep.edp == 1)
  ip.tlo[ip.on.edp] <- ip.tlo[ip.on.edp] + log(prep.adhr.edp.rr[ip.prepclass.edp[ip.on.edp]])

  # Acute-stage multipliers
  isAcute <- ip.stage %in% 1:2
  ip.tlo[isAcute] <- ip.tlo[isAcute] + log(hiv.acute.or)

  ## Multiplier for STI
  is.rGC <- which(ip.rGC == 1)
  is.rCT <- which(ip.rCT == 1)
  is.rect.dual <- intersect(is.rGC, is.rCT)
  is.rGC.sing <- setdiff(is.rGC, is.rect.dual)
  is.rCT.sing <- setdiff(is.rCT, is.rect.dual)
  ip.tlo[is.rGC.sing] <- ip.tlo[is.rGC.sing] + log(hiv.rgc.or)
  ip.tlo[is.rCT.sing] <- ip.tlo[is.rCT.sing] + log(hiv.rct.or)
  ip.tlo[is.rect.dual] <- ip.tlo[is.rect.dual] +
    max(log(hiv.rgc.or), log(hiv.rct.or)) +
    min(log(hiv.rgc.or), log(hiv.rct.or)) * hiv.dual.rr

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
  rp.circ   <- circ[disc.rp[, 1]]
  rp.prep   <- prepStat[disc.rp[, 1]]
  rp.prepcl <- prepClass[disc.rp[, 1]]
  rp.uGC    <- uGC[disc.rp[, 1]]
  rp.uCT    <- uCT[disc.rp[, 1]]

  # EDP attributes of susceptible
  rp.prep.edp   <- edp.prepStat[disc.rp[, 1]]
  rp.prepcl.edp <- edp.prepClass[disc.rp[, 1]]

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
  condom.rr <- rep(NA, nrow(disc.rp))
  races <- sort(unique(race[disc.rp[, 1]]))
  for (i in races) {
    not.UAI.race <- intersect(not.UAI, which(race[disc.rp[, 1]] == i))
    condom.rr[not.UAI.race] <- 1 - (hiv.cond.eff.or - hiv.cond.fail.or[i])
  }
  rp.tlo[not.UAI] <- rp.tlo[not.UAI] + log(condom.rr[not.UAI])

  # Daily PrEP, by adherence class
  rp.on.prep <- which(rp.prep == 1)
  rp.tlo[rp.on.prep] <- rp.tlo[rp.on.prep] + log(prep.adhr.or[rp.prepcl[rp.on.prep]])

  # EDP, by adherence class
  rp.on.edp <- which(rp.prep.edp == 1)
  rp.tlo[rp.on.edp] <- rp.tlo[rp.on.edp] + log(prep.adhr.edp.rr[rp.prepcl.edp[rp.on.edp]])

  # Acute-stage multipliers
  isAcute <- rp.stage %in% 1:2
  rp.tlo[isAcute] <- rp.tlo[isAcute] + log(hiv.acute.or)

  ## Multiplier for STI
  is.uGC <- which(rp.uGC == 1)
  is.uCT <- which(rp.uCT == 1)
  is.ureth.dual <- intersect(is.uGC, is.uCT)
  is.uGC.sing <- setdiff(is.uGC, is.ureth.dual)
  is.uCT.sing <- setdiff(is.uCT, is.ureth.dual)
  rp.tlo[is.uGC.sing] <- rp.tlo[is.uGC.sing] + log(hiv.ugc.or)
  rp.tlo[is.uCT.sing] <- rp.tlo[is.uCT.sing] + log(hiv.uct.or)
  rp.tlo[is.ureth.dual] <- rp.tlo[is.ureth.dual] +
    max(log(hiv.ugc.or), log(hiv.uct.or)) +
    min(log(hiv.ugc.or), log(hiv.uct.or)) * hiv.dual.rr

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

  # EDP attributes
  dat <- set_attr(dat, "edp.prepClass", edp.prepClass)
  dat <- set_attr(dat, "sex.edp", sex.edp)
  dat <- set_attr(dat, "timeLastSex", timeLastSex)
  dat <- set_attr(dat, "lastPrepCombo", lastPrepCombo)
  dat <- set_attr(dat, "pillCount", pillCount)

  # Set Epi Trackers
  dat <- set_epi(dat, "incid", at,  length(infected))
  dat <- set_epi(dat, "incid.B", at,  sum(race[infected] == 1))
  dat <- set_epi(dat, "incid.H", at,  sum(race[infected] == 2))
  dat <- set_epi(dat, "incid.W", at,  sum(race[infected] == 3))

  # Trackers for outcome: percent of sex acts and discordant sex acts covered
  dat <- set_epi(dat, "sex.do", at, length(which(sex.do == 1)))
  dat <- set_epi(dat, "sex.do.hi", at, length(which(sex.do.hi == 1)))
  dat <- set_epi(dat, "sex.do.disc", at, length(which(sex.do.disc == 1)))
  dat <- set_epi(dat, "sex.do.disc.hi", at, length(which(sex.do.disc.hi == 1)))

  dat <- set_epi(dat, "sex.edp", at, length(which(sex.edp == 1)))
  dat <- set_epi(dat, "sex.edp.4", at, length(which(edp.prepClass == 4)))
  dat <- set_epi(dat, "sex.edp.disc", at, length(which(sex.edp.disc == 1)))
  dat <- set_epi(dat, "sex.edp.disc.4", at, length(which(sex.edp.disc.4 == 1)))

  # EDP Epi Trackers
  ## Adherence trackers
  dat <- set_epi(dat, "edp.class.1", at, length(which(edp.prepClass == 1))/length(which(sex.edp == 1)))
  dat <- set_epi(dat, "edp.class.2", at, length(which(edp.prepClass == 2))/length(which(sex.edp == 1)))
  dat <- set_epi(dat, "edp.class.3", at, length(which(edp.prepClass == 3))/length(which(sex.edp == 1)))
  dat <- set_epi(dat, "edp.class.4", at, length(which(edp.prepClass == 4))/length(which(sex.edp == 1)))

  ## Incidence trackers by EDP regimen and adherence class
  dat <- set_epi(dat, "incid.daily", at, sum(prepStat[infected] == 1, na.rm = TRUE))
  dat <- set_epi(dat, "incid.edp", at, sum(edp.prepStat[infected] == 1, na.rm = TRUE))
  dat <- set_epi(dat, "incid.edp.1", at, sum(edp.prepClass[infected] == 1, na.rm = TRUE))
  dat <- set_epi(dat, "incid.edp.2", at, sum(edp.prepClass[infected] == 2, na.rm = TRUE))
  dat <- set_epi(dat, "incid.edp.3", at, sum(edp.prepClass[infected] == 3, na.rm = TRUE))
  dat <- set_epi(dat, "incid.edp.4", at, sum(edp.prepClass[infected] == 4, na.rm = TRUE))

  #if (dat$epi$incid.edp[at] < sum(c(dat$epi$incid.edp.1[at], dat$epi$incid.edp.2[at],
  #                                  dat$epi$incid.edp.3[at], dat$epi$incid.edp.4[at]),
  #                                na.rm = TRUE)) browser()

  ## Pill count trackers
  dat <- set_epi(dat, "pills", at, sum(pillCount, na.rm = TRUE))
  dat <- set_epi(dat, "pillspp", at, sum(pillCount, na.rm = TRUE) / length(which(edp.prepStat == 1)))

  ## Record attribute history for EDP users
  #dat <- edp_attr_logger(dat, at)

  return(dat)
}

#' @title EDP Attribute History Tracker Module
#'
#' @description Tracks time-varying EDP attributes over time for EDP users
#'
#' @inheritParams aging_msm
#'
#' @export
#'

edp_attr_logger <- function(dat, at) {

  prep.edp.start <- get_param(dat, "prep.edp.start")

  # Run after EDP is introduced
  if (at >= prep.edp.start) {

    # Attributes
    edp.prepClass    <- get_attr(dat, "edp.prepClass")
    edp.prepStat     <- get_attr(dat, "edp.prepStat")
    sex.edp          <- get_attr(dat, "sex.edp")
    lastPrepCombo    <- get_attr(dat, "lastPrepCombo")
    timeLastSex      <- get_attr(dat, "timeLastSex")
    pillCount        <- get_attr(dat, "pillCount")

    edp.user <- which(edp.prepStat == 1)

    dat <- EpiModel::record_attr_history(
      dat, at,
      "edp.prepClass",
      edp.user,
      edp.prepClass[edp.user]
    )

    dat <- EpiModel::record_attr_history(
      dat, at,
      "sex.edp",
      edp.user,
      sex.edp[edp.user]
    )

    dat <- EpiModel::record_attr_history(
      dat, at,
      "lastPrepCombo",
      edp.user,
      lastPrepCombo[edp.user]
    )

    dat <- EpiModel::record_attr_history(
      dat, at,
      "timeLastSex",
      edp.user,
      timeLastSex[edp.user]
    )

    dat <- EpiModel::record_attr_history(
      dat, at,
      "pillCount",
      edp.user,
      pillCount[edp.user]
    )

  }

  # Output
  return(dat)
}
