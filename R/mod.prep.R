
#' @title PrEP Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging_msm
#'
#' @export
#'
prep_msm <- function(dat, at) {

  # Function Selection ------------------------------------------------------

  if (at >= get_param(dat, "riskh.start")) {
    dat <- riskhist_msm(dat, at)
  } else {
    return(dat)
  }

  if (at < get_param(dat, "prep.start")) {
    return(dat)
  }

  ## Input

  # Attributes ----------------------------------------------------------------

  # Core Attributes
  active      <- get_attr(dat, "active")
  diag.status <- get_attr(dat, "diag.status")
  lnt         <- get_attr(dat, "last.neg.test")
  part.scrnd  <- get_attr(dat, "part.scrnd")
  race        <- get_attr(dat, "race")

  # Shared PrEP Attributes
  prepLastRisk       <- get_attr(dat, "prepLastRisk")
  prepStartTime      <- get_attr(dat, "prepStartTime")
  prepLastStiScreen  <- get_attr(dat, "prepLastStiScreen")
  prep.start.counter <- get_attr(dat, "prep.start.counter")

  # Daily PrEP Attributes
  prepElig           <- get_attr(dat, "prepElig")
  prepStat           <- get_attr(dat, "prepStat")
  prepClass          <- get_attr(dat, "prepClass")

  # EDP Attributes
  edp.prepElig       <- get_attr(dat, "edp.prepElig")
  edp.prepStat       <- get_attr(dat, "edp.prepStat")

  # DO PrEP Indications
  ind1 <- get_attr(dat, "prep.ind.uai.mono")
  ind2 <- get_attr(dat, "prep.ind.uai.conc")
  ind3 <- get_attr(dat, "prep.ind.sti")

  # EDP Indications
  edp.ind2 <- get_attr(dat, "edp.ind2")
  edp.ind3 <- get_attr(dat, "edp.ind3")
  edp.ind4 <- get_attr(dat, "edp.ind4")

  # Parameters ----------------------------------------------------------------

  # Shared or daily PrEP Parameters
  prep.start.prob           <- get_param(dat, "prep.start.prob")
  part.prep.start.prob      <- get_param(dat, "part.prep.start.prob")
  prep.adhr.dist            <- get_param(dat, "prep.adhr.dist")
  prep.require.lnt          <- get_param(dat, "prep.require.lnt")
  prep.risk.reassess.int    <- get_param(dat, "prep.risk.reassess.int")
  prep.risk.int             <- get_param(dat, "prep.risk.int")
  prep.discont.int          <- get_param(dat, "prep.discont.int")
  prep.switch.allow         <- get_param(dat, "prep.switch.allow")
  prep.switch.prob          <- get_param(dat, "prep.switch.prob")

  prep.discont.rate <- 1 - 2^(-1 / prep.discont.int)

  twind <- at - prep.risk.int

  # EDP Parameters
  prep.edp.start      <- get_param(dat, "prep.edp.start")
  edp.switch.prob     <- get_param(dat, "edp.switch.prob")
  edp.start.scenario  <- get_param(dat, "edp.start.scenario")
  edp.prep.start.prob <- get_param(dat, "edp.prep.start.prob")


  ## PrEP Eligibility ----------------------------------------------------------

  # No indications in window
  idsNoIndic <- which(
    (ind1 < twind | is.na(ind1)) &
    (ind2 < twind | is.na(ind2)) &
    (ind3 < twind | is.na(ind3))
  )
  base.cond.no <- which(!active | diag.status == 1)
  idsNoIndic <- union(idsNoIndic, base.cond.no)

  # Indications in window
  idsIndic <- which(ind1 >= twind | ind2 >= twind | ind3 >= twind)
  base.cond.yes <- which(active & diag.status == 0)
  idsIndic <- intersect(idsIndic, base.cond.yes)

  # Set daily PrEP eligibility to 1 if indications and 0 if no indications
  prepElig[idsIndic] <- 1
  prepElig[idsNoIndic] <- 0

  # No EDP starts
  if (at < prep.edp.start) {
    edp.prepElig <- 0
  }

  if (edp.start.scenario == 1 && (at >= prep.edp.start)) {
    # Set EDP eligibility to 1 if indications and 0 if no indications
    edp.prepElig[idsIndic] <- 1
    edp.prepElig[idsNoIndic] <- 0
  }

  if (edp.start.scenario == 2 && (at >= prep.edp.start)) {
    # Set EDP eligibility to 1 if any sexual activity and 0 if no sexual activity
    idsNoIndic.edp2 <- which(edp.ind2 < twind | is.na(edp.ind2))
    idsNoIndic.edp2 <- union(idsNoIndic.edp2, base.cond.no)

    idsIndic.edp2 <- which(edp.ind2 >= twind)
    idsIndic.edp2 <- intersect(idsIndic.edp2, base.cond.yes)

    edp.prepElig[idsIndic.edp2] <- 1
    edp.prepElig[idsNoIndic.edp2] <- 0
  }

  if (edp.start.scenario == 3 && (at >= prep.edp.start)) {
    # Set EDP eligibility to 1 if any sexual activity AND NO DO PrEP indications
    # 0 otherwise
    idsNoIndic.edp3 <- which(edp.ind3 < twind | is.na(edp.ind3))
    idsNoIndic.edp3 <- union(idsNoIndic.edp3, base.cond.no)

    idsIndic.edp3 <- which(edp.ind3 >= twind)
    idsIndic.edp3 <- intersect(idsIndic.edp3, base.cond.yes)

    edp.prepElig[idsIndic.edp3] <- 1
    edp.prepElig[idsNoIndic.edp3] <- 0
  }

  if (edp.start.scenario == 4 && (at >= prep.edp.start)) {
    # Set EDP eligibility to 1 if low adherence to DO PrEP and 0 otherwise
    idsNoIndic.edp4 <- which(edp.ind4 < twind | is.na(edp.ind4))
    idsNoIndic.edp4 <- union(idsNoIndic.edp4, base.cond.no)

    idsIndic.edp4 <- which(edp.ind4 >= twind)
    idsIndic.edp4 <- intersect(idsIndic.edp4, base.cond.yes)

    edp.prepElig[idsIndic.edp4] <- 1
    edp.prepElig[idsNoIndic.edp4] <- 0
  }


  ## Stoppage ------------------------------------------------------------------
  # Indication lapse
  idsRiskAssess <- which(
    active &
      (prepStat == 1 | edp.prepStat == 1) &
      lnt == at &
      (at - prepLastRisk) >= prep.risk.reassess.int
  )
  prepLastRisk[idsRiskAssess] <- at
  idsStpInd <- intersect(idsNoIndic, idsRiskAssess)

  # Random discontinuation
  idsEligStpRand <- which(active & (prepStat == 1 | edp.prepStat == 1))
  prob.discont <- prep.discont.rate[race[idsEligStpRand]]
  vecStpRand <- runif(length(idsEligStpRand)) < prob.discont
  idsStpRand <- idsEligStpRand[vecStpRand]

  # Diagnosis
  idsStpDx <- which(active & (prepStat == 1 | edp.prepStat == 1) & diag.status == 1)

  # Death
  idsStpDth <- which(!active & (prepStat == 1 | edp.prepStat == 1))

  # Reset PrEP status
  idsStp <- c(idsStpInd, idsStpRand, idsStpDx, idsStpDth)

  # Update attributes for stoppers
  prepStat[idsStp] <- 0
  edp.prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  prepLastStiScreen[idsStp] <- NA

  ## Switching -----------------------------------------------------------------
    if (prep.switch.allow && (at >= prep.edp.start)) {

    # Switching every 3 months
    timeSincePrepStart <- at - prepStartTime

    # Switching from daily to EDP
    idsSwitchElig <- which(prepStat == 1 & (timeSincePrepStart %% 90 == 0))
    vecSwitch <- runif(length(idsSwitchElig)) < prep.switch.prob
    idsSwitch <- idsSwitchElig[vecSwitch]

    # Switching from EDP to daily
    edp.idsSwitchElig <- which(edp.prepStat == 1 & (timeSincePrepStart %% 90 == 0))
    edp.vecSwitch <- runif(length(edp.idsSwitchElig)) < edp.switch.prob
    edp.idsSwitch <- edp.idsSwitchElig[edp.vecSwitch]

    # Update PrEP status attributes
    prepStat[idsSwitch] <- 0
    edp.prepStat[idsSwitch] <- 1

    prepStat[edp.idsSwitch] <- 1
    edp.prepStat[edp.idsSwitch] <- 0

    needPC <- which(is.na(prepClass[edp.idsSwitch]))
    prepClass[edp.idsSwitch[needPC]] <- sample(
      x = 1:3,
      size = length(needPC),
      replace = TRUE,
      prob = prep.adhr.dist
    )
  }

  ## Initiation ----------------------------------------------------------------

  ## Eligibility ##

  # Indications to start: identified partners
  idsStartPart <- numeric()
  idsEligStartPart <- numeric()

  idsEligStartPart <- which(
    part.scrnd == at &
      (prepStat == 0 & edp.prepStat == 0) &
      (prepElig == 1 | edp.prepElig == 1)
  )

  prob.start.part <- part.prep.start.prob[race[idsEligStartPart]]
  vecStartPart <- runif(length(idsEligStartPart)) < prob.start.part
  idsStartPart <- idsEligStartPart[vecStartPart]

  #NOTE: for validation only
  dat <- set_epi(dat, "prepStartPart", at, length(idsStartPart))

  # Update partner attributes
  dat <- set_attr(dat, "prep.start.part", at, posit_ids = idsStartPart)

  # Indications to start: general prep initiation
  if (prep.require.lnt) {
    idsEligStartDaily <- which(prepElig & prepStat == 0 & edp.prepStat == 0 & lnt == at)
    idsEligStartEDP   <- which(edp.prepElig & prepStat == 0 & edp.prepStat == 0 & lnt == at)
  } else {
    idsEligStartDaily <- which(prepElig & prepStat == 0 & edp.prepStat == 0)
    idsEligStartEDP <- which(edp.prepElig & prepStat == 0 & edp.prepStat == 0)
  }

  idsEligStartDaily <- setdiff(idsEligStartDaily, idsEligStartPart)

  # Daily PrEP start
  prob.start <- prep.start.prob[race[idsEligStartDaily]]
  vecStartDaily <- runif(length(idsEligStartDaily)) < prob.start
  idsStartDaily <- idsEligStartDaily[vecStartDaily]

  # Remove those starting daily PrEP from those eligible for EDP
  idsEligStartEDP <- setdiff(idsEligStartEDP, idsStartDaily)

  # EDP start
  edp.prob.start <- edp.prep.start.prob[race[idsEligStartEDP]]
  vecStartEDP <- runif(length(idsEligStartEDP)) < edp.prob.start
  idsStartEDP <- idsEligStartEDP[vecStartEDP]

  ## Attribute Update ##
  idsPrepStart <- c(idsStartPart, idsStartDaily, idsStartEDP)

  # Set attributes for starters
  if (length(idsPrepStart) > 0) {
    prepStat[idsStartDaily] <- 1
    edp.prepStat[idsStartEDP] <- 1
    prepStartTime[idsPrepStart] <- at
    prepLastRisk[idsPrepStart] <- at

    prep.start.counter[idsPrepStart] <- ifelse(
      is.na(prep.start.counter[idsPrepStart]),
      yes = 1,
      no = prep.start.counter[idsPrepStart] + 1
    )

    dat <- set_attr(dat, "prep.start.counter", prep.start.counter)

    # Daily PrEP adherence class
    needPC <- which(is.na(prepClass[idsStartDaily]))
    prepClass[idsStartDaily[needPC]] <- sample(
      x = 1:3,
      size = length(needPC),
      replace = TRUE,
      prob = prep.adhr.dist
    )
  }

  # Test IDs cannot be on both daily PrEP and EDP
  idsPrEPBoth <- which(prepStat == 1 & edp.prepStat == 1)
  stopifnot(length(idsPrEPBoth) == 0)

  ## Output --------------------------------------------------------------------

  # Random discontinuation
  dat <- set_epi(dat, "prep.rand.stop", at, length(idsStpRand))

  # Attributes
  dat <- set_attr(dat, "prepElig", prepElig)
  dat <- set_attr(dat, "prepStat", prepStat)
  dat <- set_attr(dat, "prepClass", prepClass)

  if (at >= prep.edp.start) {
    dat <- set_attr(dat, "edp.prepElig", edp.prepElig)
    dat <- set_attr(dat, "edp.prepStat", edp.prepStat)
  }

  dat <- set_attr(dat, "prepStartTime", prepStartTime)
  dat <- set_attr(dat, "prepLastRisk", prepLastRisk)
  dat <- set_attr(dat, "prepLastStiScreen", prepLastStiScreen)

  # Counters/summary stats to track daily PrEP and EDP starters over time
  dat <- set_epi(dat, "daily.starters", at, length(idsStartDaily))
  dat <- set_epi(dat, "edp.starters", at, length(idsStartEDP))

  if (prep.switch.allow && (at >= prep.edp.start)) {
    dat <- set_epi(dat, "prep.switch", at, length(idsSwitch))
    dat <- set_epi(dat, "edp.switch", at, length(edp.idsSwitch))
  }

  return(dat)
}


#' @title Risk History Sub-Module
#'
#' @description Sub-Module function to track the risk history of uninfected
#'              persons for purpose of PrEP targeting.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_msm <- function(dat, at) {

  ## Input
  # Attributes
  active <- get_attr(dat, "active")
  last.neg.test <- get_attr(dat, "last.neg.test")
  dx <- get_attr(dat, "diag.status")
  rGC.tx <- get_attr(dat, "rGC.tx")
  uGC.tx <- get_attr(dat, "uGC.tx")
  rCT.tx <- get_attr(dat, "rCT.tx")
  uCT.tx <- get_attr(dat, "uCT.tx")
  prep.ind.uai.mono <- get_attr(dat, "prep.ind.uai.mono",
                                override.null.error = TRUE)
  prep.ind.uai.conc <- get_attr(dat, "prep.ind.uai.conc",
                                override.null.error = TRUE)
  prepClass <- get_attr(dat, "prepClass")

  n <- length(active)
  since.test <- at - last.neg.test

  # Parameters
  time.unit <- get_param(dat, "time.unit")

  ## Edgelist, adds uai summation per partnership from act list
  # nolint start
  pid <- NULL # For R CMD Check
  # nolint end
  uai <- as.vector(rowsum(
    dat[["temp"]][["al"]][, "uai"],
    dat[["temp"]][["al"]][, "pid"]
  ))
  el <- as.data.frame(cbind(dat[["temp"]][["el"]], uai))

  if (max(el[, 1:2]) > n) stop("riskhist max(el) > n")

  # Remove concordant positive edges
  el2 <- el[el[["st2"]] == 0, ]

  # Initialize attributes
  if (is.null(prep.ind.uai.mono)) {
    dat <- set_attr(dat, "prep.ind.uai.mono", rep(NA, n))
    dat <- set_attr(dat, "prep.ind.uai.nmain", rep(NA, n))
    dat <- set_attr(dat, "prep.ind.sti", rep(NA, n))
    dat <- set_attr(dat, "edp.ind2", rep(NA, n))
    dat <- set_attr(dat, "edp.ind3", rep(NA, n))
    dat <- set_attr(dat, "edp.ind4", rep(NA, n))
  }
  if (is.null(prep.ind.uai.conc)) {
    dat <- set_attr(dat, "prep.ind.uai.conc", rep(NA, n))
  }

  ## Degree ##
  main.deg <- get_degree(dat[["el"]][[1]])
  casl.deg <- get_degree(dat[["el"]][[2]])
  inst.deg <- get_degree(dat[["el"]][[3]])


  ## Preconditions ##

  # Any UAI
  uai.any <- unique(c(
    el2[["p1"]][el2[["uai"]] > 0],
    el2[["p2"]][el2[["uai"]] > 0]
  ))

  # Monogamous partnerships: 1-sided
  tot.deg <- main.deg + casl.deg + inst.deg
  uai.mono1 <- intersect(which(tot.deg == 1), uai.any)

  # "Negative" partnerships
  tneg <- unique(c(
    el2[["p1"]][el2[["st1"]] == 0],
    el2[["p2"]][el2[["st1"]] == 0]
  ))
  fneg <- unique(c(
    el2[["p1"]][which(dx[el2[["p1"]]] == 0)],
    el2[["p2"]][which(dx[el2[["p1"]]] == 0)]
  ))
  all.neg <- c(tneg, fneg)

  ## Condition 1b: UAI in 1-sided "monogamous" "negative" partnership,
  ##               partner not tested in past 6 months
  uai.mono1.neg <- intersect(uai.mono1, all.neg)
  part.id1 <- union(
    el2[el2[["p1"]] %in% uai.mono1.neg, 2],
    el2[el2[["p2"]] %in% uai.mono1.neg, 1]
  )
  not.tested.6mo <- since.test[part.id1] > (182 / time.unit)
  not.tested.6mo <- which(not.tested.6mo)
  part.not.tested.6mo <- union(
    el2[el2[["p1"]] %in% part.id1[not.tested.6mo], 2],
    el2[el2[["p2"]] %in% part.id1[not.tested.6mo], 1]
  )

  dat <- set_attr(dat, "prep.ind.uai.mono", at, posit_ids = part.not.tested.6mo)

  ## Condition 2a: UAI + concurrency
  el2.uai <- el2[el2[["uai"]] > 0, ]
  vec <- c(el2.uai[, 1], el2.uai[, 2])
  uai.conc <- unique(vec[duplicated(vec)])
  dat <- set_attr(dat, "prep.ind.uai.conc", at, posit_ids = uai.conc)

  ## Condition 2b: UAI in non-main partnerships
  uai.nmain <- unique(c(
    el2[["p1"]][
      el2[["st1"]] == 0 &
      el2[["uai"]] > 0 &
      el2[["ptype"]] %in% 2:3
    ],
    el2[["p2"]][
      el2[["uai"]] > 0 &
      el2[["ptype"]] %in% 2:3
    ]
  ))
  dat <- set_attr(dat, "prep.ind.uai.nmain", at, posit_ids = uai.nmain)

  ## Condition 4, any STI diagnosis
  idsDx <- which(rGC.tx == 1 | uGC.tx == 1 | rCT.tx == 1 | uCT.tx == 1)
  dat <- set_attr(dat, "prep.ind.sti", at, posit_ids = idsDx)

  # EDP indications
  ## Any sexual activity
  idsAnySex <- which(tot.deg > 0)
  dat <- set_attr(dat, "edp.ind2", at, posit_ids = idsAnySex)

  ## Any sex but does not meet any of DO PrEP indications
  idsPrepInd <- union(idsDx, uai.conc)
  idsPrepInd <- union(idsPrepInd, part.not.tested.6mo)
  idsAnySex.noInd <- setdiff(idsAnySex, idsPrepInd)
  dat <- set_attr(dat, "edp.ind3", at, posit_ids = idsAnySex.noInd)

  ## Low DO PrEP adherence and not on PrEP
  idsLowAdhr <- which(prepClass == 1)
  dat <- set_attr(dat, "edp.ind4", at, posit_ids = idsLowAdhr)

  return(dat)
}


#' @title Proportionally Reallocate EDP PrEP Adherence Class Probability
#'
#' @description Shifts probabilities from the excellent-adherence category to the lower
#'              three adherence categories while maintaining the proportional
#'              distribution of those lower categories.
#'
#' @param in.pcp Input vector of length four for the \code{prep.adhr.edp.dist}
#'        parameter.
#' @param reall The pure percentage points to shift from the high adherence
#'        group to the lower three groups.
#'
#' @export
#'
reallocate_pcp <- function(in.pcp = c(0.11, 0.06, 0.09, 0.74), reall = 0) {

  dist <- in.pcp[1]/sum(in.pcp[1:3])
  dist[2] <- in.pcp[2]/sum(in.pcp[1:3])
  dist[3] <- in.pcp[3]/sum(in.pcp[1:3])

  out.pcp <- rep(NA, 4)
  out.pcp[1:3] <- in.pcp[1:3] - (dist * reall)
  out.pcp[4] <- 1 - sum(out.pcp[1:3])

  return(out.pcp)
}
