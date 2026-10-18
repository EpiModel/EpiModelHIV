
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
  # Core Attributes
  active      <- get_attr(dat, "active")
  status      <- get_attr(dat, "status")
  diag.status <- get_attr(dat, "diag.status")
  diag.time   <- get_attr(dat, "diag.time")      #to be removed
  lnt         <- get_attr(dat, "last.neg.test")
  part.scrnd  <- get_attr(dat, "part.scrnd")
  race        <- get_attr(dat, "race")

  # PrEP Attributes
  prepElig           <- get_attr(dat, "prepElig")
  prepStat           <- get_attr(dat, "prepStat")
  prepClass          <- get_attr(dat, "prepClass")
  prepLastRisk       <- get_attr(dat, "prepLastRisk")
  prepStartTime      <- get_attr(dat, "prepStartTime")
  prepLastStiScreen  <- get_attr(dat, "prepLastStiScreen")
  prep.start.counter <- get_attr(dat, "prep.start.counter")
  # Indications
  ind1 <- get_attr(dat, "prep.ind.uai.mono")
  ind2 <- get_attr(dat, "prep.ind.uai.conc")
  ind3 <- get_attr(dat, "prep.ind.sti")

  # Parameters
  prep.start.prob           <- get_param(dat, "prep.start.prob")
  part.prep.start.prob      <- get_param(dat, "part.prep.start.prob")
  prep.adhr.dist            <- get_param(dat, "prep.adhr.dist")
  prep.require.lnt          <- get_param(dat, "prep.require.lnt")
  prep.risk.reassess.int    <- get_param(dat, "prep.risk.reassess.int")
  prep.risk.int             <- get_param(dat, "prep.risk.int")
  prep.discont.int          <- get_param(dat, "prep.discont.int")

  prep.discont.rate <- 1 - 2^(-1 / prep.discont.int)


  twind <- at - prep.risk.int

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

  # Set eligibility to 1 if indications
  prepElig[idsIndic] <- 1

  # Set eligibility to 0 if no indications
  prepElig[idsNoIndic] <- 0


  ## Stoppage ------------------------------------------------------------------

  # Indication lapse
  # Rules = None, instant, yearly (CDC guidelines)
  # generalize for when risk reassessment should happen in general language
    idsRiskAssess <- which(
      active &
      prepStat == 1 &
      lnt == at &
      (at - prepLastRisk) >= prep.risk.reassess.int
    )
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)

  # Random discontinuation
  idsEligStpRand <- which(active & prepStat == 1)
  prob.discont <- prep.discont.rate[race[idsEligStpRand]]
  vecStpRand <- runif(length(idsEligStpRand)) < prob.discont
  idsStpRand <- idsEligStpRand[vecStpRand]

  # Diagnosis
  idsStpDx <- which(active & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(!active & prepStat == 1)

  # Reset PrEP status
  idsStp <- c(idsStpInd, idsStpRand, idsStpDx, idsStpDth)


  # Update attributes for stoppers
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  prepLastStiScreen[idsStp] <- NA

  ## Initiation ----------------------------------------------------------------

  ## Eligibility ##
  # Indications to start: identified partners
  idsStartPart <- numeric()
  idsEligStartPart <- numeric()

  idsEligStartPart <- which(
    part.scrnd == at &
      lnt        == at &
      prepStat   == 0 &
      prepElig   == 1
  )

  prob.start.part <- part.prep.start.prob[race[idsEligStartPart]]
  vecStartPart <- runif(length(idsEligStartPart)) < prob.start.part
  idsStartPart <- idsEligStartPart[vecStartPart]

  #For validation only
  dat <- set_epi(dat, "part.scrnd.prep", at, length(which(part.scrnd == at)))
  dat <- set_epi(dat, "scrnd.neg", at, length(which(part.scrnd == at
                                                    & lnt == at)))
  dat <- set_epi(dat, "scrnd.pos", at, length(which(part.scrnd == at
                                                    & diag.time == at)))
  dat <- set_epi(dat, "scrnd.prepon", at, length(which(part.scrnd == at
                                                       & lnt == at
                                                       & prepStat == 1)))
  dat <- set_epi(dat, "scrnd.noprep", at, length(which(part.scrnd == at
                                                       & lnt == at
                                                       & prepStat == 0)))
  dat <- set_epi(dat, "scrnd.noprepnorisk", at, length(which(part.scrnd == at
                                                             & lnt == at
                                                             & prepStat == 0
                                                             & prepElig == 0)))
  dat <- set_epi(dat, "elig.prepStartPart", at, length(idsEligStartPart))
  dat <- set_epi(dat, "prepStartPart", at, length(idsStartPart))

  # Update partner attributes
  dat <- set_attr(dat, "prep.start.part", at, posit_ids = idsStartPart)




  # Indications to start: general prep initiation
  if (prep.require.lnt) {
    idsEligStart <- which(prepElig & prepStat == 0 & lnt == at)
  } else {
    idsEligStart <- which(prepElig & prepStat == 0)
  }

  idsEligStart <- setdiff(idsEligStart, idsStartPart)

  prob.start <- prep.start.prob[race[idsEligStart]]
  vecStart <- runif(length(idsEligStart)) < prob.start
  idsStart <- idsEligStart[vecStart]

  dat <- set_epi(dat, "elig.prepStartGen", at, length(idsEligStart))
  dat <- set_epi(dat, "prepStartGen", at, length(idsStart))

  ## Attribute Update ##
  idsPrepStart <- c(idsStartPart, idsStart)
  dat <- set_epi(dat, "prepStartAll", at, length(idsPrepStart))

  # Set attributes for starters
  if (length(idsPrepStart) > 0) {
    prepStat[idsPrepStart] <- 1
    prepStartTime[idsPrepStart] <- at
    prepLastRisk[idsPrepStart] <- at

    prep.start.counter[idsPrepStart] <- ifelse(
      is.na(prep.start.counter[idsPrepStart]),
      yes = 1,
      no = prep.start.counter[idsPrepStart] + 1
    )

    dat <- set_attr(dat, "prep.start.counter", prep.start.counter)


    # PrEP adherence class
    needPC <- which(is.na(prepClass[idsPrepStart]))
    prepClass[idsPrepStart[needPC]] <- sample(
      x = 1:3,
      size = length(needPC),
      replace = TRUE,
      prob = prep.adhr.dist
    )
  }



  ## Output --------------------------------------------------------------------

  # Discontinuations
  dat <- set_epi(dat, "prep.rand.stop", at, length(idsStpRand))
  dat <- set_epi(dat, "prep.all.stop", at, length(idsStp))

  # Attributes
  dat <- set_attr(dat, "prepElig", prepElig)
  dat <- set_attr(dat, "prepStat", prepStat)
  dat <- set_attr(dat, "prepClass", prepClass)

  dat <- set_attr(dat, "prepStartTime", prepStartTime)
  dat <- set_attr(dat, "prepLastRisk", prepLastRisk)
  dat <- set_attr(dat, "prepLastStiScreen", prepLastStiScreen)


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

  return(dat)
}


#' @title Proportionally Reallocate PrEP Adherence Class Probability
#'
#' @description Shifts probabilities from the high-adherence category to the lower
#'              three adherence categories while maintaining the proportional
#'              distribution of those lower categories.
#'
#' @param in.pcp Input vector of length four for the \code{prep.adhr.dist}
#'        parameters.
#' @param reall The pure percentage points to shift from the high adherence
#'        group to the lower three groups.
#'
#' @export
#'
reallocate_pcp <- function(in.pcp = c(0.089, 0.127, 0.784), reall = 0) {

  dist <- in.pcp[1] / sum(in.pcp[1:2])
  dist[2] <- in.pcp[2] / sum(in.pcp[1:2])

  out.pcp <- rep(NA, 3)
  out.pcp[1:2] <- in.pcp[1:2] - (dist * reall)
  out.pcp[3] <- 1 - sum(out.pcp[1:2])

  return(out.pcp)
}
