
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
  diag.status <- get_attr(dat, "diag.status")
  lnt         <- get_attr(dat, "last.neg.test")
  part.scrnd  <- get_attr(dat, "part.scrnd")
  race        <- get_attr(dat, "race")

  # PrEP Attributes
  prepElig           <- get_attr(dat, "prepElig")
  prepStat           <- get_attr(dat, "prepStat")
  prepClass          <- get_attr(dat, "prepClass")
  prepLastRisk       <- get_attr(dat, "prepLastRisk")
  prepStartTime      <- get_attr(dat, "prepStartTime")
  prep.start.counter <- get_attr(dat, "prep.start.counter")

  # mdd/mde Attributes
  mdeCurr                 <- get_attr(dat, "mdeCurr")
  sexStigmaClass          <- get_attr(dat, "sexStigmaClass")
  mdd.txintv.status       <- get_attr(dat, "mdd.txintv.status")
  mdd.txintvinit.time     <- get_attr(dat, "mdd.txintvinit.time")
  mdd.txintv.last.endtime <- get_attr(dat, "mdd.txintv.last.endtime")
  mdd.tx.status           <- get_attr(dat, "mdd.tx.status")
  mdd.currtx.tottime      <- get_attr(dat, "mdd.currtx.tottime")

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

  prepinit.prob.mult        <- get_param(dat, "prepinit.prob.mult")
  prephalt.prob.mult        <- get_param(dat, "prephalt.prob.mult")
  minadeq.txdur             <- get_param(dat, "minadeq.txdur")
  stigcounse.efxdur         <- get_param(dat, "stigcounse.efxdur")
  mh.stigtxefx.prepinit     <- get_param(dat, "mh.stigtxefx.prepinit")
  mh.stigtxefx.prephalt     <- get_param(dat, "mh.stigtxefx.prephalt")



  prep.discont.rate <- 1 - 2^(-1 / prep.discont.int)

  #create stigma_mental distress classes
  stigmde <- rep(0, length(sexStigmaClass))
  mdeupdated <- mdeCurr
  mdeupdated[mdd.tx.status == 1 & mdd.currtx.tottime >= minadeq.txdur] <- 0

  stigmde[sexStigmaClass == 1 & mdeupdated == 1] <- 1
  stigmde[sexStigmaClass == 2 & mdeupdated == 1] <- 2
  stigmde[sexStigmaClass == 3 & mdeupdated == 1] <- 3
  stigmde[sexStigmaClass == 4 & mdeupdated == 1] <- 4
  stigmde[sexStigmaClass == 1 & mdeupdated == 0] <- 5
  stigmde[sexStigmaClass == 2 & mdeupdated == 0] <- 6
  stigmde[sexStigmaClass == 3 & mdeupdated == 0] <- 7
  stigmde[sexStigmaClass == 4 & mdeupdated == 0] <- 8



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

      #add depression effects
      if(at >= get_param(dat, "mhefx.start")) {
        mdeCurr.elig <- mdeupdated[idsEligStpRand]

        #get eligible individuals' txintv participation status (id who's on stigma-impacting
        #mdd tx), tottime (if on txintv) or time since txintv ended if done recently
        txintv.status          <- mdd.txintv.status[idsEligStpRand]
        txintv.tottime         <- at - mdd.txintvinit.time[idsEligStpRand]
        txintv.lastendtime     <- at - mdd.txintv.last.endtime[idsEligStpRand]

        #identify those eligible for multiplier adjustment because they are on stigma-impacting mdd tx
        #here this would be individual's with any stigma class + mde present
        #[i.e. stigmde == 1,2,3 or 4] and those with anticipated healthcare sigma without
        #mde present [i.e. stigmde==6], who have received minimally adequate tx or recently
        #stopped tx after remission
        txintv.multadj.elig    <- rep(0, length(idsEligStpRand))
        txintv.multadj.elig[which(txintv.tottime >= minadeq.txdur |
                                    txintv.lastendtime <= stigcounse.efxdur)] <- 1

        #get multiplier (pulled from param) based on individual's stigmde combo
        prephalt.mult <- prephalt.prob.mult[mdeCurr.elig + 1]
        #remove stigma impacts for those on txintv
        prephalt.mult[txintv.multadj.elig == 1] <- mh.stigtxefx.prephalt


        prob.discont <- exp(log(prob.discont) + log(prephalt.mult))
      }

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




  ## Initiation ----------------------------------------------------------------

  # Partner Services: identified partners
  idsStartPart <- numeric()
  idsEligStartPart <- numeric()

  idsEligStartPart <- which(
    part.scrnd == at &
      prepStat == 0 &
      prepElig == 1
  )

  prob.start.part <- part.prep.start.prob[race[idsEligStartPart]]

      #add depression effects
      if(at >= get_param(dat, "mhefx.start")) {
        #get eligible individuals' stigmaClass and mde combo (aim1b findings will be used to
        #adjust prepinit prob for each individual based on the individual's stigmde combo)
        stigmde.combo          <- stigmde[idsEligStartPart]

        #get eligible individuals' txintv participation status (id who's on stigma-impacting
        #mdd tx), tottime (if on txintv) or time since txintv ended if done recently
        txintv.status          <- mdd.txintv.status[idsEligStartPart]
        txintv.tottime         <- at - mdd.txintvinit.time[idsEligStartPart]
        txintv.lastendtime     <- at - mdd.txintv.last.endtime[idsEligStartPart]

        #identify those eligible for multiplier adjustment because they are on stigma-impacting mdd tx
        #here this would be individual's with any stigma class + mde present
        #[i.e. stigmde == 1,2,3 or 4] and those with anticipated healthcare sigma without
        #mde present [i.e. stigmde==6], who have received minimally adequate tx or recently
        #stopped tx after remission
        txintv.multadj.elig        <- rep(0, length(idsEligStartPart))
        txintv.multadj.elig[which(txintv.tottime >= minadeq.txdur |
                                     txintv.lastendtime <= stigcounse.efxdur)] <- 1

        #get multiplier (pulled from param) based on individual's stigmde combo
        prob.mult <- prepinit.prob.mult[stigmde.combo]
        #remove stigma impacts for those on txintv
        prob.mult[txintv.multadj.elig == 1] <- mh.stigtxefx.prepinit

        #apply adjusted multiplier to test rate
        prob.start.part <- plogis(qlogis(prob.start.part) + log(prob.mult))
      }

  vecStartPart <- runif(length(idsEligStartPart)) < prob.start.part
  idsStartPart <- idsEligStartPart[vecStartPart]

  # Update partner attributes
  dat <- set_attr(dat, "prep.start.part", at, posit_ids = idsStartPart)



  # Indications to start: general prep initiation
  if (prep.require.lnt) {
    idsEligStart <- which(prepElig & prepStat == 0 & lnt == at)
  } else {
    idsEligStart <- which(prepElig & prepStat == 0)
  }
  idsEligStart <- setdiff(idsEligStart, idsEligStartPart)

  prob.start <- prep.start.prob[race[idsEligStart]]
  #add depression effects
  if(at >= get_param(dat, "mhefx.start")) {
    #get eligible individuals' stigmaClass and mde combo and aim1b findings (to be used for
    #adjusting prepinit prob for each individual based on the individual's stigmde combo)
    stigmde.combo          <- stigmde[idsEligStart]

    #get eligible individuals' txintv participation status (id who's on stigma-impacting
    #mdd tx), tottime (if on txintv) or time since txintv ended if done recently
    txintv.status          <- mdd.txintv.status[idsEligStart]
    txintv.tottime         <- at - mdd.txintvinit.time[idsEligStart]
    txintv.lastendtime     <- at - mdd.txintv.last.endtime[idsEligStart]

    #identify those eligible for multiplier adjustment because they are on stigma-impacting mdd tx
    #here this would be individual's with any stigma class + mde present
    #[i.e. stigmde == 1,2,3 or 4] and those with anticipated healthcare sigma without
    #mde present [i.e. stigmde==6], who have received minimally adequate tx or recently
    #stopped tx after remission
    txintv.multelig        <- rep(0, length(idsEligStart))
    txintv.multelig[which(txintv.tottime >= minadeq.txdur |
                             txintv.lastendtime <= stigcounse.efxdur)] <- 1

    #get multiplier (pulled from param) based on individual's stigmde combo
    prob.mult <- prepinit.prob.mult[stigmde.combo]
    #remove stigma impacts for those on txintv
    prob.mult[txintv.multelig == 1] <- mh.stigtxefx.prepinit

    #apply adjusted multiplier to test rate
    prob.start <- plogis(qlogis(prob.start) + log(prob.mult))
  }

  vecStart <- runif(length(idsEligStart)) < prob.start
  idsStart <- idsEligStart[vecStart]


  ## Attribute Update ##
  idsPrepStart <- c(idsStartPart, idsStart)



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
  # Attributes
  dat <- set_attr(dat, "prepElig", prepElig)
  dat <- set_attr(dat, "prepStat", prepStat)
  dat <- set_attr(dat, "prepClass", prepClass)

  dat <- set_attr(dat, "prepStartTime", prepStartTime)
  dat <- set_attr(dat, "prepLastRisk", prepLastRisk)


  #epi trackers
  dat <- set_epi(dat, "elig.prepStartGen", at, length(idsEligStart))
  dat <- set_epi(dat, "prepStartGen", at, length(idsStart))

  dat <- set_epi(dat, "elig.prepStartPart", at, length(idsEligStartPart))
  dat <- set_epi(dat, "prepStartPart", at, length(idsStartPart))

  dat <- set_epi(dat, "prepStart.numstigma1", at,
                 length(intersect(idsPrepStart, which(sexStigmaClass == 1))))
  dat <- set_epi(dat, "prepStart.numstigma2", at,
                 length(intersect(idsPrepStart, which(sexStigmaClass == 2))))
  dat <- set_epi(dat, "prepStart.numstigma3", at,
                 length(intersect(idsPrepStart, which(sexStigmaClass == 3))))
  dat <- set_epi(dat, "prepStart.numstigma4", at,
                 length(intersect(idsPrepStart, which(sexStigmaClass == 4))))

  dat <- set_epi(dat, "prepStart.numactivemde0", at,
                 length(intersect(idsPrepStart, which(mdeCurr == 0))))
  dat <- set_epi(dat, "prepStart.numactivemde1", at,
                 length(intersect(idsPrepStart, which(mdeCurr == 1))))

  dat <- set_epi(dat, "prep.all.stop", at, length(idsStp))
  dat <- set_epi(dat, "prep.rand.stop", at, length(idsStpRand))



  return(dat)
}


#' @title Risk History Sub-Module
#'
#' @description Sub-Module function to track the behavioral risk and STI diagnosis history of
#'              HIV-negative MSM persons for purpose of PrEP indications in [`prep_msm`].
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

  gc.last.dxpos.time <- get_attr(dat, "gc.last.dxpos.time")
  ct.last.dxpos.time <- get_attr(dat, "ct.last.dxpos.time")

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
  last.dx.time <- pmax(gc.last.dxpos.time, ct.last.dxpos.time, na.rm = TRUE)
  dat <- set_attr(dat, "prep.ind.sti", last.dx.time)

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
