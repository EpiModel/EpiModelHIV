#' @title STI Testing : PAF
#'
#' @description Module function for STI screening of asymptomatic persons.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This testing module supports two testing parameterizations, input via the
#' \code{testing.pattern} parameter: memoryless for stochastic and
#' geometrically-distributed waiting times to test (constant hazard); and
#' interval for deterministic tested after defined waiting time intervals.
#' Symptomatic testing is handled in the STI treatment module.
#'
#' @return
#' This function returns the \code{dat} object with updated
#' \code{last.neg.test}, \code{diag.status} and \code{diag.time} attributes for
#' each STI.
#'
#' @keywords module msm
#'
#' @export
sti_test_msm_paf <- function(dat, at) {

  # risk history pre-module
  dat <- riskhist_stitest_msm(dat, at)

  # 1. Setup ----------------------------------------------------------------

  # Attributes
  tt.traj <- get_attr(dat, "tt.traj")
  tt.traj.gc.hivpos <- get_attr(dat, "tt.traj.gc.hivpos")
  tt.traj.gc.hivneg <- get_attr(dat, "tt.traj.gc.hivneg")
  tt.traj.ct.hivpos <- get_attr(dat, "tt.traj.ct.hivpos")
  tt.traj.ct.hivneg <- get_attr(dat, "tt.traj.ct.hivneg")
  diag.status.gc <- get_attr(dat, "diag.status.gc")
  diag.status.ct <- get_attr(dat, "diag.status.ct")
  diag.status <- get_attr(dat, "diag.status")
  status <- get_attr(dat, "status")
  rGC <- get_attr(dat, "rGC")
  uGC <- get_attr(dat, "uGC")
  rCT <- get_attr(dat, "rCT")
  uCT <- get_attr(dat, "uCT")
  role.class <- get_attr(dat, "role.class")
  last.neg.test.rgc <- get_attr(dat, "last.neg.test.rgc")
  last.neg.test.ugc <- get_attr(dat, "last.neg.test.ugc")
  last.neg.test.rct <- get_attr(dat, "last.neg.test.rct")
  last.neg.test.uct <- get_attr(dat, "last.neg.test.uct")
  last.diag.time.gc <- get_attr(dat, "last.diag.time.gc")
  last.diag.time.ct <- get_attr(dat, "last.diag.time.ct")

  tsinceltst.rgc <- 1 + get_attr(dat, "time.since.last.test.rgc")
  tsinceltst.ugc <- 1 + get_attr(dat, "time.since.last.test.ugc")
  tsinceltst.rct <- 1 + get_attr(dat, "time.since.last.test.rct")
  tsinceltst.uct <- 1 + get_attr(dat, "time.since.last.test.uct")
  tsinceltst.gc <- pmin(tsinceltst.rgc, tsinceltst.ugc)
  tsinceltst.ct <- pmin(tsinceltst.rct, tsinceltst.uct)

  prepStat <- get_attr(dat, "prepStat")
  stitestind1 <- get_attr(dat, "stitest.ind.active")
  stitestind2 <- get_attr(dat, "stitest.ind.recentpartners")

  # Parameters
  stianntest.ct.hivneg.coverage <- get_param(dat, "stianntest.ct.hivneg.coverage")
  stihighrisktest.ct.hivneg.coverage <- get_param(dat, "stihighrisktest.ct.hivneg.coverage")

      #only ct here because used for overall gcct measure below
  stianntest.ct.hivpos.coverage <- get_param(dat, "stianntest.ct.hivpos.coverage")
  stihighrisktest.ct.hivpos.coverage <- get_param(dat, "stihighrisktest.ct.hivpos.coverage")

  testing.pattern.sti <- get_param(dat, "testing.pattern.sti")
  stitest.active.int <- get_param(dat, "stitest.active.int")
  sti.highrisktest.int <- get_param(dat, "sti.highrisktest.int")
  tst.rect.sti.rr.hivneg <- get_param(dat, "tst.rect.sti.rr.hivneg")
  tst.rect.sti.rr.hivpos <- get_param(dat, "tst.rect.sti.rr.hivpos")

  # Eligibility and trajectory
  # Annual indications- sexually active in last year
  idsactive.hivpos <- which(stitestind1 == 1 & diag.status == 1)# & tt.traj %in% 3:4)
  idsactive.hivneg <- setdiff((which(stitestind1 == 1)), idsactive.hivpos)

  # STI testing higher-risk eligibility scenarios
  idshighrisk.hivpos <- which(stitestind2 == 1 & diag.status == 1)# & tt.traj %in% 3:4)
  idshighrisk.hivneg <- setdiff((which(stitestind2 == 1)), idshighrisk.hivpos)


  # 2. Indication Trajectory Stoppage ---------------------------------------

  # Reduce testing trajectory to NA if no longer indicated for more frequent high-risk testing
  idsnothighriskelig.hivpos <- which((tt.traj.gc.hivpos == 2 |
                                        tt.traj.ct.hivpos == 2) & stitestind2 != 1)
  tt.traj.gc.hivpos[idsnothighriskelig.hivpos] <- NA
  tt.traj.ct.hivpos[idsnothighriskelig.hivpos] <- NA

  idsnothighriskelig.hivneg <- which(
    (tt.traj.gc.hivneg == 2 | tt.traj.ct.hivneg == 2)
    & stitestind2 != 1
  )

  tt.traj.gc.hivneg[idsnothighriskelig.hivneg] <- NA
  tt.traj.ct.hivneg[idsnothighriskelig.hivneg] <- NA

  # Reduce testing trajectory to NA if no longer indicated for lower-risk testing
  idsnotactiveelig.hivpos <- which(
    (tt.traj.gc.hivpos == 1 | tt.traj.ct.hivpos == 1) &
    stitestind1 != 1
  )

  tt.traj.gc.hivpos[idsnotactiveelig.hivpos] <- NA
  tt.traj.ct.hivpos[idsnotactiveelig.hivpos] <- NA

  idsnotactiveelig.hivneg <- which(
    (tt.traj.gc.hivneg == 1 | tt.traj.ct.hivneg == 1) &
    stitestind1 != 1
  )
  tt.traj.gc.hivneg[idsnotactiveelig.hivneg] <- NA
  tt.traj.ct.hivneg[idsnotactiveelig.hivneg] <- NA

  # 3. Screening for non-HIV Diagnosed --------------------------------------

  ## 3a. High-Risk Testing ##

  ## Assume coverage and people are same for NG and CT - correlation
  idsEligSt <- idshighrisk.hivneg
  nEligSt <- length(idshighrisk.hivneg)

  ## Evaluate existing coverage
  stihighrisktestCov.ct <- sum(tt.traj.ct.hivneg == 2, na.rm = TRUE) / nEligSt
  stihighrisktestCov.ct <- ifelse(is.na(stihighrisktestCov.ct), 0, stihighrisktestCov.ct)
  stihighrisktestCov.gc <- sum(tt.traj.gc.hivneg == 2, na.rm = TRUE) / nEligSt
  stihighrisktestCov.gc <- ifelse(is.na(stihighrisktestCov.gc), 0, stihighrisktestCov.gc)

  ## Count how many are eligible to start a new trajectory
  nStart.gcct <- max(0, min(
    nEligSt,
    round(
      (stihighrisktest.ct.hivneg.coverage - stihighrisktestCov.ct) *
      length(idsEligSt)
    )
  ))

  ## Sample individuals
  idsStart.gcct <- numeric()
  if (nStart.gcct > 0) {
    idsStart.gcct <- ssample(idsEligSt, nStart.gcct)
  }

  ## Update testing trajectory for higher-risk
  if (length(idsStart.gcct) > 0) {
    tt.traj.gc.hivneg[idsStart.gcct] <- 2
    tt.traj.ct.hivneg[idsStart.gcct] <- 2
  }

  ## 3b. Sexually Active (non-HR) Testing ##

  ## Assume coverage and people are same for NG and CT - correlation
  idsEligSt.gcct <- setdiff(idsactive.hivneg, which(tt.traj.ct.hivneg == 2))
  nEligSt.gcct <- length(idsEligSt.gcct)

  ## Evaluate existing coverage
  stianntestCov.ct <- sum(tt.traj.ct.hivneg == 1, na.rm = TRUE) / nEligSt.gcct
  stianntestCov.ct <- ifelse(is.na(stianntestCov.ct), 0, stianntestCov.ct)
  stianntestCov.gc <- sum(tt.traj.gc.hivneg == 1, na.rm = TRUE) / nEligSt.gcct
  stianntestCov.gc <- ifelse(is.na(stianntestCov.gc), 0, stianntestCov.gc)

  ## Count how many are eligible to start a new trajectory
  nStart.gcct <- max(0, min(
    nEligSt.gcct,
    round((stianntest.ct.hivneg.coverage - stianntestCov.ct) * nEligSt.gcct)
  ))

  ## Sample individuals
  idsStart.gcct <- numeric()
  if (nStart.gcct > 0) {
    idsStart.gcct <- ssample(idsEligSt.gcct, nStart.gcct)
  }

  ## Update testing trajectory for lower-risk
  if (length(idsStart.gcct) > 0) {
    tt.traj.ct.hivneg[idsStart.gcct] <- 1
    tt.traj.gc.hivneg[idsStart.gcct] <- 1
  }

  ## 3c. Asymptomatic Screening ##

  # GC
  if (testing.pattern.sti == "interval" ) {
    tst.gc.annual.interval <- which(
      tt.traj.gc.hivneg == 1 &
      (diag.status.gc == 0 | is.na(diag.status.gc)) &
      (tsinceltst.ugc >= stitest.active.int |
      tsinceltst.uct >= stitest.active.int) &
      prepStat == 0
    )

    tst.gc.highrisk.interval <- which(
      tt.traj.gc.hivneg == 2 &
      (diag.status.gc == 0 | is.na(diag.status.gc)) &
      (tsinceltst.ugc >= sti.highrisktest.int |
      tsinceltst.uct >= sti.highrisktest.int) &
      prepStat == 0
    )

    tst.gc.nprep.hivneg <- c(tst.gc.annual.interval, tst.gc.highrisk.interval)
  }

  # CT
  if (testing.pattern.sti == "interval" ) {
    tst.ct.annual.interval <- which(
      tt.traj.ct.hivneg == 1 &
      (diag.status.ct == 0 | is.na(diag.status.ct)) &
      (tsinceltst.ugc >= stitest.active.int |
      tsinceltst.uct >= stitest.active.int) &
      prepStat == 0
    )

    tst.ct.highrisk.interval <- which(
      tt.traj.ct.hivneg == 2 &
      (diag.status.ct == 0 | is.na(diag.status.ct)) &
      (tsinceltst.ugc >= sti.highrisktest.int |
      tsinceltst.uct >= sti.highrisktest.int) &
      prepStat == 0
    )

    tst.ct.nprep.hivneg <- c(tst.ct.annual.interval, tst.ct.highrisk.interval)
  }

  # GC non-PrEP testing
  tst.rgc.hivneg <- tst.gc.nprep.hivneg
  tst.ugc.hivneg <- tst.gc.nprep.hivneg
  tst.rgc.hivneg <- sample(
    tst.rgc.hivneg,
    tst.rect.sti.rr.hivneg * length(tst.rgc.hivneg)
  )

  tst.rgc.pos.hivneg <- tst.rgc.hivneg[rGC[tst.rgc.hivneg] == 1]
  tst.ugc.pos.hivneg <- tst.ugc.hivneg[uGC[tst.ugc.hivneg] == 1]
  tst.rgc.neg.hivneg <- setdiff(tst.rgc.hivneg, tst.rgc.pos.hivneg)
  tst.ugc.neg.hivneg <- setdiff(tst.ugc.hivneg, tst.ugc.pos.hivneg)
  tst.gc.pos.hivneg <- unique(c(tst.rgc.pos.hivneg, tst.ugc.pos.hivneg))

  # CT non-PrEP testing
  tst.rct.hivneg <- tst.uct.hivneg <- tst.ct.nprep.hivneg
  tst.rct.hivneg <- sample(
    tst.rct.hivneg,
    tst.rect.sti.rr.hivneg * length(tst.rct.hivneg)
  )

  tst.rct.pos.hivneg <- tst.rct.hivneg[rCT[tst.rct.hivneg] == 1]
  tst.uct.pos.hivneg <- tst.uct.hivneg[uCT[tst.uct.hivneg] == 1]
  tst.rct.neg.hivneg <- setdiff(tst.rct.hivneg, tst.rct.pos.hivneg)
  tst.uct.neg.hivneg <- setdiff(tst.uct.hivneg, tst.uct.pos.hivneg)
  tst.ct.pos.hivneg <- unique(c(tst.rct.pos.hivneg, tst.uct.pos.hivneg))

  # GC Attributes
  last.neg.test.rgc[tst.rgc.neg.hivneg] <- at
  last.neg.test.ugc[tst.ugc.neg.hivneg] <- at
  last.neg.test.rgc[tst.rgc.pos.hivneg] <- NA
  last.neg.test.ugc[tst.ugc.pos.hivneg] <- NA
  diag.status.gc[tst.gc.pos.hivneg] <- 1
  last.diag.time.gc[tst.gc.pos.hivneg] <- at
  tsinceltst.rgc[tst.rgc.hivneg] <- 0
  tsinceltst.ugc[tst.ugc.hivneg] <- 0

  # CT Attributes
  last.neg.test.rct[tst.rct.neg.hivneg] <- at
  last.neg.test.uct[tst.uct.neg.hivneg] <- at
  last.neg.test.rct[tst.rct.pos.hivneg] <- NA
  last.neg.test.uct[tst.uct.pos.hivneg] <- NA
  diag.status.ct[tst.ct.pos.hivneg] <- 1
  last.diag.time.ct[tst.ct.pos.hivneg] <- at
  tsinceltst.rct[tst.rct.hivneg] <- 0
  tsinceltst.uct[tst.uct.hivneg] <- 0


  # 4. Screening for HIV-Diagnosed ------------------------------------------

  ## 4a. High-Risk Testing ##

  ## Assume coverage and people are same for NG and CT - correlation
  idsEligSt <- idshighrisk.hivpos
  nEligSt <- length(idshighrisk.hivpos)

  ## Evaluate existing coverage
  stihighrisktestCov.ct <- sum(tt.traj.ct.hivpos == 2, na.rm = TRUE) / nEligSt
  stihighrisktestCov.ct <- ifelse(
    is.na(stihighrisktestCov.ct),
    0, stihighrisktestCov.ct
  )
  stihighrisktestCov.gc <- sum(tt.traj.gc.hivpos == 2, na.rm = TRUE) / nEligSt
  stihighrisktestCov.gc <- ifelse(
    is.na(stihighrisktestCov.gc),
    0, stihighrisktestCov.gc
  )

  ## Count how many are eligible to start a new trajectory
  nStart.gcct <- max(0, min(
    nEligSt, round(
      (stihighrisktest.ct.hivpos.coverage - stihighrisktestCov.ct) *
      length(idshighrisk.hivpos)
    )
  ))

  ## Sample individuals
  idsStart.gcct <- numeric()
  if (nStart.gcct > 0) {
    idsStart.gcct <- ssample(idsEligSt, nStart.gcct)
  }

  ## Update testing trajectory for higher-risk
  if (length(idsStart.gcct) > 0) {
    tt.traj.ct.hivpos[idsStart.gcct] <- 2
    tt.traj.gc.hivpos[idsStart.gcct] <- 2
  }

  ## 4b. Sexually Active (non-HR) Testing ##

  ## Assume coverage and people are same for NG and CT - correlation
  idsEligSt.gcct <- setdiff(idsactive.hivpos, which(tt.traj.ct.hivpos == 2))
  nEligSt.gcct <- length(idsEligSt.gcct)

  ## Evaluate existing coverage
  stianntestCov.ct <- sum(tt.traj.ct.hivpos == 1, na.rm = TRUE) / nEligSt.gcct
  stianntestCov.ct <- ifelse(is.na(stianntestCov.ct), 0, stianntestCov.ct)
  stianntestCov.gc <- sum(tt.traj.gc.hivpos == 1, na.rm = TRUE) / nEligSt.gcct
  stianntestCov.gc <- ifelse(is.na(stianntestCov.gc), 0, stianntestCov.gc)

  ## Count how many are eligible to start a new trajectory
  nStart.gcct <- max(0, min(
    nEligSt.gcct, round(
      (stianntest.ct.hivpos.coverage - stianntestCov.ct) * nEligSt.gcct
    )
  ))

  ## Sample individuals
  idsStart.gcct <- numeric()
  if (nStart.gcct > 0) {
    idsStart.gcct <- ssample(idsEligSt.gcct, nStart.gcct)
  }

  ## Update testing trajectory for lower-risk
  if (length(idsStart.gcct) > 0) {
    tt.traj.ct.hivpos[idsStart.gcct] <- 1
    tt.traj.gc.hivpos[idsStart.gcct] <- 1
  }

  ## 4c. Asymptomatic screening ##

  ## GC
  if (testing.pattern.sti == "interval" ) {
    tst.gc.annual.interval <- which(
      tt.traj.gc.hivpos == 1 &
      (diag.status.gc == 0 | is.na(diag.status.gc)) &
      (tsinceltst.ugc >= stitest.active.int |
      tsinceltst.uct >= stitest.active.int) &
      prepStat == 0
    )

    tst.gc.highrisk.interval <- which(
      tt.traj.gc.hivpos == 2 &
      (diag.status.gc == 0 | is.na(diag.status.gc)) &
      (tsinceltst.ugc >= sti.highrisktest.int |
      tsinceltst.uct >= sti.highrisktest.int) &
      prepStat == 0
    )

    tst.gc.nprep.hivpos <- c(tst.gc.annual.interval, tst.gc.highrisk.interval)
  }

  ## CT
  if (testing.pattern.sti == "interval" ) {
    tst.ct.annual.interval <- which(
      tt.traj.ct.hivpos == 1 &
      (diag.status.ct == 0 | is.na(diag.status.ct)) &
      (tsinceltst.ugc >= stitest.active.int |
      tsinceltst.uct >= stitest.active.int) &
      prepStat == 0
    )

    tst.ct.highrisk.interval <- which(
      tt.traj.ct.hivpos == 2 &
      (diag.status.ct == 0 | is.na(diag.status.ct)) &
      (tsinceltst.ugc >= sti.highrisktest.int |
        tsinceltst.uct >= sti.highrisktest.int) &
      prepStat == 0
    )

    tst.ct.nprep.hivpos <- c(tst.ct.annual.interval, tst.ct.highrisk.interval)
  }

  # GC non-PrEP testing
  tst.rgc.hivpos <- tst.gc.nprep.hivpos
  tst.ugc.hivpos <- tst.gc.nprep.hivpos
  tst.rgc.hivpos <- sample(
    tst.rgc.hivpos,
    tst.rect.sti.rr.hivpos * length(tst.rgc.hivpos)
  )

  tst.rgc.pos.hivpos <- tst.rgc.hivpos[rGC[tst.rgc.hivpos] == 1]
  tst.ugc.pos.hivpos <- tst.ugc.hivpos[uGC[tst.ugc.hivpos] == 1]
  tst.rgc.neg.hivpos <- setdiff(tst.rgc.hivpos, tst.rgc.pos.hivpos)
  tst.ugc.neg.hivpos <- setdiff(tst.ugc.hivpos, tst.ugc.pos.hivpos)
  tst.gc.pos.hivpos <- unique(c(tst.rgc.pos.hivpos, tst.ugc.pos.hivpos))

  # CT non-PrEP testing
  tst.rct.hivpos <- tst.uct.hivpos <- tst.ct.nprep.hivpos
  tst.rct.hivpos <- sample(
    tst.rct.hivpos,
    tst.rect.sti.rr.hivpos * length(tst.rct.hivpos)
  )

  tst.rct.pos.hivpos <- tst.rct.hivpos[rCT[tst.rct.hivpos] == 1]
  tst.uct.pos.hivpos <- tst.uct.hivpos[uCT[tst.uct.hivpos] == 1]
  tst.rct.neg.hivpos <- setdiff(tst.rct.hivpos, tst.rct.pos.hivpos)
  tst.uct.neg.hivpos <- setdiff(tst.uct.hivpos, tst.uct.pos.hivpos)
  tst.ct.pos.hivpos <- unique(c(tst.rct.pos.hivpos, tst.uct.pos.hivpos))

  # GC Attributes
  last.neg.test.rgc[tst.rgc.neg.hivpos] <- at
  last.neg.test.ugc[tst.ugc.neg.hivpos] <- at
  last.neg.test.rgc[tst.rgc.pos.hivpos] <- NA
  last.neg.test.ugc[tst.ugc.pos.hivpos] <- NA
  diag.status.gc[tst.gc.pos.hivpos] <- 1
  last.diag.time.gc[tst.gc.pos.hivpos] <- at
  tsinceltst.rgc[tst.rgc.hivpos] <- 0
  tsinceltst.ugc[tst.ugc.hivpos] <- 0

  # CT Attributes
  last.neg.test.rct[tst.rct.neg.hivpos] <- at
  last.neg.test.uct[tst.uct.neg.hivpos] <- at
  last.neg.test.rct[tst.rct.pos.hivpos] <- NA
  last.neg.test.uct[tst.uct.pos.hivpos] <- NA
  diag.status.ct[tst.ct.pos.hivpos] <- 1
  last.diag.time.ct[tst.ct.pos.hivpos] <- at
  tsinceltst.rct[tst.rct.hivpos] <- 0
  tsinceltst.uct[tst.uct.hivpos] <- 0


  ## Output -----------------------------------------------------------------
  #
  # GC
  #    HIV-pos
  dat <- set_epi(dat, "rGC.tot.test.hivpos", at, length(tst.rgc.hivpos))
  dat <- set_epi(dat, "uGC.tot.test.hivpos", at, length(tst.ugc.hivpos))
  dat <- set_epi(dat, "rGC.pos.test.hivpos", at, length(tst.rgc.pos.hivpos))
  dat <- set_epi(dat, "uGC.pos.test.hivpos", at, length(tst.ugc.pos.hivpos))
  #    HIV-neg
  dat <- set_epi(dat, "rGC.tot.test.hivneg", at, length(tst.rgc.hivneg))
  dat <- set_epi(dat, "uGC.tot.test.hivneg", at, length(tst.ugc.hivneg))
  dat <- set_epi(dat, "rGC.pos.test.hivneg", at, length(tst.rgc.pos.hivneg))
  dat <- set_epi(dat, "uGC.pos.test.hivneg", at, length(tst.ugc.pos.hivneg))

  # CT
  #    HIV-pos
  dat <- set_epi(dat, "rCT.tot.test.hivpos", at, length(tst.rct.hivpos))
  dat <- set_epi(dat, "uCT.tot.test.hivpos", at, length(tst.uct.hivpos))
  dat <- set_epi(dat, "rCT.pos.test.hivpos", at, length(tst.rct.pos.hivpos))
  dat <- set_epi(dat, "uCT.pos.test.hivpos", at, length(tst.uct.pos.hivpos))
  #    HIV-neg
  dat <- set_epi(dat, "rCT.tot.test.hivneg", at, length(tst.rct.hivneg))
  dat <- set_epi(dat, "uCT.tot.test.hivneg", at, length(tst.uct.hivneg))
  dat <- set_epi(dat, "rCT.pos.test.hivneg", at, length(tst.rct.pos.hivneg))
  dat <- set_epi(dat, "uCT.pos.test.hivneg", at, length(tst.uct.pos.hivneg))

  # Screening episodes
  dat <- set_epi(dat, "sti.screening.ep.hivpos", at,
    length(unique(c(
      tst.rgc.hivpos,
      tst.ugc.hivpos,
      tst.rct.hivpos,
      tst.uct.hivpos
    )))
  )
  dat <- set_epi(dat, "sti.screening.ep.hivneg", at,
    length(unique(c(
      tst.rgc.hivneg,
      tst.ugc.hivneg,
      tst.rct.hivneg,
      tst.uct.hivneg
    )))
  )

  ## Output by true HIV status ----------------------------------------------
  # GC
  # tot test
  #  rgc
  rgc.tot.tests <- union(tst.rgc.hivpos, tst.rgc.hivneg)
  rgc.tot.tests.hivpos <- intersect(rgc.tot.tests, which(status == 1))
  rgc.tot.tests.hivneg <- intersect(rgc.tot.tests, which(status == 0))
  dat <- set_epi(dat, "rgc.tot.test.hivpos", at, length(rgc.tot.tests.hivpos))
  dat <- set_epi(dat, "rgc.tot.test.hivneg", at, length(rgc.tot.tests.hivneg))
  #  ugc
  ugc.tot.tests <- union(tst.ugc.hivpos, tst.ugc.hivneg)
  ugc.tot.tests.hivpos <- intersect(ugc.tot.tests, which(status == 1))
  ugc.tot.tests.hivneg <- intersect(ugc.tot.tests, which(status == 0))
  dat <- set_epi(dat, "ugc.tot.test.hivpos", at, length(ugc.tot.tests.hivpos))
  dat <- set_epi(dat, "ugc.tot.test.hivneg", at, length(ugc.tot.tests.hivneg))
  # pos test
  #  rgc
  rgc.pos.tests <- union(tst.rgc.pos.hivpos, tst.rgc.pos.hivneg)
  rgc.pos.tests.hivpos <- intersect(rgc.pos.tests, which(status == 1))
  rgc.pos.tests.hivneg <- intersect(rgc.pos.tests, which(status == 0))
  dat <- set_epi(dat, "rgc.pos.test.hivpos", at, length(rgc.pos.tests.hivpos))
  dat <- set_epi(dat, "rgc.pos.test.hivneg", at, length(rgc.pos.tests.hivneg))
  #  ugc
  ugc.pos.tests <- union(tst.ugc.pos.hivpos, tst.ugc.pos.hivneg)
  ugc.pos.tests.hivpos <- intersect(ugc.pos.tests, which(status == 1))
  ugc.pos.tests.hivneg <- intersect(ugc.pos.tests, which(status == 0))
  dat <- set_epi(dat, "ugc.pos.test.hivpos", at, length(ugc.pos.tests.hivpos))
  dat <- set_epi(dat, "ugc.pos.test.hivneg", at, length(ugc.pos.tests.hivneg))

  # CT
  # tot test
  #  rct
  rct.tot.tests <- union(tst.rct.hivpos, tst.rct.hivneg)
  rct.tot.tests.hivpos <- intersect(rct.tot.tests, which(status == 1))
  rct.tot.tests.hivneg <- intersect(rct.tot.tests, which(status == 0))
  dat <- set_epi(dat, "rct.tot.test.hivpos", at, length(rct.tot.tests.hivpos))
  dat <- set_epi(dat, "rct.tot.test.hivneg", at, length(rct.tot.tests.hivneg))
  #  uct
  uct.tot.tests <- union(tst.uct.hivpos, tst.uct.hivneg)
  uct.tot.tests.hivpos <- intersect(uct.tot.tests, which(status == 1))
  uct.tot.tests.hivneg <- intersect(uct.tot.tests, which(status == 0))
  dat <- set_epi(dat, "uct.tot.test.hivpos", at, length(uct.tot.tests.hivpos))
  dat <- set_epi(dat, "uct.tot.test.hivneg", at, length(uct.tot.tests.hivneg))
  # pos test
  #  rct
  rct.pos.tests <- union(tst.rct.pos.hivpos, tst.rct.pos.hivneg)
  rct.pos.tests.hivpos <- intersect(rct.pos.tests, which(status == 1))
  rct.pos.tests.hivneg <- intersect(rct.pos.tests, which(status == 0))
  dat <- set_epi(dat, "rct.pos.test.hivpos", at, length(rct.pos.tests.hivpos))
  dat <- set_epi(dat, "rct.pos.test.hivneg", at, length(rct.pos.tests.hivneg))
  #  uct
  uct.pos.tests <- union(tst.uct.pos.hivpos, tst.uct.pos.hivneg)
  uct.pos.tests.hivpos <- intersect(uct.pos.tests, which(status == 1))
  uct.pos.tests.hivneg <- intersect(uct.pos.tests, which(status == 0))
  dat <- set_epi(dat, "uct.pos.test.hivpos", at, length(uct.pos.tests.hivpos))
  dat <- set_epi(dat, "uct.pos.test.hivneg", at, length(uct.pos.tests.hivneg))


  # dat$epi$stiscreeningepisodes[at] <- length(unique(c(tst.rgc.hivneg, tst.rgc.hivpos, tst.ugc.hivneg, tst.ugc.hivpos,
  #                                                     tst.rct.hivneg, tst.rct.hivpos, tst.uct.hivneg, tst.uct.hivpos)))

  # # Number of people on each testing trajectory by serostatus
  # dat$epi$tt.traj.gc1.hivneg[at] <- length(which(tt.traj.gc.hivneg == 1))
  # dat$epi$tt.traj.ct1.hivneg[at] <- length(which(tt.traj.ct.hivneg == 1))
  # dat$epi$tt.traj.gc2.hivneg[at] <- length(which(tt.traj.gc.hivneg == 2))
  # dat$epi$tt.traj.ct2.hivneg[at] <- length(which(tt.traj.ct.hivneg == 2))
  # dat$epi$tt.traj.gc1.hivpos[at] <- length(which(tt.traj.gc.hivpos == 1)) dat$epi$tt.traj.ct1.hivpos[at] <- length(which(tt.traj.ct.hivpos == 1))
  # dat$epi$tt.traj.gc2.hivpos[at] <- length(which(tt.traj.gc.hivpos == 2))
  # dat$epi$tt.traj.ct2.hivpos[at] <- length(which(tt.traj.ct.hivpos == 2))

  # # Number of people on each testing trajectory
  # dat$epi$tt.traj.gc1[at] <- length(which(tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1))
  # dat$epi$tt.traj.ct1[at] <- length(which(tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1))
  # dat$epi$tt.traj.gc2[at] <- length(which(tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2))
  # dat$epi$tt.traj.ct2[at] <- length(which(tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2))

  # dat$epi$tt.traj.sti1[at] <- length(which(tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1 |
  #                                            tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1))
  # dat$epi$tt.traj.sti2[at] <- length(which(tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2 |
  #                                            tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2))

  # # Number of screening tests and number of positive results on screening tests
  # # GC overall
  # dat$epi$rGCasympttests[at] <- length(tst.rgc.hivpos) +
  #   length(tst.rgc.hivneg)
  # dat$epi$uGCasympttests[at] <- length(tst.ugc.hivpos) +
  #   length(tst.ugc.hivneg)
  # dat$epi$GCasympttests[at] <- length(tst.rgc.hivpos) + length(tst.ugc.hivpos) +
  #   length(tst.rgc.hivneg) + length(tst.ugc.hivneg)

  # # GC positive tests
  # dat$epi$rGCasympttests.pos[at] <- length(tst.rgc.pos.hivpos) +
  #   length(tst.rgc.pos.hivneg)
  # dat$epi$uGCasympttests.pos[at] <- length(tst.ugc.pos.hivpos) +
  #   length(tst.ugc.pos.hivneg)
  # dat$epi$GCasympttests.pos[at] <- length(tst.rgc.pos.hivpos) + length(tst.ugc.pos.hivpos) +
  #   length(tst.rgc.pos.hivneg) + length(tst.ugc.pos.hivneg)

  # # CT overall
  # dat$epi$rCTasympttests[at] <- length(tst.rct.hivpos) +
  #   length(tst.rct.hivneg)
  # dat$epi$uCTasympttests[at] <- length(tst.uct.hivpos) +
  #   length(tst.uct.hivneg)
  # dat$epi$CTasympttests[at] <- length(tst.rct.hivpos) + length(tst.uct.hivpos) +
  #   length(tst.rct.hivneg) + length(tst.uct.hivneg)

  # # CT positive tests
  # dat$epi$rCTasympttests.pos[at] <- length(tst.rct.pos.hivpos) +
  #   length(tst.rct.pos.hivneg)
  # dat$epi$uCTasympttests.pos[at] <- length(tst.uct.pos.hivpos) +
  #   length(tst.uct.pos.hivneg)
  # dat$epi$CTasympttests.pos[at] <- length(tst.rct.pos.hivpos) +
  #   length(tst.uct.pos.hivpos) + length(tst.rct.pos.hivneg) +
  #   length(tst.uct.pos.hivneg)

  # # Overall STI tests and positive STI tests
  # dat$epi$stiasympttests[at] <- length(tst.rgc.hivpos) + length(tst.ugc.hivpos) +
  #   length(tst.rct.hivpos) + length(tst.uct.hivpos) +
  #   length(tst.rgc.hivneg) + length(tst.ugc.hivneg) +
  #   length(tst.rct.hivneg) + length(tst.uct.hivneg)
  # dat$epi$stiasympttests.pos[at] <- length(tst.rgc.pos.hivpos) + length(tst.ugc.pos.hivpos) +
  #   length(tst.rct.pos.hivpos) + length(tst.uct.pos.hivpos) +
  #   length(tst.rgc.pos.hivneg) + length(tst.ugc.pos.hivneg) +
  #   length(tst.rct.pos.hivneg) + length(tst.uct.pos.hivneg)

  # # Screening events for cost-effectiveness
  # dat$epi$stiscreeningepisodes[at] <- length(unique(c(tst.rgc.hivneg, tst.rgc.hivpos, tst.ugc.hivneg, tst.ugc.hivpos,
  #                                                     tst.rct.hivneg, tst.rct.hivpos, tst.uct.hivneg, tst.uct.hivpos)))

  # # Risk group-specific test counters
  # # GC lower-risk
  # dat$epi$rGCasympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 1)) +
  #   length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 1))
  # dat$epi$uGCasympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 1)) +
  #   length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 1))
  # dat$epi$GCasympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 1)) +
  #   length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 1)) +
  #   length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 1)) +
  #   length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 1))

  # # GC higher-risk
  # dat$epi$rGCasympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 2)) +
  #   length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 2))
  # dat$epi$uGCasympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 2)) +
  #   length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 2))
  # dat$epi$GCasympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 2)) +
  #   length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 2)) +
  #   length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 2)) +
  #   length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 2))

  # # CT lower-risk
  # dat$epi$rCTasympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 1)) +
  #   length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 1))
  # dat$epi$uCTasympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 1)) +
  #   length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 1))
  # dat$epi$CTasympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 1)) +
  #   length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 1)) +
  #   length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 1)) +
  #   length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 1))

  # # CT higher-risk
  # dat$epi$rCTasympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 2)) +
  #   length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 2))
  # dat$epi$uCTasympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 2)) +
  #   length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 2))
  # dat$epi$CTasympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 2)) +
  #   length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 2)) +
  #   length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 2)) +
  #   length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 2))

  # # STI
  # dat$epi$stiasympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 1)) +
  #   length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 1)) +
  #   length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 1)) +
  #   length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 1)) +
  #   length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 1)) +
  #   length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 1)) +
  #   length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 1)) +
  #   length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 1))

  # dat$epi$stiasympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 2)) +
  #   length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 2)) +
  #   length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 2)) +
  #   length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 2)) +
  #   length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 2)) +
  #   length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 2)) +
  #   length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 2)) +
  #   length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 2))

  # Update Attributes
  # Stoppage attributes
  dat <- set_attr(dat, "stihighrisktestLastElig", at, posit_ids = idsnothighriskelig.hivneg)
  dat <- set_attr(dat, "stianntestLastElig", at, posit_ids = idsnotactiveelig.hivneg)
  dat <- set_attr(dat, "stihighrisktestLastElig", at, posit_ids = idsnothighriskelig.hivpos)
  dat <- set_attr(dat, "stianntestLastElig", at, posit_ids = idsnotactiveelig.hivpos)

  # GC Attributes
  dat <- set_attr(dat, "last.neg.test.rgc", last.neg.test.rgc)
  dat <- set_attr(dat, "last.neg.test.ugc", last.neg.test.ugc)
  dat <- set_attr(dat, "diag.status.gc", diag.status.gc)
  dat <- set_attr(dat, "last.diag.time.gc", last.diag.time.gc)
  dat <- set_attr(dat, "tt.traj.gc.hivneg", tt.traj.gc.hivneg)
  dat <- set_attr(dat, "tt.traj.gc.hivpos", tt.traj.gc.hivpos)
  dat <- set_attr(dat, "time.since.last.test.rgc", tsinceltst.rgc)
  dat <- set_attr(dat, "time.since.last.test.ugc", tsinceltst.ugc)

  # CT Attributes
  dat <- set_attr(dat, "last.neg.test.rct", last.neg.test.rct)
  dat <- set_attr(dat, "last.neg.test.uct", last.neg.test.uct)
  dat <- set_attr(dat, "diag.status.ct", diag.status.ct)
  dat <- set_attr(dat, "last.diag.time.ct", last.diag.time.ct)
  dat <- set_attr(dat, "tt.traj.ct.hivneg", tt.traj.ct.hivneg)
  dat <- set_attr(dat, "tt.traj.ct.hivpos", tt.traj.ct.hivpos)
  dat <- set_attr(dat, "time.since.last.test.rct", tsinceltst.rct)
  dat <- set_attr(dat, "time.since.last.test.uct", tsinceltst.uct)

  return(dat)
}

#' @title Risk History for STI Testing Module
#'
#' @description Module function to track the risk history of uninfected persons
#'              for purpose of STI testing prevention intervention targeting.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
riskhist_stitest_msm <- function(dat, at) {

  ## Parameters
  partnercutoff <- get_param(dat, "partnercutoff")
  stitest.active.int <- get_param(dat, "stitest.active.int")
  ept.start <- get_param(dat, "ept.start")
  stitest.start <- get_param(dat, "stitest.start")

  ## Attributes
  uid <- get_unique_ids(dat)
  active <- get_attr(dat, "active")
  stitestind1 <- get_attr(dat, "time.last.sex")

  # Indications -------------------------------------------------------------

  part.list <- dat$temp$plist

  ### Lower risk - sexually active in the last year
  idsactive <- which((at - stitestind1) <= stitest.active.int)
  idsnotactive <- setdiff(which(active == 1), idsactive)
  # these are relative ids of nodes in partner list

  ### High-risk: Have more than one sex partner in last x months
  # Reset # of partners at each time step- length of "recent" interval is drawn from interval of partner list lookback
  recentpartners <- rep(0, length(active))

  if (at >= stitest.start | at >= ept.start) {
    # For those who had partners, calculate # of occurrences in partner list
    part.count <- as.data.frame(table(part.list[, c("p1", "p2")]))

    if (nrow(part.count) > 1) {
      # Calculate # of recent partners: 0 for those not in part list, update numbers for only actives in part list
      recents <- part.count[which(part.count$Var1 %in% uid), 2]
      recentpartners[which(uid %in% part.count$Var1)] <- recents

      # Choose those who have had more than X partners in last x months
      idsrecentpartners <- which(recentpartners > partnercutoff)
      idsnotrecentpartners <- setdiff(
        which(active == 1),
        idsrecentpartners
      )

      dat <- set_attr(dat, "stitest.ind.recentpartners", 1, posit_ids = idsrecentpartners)
      dat <- set_attr(dat, "stitest.ind.recentpartners", 0, posit_ids = idsnotrecentpartners)
      dat <- set_epi(dat, "recentpartners", at, length(idsrecentpartners))

    }

  }
  ### Update STI indication attributes
  dat <- set_attr(dat, "recentpartners", recentpartners)
  dat <- set_attr(dat, "stitest.ind.active", 1, posit_ids = idsactive)
  dat <- set_attr(dat, "stitest.ind.active", 0, posit_ids = idsnotactive)

  dat <- set_epi(dat, "stiactiveind", at, length(idsactive))

  return(dat)
}

