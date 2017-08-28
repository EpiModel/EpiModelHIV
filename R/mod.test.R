
#' @title HIV Testing Module
#'
#' @description Module function for HIV diagnostic testing of infected persons.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This testing module supports two testing parameterizations, input via the
#' \code{testing.pattern} parameter: memoryless for stochastic and
#' geometrically-distributed waiting times to test (constant hazard); and
#' interval for deterministic tested after defined waiting time intervals.
#'
#' @return
#' This function returns the \code{dat} object with updated
#' \code{last.neg.test}, \code{diag.status} and \code{diag.time} attributes.
#'
#' @keywords module msm
#'
#' @export
#'
hiv_test_msm <- function(dat, at) {

  ## Variables

  # Attributes
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  tt.traj <- dat$attr$tt.traj
  status <- dat$attr$status
  inf.time <- dat$attr$inf.time

  prepStat <- dat$attr$prepStat
  prep.tst.int <- dat$param$prep.tst.int

  # Parameters
  testing.pattern <- dat$param$testing.pattern
  mean.test.B.int <- dat$param$mean.test.B.int
  mean.test.W.int <- dat$param$mean.test.W.int
  twind.int <- dat$param$test.window.int

  tsincelntst <- at - dat$attr$last.neg.test
  tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]

  ## Process

  if (testing.pattern == "memoryless") {
    elig.B <- which(race == "B" &
                      tt.traj != 1 &
                      (diag.status == 0 | is.na(diag.status)) &
                      prepStat == 0)
    rates.B <- rep(1/mean.test.B.int, length(elig.B))
    tst.B <- elig.B[rbinom(length(elig.B), 1, rates.B) == 1]

    elig.W <- which(race == "W" &
                      tt.traj != 1 &
                      (diag.status == 0 | is.na(diag.status)) &
                      prepStat == 0)
    rates.W <- rep(1/mean.test.W.int, length(elig.W))
    tst.W <- elig.W[rbinom(length(elig.W), 1, rates.W) == 1]
    tst.nprep <- c(tst.B, tst.W)
  }

  if (testing.pattern == "interval") {
    tst.B <- which(race == "B" &
                     tt.traj != 1 &
                     (diag.status == 0 | is.na(diag.status)) &
                     tsincelntst >= 2*(mean.test.B.int) &
                     prepStat == 0)

    tst.W <- which(race == "W" &
                     tt.traj != 1 &
                     (diag.status == 0 | is.na(diag.status)) &
                     tsincelntst >= 2*(mean.test.W.int) &
                     prepStat == 0)
    tst.nprep <- c(tst.B, tst.W)
  }

  # PrEP testing
  tst.prep <- which((diag.status == 0 | is.na(diag.status)) &
                      prepStat == 1 &
                      tsincelntst >= prep.tst.int)

  tst.all <- c(tst.nprep, tst.prep)

  tst.pos <- tst.all[status[tst.all] == 1 & inf.time[tst.all] <= at - twind.int]
  tst.neg <- setdiff(tst.all, tst.pos)

  # Attributes
  dat$attr$last.neg.test[tst.neg] <- at
  dat$attr$diag.status[tst.pos] <- 1
  dat$attr$diag.time[tst.pos] <- at

  # Tests
  dat$epi$hivtests.prep[at] <- length(tst.prep)
  dat$epi$hivtests.nprep[at] <- length(tst.nprep)
  dat$epi$hivtests.pos[at] <- length(tst.pos)

  return(dat)
}


#' @title STI Testing Module
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
#'
sti_test_msm <- function(dat, at) {

  ## Intervention -------------------------------------------------------------

  # Attributes
  tt.traj.ct <- dat$attr$tt.traj.ct
  tt.traj.gc <- dat$attr$tt.traj.gc
  tt.traj.syph <- dat$attr$tt.traj.syph
  diag.status.gc <- dat$attr$diag.status.gc
  diag.status.ct <- dat$attr$diag.status.ct
  diag.status.syph <- dat$attr$diag.status.syph
  syphilis <- dat$attr$syphilis
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT
  stage.syph <- dat$attr$stage.syph
  role.class <- dat$attr$role.class
  last.neg.test.rgc <- dat$attr$last.neg.test.rgc
  last.neg.test.ugc <- dat$attr$last.neg.test.ugc
  last.neg.test.rct <- dat$attr$last.neg.test.rct
  last.neg.test.uct <- dat$attr$last.neg.test.uct
  last.neg.test.syph <- dat$attr$last.neg.test.syph
  last.diag.time.gc <- dat$attr$last.diag.time.gc
  last.diag.time.ct <- dat$attr$last.diag.time.ct
  last.diag.time.syph <- dat$attr$last.diag.time.syph
  tsinceltst.syph <- dat$attr$time.since.last.test.syph + 1
  tsinceltst.rgc <- dat$attr$time.since.last.test.rgc + 1
  tsinceltst.ugc <- dat$attr$time.since.last.test.ugc + 1
  tsinceltst.rct <- dat$attr$time.since.last.test.rct + 1
  tsinceltst.uct <- dat$attr$time.since.last.test.uct + 1
  tsinceltst.gc <- pmin(tsinceltst.rgc, tsinceltst.ugc)
  tsinceltst.ct <- pmin(tsinceltst.rct, tsinceltst.uct)
  race <- dat$attr$race
  prepStat <- dat$attr$prepStat
  stitestind1 <- dat$attr$stitest.ind.active
  stitestind2 <- dat$attr$stitest.ind.recentpartners

  # Parameters
  stianntest.coverage <- dat$param$stianntest.coverage
  stianntest.cov.rate <- dat$param$stianntest.cov.rate
  stihighrisktest.coverage <- dat$param$stihighrisktest.coverage
  stihighrisktest.cov.rate <- dat$param$stihighrisktest.cov.rate
  testing.pattern.sti <- dat$param$testing.pattern.sti
  stitest.active.int <- dat$param$stitest.active.int
  sti.highrisktest.int <- dat$param$sti.highrisktest.int
  tst.rect.sti.rr <- dat$param$tst.rect.sti.rr

  # Eligibility and trajectory
  # Base eligibility
  idsEligTest <- which(race %in% c("B", "W"))

  # Annual indications- sexually active in last year
  idsactive <- which(stitestind1 == 1)

  # STI testing higher-risk eligibility scenarios
  idshighrisk <- which(stitestind2 == 1)

  ## Stoppage (tt.traj.gc/.ct/.syph <- NA ------------------------------------
  # Reduce testing trajectory to NA if no longer indicated for more frequent high-risk testing
  idsnothighriskelig <- which(tt.traj.syph == 2 & stitestind2 != 1)
  tt.traj.syph[idsnothighriskelig] <- tt.traj.gc[idsnothighriskelig] <- tt.traj.ct[idsnothighriskelig] <- NA

  # Reduce testing trajectory to NA if no longer indicated for lower-risk testing
  idsnotactiveelig <- which(tt.traj.syph == 1 & stitestind1 != 1)
  tt.traj.syph[idsnotactiveelig] <- tt.traj.gc[idsnotactiveelig] <- tt.traj.ct[idsnotactiveelig] <- NA

  ## Initiation -------------------------------------------------------------

  ### Testing coverage for high risk
  stihighrisktestCov <- sum(tt.traj.ct == 2, na.rm = TRUE) / length(idshighrisk)
  stihighrisktestCov <- ifelse(is.nan(stihighrisktestCov), 0, stihighrisktestCov)

  idsEligSt <- idshighrisk
  nEligSt <- length(idshighrisk)

  nStart <- max(0, min(nEligSt, round((stihighrisktest.coverage - stihighrisktestCov) *
                                        length(idshighrisk))))
  idsStart <- NULL
  if (nStart > 0) {
    if (stihighrisktest.cov.rate >= 1) {
      idsStart <- ssample(idsEligSt, nStart)
    } else {
      idsStart <- idsEligSt[rbinom(nStart, 1, stihighrisktest.cov.rate) == 1]
    }
  }

  ## Update testing trajectory for higher-risk
  if (length(idsStart) > 0) {
    tt.traj.syph[idsStart] <- tt.traj.gc[idsStart] <- tt.traj.ct[idsStart] <- 2
  }

  ### Testing coverage for annual - all those sexually active without high-risk indications
  stianntestCov <- sum(tt.traj.ct == 1, na.rm = TRUE) / length(setdiff(idsactive, which(tt.traj.syph == 2)))
  stianntestCov <- ifelse(is.nan(stianntestCov), 0, stianntestCov)

  idsEligSt <- setdiff(idsactive, idshighrisk)
  nEligSt <- length(setdiff(idsactive, idshighrisk))

  nStart <- max(0, min(nEligSt, round((stianntest.coverage - stianntestCov) *
                                        length(setdiff(idsactive, idshighrisk)))))
  idsStart <- NULL
  if (nStart > 0) {
    if (stianntest.cov.rate >= 1) {
      idsStart <- ssample(idsEligSt, nStart)
    } else {
      idsStart <- idsEligSt[rbinom(nStart, 1, stianntest.cov.rate) == 1]
    }
  }

  ## Update testing trajectory for lower-risk
  if (length(idsStart) > 0) {
    tt.traj.syph[idsStart] <- tt.traj.gc[idsStart] <- tt.traj.ct[idsStart] <- 1
  }

  # Testing Rates by serostatus/race?
  # All MSM with HIV infection entering care should be screened for GC and CT
  # ct appropriate anatomic sites of exposure, as well as for syphilis
  # For sexually active individuals, screen at first HIV evaluation,
  # and at least annually thereafter

  ## Process for asymptomatic syphilis screening
  if (testing.pattern.sti == "memoryless") {
    elig.syph.ann <- which(tt.traj.syph == 1 &
                             (diag.status.syph == 0 | is.na(diag.status.syph)) &
                             prepStat == 0)
    rates.syph <- rep(1/stitest.active.int, length(elig.syph.ann))
    tst.syph.nprep.ann <- elig.syph.ann[rbinom(length(elig.syph.ann), 1, rates.syph) == 1]

    elig.syph.highrisk <- which(tt.traj.syph == 2 &
                                  (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                  prepStat == 0)
    rates.syph <- rep(1/sti.highrisktest.int, length(elig.syph.highrisk))
    tst.syph.nprep.highrisk <- elig.syph.highrisk[rbinom(length(elig.syph.highrisk), 1, rates.syph) == 1]
    tst.syph.nprep <- c(tst.syph.nprep.ann, tst.syph.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.syph.annual.interval <- which(tt.traj.syph == 1 &
                                        (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                        tsinceltst.syph >= (stitest.active.int) &
                                        prepStat == 0)
    tst.syph.highrisk.interval <- which(tt.traj.syph == 2 &
                                          (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                          tsinceltst.syph >= (sti.highrisktest.int) &
                                          prepStat == 0)
    tst.syph.nprep <- c(tst.syph.annual.interval, tst.syph.highrisk.interval)
  }

  ## Process for GC asymptomatic screening
  if (testing.pattern.sti == "memoryless") {
    elig.gc.ann <- which(tt.traj.gc == 1 &
                           (diag.status.gc == 0 | is.na(diag.status.gc)) &
                           prepStat == 0)
    rates.gc <- rep(1/stitest.active.int, length(elig.gc.ann))
    tst.gc.nprep.ann <- elig.gc.ann[rbinom(length(elig.gc.ann), 1, rates.gc) == 1]

    elig.gc.highrisk <- which(tt.traj.gc == 2 &
                                (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                prepStat == 0)
    rates.gc <- rep(1/sti.highrisktest.int, length(elig.gc.highrisk))
    tst.gc.nprep.highrisk <- elig.gc.highrisk[rbinom(length(elig.gc.highrisk), 1, rates.gc) == 1]
    tst.gc.nprep <- c(tst.gc.nprep.ann, tst.gc.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.gc.annual.interval <- which(tt.traj.gc == 1 &
                                      (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                      tsinceltst.gc >= (stitest.active.int) &
                                      prepStat == 0)
    tst.gc.highrisk.interval <- which(tt.traj.gc == 2 &
                                        (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                        tsinceltst.gc >= (sti.highrisktest.int) &
                                        prepStat == 0)
    tst.gc.nprep <- c(tst.gc.annual.interval, tst.gc.highrisk.interval)
  }

  ## Process for CT asymptomatic screening
  if (testing.pattern.sti == "memoryless") {
    elig.ct.ann <- which(tt.traj.ct == 1 &
                           (diag.status.ct == 0 | is.na(diag.status.ct)) &
                           prepStat == 0)
    rates.ct <- rep(1/stitest.active.int, length(elig.ct.ann))
    tst.ct.nprep.ann <- elig.ct.ann[rbinom(length(elig.ct.ann), 1, rates.ct) == 1]

    elig.ct.highrisk <- which(tt.traj.ct == 2 &
                                (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                prepStat == 0)
    rates.ct <- rep(1/sti.highrisktest.int, length(elig.ct.highrisk))
    tst.ct.nprep.highrisk <- elig.ct.highrisk[rbinom(length(elig.ct.highrisk), 1, rates.ct) == 1]

    tst.ct.nprep <- c(tst.ct.nprep.ann, tst.ct.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.ct.annual.interval <- which(tt.traj.ct == 1 &
                                      (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                      tsinceltst.ct >= (stitest.active.int) &
                                      prepStat == 0)
    tst.ct.highrisk.interval <- which(tt.traj.ct == 2 &
                                        (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                        tsinceltst.ct >= (sti.highrisktest.int) &
                                        prepStat == 0)
    tst.ct.nprep <- c(tst.ct.annual.interval, tst.ct.highrisk.interval)
  }

  # Syphilis non-PrEP testing
  tst.syph.pos <- tst.syph.nprep[syphilis[tst.syph.nprep] == 1 &
                                   stage.syph[tst.syph.nprep] %in% c(2, 3, 4, 5, 6)]
  tst.syph.neg <- setdiff(tst.syph.nprep, tst.syph.pos)
  tst.earlysyph.pos <- tst.syph.nprep[syphilis[tst.syph.nprep] == 1 &
                                   stage.syph[tst.syph.nprep] %in% c(2, 3, 4)]
  tst.latesyph.pos <- tst.syph.nprep[syphilis[tst.syph.nprep] == 1 &
                                   stage.syph[tst.syph.nprep] %in% c(5, 6)]

  # GC non-PrEP testing
  tst.rgc <- tst.gc.nprep[role.class[tst.gc.nprep] %in% c("R", "V")]
  #tst.rgc <- sample(tst.rgc, tst.rect.sti.rr * length(tst.rgc))
  tst.ugc <- tst.gc.nprep[role.class[tst.gc.nprep] %in% c("I", "V")]
  tst.rgc.pos <- tst.rgc[rGC[tst.rgc] == 1]
  tst.ugc.pos <- tst.ugc[uGC[tst.ugc] == 1]
  tst.rgc.neg <- setdiff(tst.rgc, tst.rgc.pos)
  tst.ugc.neg <- setdiff(tst.ugc, tst.ugc.pos)
  tst.gc.pos <- unique(c(tst.rgc.pos, tst.ugc.pos))

  # CT non-PrEP testing
  tst.rct <- tst.ct.nprep[role.class[tst.ct.nprep] %in% c("R", "V")]
  #tst.rct <- sample(tst.rct, tst.rect.sti.rr * length(tst.rct))
  tst.uct <- tst.ct.nprep[role.class[tst.ct.nprep] %in% c("I", "V")]
  tst.rct.pos <- tst.rct[rCT[tst.rct] == 1]
  tst.uct.pos <- tst.uct[uCT[tst.uct] == 1]
  tst.rct.neg <- setdiff(tst.rct, tst.rct.pos)
  tst.uct.neg <- setdiff(tst.uct, tst.uct.pos)
  tst.ct.pos <- unique(c(tst.rct.pos, tst.uct.pos))

  # Syphilis Attributes
  last.neg.test.syph[tst.syph.neg] <- at
  last.neg.test.syph[tst.syph.pos] <- NA
  diag.status.syph[tst.syph.pos] <- 1
  last.diag.time.syph[tst.syph.pos] <- at
  tsinceltst.syph[tst.syph.nprep] <- 0

  # GC Attributes
  last.neg.test.rgc[tst.rgc.neg] <- at
  last.neg.test.ugc[tst.ugc.neg] <- at
  last.neg.test.rgc[tst.rgc.pos] <- NA
  last.neg.test.ugc[tst.ugc.pos] <- NA
  diag.status.gc[tst.gc.pos] <- 1
  last.diag.time.gc[tst.gc.pos] <- at
  tsinceltst.rgc[tst.rgc] <- 0
  tsinceltst.ugc[tst.ugc] <- 0

  # CT Attributes
  last.neg.test.rct[tst.rct.neg] <- at
  last.neg.test.uct[tst.uct.neg] <- at
  last.neg.test.rct[tst.rct.pos] <- NA
  last.neg.test.uct[tst.uct.pos] <- NA
  diag.status.ct[tst.ct.pos] <- 1
  last.diag.time.ct[tst.ct.pos] <- at
  tsinceltst.rct[tst.rct] <- 0
  tsinceltst.uct[tst.uct] <- 0

  # Number of tests for asymptomatic
  dat$epi$rGCasympttests[at] <- length(tst.rgc)
  dat$epi$uGCasympttests[at] <- length(tst.ugc)
  dat$epi$GCasympttests[at] <- length(tst.rgc) + length(tst.ugc)

  dat$epi$rGCasympttests.pos[at] <- length(tst.rgc.pos)
  dat$epi$uGCasympttests.pos[at] <- length(tst.ugc.pos)
  dat$epi$GCasympttests.pos[at] <- length(tst.rgc.pos) + length(tst.ugc.pos)

  dat$epi$rCTasympttests[at] <- length(tst.rct)
  dat$epi$uCTasympttests[at] <- length(tst.uct)
  dat$epi$CTasympttests[at] <- length(tst.rct) + length(tst.uct)

  dat$epi$rCTasympttests.pos[at] <- length(tst.rct.pos)
  dat$epi$uCTasympttests.pos[at] <- length(tst.uct.pos)
  dat$epi$CTasympttests.pos[at] <- length(tst.rct.pos) + length(tst.uct.pos)

  dat$epi$syphasympttests[at] <- length(tst.syph.nprep)
  dat$epi$syphasympttests.pos[at] <- length(tst.syph.pos)
  dat$epi$syphearlyasympttests.pos[at] <- length(tst.earlysyph.pos)
  dat$epi$syphlateasympttests.pos[at] <- length(tst.latesyph.pos)

  dat$epi$stiasympttests.pos[at] <- length(tst.rgc) + length(tst.ugc) +
    length(tst.rct) + length(tst.uct) + length(tst.syph.nprep)
  dat$epi$stiasympttests.pos[at] <- length(tst.rgc.pos) + length(tst.ugc.pos) +
    length(tst.rct.pos) + length(tst.uct.pos) + length(tst.syph.pos)

  ## Output -----------------------------------------------------------------

  # Attributes

  # Stoppage attributes
  dat$attr$stihighrisktestLastElig[idsnothighriskelig] <- at
  dat$attr$stianntestLastElig[idsnotactiveelig] <- at

  # Syphilis Attributes
  dat$attr$last.neg.test.syph <- last.neg.test.syph
  dat$attr$diag.status.syph <- diag.status.syph
  dat$attr$last.diag.time.syph <- last.diag.time.syph
  dat$attr$tt.traj.syph <- tt.traj.syph
  dat$attr$time.since.last.test.syph <- tsinceltst.syph

  # GC Attributes
  dat$attr$last.neg.test.rgc <- last.neg.test.rgc
  dat$attr$last.neg.test.ugc <- last.neg.test.ugc
  dat$attr$diag.status.gc <- diag.status.gc
  dat$attr$last.diag.time.gc <- last.diag.time.gc
  dat$attr$tt.traj.gc <- tt.traj.gc
  dat$attr$time.since.last.test.rgc <- tsinceltst.rgc
  dat$attr$time.since.last.test.ugc <- tsinceltst.ugc

  # CT Attributes
  dat$attr$last.neg.test.rct <- last.neg.test.rct
  dat$attr$last.neg.test.uct <- last.neg.test.uct
  dat$attr$diag.status.ct <- diag.status.ct
  dat$attr$last.diag.time.ct <- last.diag.time.ct
  dat$attr$tt.traj.ct <- tt.traj.ct
  dat$attr$time.since.last.test.rct <- tsinceltst.rct
  dat$attr$time.since.last.test.uct <- tsinceltst.uct

  return(dat)
}



#' @title HIV Diagnosis Module
#'
#' @description Module function for simulating HIV diagnosis after infection,
#'              currently based on diagnosis at treatment initiation.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
dx_het <- function(dat, at) {

  # Variables
  status <- dat$attr$status
  txCD4min <- dat$attr$txCD4min
  cd4Count <- dat$attr$cd4Count
  dxStat <- dat$attr$dxStat

  # Process
  tested <- which(status == 1 & dxStat == 0 & cd4Count <= txCD4min)


  # Results
  if (length(tested) > 0) {
    dat$attr$dxStat[tested] <- 1
    dat$attr$txStat[tested] <- 0
    dat$attr$dxTime[tested] <- at
  }

  return(dat)
}

