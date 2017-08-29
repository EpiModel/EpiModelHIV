
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
  dat$attr$tt.traj.syph.hivneg[tst.pos] <- NA
  dat$attr$tt.traj.gc.hivneg[tst.pos] <- NA
  dat$attr$tt.traj.ct.hivneg[tst.pos] <- NA
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
  tt.traj.syph.hivpos <- dat$attr$tt.traj.syph.hivpos
  tt.traj.syph.hivneg <- dat$attr$tt.traj.syph.hivneg
  tt.traj.gc.hivpos <- dat$attr$tt.traj.gc.hivpos
  tt.traj.gc.hivneg <- dat$attr$tt.traj.gc.hivneg
  tt.traj.ct.hivpos <- dat$attr$tt.traj.ct.hivpos
  tt.traj.ct.hivneg <- dat$attr$tt.traj.ct.hivneg
  diag.status.gc <- dat$attr$diag.status.gc
  diag.status.ct <- dat$attr$diag.status.ct
  diag.status.syph <- dat$attr$diag.status.syph
  diag.status <- dat$attr$diag.status
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
  stianntest.hivneg.coverage <- dat$param$stianntest.hivneg.coverage
  stihighrisktest.hivneg.coverage <- dat$param$stihighrisktest.hivneg.coverage
  stianntest.hivpos.coverage <- dat$param$stianntest.hivpos.coverage
  stihighrisktest.hivpos.coverage <- dat$param$stihighrisktest.hivpos.coverage
  stianntest.cov.rate <- dat$param$stianntest.cov.rate
  stihighrisktest.cov.rate <- dat$param$stihighrisktest.cov.rate
  testing.pattern.sti <- dat$param$testing.pattern.sti
  stitest.active.int <- dat$param$stitest.active.int
  sti.highrisktest.int <- dat$param$sti.highrisktest.int
  tst.rect.sti.rr <- dat$param$tst.rect.sti.rr

  # Eligibility and trajectory
  # Base eligibility
  idsEligTest.hivpos <- which(race %in% c("B", "W") & diag.status == 1)
  idsEligTest.hivneg <- which(race %in% c("B", "W") & (diag.status == 0 | is.na(diag.status)))

  # Annual indications- sexually active in last year
  idsactive.hivpos <- which(stitestind1 == 1 & diag.status == 1)
  idsactive.hivneg <- which(stitestind1 == 1 & (diag.status == 0 | is.na(diag.status)))

  # STI testing higher-risk eligibility scenarios
  idshighrisk.hivpos <- which(stitestind2 == 1 & diag.status == 1)
  idshighrisk.hivneg <- which(stitestind2 == 1 & (diag.status == 0 | is.na(diag.status)))

  ## Stoppage (tt.traj.gc/.ct/.syph <- NA ------------------------------------
  # Reduce testing trajectory to NA if no longer indicated for more frequent high-risk testing
  idsnothighriskelig.hivpos <- which(tt.traj.syph.hivpos == 2 & stitestind2 != 1)
  tt.traj.syph.hivpos[idsnothighriskelig.hivpos] <-
    tt.traj.gc.hivpos[idsnothighriskelig.hivpos] <-
    tt.traj.ct.hivpos[idsnothighriskelig.hivpos] <-
    NA
  idsnothighriskelig.hivneg <- which(tt.traj.syph.hivneg == 2 & stitestind2 != 1)
  tt.traj.syph.hivneg[idsnothighriskelig.hivneg] <-
    tt.traj.gc.hivneg[idsnothighriskelig.hivneg] <-
    tt.traj.ct.hivneg[idsnothighriskelig.hivneg] <-
    NA

  # Reduce testing trajectory to NA if no longer indicated for lower-risk testing
  idsnotactiveelig.hivpos <- which(tt.traj.syph.hivpos == 2 & stitestind1 != 1)
  tt.traj.syph.hivpos[idsnotactiveelig.hivpos] <-
    tt.traj.gc.hivpos[idsnotactiveelig.hivpos] <-
    tt.traj.ct.hivpos[idsnotactiveelig.hivpos] <-
    NA
  idsnotactiveelig.hivneg <- which(tt.traj.syph.hivneg == 2 & stitestind1 != 1)
  tt.traj.syph.hivneg[idsnotactiveelig.hivneg] <-
    tt.traj.gc.hivneg[idsnotactiveelig.hivneg] <-
    tt.traj.ct.hivneg[idsnotactiveelig.hivneg] <-
    NA

  ## Initiation (non-HIV diagnosed) --------------------------------------------
  ### Testing coverage for high risk
  stihighrisktestCov <- sum(tt.traj.ct.hivneg == 2, na.rm = TRUE) / length(idshighrisk.hivneg)
  stihighrisktestCov <- ifelse(is.nan(stihighrisktestCov), 0, stihighrisktestCov)

  idsEligSt <- idshighrisk.hivneg
  nEligSt <- length(idshighrisk.hivneg)

  nStart <- max(0, min(nEligSt, round((stihighrisktest.hivneg.coverage - stihighrisktestCov) *
                                        length(idshighrisk.hivneg))))
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
    tt.traj.syph.hivneg[idsStart] <- tt.traj.gc.hivneg[idsStart] <- tt.traj.ct.hivneg[idsStart] <- 2
  }

  ### Testing coverage for annual - all those sexually active without high-risk indications
  stianntestCov <- sum(tt.traj.ct.hivneg == 1, na.rm = TRUE) / length(setdiff(idsactive.hivneg, which(tt.traj.syph.hivneg == 2)))
  stianntestCov <- ifelse(is.nan(stianntestCov), 0, stianntestCov)

  idsEligSt <- setdiff(idsactive.hivneg, idshighrisk.hivneg)
  nEligSt <- length(setdiff(idsactive.hivneg, idshighrisk.hivneg))

  nStart <- max(0, min(nEligSt, round((stianntest.hivneg.coverage - stianntestCov) *
                                        length(setdiff(idsactive.hivneg, idshighrisk.hivneg)))))
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
    tt.traj.syph.hivneg[idsStart] <- tt.traj.gc.hivneg[idsStart] <- tt.traj.ct.hivneg[idsStart] <- 1
  }

  ## Process for asymptomatic syphilis screening
  if (testing.pattern.sti == "memoryless") {
    elig.syph.ann <- which(tt.traj.syph.hivneg == 1 &
                             (diag.status.syph == 0 | is.na(diag.status.syph)) &
                             prepStat == 0)
    rates.syph <- rep(1/stitest.active.int, length(elig.syph.ann))
    tst.syph.nprep.ann <- elig.syph.ann[rbinom(length(elig.syph.ann), 1, rates.syph) == 1]

    elig.syph.highrisk <- which(tt.traj.syph.hivneg == 2 &
                                  (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                  prepStat == 0)
    rates.syph <- rep(1/sti.highrisktest.int, length(elig.syph.highrisk))
    tst.syph.nprep.highrisk <- elig.syph.highrisk[rbinom(length(elig.syph.highrisk), 1, rates.syph) == 1]
    tst.syph.nprep.hivneg <- c(tst.syph.nprep.ann, tst.syph.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.syph.annual.interval <- which(tt.traj.syph.hivneg == 1 &
                                        (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                        tsinceltst.syph >= (stitest.active.int) &
                                        prepStat == 0)
    tst.syph.highrisk.interval <- which(tt.traj.syph.hivneg == 2 &
                                          (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                          tsinceltst.syph >= (sti.highrisktest.int) &
                                          prepStat == 0)
    tst.syph.nprep.hivneg <- c(tst.syph.annual.interval, tst.syph.highrisk.interval)
  }

  ## Process for GC asymptomatic screening
  if (testing.pattern.sti == "memoryless") {
    elig.gc.ann <- which(tt.traj.gc.hivneg == 1 &
                           (diag.status.gc == 0 | is.na(diag.status.gc)) &
                           prepStat == 0)
    rates.gc <- rep(1/stitest.active.int, length(elig.gc.ann))
    tst.gc.nprep.ann <- elig.gc.ann[rbinom(length(elig.gc.ann), 1, rates.gc) == 1]

    elig.gc.highrisk <- which(tt.traj.gc.hivneg == 2 &
                                (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                prepStat == 0)
    rates.gc <- rep(1/sti.highrisktest.int, length(elig.gc.highrisk))
    tst.gc.nprep.highrisk <- elig.gc.highrisk[rbinom(length(elig.gc.highrisk), 1, rates.gc) == 1]
    tst.gc.nprep.hivneg <- c(tst.gc.nprep.ann, tst.gc.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.gc.annual.interval <- which(tt.traj.gc.hivneg == 1 &
                                      (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                      tsinceltst.gc >= (stitest.active.int) &
                                      prepStat == 0)
    tst.gc.highrisk.interval <- which(tt.traj.gc.hivneg == 2 &
                                        (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                        tsinceltst.gc >= (sti.highrisktest.int) &
                                        prepStat == 0)
    tst.gc.nprep.hivneg <- c(tst.gc.annual.interval, tst.gc.highrisk.interval)
  }

  ## Process for CT asymptomatic screening
  if (testing.pattern.sti == "memoryless") {
    elig.ct.ann <- which(tt.traj.ct.hivneg == 1 &
                           (diag.status.ct == 0 | is.na(diag.status.ct)) &
                           prepStat == 0)
    rates.ct <- rep(1/stitest.active.int, length(elig.ct.ann))
    tst.ct.nprep.ann <- elig.ct.ann[rbinom(length(elig.ct.ann), 1, rates.ct) == 1]

    elig.ct.highrisk <- which(tt.traj.ct.hivneg == 2 &
                                (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                prepStat == 0)
    rates.ct <- rep(1/sti.highrisktest.int, length(elig.ct.highrisk))
    tst.ct.nprep.highrisk <- elig.ct.highrisk[rbinom(length(elig.ct.highrisk), 1, rates.ct) == 1]

    tst.ct.nprep.hivneg <- c(tst.ct.nprep.ann, tst.ct.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.ct.annual.interval <- which(tt.traj.ct.hivneg == 1 &
                                      (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                      tsinceltst.ct >= (stitest.active.int) &
                                      prepStat == 0)
    tst.ct.highrisk.interval <- which(tt.traj.ct.hivneg == 2 &
                                        (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                        tsinceltst.ct >= (sti.highrisktest.int) &
                                        prepStat == 0)
    tst.ct.nprep.hivneg <- c(tst.ct.annual.interval, tst.ct.highrisk.interval)
  }

  # Syphilis non-PrEP testing
  tst.syph.pos.hivneg <- tst.syph.nprep.hivneg[syphilis[tst.syph.nprep.hivneg] == 1 &
                                   stage.syph[tst.syph.nprep.hivneg] %in% c(2, 3, 4, 5, 6)]
  tst.syph.neg.hivneg <- setdiff(tst.syph.nprep.hivneg, tst.syph.pos.hivneg)
  tst.earlysyph.pos.hivneg <- tst.syph.nprep.hivneg[syphilis[tst.syph.nprep.hivneg] == 1 &
                                   stage.syph[tst.syph.nprep.hivneg] %in% c(2, 3, 4)]
  tst.latesyph.pos.hivneg <- tst.syph.nprep.hivneg[syphilis[tst.syph.nprep.hivneg] == 1 &
                                   stage.syph[tst.syph.nprep.hivneg] %in% c(5, 6)]

  # GC non-PrEP testing
  tst.rgc.hivneg <- tst.gc.nprep.hivneg[role.class[tst.gc.nprep.hivneg] %in% c("R", "V")]
  #tst.rgc <- sample(tst.rgc.hivneg, tst.rect.sti.rr * length(tst.rgc.hivneg))
  tst.ugc.hivneg <- tst.gc.nprep.hivneg[role.class[tst.gc.nprep.hivneg] %in% c("I", "V")]
  tst.rgc.pos.hivneg <- tst.rgc.hivneg[rGC[tst.rgc.hivneg] == 1]
  tst.ugc.pos.hivneg <- tst.ugc.hivneg[uGC[tst.ugc.hivneg] == 1]
  tst.rgc.neg.hivneg <- setdiff(tst.rgc.hivneg, tst.rgc.pos.hivneg)
  tst.ugc.neg.hivneg <- setdiff(tst.ugc.hivneg, tst.ugc.pos.hivneg)
  tst.gc.pos.hivneg <- unique(c(tst.rgc.pos.hivneg, tst.ugc.pos.hivneg))

  # CT non-PrEP testing
  tst.rct.hivneg <- tst.ct.nprep.hivneg[role.class[tst.ct.nprep.hivneg] %in% c("R", "V")]
  #tst.rct <- sample(tst.rct.hivneg, tst.rect.sti.rr * length(tst.rct.hivneg))
  tst.uct.hivneg <- tst.ct.nprep.hivneg[role.class[tst.ct.nprep.hivneg] %in% c("I", "V")]
  tst.rct.pos.hivneg <- tst.rct.hivneg[rCT[tst.rct.hivneg] == 1]
  tst.uct.pos.hivneg <- tst.uct.hivneg[uCT[tst.uct.hivneg] == 1]
  tst.rct.neg.hivneg <- setdiff(tst.rct.hivneg, tst.rct.pos.hivneg)
  tst.uct.neg.hivneg <- setdiff(tst.uct.hivneg, tst.uct.pos.hivneg)
  tst.ct.pos.hivneg <- unique(c(tst.rct.pos.hivneg, tst.uct.pos.hivneg))

  # Syphilis Attributes
  last.neg.test.syph[tst.syph.neg.hivneg] <- at
  last.neg.test.syph[tst.syph.pos.hivneg] <- NA
  diag.status.syph[tst.syph.pos.hivneg] <- 1
  last.diag.time.syph[tst.syph.pos.hivneg] <- at
  tsinceltst.syph[tst.syph.nprep.hivneg] <- 0

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

  ## Initiation (HIV diagnosed) --------------------------------------------
  ### Testing coverage for high risk
  stihighrisktestCov <- sum(tt.traj.ct.hivpos == 2, na.rm = TRUE) / length(idshighrisk.hivpos)
  stihighrisktestCov <- ifelse(is.nan(stihighrisktestCov), 0, stihighrisktestCov)

  idsEligSt <- idshighrisk.hivpos
  nEligSt <- length(idshighrisk.hivpos)

  nStart <- max(0, min(nEligSt, round((stihighrisktest.hivpos.coverage - stihighrisktestCov) *
                                        length(idshighrisk.hivpos))))
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
    tt.traj.syph.hivpos[idsStart] <- tt.traj.gc.hivpos[idsStart] <- tt.traj.ct.hivpos[idsStart] <- 2
  }

  ### Testing coverage for annual - all those sexually active without high-risk indications
  stianntestCov <- sum(tt.traj.ct.hivpos == 1, na.rm = TRUE) / length(setdiff(idsactive.hivpos, which(tt.traj.syph.hivpos == 2)))
  stianntestCov <- ifelse(is.nan(stianntestCov), 0, stianntestCov)

  idsEligSt <- setdiff(idsactive.hivpos, idshighrisk.hivpos)
  nEligSt <- length(setdiff(idsactive.hivpos, idshighrisk.hivpos))

  nStart <- max(0, min(nEligSt, round((stianntest.hivpos.coverage - stianntestCov) *
                                        length(setdiff(idsactive.hivpos, idshighrisk.hivpos)))))
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
    tt.traj.syph.hivpos[idsStart] <- tt.traj.gc.hivpos[idsStart] <- tt.traj.ct.hivpos[idsStart] <- 1
  }

  ## Process for asymptomatic syphilis screening
  if (testing.pattern.sti == "memoryless") {
    elig.syph.ann <- which(tt.traj.syph.hivpos == 1 &
                             (diag.status.syph == 0 | is.na(diag.status.syph)) &
                             prepStat == 0)
    rates.syph <- rep(1/stitest.active.int, length(elig.syph.ann))
    tst.syph.nprep.ann <- elig.syph.ann[rbinom(length(elig.syph.ann), 1, rates.syph) == 1]

    elig.syph.highrisk <- which(tt.traj.syph.hivpos == 2 &
                                  (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                  prepStat == 0)
    rates.syph <- rep(1/sti.highrisktest.int, length(elig.syph.highrisk))
    tst.syph.nprep.highrisk <- elig.syph.highrisk[rbinom(length(elig.syph.highrisk), 1, rates.syph) == 1]
    tst.syph.nprep.hivpos <- c(tst.syph.nprep.ann, tst.syph.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.syph.annual.interval <- which(tt.traj.syph.hivpos == 1 &
                                        (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                        tsinceltst.syph >= (stitest.active.int) &
                                        prepStat == 0)
    tst.syph.highrisk.interval <- which(tt.traj.syph.hivpos == 2 &
                                          (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                          tsinceltst.syph >= (sti.highrisktest.int) &
                                          prepStat == 0)
    tst.syph.nprep.hivpos <- c(tst.syph.annual.interval, tst.syph.highrisk.interval)
  }

  ## Process for GC asymptomatic screening
  if (testing.pattern.sti == "memoryless") {
    elig.gc.ann <- which(tt.traj.gc.hivpos == 1 &
                           (diag.status.gc == 0 | is.na(diag.status.gc)) &
                           prepStat == 0)
    rates.gc <- rep(1/stitest.active.int, length(elig.gc.ann))
    tst.gc.nprep.ann <- elig.gc.ann[rbinom(length(elig.gc.ann), 1, rates.gc) == 1]

    elig.gc.highrisk <- which(tt.traj.gc.hivpos == 2 &
                                (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                prepStat == 0)
    rates.gc <- rep(1/sti.highrisktest.int, length(elig.gc.highrisk))
    tst.gc.nprep.highrisk <- elig.gc.highrisk[rbinom(length(elig.gc.highrisk), 1, rates.gc) == 1]
    tst.gc.nprep.hivpos <- c(tst.gc.nprep.ann, tst.gc.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.gc.annual.interval <- which(tt.traj.gc.hivpos == 1 &
                                      (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                      tsinceltst.gc >= (stitest.active.int) &
                                      prepStat == 0)
    tst.gc.highrisk.interval <- which(tt.traj.gc.hivpos == 2 &
                                        (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                        tsinceltst.gc >= (sti.highrisktest.int) &
                                        prepStat == 0)
    tst.gc.nprep.hivpos <- c(tst.gc.annual.interval, tst.gc.highrisk.interval)
  }

  ## Process for CT asymptomatic screening
  if (testing.pattern.sti == "memoryless") {
    elig.ct.ann <- which(tt.traj.ct.hivpos == 1 &
                           (diag.status.ct == 0 | is.na(diag.status.ct)) &
                           prepStat == 0)
    rates.ct <- rep(1/stitest.active.int, length(elig.ct.ann))
    tst.ct.nprep.ann <- elig.ct.ann[rbinom(length(elig.ct.ann), 1, rates.ct) == 1]

    elig.ct.highrisk <- which(tt.traj.ct.hivpos == 2 &
                                (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                prepStat == 0)
    rates.ct <- rep(1/sti.highrisktest.int, length(elig.ct.highrisk))
    tst.ct.nprep.highrisk <- elig.ct.highrisk[rbinom(length(elig.ct.highrisk), 1, rates.ct) == 1]

    tst.ct.nprep.hivpos <- c(tst.ct.nprep.ann, tst.ct.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.ct.annual.interval <- which(tt.traj.ct.hivpos == 1 &
                                      (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                      tsinceltst.ct >= (stitest.active.int) &
                                      prepStat == 0)
    tst.ct.highrisk.interval <- which(tt.traj.ct.hivpos == 2 &
                                        (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                        tsinceltst.ct >= (sti.highrisktest.int) &
                                        prepStat == 0)
    tst.ct.nprep.hivpos <- c(tst.ct.annual.interval, tst.ct.highrisk.interval)
  }

  # Syphilis non-PrEP testing
  tst.syph.pos.hivpos <- tst.syph.nprep.hivpos[syphilis[tst.syph.nprep.hivpos] == 1 &
                                                 stage.syph[tst.syph.nprep.hivpos] %in% c(2, 3, 4, 5, 6)]
  tst.syph.neg.hivpos <- setdiff(tst.syph.nprep.hivpos, tst.syph.pos.hivpos)
  tst.earlysyph.pos.hivpos <- tst.syph.nprep.hivpos[syphilis[tst.syph.nprep.hivpos] == 1 &
                                                      stage.syph[tst.syph.nprep.hivpos] %in% c(2, 3, 4)]
  tst.latesyph.pos.hivpos <- tst.syph.nprep.hivpos[syphilis[tst.syph.nprep.hivpos] == 1 &
                                                     stage.syph[tst.syph.nprep.hivpos] %in% c(5, 6)]

  # GC non-PrEP testing
  tst.rgc.hivpos <- tst.gc.nprep.hivpos[role.class[tst.gc.nprep.hivpos] %in% c("R", "V")]
  #tst.rgc <- sample(tst.rgc.hivpos, tst.rect.sti.rr * length(tst.rgc.hivpos))
  tst.ugc.hivpos <- tst.gc.nprep.hivpos[role.class[tst.gc.nprep.hivpos] %in% c("I", "V")]
  tst.rgc.pos.hivpos <- tst.rgc.hivpos[rGC[tst.rgc.hivpos] == 1]
  tst.ugc.pos.hivpos <- tst.ugc.hivpos[uGC[tst.ugc.hivpos] == 1]
  tst.rgc.neg.hivpos <- setdiff(tst.rgc.hivpos, tst.rgc.pos.hivpos)
  tst.ugc.neg.hivpos <- setdiff(tst.ugc.hivpos, tst.ugc.pos.hivpos)
  tst.gc.pos.hivpos <- unique(c(tst.rgc.pos.hivpos, tst.ugc.pos.hivpos))

  # CT non-PrEP testing
  tst.rct.hivpos <- tst.ct.nprep.hivpos[role.class[tst.ct.nprep.hivpos] %in% c("R", "V")]
  #tst.rct <- sample(tst.rct.hivpos, tst.rect.sti.rr * length(tst.rct.hivpos))
  tst.uct.hivpos <- tst.ct.nprep.hivpos[role.class[tst.ct.nprep.hivpos] %in% c("I", "V")]
  tst.rct.pos.hivpos <- tst.rct.hivpos[rCT[tst.rct.hivpos] == 1]
  tst.uct.pos.hivpos <- tst.uct.hivpos[uCT[tst.uct.hivpos] == 1]
  tst.rct.neg.hivpos <- setdiff(tst.rct.hivpos, tst.rct.pos.hivpos)
  tst.uct.neg.hivpos <- setdiff(tst.uct.hivpos, tst.uct.pos.hivpos)
  tst.ct.pos.hivpos <- unique(c(tst.rct.pos.hivpos, tst.uct.pos.hivpos))

  # Syphilis Attributes
  last.neg.test.syph[tst.syph.neg.hivpos] <- at
  last.neg.test.syph[tst.syph.pos.hivpos] <- NA
  diag.status.syph[tst.syph.pos.hivpos] <- 1
  last.diag.time.syph[tst.syph.pos.hivpos] <- at
  tsinceltst.syph[tst.syph.nprep.hivpos] <- 0

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

  # Number of tests for asymptomatic (non-HIV diagnosed)
  dat$epi$rGCasympttests.hivneg[at] <- length(tst.rgc.hivneg)
  dat$epi$uGCasympttests.hivneg[at] <- length(tst.ugc.hivneg)
  dat$epi$GCasympttests.hivneg[at] <- length(tst.rgc.hivneg) + length(tst.ugc.hivneg)

  dat$epi$rGCasympttests.pos.hivneg[at] <- length(tst.rgc.pos.hivneg)
  dat$epi$uGCasympttests.pos.hivneg[at] <- length(tst.ugc.pos.hivneg)
  dat$epi$GCasympttests.pos.hivneg[at] <- length(tst.rgc.pos.hivneg) + length(tst.ugc.pos.hivneg)

  dat$epi$rCTasympttests.hivneg[at] <- length(tst.rct.hivneg)
  dat$epi$uCTasympttests.hivneg[at] <- length(tst.uct.hivneg)
  dat$epi$CTasympttests.hivneg[at] <- length(tst.rct.hivneg) + length(tst.uct.hivneg)

  dat$epi$rCTasympttests.pos.hivneg[at] <- length(tst.rct.pos.hivneg)
  dat$epi$uCTasympttests.pos.hivneg[at] <- length(tst.uct.pos.hivneg)
  dat$epi$CTasympttests.pos.hivneg[at] <- length(tst.rct.pos.hivneg) + length(tst.uct.pos.hivneg)

  dat$epi$syphasympttests.hivneg[at] <- length(tst.syph.nprep.hivneg)
  dat$epi$syphasympttests.pos.hivneg[at] <- length(tst.syph.pos.hivneg)
  dat$epi$syphearlyasympttests.pos.hivneg[at] <- length(tst.earlysyph.pos.hivneg)
  dat$epi$syphlateasympttests.pos.hivneg[at] <- length(tst.latesyph.pos.hivneg)

  dat$epi$stiasympttests.pos.hivneg[at] <- length(tst.rgc.hivneg) + length(tst.ugc.hivneg) +
    length(tst.rct.hivneg) + length(tst.uct.hivneg) + length(tst.syph.nprep.hivneg)
  dat$epi$stiasympttests.pos.hivneg[at] <- length(tst.rgc.pos.hivneg) + length(tst.ugc.pos.hivneg) +
    length(tst.rct.pos.hivneg) + length(tst.uct.pos.hivneg) + length(tst.syph.pos.hivneg)

  # Number of tests for asymptomatic (HIV-diagnosed)
  dat$epi$rGCasympttests.hivpos[at] <- length(tst.rgc.hivpos)
  dat$epi$uGCasympttests.hivpos[at] <- length(tst.ugc.hivpos)
  dat$epi$GCasympttests.hivpos[at] <- length(tst.rgc.hivpos) + length(tst.ugc.hivpos)

  dat$epi$rGCasympttests.pos.hivpos[at] <- length(tst.rgc.pos.hivpos)
  dat$epi$uGCasympttests.pos.hivpos[at] <- length(tst.ugc.pos.hivpos)
  dat$epi$GCasympttests.pos.hivpos[at] <- length(tst.rgc.pos.hivpos) + length(tst.ugc.pos.hivpos)

  dat$epi$rCTasympttests.hivpos[at] <- length(tst.rct.hivpos)
  dat$epi$uCTasympttests.hivpos[at] <- length(tst.uct.hivpos)
  dat$epi$CTasympttests.hivpos[at] <- length(tst.rct.hivpos) + length(tst.uct.hivpos)

  dat$epi$rCTasympttests.pos.hivpos[at] <- length(tst.rct.pos.hivpos)
  dat$epi$uCTasympttests.pos.hivpos[at] <- length(tst.uct.pos.hivpos)
  dat$epi$CTasympttests.pos.hivpos[at] <- length(tst.rct.pos.hivpos) + length(tst.uct.pos.hivpos)

  dat$epi$syphasympttests.hivpos[at] <- length(tst.syph.nprep.hivpos)
  dat$epi$syphasympttests.pos.hivpos[at] <- length(tst.syph.pos.hivpos)
  dat$epi$syphearlyasympttests.pos.hivpos[at] <- length(tst.earlysyph.pos.hivpos)
  dat$epi$syphlateasympttests.pos.hivpos[at] <- length(tst.latesyph.pos.hivpos)

  dat$epi$stiasympttests.pos.hivpos[at] <- length(tst.rgc.hivpos) + length(tst.ugc.hivpos) +
    length(tst.rct.hivpos) + length(tst.uct.hivpos) + length(tst.syph.nprep.hivpos)
  dat$epi$stiasympttests.pos.hivpos[at] <- length(tst.rgc.pos.hivpos) + length(tst.ugc.pos.hivpos) +
    length(tst.rct.pos.hivpos) + length(tst.uct.pos.hivpos) + length(tst.syph.pos.hivpos)

  # Number of tests for asymptomatic (total)
  dat$epi$rGCasympttests[at] <- length(tst.rgc.hivpos) + length(tst.rgc.hivneg)
  dat$epi$uGCasympttests[at] <- length(tst.ugc.hivpos) + length(tst.ugc.hivneg)
  dat$epi$GCasympttests[at] <- length(tst.rgc.hivpos) + length(tst.ugc.hivpos) +
    length(tst.rgc.hivneg) + length(tst.ugc.hivneg)

  dat$epi$rGCasympttests.pos[at] <- length(tst.rgc.pos.hivpos) + length(tst.rgc.pos.hivneg)
  dat$epi$uGCasympttests.pos[at] <- length(tst.ugc.pos.hivpos) + length(tst.ugc.pos.hivneg)
  dat$epi$GCasympttests.pos[at] <- length(tst.rgc.pos.hivpos) + length(tst.ugc.pos.hivpos) +
    length(tst.rgc.pos.hivneg) + length(tst.ugc.pos.hivneg)

  dat$epi$rCTasympttests[at] <- length(tst.rct.hivpos) + length(tst.rct.hivneg)
  dat$epi$uCTasympttests[at] <- length(tst.uct.hivpos) + length(tst.uct.hivneg)
  dat$epi$CTasympttests[at] <- length(tst.rct.hivpos) + length(tst.uct.hivpos) +
    length(tst.rct.hivneg) + length(tst.uct.hivneg)

  dat$epi$rCTasympttests.pos[at] <- length(tst.rct.pos.hivpos) + length(tst.rct.pos.hivneg)
  dat$epi$uCTasympttests.pos[at] <- length(tst.uct.pos.hivpos) + length(tst.uct.pos.hivneg)
  dat$epi$CTasympttests.pos[at] <- length(tst.rct.pos.hivpos) + length(tst.uct.pos.hivpos) +
    length(tst.rct.pos.hivneg) + length(tst.uct.pos.hivneg)

  dat$epi$syphasympttests[at] <- length(tst.syph.nprep.hivpos) + length(tst.syph.nprep.hivneg)
  dat$epi$syphasympttests.pos[at] <- length(tst.syph.pos.hivpos) + length(tst.syph.pos.hivneg)
  dat$epi$syphearlyasympttests.pos[at] <- length(tst.earlysyph.pos.hivpos) + length(tst.earlysyph.pos.hivneg)
  dat$epi$syphlateasympttests.pos[at] <- length(tst.latesyph.pos.hivpos) + length(tst.latesyph.pos.hivneg)

  dat$epi$stiasympttests.pos[at] <- length(tst.rgc.hivpos) + length(tst.ugc.hivpos) +
    length(tst.rct.hivpos) + length(tst.uct.hivpos) + length(tst.syph.nprep.hivpos) +
    length(tst.rgc.hivneg) + length(tst.ugc.hivneg) +
    length(tst.rct.hivneg) + length(tst.uct.hivneg) + length(tst.syph.nprep.hivneg)
  dat$epi$stiasympttests.pos[at] <- length(tst.rgc.pos.hivpos) + length(tst.ugc.pos.hivpos) +
    length(tst.rct.pos.hivpos) + length(tst.uct.pos.hivpos) + length(tst.syph.pos.hivpos) +
    length(tst.rgc.pos.hivneg) + length(tst.ugc.pos.hivneg) +
    length(tst.rct.pos.hivneg) + length(tst.uct.pos.hivneg) + length(tst.syph.pos.hivneg)

  # Attributes
  # Stoppage attributes
  dat$attr$stihighrisktestLastElig[idsnothighriskelig.hivneg] <- at
  dat$attr$stianntestLastElig[idsnotactiveelig.hivneg] <- at
  dat$attr$stihighrisktestLastElig[idsnothighriskelig.hivpos] <- at
  dat$attr$stianntestLastElig[idsnotactiveelig.hivpos] <- at

  # Syphilis Attributes
  dat$attr$last.neg.test.syph <- last.neg.test.syph
  dat$attr$diag.status.syph <- diag.status.syph
  dat$attr$last.diag.time.syph <- last.diag.time.syph
  dat$attr$tt.traj.syph.hivneg <- tt.traj.syph.hivneg
  dat$attr$time.since.last.test.syph <- tsinceltst.syph

  # GC Attributes
  dat$attr$last.neg.test.rgc <- last.neg.test.rgc
  dat$attr$last.neg.test.ugc <- last.neg.test.ugc
  dat$attr$diag.status.gc <- diag.status.gc
  dat$attr$last.diag.time.gc <- last.diag.time.gc
  dat$attr$tt.traj.gc.hivneg <- tt.traj.gc.hivneg
  dat$attr$time.since.last.test.rgc <- tsinceltst.rgc
  dat$attr$time.since.last.test.ugc <- tsinceltst.ugc

  # CT Attributes
  dat$attr$last.neg.test.rct <- last.neg.test.rct
  dat$attr$last.neg.test.uct <- last.neg.test.uct
  dat$attr$diag.status.ct <- diag.status.ct
  dat$attr$last.diag.time.ct <- last.diag.time.ct
  dat$attr$tt.traj.ct.hivneg <- tt.traj.ct.hivneg
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

