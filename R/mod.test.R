
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

  if (at < dat$param$stitest.start) {

    ## Background screening ----------------------------------------------
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
    lastdiag.time.gc <- dat$attr$lastdiag.time.gc
    lastdiag.time.ct <- dat$attr$lastdiag.time.ct
    lastdiag.time.syph <- dat$attr$lastdiag.time.syph

    # Parameters
    tst.rect.sti.rr <- dat$param$tst.rect.sti.rr
    asympt.screen.prob <- dat$param$asympt.screen.prob

    # Eligibility for diagnosis of asymptomatic infection (pre-intervention)
    # Syphilis
    screen.elig.syph <- which((diag.status.syph == 0 | is.na(diag.status.syph)))
    screen.rates.syph <- rep(asympt.screen.prob, length(screen.elig.syph))
    screen.syph <- screen.elig.syph[rbinom(length(screen.elig.syph), 1, screen.rates.syph) == 1]

    # CT
    screen.elig.ct <- which((diag.status.ct == 0 | is.na(diag.status.ct)))
    screen.rates.ct <- rep(asympt.screen.prob, length(screen.elig.ct))
    screen.ct <- screen.elig.ct[rbinom(length(screen.elig.ct), 1, screen.rates.ct) == 1]

    # GC
    screen.elig.gc <- which((diag.status.gc == 0 | is.na(diag.status.gc)))
    screen.rates.gc <- rep(asympt.screen.prob, length(screen.elig.gc))
    screen.gc <- screen.elig.gc[rbinom(length(screen.elig.gc), 1, screen.rates.gc) == 1]

    # Syphilis screening
    screen.syph.pos <- screen.syph[syphilis[screen.syph] == 1 &
                                   stage.syph[screen.syph] %in% c(2, 3, 4, 5, 6, 7)]
    screen.syph.neg <- setdiff(screen.syph, screen.syph.pos)

    # GC screening
    screen.rgc <- screen.gc[role.class[screen.gc] %in% c("R", "V")]
    screen.rgc <- sample(screen.rgc, tst.rect.sti.rr * length(screen.rgc))
    screen.ugc <- screen.gc[role.class[screen.gc] %in% c("I", "V")]
    screen.rgc.pos <- screen.rgc[rGC[screen.rgc] == 1]
    screen.ugc.pos <- screen.ugc[uGC[screen.ugc] == 1]
    screen.rgc.neg <- setdiff(screen.rgc, screen.rgc.pos)
    screen.ugc.neg <- setdiff(screen.ugc, screen.ugc.pos)
    screen.gc.pos <- unique(c(screen.rgc.pos, screen.ugc.pos))

    # CT screening
    screen.rct <- screen.ct[role.class[screen.gc] %in% c("R", "V")]
    screen.rct <- sample(screen.rct, tst.rect.sti.rr * length(screen.rct))
    screen.uct <- screen.ct[role.class[screen.ct] %in% c("I", "V")]
    screen.rct.pos <- screen.rct[rCT[screen.rct] == 1]
    screen.uct.pos <- screen.uct[uCT[screen.uct] == 1]
    screen.rct.neg <- setdiff(screen.rct, screen.rct.pos)
    screen.uct.neg <- setdiff(screen.uct, screen.uct.pos)
    screen.ct.pos <- unique(c(screen.rct.pos, screen.uct.pos))

    # Syphilis Attributes
    last.neg.test.syph[screen.syph.neg] <- at
    last.neg.test.syph[screen.syph.pos] <- NA
    diag.status.syph[screen.syph.pos] <- 1
    lastdiag.time.syph[screen.syph.pos] <- at

    # GC Attributes
    last.neg.test.rgc[screen.rgc.neg] <- at
    last.neg.test.ugc[screen.ugc.neg] <- at
    last.neg.test.rgc[screen.rgc.pos] <- NA
    last.neg.test.ugc[screen.ugc.pos] <- NA
    diag.status.gc[screen.gc.pos] <- 1
    lastdiag.time.gc[screen.gc.pos] <- at

    # CT Attributes
    last.neg.test.rct[screen.rct.neg] <- at
    last.neg.test.uct[screen.uct.neg] <- at
    last.neg.test.rct[screen.rct.pos] <- NA
    last.neg.test.uct[screen.uct.pos] <- NA
    diag.status.ct[screen.ct.pos] <- 1
    lastdiag.time.ct[screen.ct.pos] <- at

    ## Output
    # Syphilis Attributes
    dat$attr$last.neg.test.syph <- last.neg.test.syph
    dat$attr$diag.status.syph <- diag.status.syph
    dat$attr$lastdiag.time.syph <- lastdiag.time.syph

    # GC Attributes
    dat$attr$last.neg.test.rgc <- last.neg.test.rgc
    dat$attr$last.neg.test.ugc <- last.neg.test.ugc
    dat$attr$diag.status.gc <- diag.status.gc
    dat$attr$lastdiag.time.gc <- lastdiag.time.gc

    # CT Attributes
    dat$attr$last.neg.test.rct <- last.neg.test.rct
    dat$attr$last.neg.test.uct <- last.neg.test.uct
    dat$attr$diag.status.ct <- diag.status.ct
    dat$attr$lastdiag.time.ct <- lastdiag.time.ct

    # Number of tests for asymptomatic
    dat$epi$rGCasympttests[at] <- length(screen.rgc)
    dat$epi$uGCasympttests[at] <- length(screen.ugc)
    dat$epi$GCasympttests[at] <- length(screen.rgc) + length(screen.ugc)

    dat$epi$rGCasympttests.pos[at] <- length(screen.rgc.pos)
    dat$epi$uGCasympttests.pos[at] <- length(screen.ugc.pos)
    dat$epi$GCasympttests.pos[at] <- length(screen.rgc.pos) + length(screen.ugc.pos)

    dat$epi$rCTasympttests[at] <- length(screen.rct)
    dat$epi$uCTasympttests[at] <- length(screen.uct)
    dat$epi$CTasympttests[at] <- length(screen.rct) + length(screen.uct)

    dat$epi$rCTasympttests.pos[at] <- length(screen.rct.pos)
    dat$epi$uCTasympttests.pos[at] <- length(screen.uct.pos)
    dat$epi$CTasympttests.pos[at] <- length(screen.rct.pos) + length(screen.uct.pos)

    dat$epi$syphasympttests[at] <- length(screen.syph)
    dat$epi$syphasympttests.pos[at] <- length(screen.syph.pos)

    return(dat)
  }

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
  lastdiag.time.gc <- dat$attr$lastdiag.time.gc
  lastdiag.time.ct <- dat$attr$lastdiag.time.ct
  lastdiag.time.syph <- dat$attr$lastdiag.time.syph
  race <- dat$attr$race
  prepStat <- dat$attr$prepStat

  # Parameters
  stianntest.coverage <- dat$param$stianntest.coverage
  stianntest.cov.rate <- dat$param$stianntest.cov.rate
  stihighrisktest.coverage <- dat$param$stihighrisktest.coverage
  stihighrisktest.cov.rate <- dat$param$stihighrisktest.cov.rate
  testing.pattern.sti <- dat$param$testing.pattern.sti
  stitest.active.int <- dat$param$stitest.active.int
  sti.highrisktest.int <- dat$param$sti.highrisktest.int

  # Eligibility and trajectory
  # Base eligibility
  idsEligTest <- which(race %in% c("B", "W"))

  # Annual indications- sexually active in last year
  stitestind1 <- dat$attr$stitest.ind.active
  idsactive <- intersect(which(stitestind1 == 1), idsEligTest)

  # STI testing higher-risk eligibility scenarios
  stitestind2 <- dat$attr$stitest.ind.recentpartners
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
  stianntestCov <- sum(tt.traj.ct == 1, na.rm = TRUE) / length(setdiff(idsactive, idshighrisk))
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

  ## Testing
  tsincelntst.syph <- at - dat$attr$last.neg.test.syph
  tsincelntst.syph[is.na(tsincelntst.syph)] <- at - dat$attr$arrival.time[is.na(tsincelntst.syph)]

  tsincelntst.rgc <- at - dat$attr$last.neg.test.rgc
  tsincelntst.ugc <- at - dat$attr$last.neg.test.ugc
  tsincelntst.rgc[is.na(tsincelntst.rgc)] <- at - dat$attr$arrival.time[is.na(tsincelntst.rgc)]
  tsincelntst.ugc[is.na(tsincelntst.ugc)] <- at - dat$attr$arrival.time[is.na(tsincelntst.ugc)]
  tsincelntst.gc <- pmin(tsincelntst.rgc, tsincelntst.ugc)

  tsincelntst.rct <- at - dat$attr$last.neg.test.rct
  tsincelntst.uct <- at - dat$attr$last.neg.test.uct
  tsincelntst.rct[is.na(tsincelntst.rct)] <- at - dat$attr$arrival.time[is.na(tsincelntst.rct)]
  tsincelntst.uct[is.na(tsincelntst.uct)] <- at - dat$attr$arrival.time[is.na(tsincelntst.uct)]
  tsincelntst.ct <- pmin(tsincelntst.rct, tsincelntst.uct)

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
                                        tsincelntst.syph >= 2*(stitest.active.int) &
                                        prepStat == 0)
    tst.syph.highrisk.interval <- which(tt.traj.syph == 2 &
                                          (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                          tsincelntst.syph >= 2*(sti.highrisktest.int) &
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
                                      tsincelntst.gc >= 2*(stitest.active.int) &
                                      prepStat == 0)
    tst.gc.highrisk.interval <- which(tt.traj.gc == 2 &
                                        (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                        tsincelntst.gc >= 2*(sti.highrisktest.int) &
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
                                      tsincelntst.ct >= 2*(stitest.active.int) &
                                      prepStat == 0)
    tst.ct.highrisk.interval <- which(tt.traj.ct == 2 &
                                        (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                        tsincelntst.ct >= 2*(sti.highrisktest.int) &
                                        prepStat == 0)
    tst.ct.nprep <- c(tst.ct.annual.interval, tst.ct.highrisk.interval)
  }

  # Syphilis non-PrEP testing
  tst.syph.pos <- tst.syph.nprep[syphilis[tst.syph.nprep] == 1 &
                                   stage.syph[tst.syph.nprep] %in% c(2, 3, 4, 5, 6, 7)]
  tst.syph.neg <- setdiff(tst.syph.nprep, tst.syph.pos)

  # GC non-PrEP testing
  tst.rgc <- tst.gc.nprep[role.class[tst.gc.nprep] %in% c("R", "V")]
  tst.rgc <- sample(tst.rgc, tst.rect.sti.rr * length(tst.rgc))
  tst.ugc <- tst.gc.nprep[role.class[tst.gc.nprep] %in% c("I", "V")]
  tst.rgc.pos <- tst.rgc[rGC[tst.rgc] == 1]
  tst.ugc.pos <- tst.ugc[uGC[tst.ugc] == 1]
  tst.rgc.neg <- setdiff(tst.rgc, tst.rgc.pos)
  tst.ugc.neg <- setdiff(tst.ugc, tst.ugc.pos)
  tst.gc.pos <- unique(c(tst.rgc.pos, tst.ugc.pos))

  # CT non-PrEP testing
  tst.rct <- tst.ct.nprep[role.class[tst.ct.nprep] %in% c("R", "V")]
  tst.rct <- sample(tst.rct, tst.rect.sti.rr * length(tst.rct))
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
  lastdiag.time.syph[tst.syph.pos] <- at

  # GC Attributes
  last.neg.test.rgc[tst.rgc.neg] <- at
  last.neg.test.ugc[tst.ugc.neg] <- at
  last.neg.test.rgc[tst.rgc.pos] <- NA
  last.neg.test.ugc[tst.ugc.pos] <- NA
  diag.status.gc[tst.gc.pos] <- 1
  lastdiag.time.gc[tst.gc.pos] <- at

  # CT Attributes
  last.neg.test.rct[tst.rct.neg] <- at
  last.neg.test.uct[tst.uct.neg] <- at
  last.neg.test.rct[tst.rct.pos] <- NA
  last.neg.test.uct[tst.uct.pos] <- NA
  diag.status.ct[tst.ct.pos] <- 1
  lastdiag.time.ct[tst.ct.pos] <- at

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

  ## Output -----------------------------------------------------------------

  # Attributes

  # Stoppage attributes
  dat$attr$stihighrisktestLastElig[idsnothighriskelig] <- at
  dat$attr$stianntestLastElig[idsnotactiveelig] <- at

  # Syphilis Attributes
  dat$attr$last.neg.test.syph <- last.neg.test.syph
  dat$attr$diag.status.syph <- diag.status.syph
  dat$attr$lastdiag.time.syph <- lastdiag.time.syph
  dat$attr$tt.traj.syph <- tt.traj.syph

  # GC Attributes
  dat$attr$last.neg.test.rgc <- last.neg.test.rgc
  dat$attr$last.neg.test.ugc <- last.neg.test.ugc
  dat$attr$diag.status.gc <- diag.status.gc
  dat$attr$lastdiag.time.gc <- lastdiag.time.gc
  dat$attr$tt.traj.gc <- tt.traj.gc

  # CT Attributes
  dat$attr$last.neg.test.rct <- last.neg.test.rct
  dat$attr$last.neg.test.uct <- last.neg.test.uct
  dat$attr$diag.status.ct <- diag.status.ct
  dat$attr$lastdiag.time.ct <- lastdiag.time.ct
  dat$attr$tt.traj.ct <- tt.traj.ct


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

