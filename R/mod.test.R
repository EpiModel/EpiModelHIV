
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
  role.class <- dat$attr$role.class
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

  # Assign new STI treatment trajectory (HIV-pos) for diagnosed who will be on ART
  tt.traj.new <- tst.pos[which(tt.traj[tst.pos] %in% c(3,4))]

  # Attributes
  dat$attr$last.neg.test[tst.neg] <- at
  dat$attr$diag.status[tst.pos] <- 1
  dat$attr$tt.traj.syph.hivneg[tt.traj.new] <- NA
  dat$attr$tt.traj.gc.hivneg[tt.traj.new] <- NA
  dat$attr$tt.traj.ct.hivneg[tt.traj.new] <- NA
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

  # 1. Setup ----------------------------------------------------------------

  # Attributes
  tt.traj <- dat$attr$tt.traj
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

  prepStat <- dat$attr$prepStat
  stitestind1 <- dat$attr$stitest.ind.active
  stitestind2 <- dat$attr$stitest.ind.recentpartners

  # Parameters
  stianntest.ct.hivneg.coverage <- dat$param$stianntest.ct.hivneg.coverage
  stianntest.syph.hivneg.coverage <- dat$param$stianntest.syph.hivneg.coverage
  stihighrisktest.ct.hivneg.coverage <- dat$param$stihighrisktest.ct.hivneg.coverage
  stihighrisktest.syph.hivneg.coverage <- dat$param$stihighrisktest.syph.hivneg.coverage

  stianntest.ct.hivpos.coverage <- dat$param$stianntest.ct.hivpos.coverage
  stianntest.syph.hivpos.coverage <- dat$param$stianntest.syph.hivpos.coverage
  stihighrisktest.ct.hivpos.coverage <- dat$param$stihighrisktest.ct.hivpos.coverage
  stihighrisktest.syph.hivpos.coverage <- dat$param$stihighrisktest.syph.hivpos.coverage

  testing.pattern.sti <- dat$param$testing.pattern.sti
  stitest.active.int <- dat$param$stitest.active.int
  sti.highrisktest.int <- dat$param$sti.highrisktest.int
  tst.rect.sti.rr <- dat$param$tst.rect.sti.rr

  # Eligibility and trajectory
  # Annual indications- sexually active in last year
  idsactive.hivpos <- which(stitestind1 == 1 & diag.status == 1 & tt.traj %in% 3:4)
  idsactive.hivneg <- setdiff((which(stitestind1 == 1)), idsactive.hivpos)

  # STI testing higher-risk eligibility scenarios
  idshighrisk.hivpos <- which(stitestind2 == 1 & diag.status == 1 & tt.traj %in% 3:4)
  idshighrisk.hivneg <- setdiff((which(stitestind2 == 1)), idshighrisk.hivpos)


  # 2. Indication Trajectory Stoppage ---------------------------------------

  # Reduce testing trajectory to NA if no longer indicated for more frequent high-risk testing
  idsnothighriskelig.hivpos <- which((tt.traj.syph.hivpos == 2 |
                                        tt.traj.gc.hivpos == 2 |
                                        tt.traj.ct.hivpos == 2) & stitestind2 != 1)
  tt.traj.syph.hivpos[idsnothighriskelig.hivpos] <-
    tt.traj.gc.hivpos[idsnothighriskelig.hivpos] <-
    tt.traj.ct.hivpos[idsnothighriskelig.hivpos] <-
    NA
  idsnothighriskelig.hivneg <- which((tt.traj.syph.hivneg == 2 |
                                        tt.traj.gc.hivneg == 2 |
                                        tt.traj.ct.hivneg == 2) & stitestind2 != 1)
  tt.traj.syph.hivneg[idsnothighriskelig.hivneg] <-
    tt.traj.gc.hivneg[idsnothighriskelig.hivneg] <-
    tt.traj.ct.hivneg[idsnothighriskelig.hivneg] <-
    NA

  # Reduce testing trajectory to NA if no longer indicated for lower-risk testing
  idsnotactiveelig.hivpos <- which((tt.traj.syph.hivpos == 1 |
                                      tt.traj.gc.hivpos == 1 |
                                      tt.traj.ct.hivpos == 1) & stitestind1 != 1)
  tt.traj.syph.hivpos[idsnotactiveelig.hivpos] <-
    tt.traj.gc.hivpos[idsnotactiveelig.hivpos] <-
    tt.traj.ct.hivpos[idsnotactiveelig.hivpos] <-
    NA
  idsnotactiveelig.hivneg <- which((tt.traj.syph.hivneg == 1 |
                                      tt.traj.gc.hivneg == 1 |
                                      tt.traj.ct.hivneg == 1) & stitestind1 != 1)
  tt.traj.syph.hivneg[idsnotactiveelig.hivneg] <-
    tt.traj.gc.hivneg[idsnotactiveelig.hivneg] <-
    tt.traj.ct.hivneg[idsnotactiveelig.hivneg] <-
    NA



  # 3. Screening for non-HIV Diagnosed --------------------------------------

  ## 3a. High-Risk Testing ##

  ## Assume coverage and people are same for NG and CT - correlation
  idsEligSt <- idshighrisk.hivneg
  nEligSt <- length(idshighrisk.hivneg)

  ## Evaluate existing coverage
  stihighrisktestCov.ct <- sum(tt.traj.ct.hivneg == 2, na.rm = TRUE) / nEligSt
  stihighrisktestCov.ct <- ifelse(is.nan(stihighrisktestCov.ct), 0, stihighrisktestCov.ct)
  stihighrisktestCov.gc <- sum(tt.traj.gc.hivneg == 2, na.rm = TRUE) / nEligSt
  stihighrisktestCov.gc <- ifelse(is.nan(stihighrisktestCov.gc), 0, stihighrisktestCov.gc)
  stihighrisktestCov.syph <- sum(tt.traj.syph.hivneg == 2, na.rm = TRUE) / nEligSt
  stihighrisktestCov.syph <- ifelse(is.nan(stihighrisktestCov.syph), 0, stihighrisktestCov.syph)

  ## Count how many are eligible to start a new trajectory
  nStart.gcct <- max(0, min(nEligSt, round((stihighrisktest.ct.hivneg.coverage - stihighrisktestCov.ct) *
                                           length(idsEligSt))))
  nStart.syph <- max(0, min(nEligSt, round((stihighrisktest.syph.hivneg.coverage - stihighrisktestCov.syph) *
                                             length(idsEligSt))))

  ## Sample individuals
  idsStart.gcct <- idsStart.syph <- NULL
  if (nStart.gcct > 0) {
    idsStart.gcct <- ssample(idsEligSt, nStart.gcct)
  }
  if (nStart.syph > 0) {
    idsStart.syph <- ssample(idsEligSt, nStart.syph)
  }

  ## Update testing trajectory for higher-risk
  if (length(idsStart.gcct) > 0) {
    tt.traj.gc.hivneg[idsStart.gcct] <- 2
    tt.traj.ct.hivneg[idsStart.gcct] <- 2
  }
  if (length(idsStart.syph) > 0) {
    tt.traj.syph.hivneg[idsStart.syph] <- 2
  }

  ## 3b. Sexually Active (non-HR) Testing ##

  ## Assume coverage and people are same for NG and CT - correlation
  idsEligSt.gcct <- setdiff(idsactive.hivneg, which(tt.traj.ct.hivneg == 2))
  idsEligSt.syph <- setdiff(idsactive.hivneg, which(tt.traj.syph.hivneg == 2))
  nEligSt.gcct <- length(idsEligSt.gcct)
  nEligSt.syph <- length(idsEligSt.syph)

  ## Evaluate existing coverage
  stianntestCov.ct <- sum(tt.traj.ct.hivneg == 1, na.rm = TRUE) / nEligSt.gcct
  stianntestCov.ct <- ifelse(is.nan(stianntestCov.ct), 0, stianntestCov.ct)
  stianntestCov.gc <- sum(tt.traj.gc.hivneg == 1, na.rm = TRUE) / nEligSt.gcct
  stianntestCov.gc <- ifelse(is.nan(stianntestCov.gc), 0, stianntestCov.gc)
  stianntestCov.syph <- sum(tt.traj.syph.hivneg == 1, na.rm = TRUE) / nEligSt.syph
  stianntestCov.syph <- ifelse(is.nan(stianntestCov.syph), 0, stianntestCov.syph)

  ## Count how many are eligible to start a new trajectory
  nStart.gcct <- max(0, min(nEligSt.gcct,
                          round((stianntest.ct.hivneg.coverage - stianntestCov.ct) * nEligSt.gcct)))
  nStart.syph <- max(0, min(nEligSt.syph,
                            round((stianntest.syph.hivneg.coverage - stianntestCov.syph) * nEligSt.syph)))

  ## Sample individuals
  idsStart.gcct <- idsStart.syph <- NULL
  if (nStart.gcct > 0) {
    idsStart.gcct <- ssample(idsEligSt.gcct, nStart.gcct)
  }
  if (nStart.syph > 0) {
    idsStart.syph <- ssample(idsEligSt.syph, nStart.syph)
  }

  ## Update testing trajectory for lower-risk
  if (length(idsStart.gcct) > 0) {
    tt.traj.ct.hivneg[idsStart.gcct] <- 1
    tt.traj.gc.hivneg[idsStart.gcct] <- 1
  }
  if (length(idsStart.syph) > 0) {
    tt.traj.syph.hivneg[idsStart.syph] <- 1
  }


  ## 3c. Asymptomatic Screening ##

  # Syphilis
  if (testing.pattern.sti == "interval" ) {
    tst.syph.annual.interval <- which(tt.traj.syph.hivneg == 1 &
                                        (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                        (tsinceltst.syph >= stitest.active.int) &
                                        prepStat == 0)
    tst.syph.highrisk.interval <- which(tt.traj.syph.hivneg == 2 &
                                          (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                          (tsinceltst.syph >= sti.highrisktest.int) &
                                          prepStat == 0)
    tst.syph.nprep.hivneg <- c(tst.syph.annual.interval, tst.syph.highrisk.interval)
  }

  # GC
  if (testing.pattern.sti == "interval" ) {
    tst.gc.annual.interval <- which(tt.traj.gc.hivneg == 1 &
                                      (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                      (tsinceltst.gc >= stitest.active.int |
                                      tsinceltst.ct >= stitest.active.int) &
                                      prepStat == 0)
    tst.gc.highrisk.interval <- which(tt.traj.gc.hivneg == 2 &
                                        (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                        (tsinceltst.gc >= sti.highrisktest.int |
                                        tsinceltst.ct >= sti.highrisktest.int) &
                                        prepStat == 0)
    tst.gc.nprep.hivneg <- c(tst.gc.annual.interval, tst.gc.highrisk.interval)
  }

  # CT
  if (testing.pattern.sti == "interval" ) {
    tst.ct.annual.interval <- which(tt.traj.ct.hivneg == 1 &
                                      (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                      (tsinceltst.gc >= stitest.active.int |
                                      tsinceltst.ct >= stitest.active.int) &
                                      prepStat == 0)
    tst.ct.highrisk.interval <- which(tt.traj.ct.hivneg == 2 &
                                        (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                        (tsinceltst.gc >= sti.highrisktest.int |
                                        tsinceltst.ct >= sti.highrisktest.int) &
                                        prepStat == 0)
    tst.ct.nprep.hivneg <- c(tst.ct.annual.interval, tst.ct.highrisk.interval)
  }

  # Syphilis non-PrEP testing
  tst.syph.pos.hivneg <- tst.syph.nprep.hivneg[which(syphilis[tst.syph.nprep.hivneg] == 1 &
                                   stage.syph[tst.syph.nprep.hivneg] %in% 2:6)]
  tst.syph.neg.hivneg <- setdiff(tst.syph.nprep.hivneg, tst.syph.pos.hivneg)
  tst.earlysyph.pos.hivneg <- tst.syph.nprep.hivneg[which(syphilis[tst.syph.nprep.hivneg] == 1 &
                                   stage.syph[tst.syph.nprep.hivneg] %in% 2:3)]
  tst.latesyph.pos.hivneg <- tst.syph.nprep.hivneg[which(syphilis[tst.syph.nprep.hivneg] == 1 &
                                   stage.syph[tst.syph.nprep.hivneg] %in% 4:6)]

  # GC non-PrEP testing
  tst.rgc.hivneg <- tst.gc.nprep.hivneg[which(role.class[tst.gc.nprep.hivneg] %in% c("R", "V"))]
  tst.rgc.hivneg <- sample(tst.rgc.hivneg, tst.rect.sti.rr * length(tst.rgc.hivneg))
  tst.ugc.hivneg <- tst.gc.nprep.hivneg[which(role.class[tst.gc.nprep.hivneg] %in% c("I", "V"))]
  tst.rgc.pos.hivneg <- tst.rgc.hivneg[which(rGC[tst.rgc.hivneg] == 1)]
  tst.ugc.pos.hivneg <- tst.ugc.hivneg[which(uGC[tst.ugc.hivneg] == 1)]
  tst.rgc.neg.hivneg <- setdiff(tst.rgc.hivneg, tst.rgc.pos.hivneg)
  tst.ugc.neg.hivneg <- setdiff(tst.ugc.hivneg, tst.ugc.pos.hivneg)
  tst.gc.pos.hivneg <- unique(c(tst.rgc.pos.hivneg, tst.ugc.pos.hivneg))

  # CT non-PrEP testing
  tst.rct.hivneg <- tst.ct.nprep.hivneg[which(role.class[tst.ct.nprep.hivneg] %in% c("R", "V"))]
  tst.rct.hivneg <- sample(tst.rct.hivneg, tst.rect.sti.rr * length(tst.rct.hivneg))
  tst.uct.hivneg <- tst.ct.nprep.hivneg[which(role.class[tst.ct.nprep.hivneg] %in% c("I", "V"))]
  tst.rct.pos.hivneg <- tst.rct.hivneg[which(rCT[tst.rct.hivneg] == 1)]
  tst.uct.pos.hivneg <- tst.uct.hivneg[which(uCT[tst.uct.hivneg] == 1)]
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


  # 4. Screening for HIV-Diagnosed ------------------------------------------

  ## 4a. High-Risk Testing ##

  ## Assume coverage and people are same for NG and CT - correlation
  idsEligSt <- idshighrisk.hivpos
  nEligSt <- length(idshighrisk.hivpos)

  ## Evaluate existing coverage
  stihighrisktestCov.ct <- sum(tt.traj.ct.hivpos == 2, na.rm = TRUE) / nEligSt
  stihighrisktestCov.ct <- ifelse(is.nan(stihighrisktestCov.ct), 0, stihighrisktestCov.ct)
  stihighrisktestCov.gc <- sum(tt.traj.gc.hivpos == 2, na.rm = TRUE) / nEligSt
  stihighrisktestCov.gc <- ifelse(is.nan(stihighrisktestCov.gc), 0, stihighrisktestCov.gc)
  stihighrisktestCov.syph <- sum(tt.traj.syph.hivpos == 2, na.rm = TRUE) / nEligSt
  stihighrisktestCov.syph <- ifelse(is.nan(stihighrisktestCov.syph), 0, stihighrisktestCov.syph)

  ## Count how many are eligible to start a new trajectory
  nStart.gcct <- max(0, min(nEligSt, round((stihighrisktest.ct.hivpos.coverage - stihighrisktestCov.ct) *
                                             length(idshighrisk.hivpos))))
  nStart.syph <- max(0, min(nEligSt, round((stihighrisktest.syph.hivpos.coverage - stihighrisktestCov.syph) *
                                             length(idshighrisk.hivpos))))

  ## Sample individuals
  idsStart.gcct <- idsStart.syph <- NULL
  if (nStart.gcct > 0) {
    idsStart.gcct <- ssample(idsEligSt, nStart.gcct)
  }
  if (nStart.syph > 0) {
    idsStart.syph <- ssample(idsEligSt, nStart.syph)
  }

  ## Update testing trajectory for higher-risk
  if (length(idsStart.gcct) > 0) {
    tt.traj.ct.hivpos[idsStart.gcct] <- 2
    tt.traj.gc.hivpos[idsStart.gcct] <- 2
  }
  if (length(idsStart.syph) > 0) {
    tt.traj.syph.hivpos[idsStart.syph] <- 2
  }

  ## 4b. Sexually Active (non-HR) Testing ##

  ## Assume coverage and people are same for NG and CT - correlation
  idsEligSt.gcct <- setdiff(idsactive.hivpos, which(tt.traj.ct.hivpos == 2))
  idsEligSt.syph <- setdiff(idsactive.hivpos, which(tt.traj.syph.hivpos == 2))
  nEligSt.gcct <- length(idsEligSt.gcct)
  nEligSt.syph <- length(idsEligSt.syph)

  ## Evaluate existing coverage
  stianntestCov.ct <- sum(tt.traj.ct.hivpos == 1, na.rm = TRUE) / nEligSt.gcct
  stianntestCov.ct <- ifelse(is.nan(stianntestCov.ct), 0, stianntestCov.ct)
  stianntestCov.gc <- sum(tt.traj.gc.hivpos == 1, na.rm = TRUE) / nEligSt.gcct
  stianntestCov.gc <- ifelse(is.nan(stianntestCov.gc), 0, stianntestCov.gc)
  stianntestCov.syph <- sum(tt.traj.syph.hivpos == 1, na.rm = TRUE) / nEligSt.syph
  stianntestCov.syph <- ifelse(is.nan(stianntestCov.syph), 0, stianntestCov.syph)

  ## Count how many are eligible to start a new trajectory
  nStart.gcct <- max(0, min(nEligSt.gcct,
                            round((stianntest.ct.hivpos.coverage - stianntestCov.ct) * nEligSt.gcct)))
  nStart.syph <- max(0, min(nEligSt.syph,
                            round((stianntest.syph.hivpos.coverage - stianntestCov.syph) * nEligSt.syph)))

  ## Sample individuals
  idsStart.gcct <- idsStart.syph <- NULL
  if (nStart.gcct > 0) {
    idsStart.gcct <- ssample(idsEligSt.gcct, nStart.gcct)
  }
  if (nStart.syph > 0) {
    idsStart.syph <- ssample(idsEligSt.syph, nStart.syph)
  }

  ## Update testing trajectory for lower-risk
  if (length(idsStart.gcct) > 0) {
    tt.traj.ct.hivpos[idsStart.gcct] <- 1
    tt.traj.gc.hivpos[idsStart.gcct] <- 1
  }
  if (length(idsStart.syph) > 0) {
    tt.traj.syph.hivpos[idsStart.syph] <- 1
  }

  ## 4c. Asymptomatic screening ##

  ## Syphilis
  if (testing.pattern.sti == "interval" ) {
    tst.syph.annual.interval <- which((tt.traj.syph.hivpos == 1 &
                                        (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                         (tsinceltst.syph >= stitest.active.int) &
                                        prepStat == 0))
    tst.syph.highrisk.interval <- which((tt.traj.syph.hivpos == 2 &
                                          (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                          (tsinceltst.syph >= sti.highrisktest.int) &
                                          prepStat == 0))
    tst.syph.nprep.hivpos <- c(tst.syph.annual.interval, tst.syph.highrisk.interval)
  }

  ## GC
  if (testing.pattern.sti == "interval" ) {
    tst.gc.annual.interval <- which((tt.traj.gc.hivpos == 1 &
                                      (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                       (tsinceltst.gc >= stitest.active.int |
                                       tsinceltst.ct >= stitest.active.int) &
                                      prepStat == 0))
    tst.gc.highrisk.interval <- which((tt.traj.gc.hivpos == 2 &
                                        (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                         (tsinceltst.gc >= sti.highrisktest.int |
                                         tsinceltst.ct >= sti.highrisktest.int) &
                                        prepStat == 0))
    tst.gc.nprep.hivpos <- c(tst.gc.annual.interval, tst.gc.highrisk.interval)
  }

  ## CT
  if (testing.pattern.sti == "interval" ) {
    tst.ct.annual.interval <- which((tt.traj.ct.hivpos == 1 &
                                      (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                       (tsinceltst.gc >= stitest.active.int |
                                       tsinceltst.ct >= stitest.active.int) &
                                      prepStat == 0))
    tst.ct.highrisk.interval <- which((tt.traj.ct.hivpos == 2 &
                                         (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                         (tsinceltst.gc >= sti.highrisktest.int |
                                         tsinceltst.ct >= sti.highrisktest.int) &
                                         prepStat == 0))
    tst.ct.nprep.hivpos <- c(tst.ct.annual.interval, tst.ct.highrisk.interval)
  }

  # Syphilis non-PrEP testing
  tst.syph.pos.hivpos <- tst.syph.nprep.hivpos[which(syphilis[tst.syph.nprep.hivpos] == 1 &
                                                 stage.syph[tst.syph.nprep.hivpos] %in% c(2, 3, 4, 5, 6))]
  tst.syph.neg.hivpos <- setdiff(tst.syph.nprep.hivpos, tst.syph.pos.hivpos)
  tst.earlysyph.pos.hivpos <- tst.syph.nprep.hivpos[which(syphilis[tst.syph.nprep.hivpos] == 1 &
                                                      stage.syph[tst.syph.nprep.hivpos] %in% c(2, 3))]
  tst.latesyph.pos.hivpos <- tst.syph.nprep.hivpos[which(syphilis[tst.syph.nprep.hivpos] == 1 &
                                                     stage.syph[tst.syph.nprep.hivpos] %in% c(4, 5, 6))]

  # GC non-PrEP testing
  tst.rgc.hivpos <- tst.gc.nprep.hivpos[which(role.class[tst.gc.nprep.hivpos] %in% c("R", "V"))]
  tst.rgc.hivpos <- sample(tst.rgc.hivpos, tst.rect.sti.rr * length(tst.rgc.hivpos))
  tst.ugc.hivpos <- tst.gc.nprep.hivpos[which(role.class[tst.gc.nprep.hivpos] %in% c("I", "V"))]
  tst.rgc.pos.hivpos <- tst.rgc.hivpos[which(rGC[tst.rgc.hivpos] == 1)]
  tst.ugc.pos.hivpos <- tst.ugc.hivpos[which(uGC[tst.ugc.hivpos] == 1)]
  tst.rgc.neg.hivpos <- setdiff(tst.rgc.hivpos, tst.rgc.pos.hivpos)
  tst.ugc.neg.hivpos <- setdiff(tst.ugc.hivpos, tst.ugc.pos.hivpos)
  tst.gc.pos.hivpos <- unique(c(tst.rgc.pos.hivpos, tst.ugc.pos.hivpos))

  # CT non-PrEP testing
  tst.rct.hivpos <- tst.ct.nprep.hivpos[which(role.class[tst.ct.nprep.hivpos] %in% c("R", "V"))]
  tst.rct.hivpos <- sample(tst.rct.hivpos, tst.rect.sti.rr * length(tst.rct.hivpos))
  tst.uct.hivpos <- tst.ct.nprep.hivpos[which(role.class[tst.ct.nprep.hivpos] %in% c("I", "V"))]
  tst.rct.pos.hivpos <- tst.rct.hivpos[which(rCT[tst.rct.hivpos] == 1)]
  tst.uct.pos.hivpos <- tst.uct.hivpos[which(uCT[tst.uct.hivpos] == 1)]
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

  # Number of people on each testing trajectory by serostatus
  dat$epi$tt.traj.syph1.hivneg[at] <- length(which(tt.traj.syph.hivneg == 1))
  dat$epi$tt.traj.gc1.hivneg[at] <- length(which(tt.traj.gc.hivneg == 1))
  dat$epi$tt.traj.ct1.hivneg[at] <- length(which(tt.traj.ct.hivneg == 1))
  dat$epi$tt.traj.syph2.hivneg[at] <- length(which(tt.traj.syph.hivneg == 2))
  dat$epi$tt.traj.gc2.hivneg[at] <- length(which(tt.traj.gc.hivneg == 2))
  dat$epi$tt.traj.ct2.hivneg[at] <- length(which(tt.traj.ct.hivneg == 2))
  dat$epi$tt.traj.syph1.hivpos[at] <- length(which(tt.traj.syph.hivpos == 1))
  dat$epi$tt.traj.gc1.hivpos[at] <- length(which(tt.traj.gc.hivpos == 1))
  dat$epi$tt.traj.ct1.hivpos[at] <- length(which(tt.traj.ct.hivpos == 1))
  dat$epi$tt.traj.syph2.hivpos[at] <- length(which(tt.traj.syph.hivpos == 2))
  dat$epi$tt.traj.gc2.hivpos[at] <- length(which(tt.traj.gc.hivpos == 2))
  dat$epi$tt.traj.ct2.hivpos[at] <- length(which(tt.traj.ct.hivpos == 2))

  # Number of people on each testing trajectory
  dat$epi$tt.traj.gc1[at] <- length(which(tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1))
  dat$epi$tt.traj.ct1[at] <- length(which(tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1))
  dat$epi$tt.traj.syph1[at] <- length(which(tt.traj.syph.hivneg == 1 | tt.traj.syph.hivpos == 1))
  dat$epi$tt.traj.gc2[at] <- length(which(tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2))
  dat$epi$tt.traj.ct2[at] <- length(which(tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2))
  dat$epi$tt.traj.syph2[at] <- length(which(tt.traj.syph.hivneg == 2 | tt.traj.syph.hivpos == 2))

  dat$epi$tt.traj.sti1[at] <- length(which(tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1 |
                                            tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1 |
                                            tt.traj.syph.hivneg == 1 | tt.traj.syph.hivpos == 1))
  dat$epi$tt.traj.sti2[at] <- length(which(tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2 |
                                             tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2 |
                                             tt.traj.syph.hivneg == 2 | tt.traj.syph.hivpos == 2))

  # Number of screening tests and number of positive results on screening tests
  # GC overall
  dat$epi$rGCasympttests[at] <- length(tst.rgc.hivpos) +
    length(tst.rgc.hivneg)
  dat$epi$uGCasympttests[at] <- length(tst.ugc.hivpos) +
    length(tst.ugc.hivneg)
  dat$epi$GCasympttests[at] <- length(tst.rgc.hivpos) + length(tst.ugc.hivpos) +
    length(tst.rgc.hivneg) + length(tst.ugc.hivneg)

  # GC positive tests
  dat$epi$rGCasympttests.pos[at] <- length(tst.rgc.pos.hivpos) +
    length(tst.rgc.pos.hivneg)
  dat$epi$uGCasympttests.pos[at] <- length(tst.ugc.pos.hivpos) +
    length(tst.ugc.pos.hivneg)
  dat$epi$GCasympttests.pos[at] <- length(tst.rgc.pos.hivpos) + length(tst.ugc.pos.hivpos) +
    length(tst.rgc.pos.hivneg) + length(tst.ugc.pos.hivneg)

  # CT overall
  dat$epi$rCTasympttests[at] <- length(tst.rct.hivpos) +
    length(tst.rct.hivneg)
  dat$epi$uCTasympttests[at] <- length(tst.uct.hivpos) +
    length(tst.uct.hivneg)
  dat$epi$CTasympttests[at] <- length(tst.rct.hivpos) + length(tst.uct.hivpos) +
    length(tst.rct.hivneg) + length(tst.uct.hivneg)

  # CT positive tests
  dat$epi$rCTasympttests.pos[at] <- length(tst.rct.pos.hivpos) +
    length(tst.rct.pos.hivneg)
  dat$epi$uCTasympttests.pos[at] <- length(tst.uct.pos.hivpos) +
    length(tst.uct.pos.hivneg)
  dat$epi$CTasympttests.pos[at] <- length(tst.rct.pos.hivpos) +
    length(tst.uct.pos.hivpos) + length(tst.rct.pos.hivneg) +
    length(tst.uct.pos.hivneg)

  # Syph (overall, positive, early-stage positive, late-stage positive)
  dat$epi$syphasympttests[at] <- length(tst.syph.nprep.hivpos) +
    length(tst.syph.nprep.hivneg)
  dat$epi$syphasympttests.pos[at] <- length(tst.syph.pos.hivpos) +
    length(tst.syph.pos.hivneg)
  dat$epi$syphearlyasympttests.pos[at] <- length(tst.earlysyph.pos.hivpos) +
    length(tst.earlysyph.pos.hivneg)
  dat$epi$syphlateasympttests.pos[at] <- length(tst.latesyph.pos.hivpos) +
    length(tst.latesyph.pos.hivneg)

  # Overall STI tests and positive STI tests
  dat$epi$stiasympttests[at] <- length(tst.rgc.hivpos) + length(tst.ugc.hivpos) +
    length(tst.rct.hivpos) + length(tst.uct.hivpos) + length(tst.syph.nprep.hivpos) +
    length(tst.rgc.hivneg) + length(tst.ugc.hivneg) +
    length(tst.rct.hivneg) + length(tst.uct.hivneg) + length(tst.syph.nprep.hivneg)
  dat$epi$stiasympttests.pos[at] <- length(tst.rgc.pos.hivpos) + length(tst.ugc.pos.hivpos) +
    length(tst.rct.pos.hivpos) + length(tst.uct.pos.hivpos) + length(tst.syph.pos.hivpos) +
    length(tst.rgc.pos.hivneg) + length(tst.ugc.pos.hivneg) +
    length(tst.rct.pos.hivneg) + length(tst.uct.pos.hivneg) + length(tst.syph.pos.hivneg)

  # Risk group-specific test counters
  # GC lower-risk
  dat$epi$rGCasympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 1)) +
                                        length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 1))
  dat$epi$uGCasympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 1)) +
                                        length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 1))
  dat$epi$GCasympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 1)) +
                                        length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 1)) +
                                        length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 1)) +
                                        length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 1))

  # GC higher-risk
  dat$epi$rGCasympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 2)) +
                                        length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 2))
  dat$epi$uGCasympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 2)) +
                                        length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 2))
  dat$epi$GCasympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 2)) +
                                        length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 2)) +
                                        length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 2)) +
                                        length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 2))

  # CT lower-risk
  dat$epi$rCTasympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 1)) +
                                        length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 1))
  dat$epi$uCTasympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 1)) +
                                        length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 1))
  dat$epi$CTasympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 1)) +
                                        length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 1)) +
                                        length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 1)) +
                                        length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 1))

  # CT higher-risk
  dat$epi$rCTasympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 2)) +
                                        length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 2))
  dat$epi$uCTasympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 2)) +
                                        length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 2))
  dat$epi$CTasympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 2)) +
                                        length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 2)) +
                                        length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 2)) +
                                        length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 2))

  # Syph
  dat$epi$syphasympttests.tttraj1[at] <- length(which(tt.traj.syph.hivpos[tst.syph.nprep.hivpos] == 1)) +
                                          length(which(tt.traj.syph.hivneg[tst.syph.nprep.hivneg] == 1))

  dat$epi$syphasympttests.tttraj2[at] <- length(which(tt.traj.syph.hivpos[tst.syph.nprep.hivpos] == 2)) +
                                          length(which(tt.traj.syph.hivneg[tst.syph.nprep.hivneg] == 2))

  # STI
  dat$epi$stiasympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 1)) +
                                        length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 1)) +
                                        length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 1)) +
                                        length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 1)) +
                                        length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 1)) +
                                        length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 1)) +
                                        length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 1)) +
                                        length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 1)) +
                                        length(which(tt.traj.syph.hivpos[tst.syph.nprep.hivpos] == 1)) +
                                        length(which(tt.traj.syph.hivneg[tst.syph.nprep.hivneg] == 1))

  dat$epi$stiasympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[tst.rgc.hivpos] == 2)) +
                                        length(which(tt.traj.gc.hivneg[tst.rgc.hivneg] == 2)) +
                                        length(which(tt.traj.gc.hivpos[tst.ugc.hivpos] == 2)) +
                                        length(which(tt.traj.gc.hivneg[tst.ugc.hivneg] == 2)) +
                                        length(which(tt.traj.ct.hivpos[tst.rct.hivpos] == 2)) +
                                        length(which(tt.traj.ct.hivneg[tst.rct.hivneg] == 2)) +
                                        length(which(tt.traj.ct.hivpos[tst.uct.hivpos] == 2)) +
                                        length(which(tt.traj.ct.hivneg[tst.uct.hivneg] == 2)) +
                                        length(which(tt.traj.syph.hivpos[tst.syph.nprep.hivpos] == 2)) +
                                        length(which(tt.traj.syph.hivneg[tst.syph.nprep.hivneg] == 2))

  # Update Attributes
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
  dat$attr$tt.traj.syph.hivpos <- tt.traj.syph.hivpos
  dat$attr$time.since.last.test.syph <- tsinceltst.syph

  # GC Attributes
  dat$attr$last.neg.test.rgc <- last.neg.test.rgc
  dat$attr$last.neg.test.ugc <- last.neg.test.ugc
  dat$attr$diag.status.gc <- diag.status.gc
  dat$attr$last.diag.time.gc <- last.diag.time.gc
  dat$attr$tt.traj.gc.hivneg <- tt.traj.gc.hivneg
  dat$attr$tt.traj.gc.hivpos <- tt.traj.gc.hivpos
  dat$attr$time.since.last.test.rgc <- tsinceltst.rgc
  dat$attr$time.since.last.test.ugc <- tsinceltst.ugc

  # CT Attributes
  dat$attr$last.neg.test.rct <- last.neg.test.rct
  dat$attr$last.neg.test.uct <- last.neg.test.uct
  dat$attr$diag.status.ct <- diag.status.ct
  dat$attr$last.diag.time.ct <- last.diag.time.ct
  dat$attr$tt.traj.ct.hivneg <- tt.traj.ct.hivneg
  dat$attr$tt.traj.ct.hivpos <- tt.traj.ct.hivpos
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

