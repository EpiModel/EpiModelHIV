
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

  # STI
  tsinceltst.syph <- dat$attr$time.since.last.test.syph
  tsinceltst.rgc <- dat$attr$time.since.last.test.rgc
  tsinceltst.ugc <- dat$attr$time.since.last.test.ugc
  tsinceltst.rct <- dat$attr$time.since.last.test.rct
  tsinceltst.uct <- dat$attr$time.since.last.test.uct
  tsinceltst.gc <- pmin(tsinceltst.rgc, tsinceltst.ugc)
  tsinceltst.ct <- pmin(tsinceltst.rct, tsinceltst.uct)
  diag.status.gc <- dat$attr$diag.status.gc
  diag.status.ct <- dat$attr$diag.status.ct
  diag.status.syph <- dat$attr$diag.status.syph
  syphilis <- dat$attr$syphilis
  stage.syph <- dat$attr$stage.syph
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT

  # Parameters
  testing.pattern <- dat$param$testing.pattern
  mean.test.B.int <- dat$param$mean.test.B.int
  mean.test.W.int <- dat$param$mean.test.W.int
  twind.int <- dat$param$test.window.int
  sti.correlation.time <- dat$param$sti.correlation.time

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

  # STI testing if not tested for that in last X weeks ------------------------
  if (dat$param$sti.hivdx.correlation == "true") {

    tst.rgc <- tst.all[which(tsinceltst.gc[tst.all] > sti.correlation.time &
                               (is.na(diag.status.gc[tst.all]) | diag.status.gc[tst.all] == 0) &
                         role.class[tst.all] %in% c("R", "V"))]
    tst.ugc <- tst.all[which(tsinceltst.gc[tst.all] > sti.correlation.time &
                                 (is.na(diag.status.gc[tst.all]) | diag.status.gc[tst.all] == 0) &
                           role.class[tst.all] %in% c("I", "V"))]
    tst.rct <- tst.all[which(tsinceltst.ct[tst.all] > sti.correlation.time &
                                 (is.na(diag.status.ct[tst.all]) | diag.status.ct[tst.all] == 0) &
                           role.class[tst.all] %in% c("R", "V"))]
    tst.uct <- tst.all[which(tsinceltst.ct[tst.all] > sti.correlation.time &
                                 (is.na(diag.status.ct[tst.all]) | diag.status.ct[tst.all] == 0) &
                           role.class[tst.all] %in% c("I", "V"))]
    tst.syph <- tst.all[which(tsinceltst.syph[tst.all] > sti.correlation.time &
                          (is.na(diag.status.syph[tst.all]) | diag.status.syph[tst.all] == 0))]

    tst.syph.pos <- tst.syph[which(syphilis[tst.syph] == 1 &
                               stage.syph[tst.syph] %in% c(2, 3, 4, 5, 6))]
    tst.syph.neg <- setdiff(tst.syph, tst.syph.pos)
    tst.earlysyph.pos <- tst.syph[which(syphilis[tst.syph] == 1 &
                                    stage.syph[tst.syph] %in% c(2, 3, 4))]
    tst.latesyph.pos <- tst.syph[which(syphilis[tst.syph] == 1 &
                                  stage.syph[tst.syph] %in% c(5, 6))]

    tst.rgc.pos <- tst.rgc[which(rGC[tst.rgc] == 1)]
    tst.ugc.pos <- tst.ugc[which(uGC[tst.ugc] == 1)]
    tst.rgc.neg <- setdiff(tst.rgc, tst.ugc.pos)
    tst.ugc.neg <- setdiff(tst.ugc, tst.ugc.pos)
    tst.gc.pos <- c(tst.rgc.pos, tst.ugc.pos)

    tst.rct.pos <- tst.rct[which(rCT[tst.rct] == 1)]
    tst.uct.pos <- tst.uct[which(uCT[tst.uct] == 1)]
    tst.rct.neg <- setdiff(tst.rct, tst.uct.pos)
    tst.uct.neg <- setdiff(tst.uct, tst.uct.pos)
    tst.ct.pos <- c(tst.rct.pos, tst.uct.pos)

    # Syphilis Attributes
    dat$attr$last.neg.test.syph[tst.syph.neg] <- at
    dat$attr$last.neg.test.syph[tst.syph.pos] <- NA
    dat$attr$diag.status.syph[tst.syph.pos] <- 1
    dat$attr$last.diag.time.syph[tst.syph.pos] <- at
    dat$attr$time.since.last.test.syph[tst.syph] <- 0

    # GC Attributes
    dat$attr$last.neg.test.rgc[tst.rgc.neg] <- at
    dat$attr$last.neg.test.ugc[tst.ugc.neg] <- at
    dat$attr$last.neg.test.rgc[tst.rgc.pos] <- NA
    dat$attr$last.neg.test.ugc[tst.ugc.pos] <- NA
    dat$attr$diag.status.gc[tst.gc.pos] <- 1
    dat$attr$last.diag.time.gc[tst.gc.pos] <- at
    dat$attr$time.since.last.test.rgc[tst.rgc] <- 0
    dat$attr$time.since.last.test.ugc[tst.ugc] <- 0

    # CT Attributes
    dat$attr$last.neg.test.rct[tst.rct.neg] <- at
    dat$attr$last.neg.test.uct[tst.uct.neg] <- at
    dat$attr$last.neg.test.rct[tst.rct.pos] <- NA
    dat$attr$last.neg.test.uct[tst.uct.pos] <- NA
    dat$attr$diag.status.ct[tst.ct.pos] <- 1
    dat$attr$last.diag.time.ct[tst.ct.pos] <- at
    dat$attr$time.since.last.test.rct[tst.rct] <- 0
    dat$attr$time.since.last.test.uct[tst.uct] <- 0

    # Number of HIV-correlated STI tests
    dat$epi$rGC_hivdxtime[at] <- length(tst.rgc)
    dat$epi$uGC_hivdxtime[at] <- length(tst.ugc)
    dat$epi$rCT_hivdxtime[at] <- length(tst.rct)
    dat$epi$uCT_hivdxtime[at] <- length(tst.uct)
    dat$epi$syph_hivdxtime[at] <- length(tst.syph)

    dat$epi$rGC_pos_hivdxtime[at] <- length(tst.rgc.pos)
    dat$epi$uGC_pos_hivdxtime[at] <- length(tst.ugc.pos)
    dat$epi$rCT_pos_hivdxtime[at] <- length(tst.rct.pos)
    dat$epi$uCT_pos_hivdxtime[at] <- length(tst.uct.pos)
    dat$epi$syph_pos_hivdxtime[at] <- length(tst.syph.pos)
    dat$epi$syph_earlypos_hivdxtime[at] <- length(tst.earlysyph.pos)
    dat$epi$syph_latepos_hivdxtime[at] <- length(tst.latesyph.pos)

  }

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

  #Avoid increasing time for those with HIV-correlated testing (would have just tested for STIs)
  if (dat$param$sti.hivdx.correlation == "true") {

    tsinceltst.syph <- dat$attr$time.since.last.test.syph[which(is.na(dat$attr$diag.time) | dat$attr$diag.time != at)] + 1
    tsinceltst.rgc <- dat$attr$time.since.last.test.rgc[which(is.na(dat$attr$diag.time) | dat$attr$diag.time != at)] + 1
    tsinceltst.ugc <- dat$attr$time.since.last.test.ugc[which(is.na(dat$attr$diag.time) | dat$attr$diag.time != at)] + 1
    tsinceltst.rct <- dat$attr$time.since.last.test.rct[which(is.na(dat$attr$diag.time) | dat$attr$diag.time != at)] + 1
    tsinceltst.uct <- dat$attr$time.since.last.test.uct[which(is.na(dat$attr$diag.time) | dat$attr$diag.time != at)] + 1
    tsinceltst.gc <- pmin(tsinceltst.rgc, tsinceltst.ugc)
    tsinceltst.ct <- pmin(tsinceltst.rct, tsinceltst.uct)

  }

  if (dat$param$sti.hivdx.correlation == "false") {

    tsinceltst.syph <- dat$attr$time.since.last.test.syph + 1
    tsinceltst.rgc <- dat$attr$time.since.last.test.rgc + 1
    tsinceltst.ugc <- dat$attr$time.since.last.test.ugc + 1
    tsinceltst.rct <- dat$attr$time.since.last.test.rct + 1
    tsinceltst.uct <- dat$attr$time.since.last.test.uct + 1
    tsinceltst.gc <- pmin(tsinceltst.rgc, tsinceltst.ugc)
    tsinceltst.ct <- pmin(tsinceltst.rct, tsinceltst.uct)
  }

  prepStat <- dat$attr$prepStat
  stitestind1 <- dat$attr$stitest.ind.active
  stitestind2 <- dat$attr$stitest.ind.recentpartners

  # Parameters
  stianntest.gc.hivneg.coverage <- dat$param$stianntest.gc.hivneg.coverage
  stianntest.ct.hivneg.coverage <- dat$param$stianntest.ct.hivneg.coverage
  stianntest.syph.hivneg.coverage <- dat$param$stianntest.syph.hivneg.coverage
  stihighrisktest.gc.hivneg.coverage <- dat$param$stihighrisktest.gc.hivneg.coverage
  stihighrisktest.ct.hivneg.coverage <- dat$param$stihighrisktest.ct.hivneg.coverage
  stihighrisktest.syph.hivneg.coverage <- dat$param$stihighrisktest.syph.hivneg.coverage
  stianntest.gc.hivpos.coverage <- dat$param$stianntest.gc.hivpos.coverage
  stianntest.ct.hivpos.coverage <- dat$param$stianntest.ct.hivpos.coverage
  stianntest.syph.hivpos.coverage <- dat$param$stianntest.syph.hivpos.coverage
  stihighrisktest.gc.hivpos.coverage <- dat$param$stihighrisktest.gc.hivpos.coverage
  stihighrisktest.ct.hivpos.coverage <- dat$param$stihighrisktest.ct.hivpos.coverage
  stihighrisktest.syph.hivpos.coverage <- dat$param$stihighrisktest.syph.hivpos.coverage

  stianntest.cov.rate <- dat$param$stianntest.cov.rate
  stihighrisktest.cov.rate <- dat$param$stihighrisktest.cov.rate
  testing.pattern.sti <- dat$param$testing.pattern.sti
  stitest.active.int <- dat$param$stitest.active.int
  sti.highrisktest.int <- dat$param$sti.highrisktest.int
  #tst.rect.sti.rr <- dat$param$tst.rect.sti.rr

  # Eligibility and trajectory
  # Annual indications- sexually active in last year
  idsactive.hivpos <- which(stitestind1 == 1 & diag.status == 1 & tt.traj %in% c(3, 4))
  idsactive.hivneg <- setdiff((which(stitestind1 == 1)), idsactive.hivpos)

  # STI testing higher-risk eligibility scenarios
  idshighrisk.hivpos <- which(stitestind2 == 1 & diag.status == 1 & tt.traj %in% c(3, 4))
  idshighrisk.hivneg <- setdiff((which(stitestind2 == 1)), idshighrisk.hivpos)

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
  idsnotactiveelig.hivpos <- which(tt.traj.syph.hivpos == 1 & stitestind1 != 1)
  tt.traj.syph.hivpos[idsnotactiveelig.hivpos] <-
    tt.traj.gc.hivpos[idsnotactiveelig.hivpos] <-
    tt.traj.ct.hivpos[idsnotactiveelig.hivpos] <-
    NA
  idsnotactiveelig.hivneg <- which(tt.traj.syph.hivneg == 1 & stitestind1 != 1)
  tt.traj.syph.hivneg[idsnotactiveelig.hivneg] <-
    tt.traj.gc.hivneg[idsnotactiveelig.hivneg] <-
    tt.traj.ct.hivneg[idsnotactiveelig.hivneg] <-
    NA

  ## Initiation (non-HIV diagnosed) --------------------------------------------
  ### Testing coverage for high risk
  stihighrisktestCov.ct <- sum(tt.traj.ct.hivneg == 2, na.rm = TRUE) / length(idshighrisk.hivneg)
  stihighrisktestCov.ct <- ifelse(is.nan(stihighrisktestCov.ct), 0, stihighrisktestCov.ct)
  stihighrisktestCov.gc <- sum(tt.traj.gc.hivneg == 2, na.rm = TRUE) / length(idshighrisk.hivneg)
  stihighrisktestCov.gc <- ifelse(is.nan(stihighrisktestCov.gc), 0, stihighrisktestCov.gc)
  stihighrisktestCov.syph <- sum(tt.traj.syph.hivneg == 2, na.rm = TRUE) / length(idshighrisk.hivneg)
  stihighrisktestCov.syph <- ifelse(is.nan(stihighrisktestCov.syph), 0, stihighrisktestCov.syph)

  idsEligSt <- idshighrisk.hivneg
  nEligSt <- length(idshighrisk.hivneg)

  nStart.ct <- max(0, min(nEligSt, round((stihighrisktest.gc.hivneg.coverage - stihighrisktestCov.ct) *
                                           length(idshighrisk.hivneg))))
  nStart.gc <- max(0, min(nEligSt, round((stihighrisktest.ct.hivneg.coverage - stihighrisktestCov.gc) *
                                           length(idshighrisk.hivneg))))
  nStart.syph <- max(0, min(nEligSt, round((stihighrisktest.syph.hivneg.coverage - stihighrisktestCov.syph) *
                                             length(idshighrisk.hivneg))))
  idsStart.ct <- idsStart.gc <- idsStart.syph <- NULL
  if (nStart.ct > 0) {
    if (stihighrisktest.cov.rate >= 1) {
      idsStart.ct <- ssample(idsEligSt, nStart.ct)
    } else {
      idsStart.ct <- idsEligSt[rbinom(nStart.ct, 1, stihighrisktest.cov.rate) == 1]
    }
  }
  if (nStart.gc > 0) {
    if (stihighrisktest.cov.rate >= 1) {
      idsStart.gc <- ssample(idsEligSt, nStart.gc)
    } else {
      idsStart.gc <- idsEligSt[rbinom(nStart.gc, 1, stihighrisktest.cov.rate) == 1]
    }
  }
  if (nStart.syph > 0) {
    if (stihighrisktest.cov.rate >= 1) {
      idsStart.syph <- ssample(idsEligSt, nStart.syph)
    } else {
      idsStart.syph <- idsEligSt[rbinom(nStart.syph, 1, stihighrisktest.cov.rate) == 1]
    }
  }

  ## Update testing trajectory for higher-risk
  if (length(idsStart.ct) > 0) {
    tt.traj.ct.hivneg[idsStart.ct] <- 2
  }
  if (length(idsStart.gc) > 0) {
    tt.traj.gc.hivneg[idsStart.gc] <- 2
  }
  if (length(idsStart.syph) > 0) {
    tt.traj.syph.hivneg[idsStart.syph] <- 2
  }

  ### Testing coverage for annual - all those sexually active without high-risk indications
  stianntestCov.ct <- sum(tt.traj.ct.hivneg == 1, na.rm = TRUE) / length(setdiff(idsactive.hivneg, which(tt.traj.ct.hivneg == 2)))
  stianntestCov.ct <- ifelse(is.nan(stianntestCov.ct), 0, stianntestCov.ct)
  stianntestCov.gc <- sum(tt.traj.gc.hivneg == 1, na.rm = TRUE) / length(setdiff(idsactive.hivneg, which(tt.traj.gc.hivneg == 2)))
  stianntestCov.gc <- ifelse(is.nan(stianntestCov.gc), 0, stianntestCov.gc)
  stianntestCov.syph <- sum(tt.traj.syph.hivneg == 1, na.rm = TRUE) / length(setdiff(idsactive.hivneg, which(tt.traj.syph.hivneg == 2)))
  stianntestCov.syph <- ifelse(is.nan(stianntestCov.syph), 0, stianntestCov.syph)

  idsEligSt.ct <- setdiff(idsactive.hivneg, which(tt.traj.ct.hivneg == 2))
  idsEligSt.gc <- setdiff(idsactive.hivneg, which(tt.traj.ct.hivneg == 2))
  idsEligSt.syph <- setdiff(idsactive.hivneg, which(tt.traj.ct.hivneg == 2))
  nEligSt.ct <- length(setdiff(idsactive.hivneg, which(tt.traj.ct.hivneg == 2)))
  nEligSt.gc <- length(setdiff(idsactive.hivneg, which(tt.traj.ct.hivneg == 2)))
  nEligSt.syph <- length(setdiff(idsactive.hivneg, which(tt.traj.ct.hivneg == 2)))

  nStart.ct <- max(0, min(nEligSt.ct, round((stianntest.ct.hivneg.coverage - stianntestCov.ct) *
                                           length(setdiff(idsactive.hivneg, which(tt.traj.ct.hivneg == 2))))))
  nStart.gc <- max(0, min(nEligSt.gc, round((stianntest.gc.hivneg.coverage - stianntestCov.gc) *
                                           length(setdiff(idsactive.hivneg, which(tt.traj.gc.hivneg == 2))))))
  nStart.syph <- max(0, min(nEligSt.syph, round((stianntest.syph.hivneg.coverage - stianntestCov.syph) *
                                             length(setdiff(idsactive.hivneg, which(tt.traj.syph.hivneg == 2))))))
  idsStart.ct <- idsStart.gc <- idsStart.syph <- NULL

  if (nStart.ct > 0) {
    if (stianntest.cov.rate >= 1) {
      idsStart.ct <- ssample(idsEligSt.ct, nStart.ct)
    } else {
      idsStart.ct <- idsEligSt.ct[rbinom(nStart.ct, 1, stianntest.cov.rate) == 1]
    }
  }
  if (nStart.gc > 0) {
    if (stianntest.cov.rate >= 1) {
      idsStart.gc <- ssample(idsEligSt.gc, nStart.gc)
    } else {
      idsStart.gc <- idsEligSt.gc[rbinom(nStart.gc, 1, stianntest.cov.rate) == 1]
    }
  }
  if (nStart.syph > 0) {
    if (stianntest.cov.rate >= 1) {
      idsStart.syph <- ssample(idsEligSt.syph, nStart.syph)
    } else {
      idsStart.syph <- idsEligSt.syph[rbinom(nStart.syph, 1, stianntest.cov.rate) == 1]
    }
  }

  ## Update testing trajectory for lower-risk
  if (length(idsStart.ct) > 0) {
    tt.traj.ct.hivneg[idsStart.ct] <- 1
  }
  if (length(idsStart.gc) > 0) {
    tt.traj.gc.hivneg[idsStart.gc] <- 1
  }
  if (length(idsStart.syph) > 0) {
    tt.traj.syph.hivneg[idsStart.syph] <- 1
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
  tst.syph.pos.hivneg <- tst.syph.nprep.hivneg[which(syphilis[tst.syph.nprep.hivneg] == 1 &
                                   stage.syph[tst.syph.nprep.hivneg] %in% c(2, 3, 4, 5, 6))]
  tst.syph.neg.hivneg <- setdiff(tst.syph.nprep.hivneg, tst.syph.pos.hivneg)
  tst.earlysyph.pos.hivneg <- tst.syph.nprep.hivneg[which(syphilis[tst.syph.nprep.hivneg] == 1 &
                                   stage.syph[tst.syph.nprep.hivneg] %in% c(2, 3, 4))]
  tst.latesyph.pos.hivneg <- tst.syph.nprep.hivneg[which(syphilis[tst.syph.nprep.hivneg] == 1 &
                                   stage.syph[tst.syph.nprep.hivneg] %in% c(5, 6))]

  # GC non-PrEP testing
  tst.rgc.hivneg <- tst.gc.nprep.hivneg[which(role.class[tst.gc.nprep.hivneg] %in% c("R", "V"))]
  #tst.rgc.hivneg <- sample(tst.rgc.hivneg, tst.rect.sti.rr * length(tst.rgc.hivneg))
  tst.ugc.hivneg <- tst.gc.nprep.hivneg[which(role.class[tst.gc.nprep.hivneg] %in% c("I", "V"))]
  tst.rgc.pos.hivneg <- tst.rgc.hivneg[which(rGC[tst.rgc.hivneg] == 1)]
  tst.ugc.pos.hivneg <- tst.ugc.hivneg[which(uGC[tst.ugc.hivneg] == 1)]
  tst.rgc.neg.hivneg <- setdiff(tst.rgc.hivneg, tst.rgc.pos.hivneg)
  tst.ugc.neg.hivneg <- setdiff(tst.ugc.hivneg, tst.ugc.pos.hivneg)
  tst.gc.pos.hivneg <- unique(c(tst.rgc.pos.hivneg, tst.ugc.pos.hivneg))

  # CT non-PrEP testing
  tst.rct.hivneg <- tst.ct.nprep.hivneg[which(role.class[tst.ct.nprep.hivneg] %in% c("R", "V"))]
  #tst.rct.hivneg <- sample(tst.rct.hivneg, tst.rect.sti.rr * length(tst.rct.hivneg))
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

  ## Initiation (HIV diagnosed) --------------------------------------------
  ### Testing coverage for high risk
  stihighrisktestCov.ct <- sum(tt.traj.ct.hivpos == 2, na.rm = TRUE) / length(idshighrisk.hivpos)
  stihighrisktestCov.ct <- ifelse(is.nan(stihighrisktestCov.ct), 0, stihighrisktestCov.ct)
  stihighrisktestCov.gc <- sum(tt.traj.gc.hivpos == 2, na.rm = TRUE) / length(idshighrisk.hivpos)
  stihighrisktestCov.gc <- ifelse(is.nan(stihighrisktestCov.gc), 0, stihighrisktestCov.gc)
  stihighrisktestCov.syph <- sum(tt.traj.syph.hivpos == 2, na.rm = TRUE) / length(idshighrisk.hivpos)
  stihighrisktestCov.syph <- ifelse(is.nan(stihighrisktestCov.syph), 0, stihighrisktestCov.syph)

  idsEligSt <- idshighrisk.hivpos
  nEligSt <- length(idshighrisk.hivpos)

  nStart.ct <- max(0, min(nEligSt, round((stihighrisktest.ct.hivpos.coverage - stihighrisktestCov.ct) *
                                        length(idshighrisk.hivpos))))
  nStart.gc <- max(0, min(nEligSt, round((stihighrisktest.gc.hivpos.coverage - stihighrisktestCov.gc) *
                                           length(idshighrisk.hivpos))))
  nStart.syph <- max(0, min(nEligSt, round((stihighrisktest.syph.hivpos.coverage - stihighrisktestCov.syph) *
                                           length(idshighrisk.hivpos))))
  idsStart.ct <- idsStart.gc <- idsStart.syph <- NULL
  if (nStart.ct > 0) {
    if (stihighrisktest.cov.rate >= 1) {
      idsStart.ct <- ssample(idsEligSt, nStart.ct)
    } else {
      idsStart.ct <- idsEligSt[rbinom(nStart.ct, 1, stihighrisktest.cov.rate) == 1]
    }
  }
  if (nStart.gc > 0) {
    if (stihighrisktest.cov.rate >= 1) {
      idsStart.gc <- ssample(idsEligSt, nStart.gc)
    } else {
      idsStart.gc <- idsEligSt[rbinom(nStart.gc, 1, stihighrisktest.cov.rate) == 1]
    }
  }
  if (nStart.syph > 0) {
    if (stihighrisktest.cov.rate >= 1) {
      idsStart.syph <- ssample(idsEligSt, nStart.syph)
    } else {
      idsStart.syph <- idsEligSt[rbinom(nStart.syph, 1, stihighrisktest.cov.rate) == 1]
    }
  }

  ## Update testing trajectory for higher-risk
  if (length(idsStart.ct) > 0) {
    tt.traj.ct.hivpos[idsStart.ct] <- 2
  }
  if (length(idsStart.gc) > 0) {
    tt.traj.gc.hivpos[idsStart.gc] <- 2
  }
  if (length(idsStart.syph) > 0) {
    tt.traj.syph.hivpos[idsStart.syph] <- 2
  }

  ### Testing coverage for annual - all those sexually active without high-risk indications
  stianntestCov.ct <- sum(tt.traj.ct.hivpos == 1, na.rm = TRUE) / length(setdiff(idsactive.hivpos, which(tt.traj.ct.hivpos == 2)))
  stianntestCov.ct <- ifelse(is.nan(stianntestCov.ct), 0, stianntestCov.ct)
  stianntestCov.gc <- sum(tt.traj.gc.hivpos == 1, na.rm = TRUE) / length(setdiff(idsactive.hivpos, which(tt.traj.gc.hivpos == 2)))
  stianntestCov.gc <- ifelse(is.nan(stianntestCov.gc), 0, stianntestCov.gc)
  stianntestCov.syph <- sum(tt.traj.syph.hivpos == 1, na.rm = TRUE) / length(setdiff(idsactive.hivpos, which(tt.traj.syph.hivpos == 2)))
  stianntestCov.syph <- ifelse(is.nan(stianntestCov.syph), 0, stianntestCov.syph)

  idsEligSt.ct <- setdiff(idsactive.hivpos, which(tt.traj.ct.hivpos == 2))
  idsEligSt.gc <- setdiff(idsactive.hivpos, which(tt.traj.ct.hivpos == 2))
  idsEligSt.syph <- setdiff(idsactive.hivpos, which(tt.traj.ct.hivpos == 2))
  nEligSt.ct <- length(setdiff(idsactive.hivpos, which(tt.traj.ct.hivpos == 2)))
  nEligSt.gc <- length(setdiff(idsactive.hivpos, which(tt.traj.ct.hivpos == 2)))
  nEligSt.syph <- length(setdiff(idsactive.hivpos, which(tt.traj.ct.hivpos == 2)))

  nStart.ct <- max(0, min(nEligSt.ct, round((stianntest.ct.hivpos.coverage - stianntestCov.ct) *
                                              length(setdiff(idsactive.hivpos, which(tt.traj.ct.hivpos == 2))))))
  nStart.gc <- max(0, min(nEligSt.gc, round((stianntest.gc.hivpos.coverage - stianntestCov.gc) *
                                              length(setdiff(idsactive.hivpos, which(tt.traj.gc.hivpos == 2))))))
  nStart.syph <- max(0, min(nEligSt.syph, round((stianntest.syph.hivpos.coverage - stianntestCov.syph) *
                                                  length(setdiff(idsactive.hivpos, which(tt.traj.syph.hivpos == 2))))))
  idsStart.ct <- idsStart.gc <- idsStart.syph <- NULL

  if (nStart.ct > 0) {
    if (stianntest.cov.rate >= 1) {
      idsStart.ct <- ssample(idsEligSt.ct, nStart.ct)
    } else {
      idsStart.ct <- idsEligSt.ct[rbinom(nStart.ct, 1, stianntest.cov.rate) == 1]
    }
  }
  if (nStart.gc > 0) {
    if (stianntest.cov.rate >= 1) {
      idsStart.gc <- ssample(idsEligSt.gc, nStart.gc)
    } else {
      idsStart.gc <- idsEligSt.gc[rbinom(nStart.gc, 1, stianntest.cov.rate) == 1]
    }
  }
  if (nStart.syph > 0) {
    if (stianntest.cov.rate >= 1) {
      idsStart.syph <- ssample(idsEligSt.syph, nStart.syph)
    } else {
      idsStart.syph <- idsEligSt.syph[rbinom(nStart.syph, 1, stianntest.cov.rate) == 1]
    }
  }

  ## Update testing trajectory for lower-risk
  if (length(idsStart.ct) > 0) {
    tt.traj.ct.hivpos[idsStart.ct] <- 1
  }
  if (length(idsStart.gc) > 0) {
    tt.traj.gc.hivpos[idsStart.gc] <- 1
  }
  if (length(idsStart.syph) > 0) {
    tt.traj.syph.hivpos[idsStart.syph] <- 1
  }

  ## Process for asymptomatic syphilis screening
  if (testing.pattern.sti == "memoryless") {
    elig.syph.ann <- which((tt.traj.syph.hivpos == 1 &
                             (diag.status.syph == 0 | is.na(diag.status.syph)) &
                             prepStat == 0))
    rates.syph <- rep(1/stitest.active.int, length(elig.syph.ann))
    tst.syph.nprep.ann <- elig.syph.ann[rbinom(length(elig.syph.ann), 1, rates.syph) == 1]

    elig.syph.highrisk <- which((tt.traj.syph.hivpos == 2 &
                                  (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                  prepStat == 0))
    rates.syph <- rep(1/sti.highrisktest.int, length(elig.syph.highrisk))
    tst.syph.nprep.highrisk <- elig.syph.highrisk[rbinom(length(elig.syph.highrisk), 1, rates.syph) == 1]
    tst.syph.nprep.hivpos <- c(tst.syph.nprep.ann, tst.syph.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.syph.annual.interval <- which((tt.traj.syph.hivpos == 1 &
                                        (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                        tsinceltst.syph >= (stitest.active.int) &
                                        prepStat == 0))
    tst.syph.highrisk.interval <- which((tt.traj.syph.hivpos == 2 &
                                          (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                          tsinceltst.syph >= (sti.highrisktest.int) &
                                          prepStat == 0))
    tst.syph.nprep.hivpos <- c(tst.syph.annual.interval, tst.syph.highrisk.interval)
  }

  ## Process for GC asymptomatic screening
  if (testing.pattern.sti == "memoryless") {
    elig.gc.ann <- which((tt.traj.gc.hivpos == 1 &
                           (diag.status.gc == 0 | is.na(diag.status.gc)) &
                           prepStat == 0))
    rates.gc <- rep(1/stitest.active.int, length(elig.gc.ann))
    tst.gc.nprep.ann <- elig.gc.ann[rbinom(length(elig.gc.ann), 1, rates.gc) == 1]

    elig.gc.highrisk <- which((tt.traj.gc.hivpos == 2 &
                                (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                prepStat == 0))
    rates.gc <- rep(1/sti.highrisktest.int, length(elig.gc.highrisk))
    tst.gc.nprep.highrisk <- elig.gc.highrisk[rbinom(length(elig.gc.highrisk), 1, rates.gc) == 1]
    tst.gc.nprep.hivpos <- c(tst.gc.nprep.ann, tst.gc.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.gc.annual.interval <- which((tt.traj.gc.hivpos == 1 &
                                      (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                      tsinceltst.gc >= (stitest.active.int) &
                                      prepStat == 0))
    tst.gc.highrisk.interval <- which((tt.traj.gc.hivpos == 2 &
                                        (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                        tsinceltst.gc >= (sti.highrisktest.int) &
                                        prepStat == 0))
    tst.gc.nprep.hivpos <- c(tst.gc.annual.interval, tst.gc.highrisk.interval)
  }

  ## Process for CT asymptomatic screening
  if (testing.pattern.sti == "memoryless") {
    elig.ct.ann <- which((tt.traj.ct.hivpos == 1 &
                           (diag.status.ct == 0 | is.na(diag.status.ct)) &
                           prepStat == 0))
    rates.ct <- rep(1/stitest.active.int, length(elig.ct.ann))
    tst.ct.nprep.ann <- elig.ct.ann[rbinom(length(elig.ct.ann), 1, rates.ct) == 1]

    elig.ct.highrisk <- which((tt.traj.ct.hivpos == 2 &
                                (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                prepStat == 0))
    rates.ct <- rep(1/sti.highrisktest.int, length(elig.ct.highrisk))
    tst.ct.nprep.highrisk <- elig.ct.highrisk[rbinom(length(elig.ct.highrisk), 1, rates.ct) == 1]

    tst.ct.nprep.hivpos <- c(tst.ct.nprep.ann, tst.ct.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.ct.annual.interval <- which((tt.traj.ct.hivpos == 1 &
                                      (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                      tsinceltst.ct >= (stitest.active.int) &
                                      prepStat == 0))
    tst.ct.highrisk.interval <- which((tt.traj.ct.hivpos == 2 &
                                         (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                         tsinceltst.ct >= (stitest.active.int) &
                                         prepStat == 0))
    tst.ct.nprep.hivpos <- c(tst.ct.annual.interval, tst.ct.highrisk.interval)
  }

  # Syphilis non-PrEP testing
  tst.syph.pos.hivpos <- tst.syph.nprep.hivpos[which(syphilis[tst.syph.nprep.hivpos] == 1 &
                                                 stage.syph[tst.syph.nprep.hivpos] %in% c(2, 3, 4, 5, 6))]
  tst.syph.neg.hivpos <- setdiff(tst.syph.nprep.hivpos, tst.syph.pos.hivpos)
  tst.earlysyph.pos.hivpos <- tst.syph.nprep.hivpos[which(syphilis[tst.syph.nprep.hivpos] == 1 &
                                                      stage.syph[tst.syph.nprep.hivpos] %in% c(2, 3, 4))]
  tst.latesyph.pos.hivpos <- tst.syph.nprep.hivpos[which(syphilis[tst.syph.nprep.hivpos] == 1 &
                                                     stage.syph[tst.syph.nprep.hivpos] %in% c(5, 6))]

  # GC non-PrEP testing
  tst.rgc.hivpos <- tst.gc.nprep.hivpos[which(role.class[tst.gc.nprep.hivpos] %in% c("R", "V"))]
  #tst.rgc.hivpos <- sample(tst.rgc.hivpos, tst.rect.sti.rr * length(tst.rgc.hivpos))
  tst.ugc.hivpos <- tst.gc.nprep.hivpos[which(role.class[tst.gc.nprep.hivpos] %in% c("I", "V"))]
  tst.rgc.pos.hivpos <- tst.rgc.hivpos[which(rGC[tst.rgc.hivpos] == 1)]
  tst.ugc.pos.hivpos <- tst.ugc.hivpos[which(uGC[tst.ugc.hivpos] == 1)]
  tst.rgc.neg.hivpos <- setdiff(tst.rgc.hivpos, tst.rgc.pos.hivpos)
  tst.ugc.neg.hivpos <- setdiff(tst.ugc.hivpos, tst.ugc.pos.hivpos)
  tst.gc.pos.hivpos <- unique(c(tst.rgc.pos.hivpos, tst.ugc.pos.hivpos))

  # CT non-PrEP testing
  tst.rct.hivpos <- tst.ct.nprep.hivpos[which(role.class[tst.ct.nprep.hivpos] %in% c("R", "V"))]
  #tst.rct.hivpos <- sample(tst.rct.hivpos, tst.rect.sti.rr * length(tst.rct.hivpos))
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

  # Number of people on each testing trajectory
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

  dat$epi$tt.traj.gc1[at] <- length(which(tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1))
  dat$epi$tt.traj.ct1[at] <- length(which(tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1))
  dat$epi$tt.traj.syph1[at] <- length(which(tt.traj.syph.hivneg == 1 | tt.traj.syph.hivpos == 1))
  dat$epi$tt.traj.gc2[at] <- length(which(tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2))
  dat$epi$tt.traj.ct2[at] <- length(which(tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2))
  dat$epi$tt.traj.syph2[at] <- length(which(tt.traj.syph.hivneg == 2 | tt.traj.syph.hivpos == 2))

  # Number of STIs tests for asymptomatic background (non-HIV diagnosed)
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

  dat$epi$stiasympttests.hivneg[at] <- length(tst.rgc.hivneg) + length(tst.ugc.hivneg) +
    length(tst.rct.hivneg) + length(tst.uct.hivneg) + length(tst.syph.nprep.hivneg)
  dat$epi$stiasympttests.pos.hivneg[at] <- length(tst.rgc.pos.hivneg) + length(tst.ugc.pos.hivneg) +
    length(tst.rct.pos.hivneg) + length(tst.uct.pos.hivneg) + length(tst.syph.pos.hivneg)

  # Number of STI tests for asymptomatic background (HIV-diagnosed)
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

  dat$epi$stiasympttests.hivpos[at] <- length(tst.rgc.hivpos) + length(tst.ugc.hivpos) +
    length(tst.rct.hivpos) + length(tst.uct.hivpos) + length(tst.syph.nprep.hivpos)
  dat$epi$stiasympttests.pos.hivpos[at] <- length(tst.rgc.pos.hivpos) + length(tst.ugc.pos.hivpos) +
    length(tst.rct.pos.hivpos) + length(tst.uct.pos.hivpos) + length(tst.syph.pos.hivpos)

  # Number of tests for asymptomatic (total - include those from HIVdx and sympt STI dx)
  dat$epi$rGCasympttests[at] <- length(tst.rgc.hivpos) +
    length(tst.rgc.hivneg) + dat$epi$rGC_hivdxtime[at]
  dat$epi$uGCasympttests[at] <- length(tst.ugc.hivpos) +
    length(tst.ugc.hivneg) + dat$epi$uGC_hivdxtime[at]
  dat$epi$GCasympttests[at] <- length(tst.rgc.hivpos) + length(tst.ugc.hivpos) +
    length(tst.rgc.hivneg) + length(tst.ugc.hivneg) + dat$epi$rGC_hivdxtime[at] +
    dat$epi$uGC_hivdxtime[at]

  dat$epi$rGCasympttests.pos[at] <- length(tst.rgc.pos.hivpos) +
    length(tst.rgc.pos.hivneg) + dat$epi$rGC_pos_hivdxtime[at]
  dat$epi$uGCasympttests.pos[at] <- length(tst.ugc.pos.hivpos) +
    length(tst.ugc.pos.hivneg) + dat$epi$uGC_pos_hivdxtime[at]
  dat$epi$GCasympttests.pos[at] <- length(tst.rgc.pos.hivpos) + length(tst.ugc.pos.hivpos) +
    length(tst.rgc.pos.hivneg) + length(tst.ugc.pos.hivneg) + dat$epi$rGC_pos_hivdxtime[at] +
    dat$epi$uGC_pos_hivdxtime[at]

  dat$epi$rCTasympttests[at] <- length(tst.rct.hivpos) +
    length(tst.rct.hivneg) + dat$epi$rCT_hivdxtime[at]
  dat$epi$uCTasympttests[at] <- length(tst.uct.hivpos) +
    length(tst.uct.hivneg) + dat$epi$uCT_hivdxtime[at]
  dat$epi$CTasympttests[at] <- length(tst.rct.hivpos) + length(tst.uct.hivpos) +
    length(tst.rct.hivneg) + length(tst.uct.hivneg) + dat$epi$rCT_hivdxtime[at] +
    dat$epi$uCT_hivdxtime[at]

  dat$epi$rCTasympttests.pos[at] <- length(tst.rct.pos.hivpos) +
    length(tst.rct.pos.hivneg) + dat$epi$rCT_pos_hivdxtime[at]
  dat$epi$uCTasympttests.pos[at] <- length(tst.uct.pos.hivpos) +
    length(tst.uct.pos.hivneg) + dat$epi$uCT_pos_hivdxtime[at]
  dat$epi$CTasympttests.pos[at] <- length(tst.rct.pos.hivpos) +
    length(tst.uct.pos.hivpos) + length(tst.rct.pos.hivneg) +
    length(tst.uct.pos.hivneg) + dat$epi$rCT_pos_hivdxtime[at] +
    dat$epi$uCT_pos_hivdxtime[at]

  dat$epi$syphasympttests[at] <- length(tst.syph.nprep.hivpos) +
    length(tst.syph.nprep.hivneg) + dat$epi$syph_hivdxtime[at]
  dat$epi$syphasympttests.pos[at] <- length(tst.syph.pos.hivpos) +
    length(tst.syph.pos.hivneg) + dat$epi$syph_pos_hivdxtime[at]
  dat$epi$syphearlyasympttests.pos[at] <- length(tst.earlysyph.pos.hivpos) +
    length(tst.earlysyph.pos.hivneg) + dat$epi$syph_earlypos_hivdxtime[at]
  dat$epi$syphlateasympttests.pos[at] <- length(tst.latesyph.pos.hivpos) +
    length(tst.latesyph.pos.hivneg) + dat$epi$syph_latepos_hivdxtime[at]

  dat$epi$stiasympttests[at] <- length(tst.rgc.hivpos) + length(tst.ugc.hivpos) +
    length(tst.rct.hivpos) + length(tst.uct.hivpos) + length(tst.syph.nprep.hivpos) +
    length(tst.rgc.hivneg) + length(tst.ugc.hivneg) +
    length(tst.rct.hivneg) + length(tst.uct.hivneg) + length(tst.syph.nprep.hivneg) +
    dat$epi$rGC_hivdxtime[at] + dat$epi$uGC_hivdxtime[at] + dat$epi$rCT_hivdxtime[at] +
    dat$epi$uCT_hivdxtime[at] + dat$epi$syph_hivdxtime[at]
  dat$epi$stiasympttests.pos[at] <- length(tst.rgc.pos.hivpos) + length(tst.ugc.pos.hivpos) +
    length(tst.rct.pos.hivpos) + length(tst.uct.pos.hivpos) + length(tst.syph.pos.hivpos) +
    length(tst.rgc.pos.hivneg) + length(tst.ugc.pos.hivneg) +
    length(tst.rct.pos.hivneg) + length(tst.uct.pos.hivneg) + length(tst.syph.pos.hivneg) +
    dat$epi$rGC_pos_hivdxtime[at] + dat$epi$uGC_pos_hivdxtime[at] +
    dat$epi$rCT_pos_hivdxtime[at] + dat$epi$uCT_pos_hivdxtime[at] +
    dat$epi$syph_pos_hivdxtime[at]

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

