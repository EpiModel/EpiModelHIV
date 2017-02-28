
#' @title HIV Testing Module
#'
#' @description Module function for HIV diagnostic testing of infected persons.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This testing module supports two testing parameterizations, input via the
#' \code{testing.pattern} parameter: memoryless for stochastic and
#' geometrically-distributed waiting times to test (constant hazard); and interval
#' for deterministic tested after defined waiting time intervals.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{last.neg.test},
#' \code{diag.status} and \code{diag.time} attributes.
#'
#' @keywords module msm
#'
#' @export
#'
test_msm <- function(dat, at) {

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
  dat$epi$hivtests[at] <- length(tst.all)
  
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
#' geometrically-distributed waiting times to test (constant hazard); and interval
#' for deterministic tested after defined waiting time intervals.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{last.neg.test},
#' \code{diag.status} and \code{diag.time} attributes for each STI.
#'
#' @keywords module msm
#'
#' @export
#'
test_sti_msm <- function(dat, at) {
    
    if (at < dat$param$stitest.start) {
        return(dat)
    }
    
    
    ## Variables
    
    # Attributes
    active <- dat$attr$active
    uid <- dat$attr$uid
    diag.status.syph <- dat$attr$diag.status.syph
    diag.status.gc <- dat$attr$diag.status.gc
    diag.status.ct <- dat$attr$diag.status.ct
    syphilis <- dat$attr$syphilis
    rGC <- dat$attr$rGC
    uGC <- dat$attr$uGC
    rCT <- dat$attr$rCT
    uCT <- dat$attr$uCT
    #inf.time.syph <- dat$attr$inf.time.syph
    
    role.class <- dat$attr$role.class
    stage.syph <- dat$attr$stage.syph
     
    prepStat <- dat$attr$prepStat
    
    tt.traj.syph <- dat$attr$tt.traj.syph
    tt.traj.gc <- dat$attr$tt.traj.gc
    tt.traj.ct <- dat$attr$tt.traj.ct
    
    # eptElig <- dat$attr$eptElig
    # eptStat <- dat$attr$eptStat
    # eptEligdate <- dat$attr$eptEligdate
    # eptLastRisk <- dat$attr$eptLastRisk
    # eptStartTime <- dat$attr$eptStartTime

    # Parameters
    stianntest.coverage <- dat$param$stianntest.coverage
    stianntest.cov.rate <- dat$param$stianntest.cov.rate
    stihighrisktest.coverage <- dat$param$stihighrisktest.coverage
    stihighrisktest.cov.rate <- dat$param$stihighrisktest.cov.rate
    testing.pattern.sti <- dat$param$testing.pattern.sti
    prep.tst.int <- dat$param$prep.tst.int
    stitest.active.int <- dat$param$stitest.active.int
    sti.highrisktest.int <- dat$param$sti.highrisktest.int
    tst.rect.sti.rr <- dat$param$tst.rect.sti.rr
    
    # Eligibility and trajectory
    # Base eligibility
    idsEligTest <- which(active == 1)
    
    # Annual indications- sexually active in last year
    stitestind1 <- dat$attr$stitest.ind.active
    
    # High-risk indications
    stitestind2 <- dat$attr$stitest.ind.sti
    stitestind3 <- dat$attr$stitest.ind.recentpartners
    stitestind4 <- dat$attr$stitest.ind.newpartners
    stitestind5 <- dat$attr$stitest.ind.concurrpartner
    stitestind6 <- dat$attr$stitest.ind.partnersti
    stitestind7 <- dat$attr$stitest.ind.uai.nmain
    
    # Annual - testing trajectory update
    activeindwindow <- at - stitest.active.int
    idsactive <- intersect(which(at - stitestind1 <= activeindwindow), idsEligTest)
    
    # High-risk - testing trajectory update
    highriskindwindow <- at - sti.highrisktest.int
    idshighrisk <- which((at - stitestind2 <= highriskindwindow) | (at - stitestind3 <= highriskindwindow) |
                             (at - stitestind4 <= highriskindwindow) | (at - stitestind5 <= highriskindwindow) |
                             (at - stitestind6 <= highriskindwindow) | (at - stitestind7 <= highriskindwindow))
    
    # Separating tt traj 1 vs 2 needs work
    ## Testing coverage for high risk ----------------------------------------------------------------
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
    
    ## Update testing trajectory
    if (length(idsStart) > 0) {
        tt.traj.syph[idsStart] <- tt.traj.gc[idsStart] <- tt.traj.ct[idsStart] <- 2
        #eptStartTime[idsStart] <- at
        #eptLastRisk[idsStart] <- at
    }   
    
    
    ## Testing coverage for annual ----------------------------------------------------------------
    
    stianntestCov <- sum(tt.traj.ct == 1, na.rm = TRUE) / length(idsactive)
    stianntestCov <- ifelse(is.nan(stianntestCov), 0, stianntestCov)
    
    idsEligSt <- idsactive
    nEligSt <- length(idsactive)
    
    nStart <- max(0, min(nEligSt, round((stianntest.coverage - stianntestCov) *
                                            length(idsactive))))
    idsStart <- NULL
    if (nStart > 0) {
        if (stianntest.cov.rate >= 1) {
            idsStart <- ssample(idsEligSt, nStart)
        } else {
            idsStart <- idsEligSt[rbinom(nStart, 1, stianntest.cov.rate) == 1]
        }
    }
    
    ## Update testing trajectory
    if (length(idsStart) > 0) {
        tt.traj.syph[idsStart] <- tt.traj.gc[idsStart] <- tt.traj.ct[idsStart] <- 1
        #eptStartTime[idsStart] <- at
        #eptLastRisk[idsStart] <- at
    }
    
    
    ## Stoppage (tt.traj.gc/.ct/.syph <- NA------------------------------------------------------------------
    
    # Remove testing trajectory if no longer indicated (idsannual includes high-risk)
    idsnottestelig <- which(active == 1 & tt.traj.syph %in% c(1, 2) & (at - stitestind1 >= activeindwindow))
    dat$attr$stitestLastElig[idsnottestelig] <- at
    tt.traj.syph[idsnottestelig] <- tt.traj.gc[idsnottestelig] <- tt.traj.ct[idsnottestelig] <- NA

    ## Testing
    
    tsincelntst.syph <- at - dat$attr$last.neg.test.syph
    tsincelntst.syph[is.na(tsincelntst.syph)] <- at - dat$attr$arrival.time[is.na(tsincelntst.syph)]
    
    tsincelntst.rgc <- at - dat$attr$last.neg.test.rgc
    tsincelntst.ugc <- at - dat$attr$last.neg.test.ugc
    tsincelntst.rgc[is.na(tsincelntst.rgc)] <- at - dat$attr$arrival.time[is.na(tsincelntst.rgc)]
    tsincelntst.ugc[is.na(tsincelntst.ugc)] <- at - dat$attr$arrival.time[is.na(tsincelntst.ugc)]
    tsincelntst.gc <- min(tsincelntst.rgc, tsincelntst.ugc)
     
    tsincelntst.rct <- at - dat$attr$last.neg.test.rct
    tsincelntst.uct <- at - dat$attr$last.neg.test.uct
    tsincelntst.rct[is.na(tsincelntst.rct)] <- at - dat$attr$arrival.time[is.na(tsincelntst.rct)]
    tsincelntst.uct[is.na(tsincelntst.uct)] <- at - dat$attr$arrival.time[is.na(tsincelntst.uct)]
    tsincelntst.ct <- min(tsincelntst.rct, tsincelntst.uct)
    
    # Testing Rates by serostatus/race?
    # All MSM with HIV infection entering care should be screened for gonorrhea and 
    # chlamydia at appropriate anatomic sites of exposure, as well as for syphilis
    # For sexually active individuals, screen at first HIV evaluation, and at least annually thereafter
    
    # More frequent STD screening (i.e., for syphilis, gonorrhea, and chlamydia) 
    # at 3–6-month intervals is indicated for MSM, including those with HIV infection 
    # if risk behaviors persist or if they or their sexual partners have multiple partners.
    
    # Mostly asymptomatic testing handled here - symptomatic testing is equated to probability of symptomatic treatment

    ## Process for syphilis
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

    ## Process for GC
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
     
    ## Process for CT
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
     
    # PrEP testing
    tst.syph.prep <- which((diag.status.syph == 0 | is.na(diag.status.syph)) &
                         prepStat == 1 &
                         tsincelntst.syph >= prep.tst.int)
    tst.gc.prep <- which((diag.status.gc == 0 | is.na(diag.status.gc)) &
                           prepStat == 1 &
                           tsincelntst.gc >= prep.tst.int)
    tst.ct.prep <- which((diag.status.ct == 0 | is.na(diag.status.ct)) &
                           prepStat == 1 &
                           tsincelntst.ct >= prep.tst.int)
    
    # Syphilis testing
    tst.syph.all <- c(tst.syph.nprep, tst.syph.prep)
    tst.syph.pos <- tst.syph.all[syphilis[tst.syph.all] == 1 & stage.syph[tst.syph.all] %in% c(2, 3, 4, 5, 6, 7)]
    tst.syph.neg <- setdiff(tst.syph.all, tst.syph.pos)
    
    # GC testing
    tst.gc.all <- c(tst.gc.nprep, tst.gc.prep)
    tst.rgc <- tst.gc.all[dat$attr$role.class %in% c("R", "V")]
    tst.rgc <- sample(tst.rgc, tst.rect.sti.rr * length(tst.rgc))
    tst.ugc <- tst.gc.all[dat$attr$role.class %in% c("I", "V")]
    tst.rgc.pos <- tst.rgc[rGC == 1]
    tst.ugc.pos <- tst.ugc[uGC == 1]
    tst.rgc.neg <- setdiff(tst.rgc, tst.rgc.pos)
    tst.ugc.neg <- setdiff(tst.ugc, tst.ugc.pos)
    tst.gc.pos <- unique(c(tst.rgc.pos, tst.ugc.pos))
    tst.gc.neg <- unique(c(tst.rgc.neg, tst.ugc.neg))
    
    # CT testing
    tst.ct.all <- c(tst.ct.nprep, tst.ct.prep)
    tst.rct <- tst.ct.all[dat$attr$role.class %in% c("R", "V")]
    tst.rct <- sample(tst.rct, tst.rect.sti.rr * length(tst.rct))
    tst.uct <- tst.ct.all[dat$attr$role.class %in% c("I", "V")]
    tst.rct.pos <- tst.rct[rCT == 1]
    tst.uct.pos <- tst.uct[uCT == 1]
    tst.rct.neg <- setdiff(tst.rct, tst.rct.pos)
    tst.uct.neg <- setdiff(tst.uct, tst.uct.pos)
    tst.ct.pos <- unique(c(tst.rct.pos, tst.uct.pos))
    tst.ct.neg <- unique(c(tst.rct.neg, tst.uct.neg))
    
    # Syphilis Attributes
    dat$attr$last.neg.test.syph[tst.syph.neg] <- at
    dat$attr$last.neg.test.syph[tst.syph.pos] <- NA
    dat$attr$diag.status.syph[tst.syph.pos] <- 1
    dat$attr$lastdiag.time.syph[tst.syph.pos] <- at
    
    # GC Attributes
    dat$attr$last.neg.test.rgc[tst.rgc.neg] <- at
    dat$attr$last.neg.test.ugc[tst.ugc.neg] <- at
    dat$attr$last.neg.test.rgc[tst.rgc.pos] <- NA
    dat$attr$last.neg.test.ugc[tst.ugc.pos] <- NA
    dat$attr$diag.status.gc[tst.gc.pos] <- 1
    dat$attr$lastdiag.time.gc[tst.gc.pos] <- at
    
    # CT Attributes
    dat$attr$last.neg.test.rct[tst.rct.neg] <- at
    dat$attr$last.neg.test.uct[tst.uct.neg] <- at
    dat$attr$last.neg.test.rct[tst.rct.pos] <- NA
    dat$attr$last.neg.test.uct[tst.uct.pos] <- NA
    dat$attr$diag.status.ct[tst.ct.pos] <- 1
    dat$attr$lastdiag.time.ct[tst.ct.pos] <- at
    
    if (is.null(dat$epi$num.asympt.tx)) {
        dat$epi$rGCasympttests <- rep(0, length(dat$control$nsteps))
        dat$epi$uGCasympttests <- rep(0, length(dat$control$nsteps))
        dat$epi$GCasympttests <- rep(0, length(dat$control$nsteps))
        dat$epi$rCTasympttests <- rep(0, length(dat$control$nsteps))
        dat$epi$rCTasympttests <- rep(0, length(dat$control$nsteps))
        dat$epi$syphasympttests <- rep(0, length(dat$control$nsteps))
    }
    
    # Number of tests for asymptomatic
    dat$epi$rGCasympttests[at] <- length(tst.rgc)
    dat$epi$uGCasympttests[at] <- length(tst.ugc)
    dat$epi$GCasympttests[at] <- length(c(tst.rgc, tst.ugc))
    
    dat$epi$rCTasympttests[at] <- length(tst.rct)
    dat$epi$uCTasympttests[at] <- length(tst.uct)
    dat$epi$CTasympttests[at] <- length(c(tst.rct, tst.uct))
    
    dat$epi$syphasympttests[at] <- length(c(tst.syph.all))
    
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

