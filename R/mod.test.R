
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

  return(dat)
}


#' @title STI Testing Module
#'
#' @description Module function for STI diagnostic testing of infected persons.
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
test_sti_msm <- function(dat, at) {
    
    ## Variables
    
    # Attributes
    # diag.status.syph <- dat$attr$diag.status.syph
    # diag.status.gc <- dat$attr$diag.status.gc
    # diag.status.ct <- dat$attr$diag.status.ct
    # race <- dat$attr$race
    # syphstatus <- dat$attr$syphstatus
    # rGC <- dat$attr$rGC
    # uGC <- dat$attr$uGC
    # rCT <- dat$attr$rCT
    # uCT <- dat$attr$uCT
    # inf.time.syph <- dat$attr$inf.time.syph
    # ttntest.syph <- dat$attr$ttntest.syph
    # stage.syph <- dat$attr$stage.syph
    # diag.time.syph <- dat$attr$diag.time.syph
    # sexactive <- dat$attr$sexactive
    # 
    # prepStat <- dat$attr$prepStat
    # prep.tst.int <- dat$param$prep.tst.int
    # 
    # # Assign new births and diagnosis
    # newBirths <- which(dat$attr$arrival.time == at)
    # diag.status.syph[newBirths] <- NA
    # dat$attr$diag.time.syph[newBirths] <- NA
    # 
    # # Parameters
    # testing.pattern.sti <- dat$param$testing.pattern.sti
    # sti.annualtest.int <- dat$param$sti.annualtest.int
    # sti.6motest.int <- dat$param$sti.6motest.int
    # sti.3motest.int <- dat$param$sti.3motest.int
    # tt.traj.syph <- dat$param$tt.traj.syph
    # tt.traj.gc <- dat$attr$tt.traj.gc
    # tt.traj.ct <- dat$attr$tt.traj.ct
    # 
    # tsincelntst.syph <- at - dat$attr$last.neg.test.syph
    # tsincelntst.syph[is.na(tsincelntst.syph)] <- at - dat$attr$arrival.time[is.na(tsincelntst.syph)]
    # 
    # tsincelntst.gc <- at - dat$attr$last.neg.test.gc
    # tsincelntst.gc[is.na(tsincelntst.gc)] <- at - dat$attr$arrival.time[is.na(tsincelntst.gc)]
    # 
    # tsincelntst.ct <- at - dat$attr$last.neg.test.ct
    # tsincelntst.ct[is.na(tsincelntst.ct)] <- at - dat$attr$arrival.time[is.na(tsincelntst.ct)]

    # Debit one unit from time until next test
    # ttntest.syph <- ttntest.syph - time.unit
    
    ## Process for syphilis
    # if (testing.pattern.sti == "memoryless") {
    #     elig.syph.ann <- which(tt.traj.syph == 3 &
    #                            (diag.status.syph == 0 | is.na(diag.status.syph)) &
    #                             prepStat == 0)
    #     rates.syph <- rep(1/sti.annualtest.int, length(elig.syph.ann))
    #     tst.syph.nprep.ann <- elig.syph.ann[rbinom(length(elig.syph.ann), 1, rates.syph) == 1]
    #     
    #     
    #     elig.syph.6mo <- which(tt.traj.syph == 4 &
    #                            (diag.status.syph == 0 | is.na(diag.status.syph)) &
    #                            prepStat == 0)
    #     rates.syph <- rep(1/sti.6motest.int, length(elig.syph.6mo))
    #     tst.syph.nprep.6mo <- elig.syph.6mo[rbinom(length(elig.syph.6mo), 1, rates.syph) == 1]
    #     
    #     
    #     elig.syph.3mo <- which(tt.traj.syph == 5 &
    #                                (diag.status.syph == 0 | is.na(diag.status.syph)) &
    #                                prepStat == 0)
    #     rates.syph <- rep(1/sti.3motest.int, length(elig.syph.3mo))
    #     tst.syph.nprep.3mo <- elig.syph.3mo[rbinom(length(elig.syph.3mo), 1, rates.syph) == 1]
    #     
    #     
    #     tst.syph.nprep <- c(tst.syph.nprep.ann, tst.syph.nprep.6mo, tst.syph.nprep.3mo)
    # }
    # 
    # if (testing.pattern.sti == "interval" ) {
    #     tst.syph.annual.interval <- which(tt.traj.syph == 3 &
    #                                       (diag.status.syph == 0 | is.na(diag.status.syph)) &
    #                                       tsincelntst.syph >= 2*(sti.annualtest.int) & 
    #                                       prepStat == 0)
    # 
    # 
    #     tst.syph.6mo.interval <- which(tt.traj.syph == 4 &
    #                                    (diag.status.syph == 0 | is.na(diag.status.syph)) &
    #                                    tsincelntst.syph >= 2*(sti.6motest.int) &
    #                                    prepStat == 0)
    #     tst.syph.nprep <- c(tst.syph.annual.interval, tst.syph.6mo.interval)
    # 
    #     tst.syph.3mo.interval <- which(tt.traj.syph == 5 &
    #                                        (diag.status.syph == 0 | is.na(diag.status.syph)) &
    #                                        tsincelntst.syph >= 2*(sti.3motest.int) &
    #                                        prepStat == 0)
    #     tst.syph.nprep <- c(tst.syph.annual.interval, tst.syph.6mo.interval, tst.syph.3mo.interval)
    # }
    
        # ## Process for GC
    # if (testing.pattern.sti == "memoryless") {
    #     elig.gc.ann <- which(tt.traj.gc == 3 &
    #                                (diag.status.gc == 0 | is.na(diag.status.gc)) &
    #                                prepStat == 0)
    #     rates.gc <- rep(1/sti.annualtest.int, length(elig.gc.ann))
    #     tst.gc.nprep.ann <- elig.gc.ann[rbinom(length(elig.gc.ann), 1, rates.gc) == 1]
    #     
    #     
    #     elig.gc.6mo <- which(tt.traj.gc == 4 &
    #                                (diag.status.gc == 0 | is.na(diag.status.gc)) &
    #                                prepStat == 0)
    #     rates.gc <- rep(1/sti.6motest.int, length(elig.gc.6mo))
    #     tst.gc.nprep.6mo <- elig.gc.6mo[rbinom(length(elig.gc.6mo), 1, rates.gc) == 1]
    #
    #
    #     elig.gc.3mo <- which(tt.traj.gc == 5 &
    #                                (diag.status.gc == 0 | is.na(diag.status.gc)) &
    #                                prepStat == 0)
    #     rates.gc <- rep(1/sti.3motest.int, length(elig.gc.3mo))
    #     tst.gc.nprep.3mo <- elig.gc.3mo[rbinom(length(elig.gc.6mo), 1, rates.gc) == 1]
    #
    #
    #     tst.gc.nprep <- c(tst.gc.nprep.ann, tst.gc.nprep.6mo, tst.gc.nprep.3mo)
    # }
    # 
    # if (testing.pattern.sti == "interval" ) {
    #     tst.gc.annual.interval <- which(tt.traj.gc == 3 &
    #                                           (diag.status.gc == 0 | is.na(diag.status.gc)) &
    #                                           tsincelntst.gc >= #2*(gc.annualtest.int) & 
    #                                           prepStat == 0)
    #     
    #     
    #     tst.gc.6mo.interval <- which(tt.traj.gc == 4 &
    #                                        (diag.status.gc == 0 | is.na(diag.status.gc)) &
    #                                        tsincelntst.gc >= #2*(gc.6motest.int) &
    #                                        prepStat == 0)
    #
    #
    #     tst.gc.3mo.interval <- which(tt.traj.gc == 5 &
    #                                        (diag.status.gc == 0 | is.na(diag.status.gc)) &
    #                                        tsincelntst.gc >= #2*(gc.3motest.int) &
    #                                        prepStat == 0)
    #
    #
    #     tst.gc.nprep <- c(tst.gc.annual.interval, tst.gc.6mo.interval, tst.gc.3mo.interval)
    # }
    # 
    # ## Process for CT
    # if (testing.pattern.sti == "memoryless") {
    #     elig.ct.ann <- which(tt.traj.ct == 3 &
    #                                (diag.status.ct == 0 | is.na(diag.status.ct)) &
    #                                prepStat == 0)
    #     rates.ct <- rep(1/ct.annualtest.int, length(elig.ct.ann))
    #     tst.ct.nprep.ann <- elig.ct.ann[rbinom(length(elig.ct.ann), 1, rates.ct) == 1]
    #     
    #     
    #     elig.ct.6mo <- which(tt.traj.ct == 4 &
    #                                (diag.status.ct == 0 | is.na(diag.status.ct)) &
    #                                prepStat == 0)
    #     rates.ct <- rep(1/ct.6motest.int, length(elig.ct.6mo))
    #     tst.ct.nprep.6mo <- elig.ct.6mo[rbinom(length(elig.ct.6mo), 1, rates.ct) == 1]
    #
    #     elig.ct.6mo <- which(tt.traj.ct == 5 &
    #                                (diag.status.ct == 0 | is.na(diag.status.ct)) &
    #                                prepStat == 0)
    #     rates.ct <- rep(1/ct.3motest.int, length(elig.ct.3mo))
    #     tst.ct.nprep.3mo <- elig.ct.3mo[rbinom(length(elig.ct.3mo), 1, rates.ct) == 1]
    #     tst.ct.nprep <- c(tst.ct.nprep.ann, tst.ct.nprep.6mo, tst.ct.nprep.3mo)
    # }
    # 
    # if (testing.pattern.sti == "interval" ) {
    #     tst.ct.annual.interval <- which(tt.traj.ct == 3 &
    #                                           (diag.status.ct == 0 | is.na(diag.status.ct)) &
    #                                           tsincelntst.ct >= 2*(ct.annualtest.int) & 
    #                                           prepStat == 0)
    #     
    #     
    #     tst.ct.6mo.interval <- which(tt.traj.ct == 4 &
    #                                        (diag.status.ct == 0 | is.na(diag.status.ct)) &
    #                                        tsincelntst.ct >= 2*(ct.6motest.int) &
    #                                        prepStat == 0)
    #
    #
    #     tst.ct.3mo.interval <- which(tt.traj.ct == 5 &
    #                                        (diag.status.ct == 0 | is.na(diag.status.ct)) &
    #                                        tsincelntst.ct >= 2*(ct.3motest.int) &
    #                                        prepStat == 0)
    #
    #     tst.ct.nprep <- c(tst.ct.annual.interval, tst.ct.6mo.interval, tst.ct.3mo.interval)
    # }
    #
    # # PrEP testing
    # tst.syph.prep <- which((diag.status.syph == 0 | is.na(diag.status.syph)) &
    #                        prepStat == 1 &
    #                        tsincelntst.syph >= prep.tst.int)
    # tst.gc.prep <- which((diag.status.gc == 0 | is.na(diag.status.gc)) &
    #                            prepStat == 1 &
    #                            tsincelntst.gc >= prep.tst.int)
    # tst.ct.prep <- which((diag.status.ct == 0 | is.na(diag.status.ct)) &
    #                            prepStat == 1 &
    #                            tsincelntst.ct >= prep.tst.int)
    # 
    # # Syphilis testing
    # tst.syph.all <- c(tst.syph.nprep, tst.syph.prep)
    # tst.syph.pos <- tst.syph.all[syphstatus[tst.syph.all] == 1 & stage.syph[tst.syph.all] %in% c(2, 3, 4, 5, 6, 7)]
    # tst.syph.neg <- setdiff(tst.syph.all, tst.syph.pos)
    # 
    # # GC testing
    # tst.gc.all <- c(tst.gc.nprep, tst.gc.prep)
    # tst.gc.pos <- tst.gc.all[rGC == 1 | uGC == 1]
    # tst.gc.neg <- setdiff(tst.gc.all, tst.gc.pos)
    # 
    # # CT testing
    # tst.ct.all <- c(tst.ct.nprep, tst.ct.prep)
    # tst.ct.pos <- tst.ct.all[rCT == 1 | uCT == 1]
    # tst.ct.neg <- setdiff(tst.ct.all, tst.ct.pos)
    # 
    # 
    # # Syphilis Attributes
    # dat$attr$last.neg.test.syph[tst.syph.neg] <- at
    # dat$attr$diag.status.syph[tst.syph.pos] <- 1
    # dat$attr$diag.time.syph[tst.syph.pos] <- at
    # #dat$attr$ttntest.syph <- ttntest.syph
    # 
    # # GC Attributes
    # dat$attr$last.neg.test.gc[tst.gc.neg] <- at
    # dat$attr$diag.status.gc[tst.gc.pos] <- 1
    # dat$attr$diag.time.gc[tst.gc.pos] <- at
    # #dat$attr$ttntest.gc <- ttntest.gc
    # 
    # # CT Attributes
    # dat$attr$last.neg.test.ct[tst.ct.neg] <- at
    # dat$attr$diag.status.ct[tst.ct.pos] <- 1
    # dat$attr$diag.time.ct[tst.ct.pos] <- at
    # #dat$attr$ttntest.ct <- ttntest.ct
    
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

