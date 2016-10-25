
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


#' @title Syphilis Testing Module
#'
#' @description Module function for syphilis diagnostic testing of infected persons.
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
test_syph_msm <- function(dat, at) {
    
    ## Variables
    
    # Attributes
    diag.status.syph <- dat$attr$diag.status.syph
    race <- dat$attr$race
    tt.traj.syph <- dat$attr$tt.traj.syph
    syphstatus <- dat$attr$syphstatus
    inf.time.syph <- dat$attr$inf.time.syph
    ttntest.syph <- dat$attr$ttntest.syph
    stage.syph <- dat$attr$stage.syph
    
    prepStat <- dat$attr$prepStat
    prep.tst.int <- dat$param$prep.tst.int
    
    # Parameters
    testing.pattern.syph <- dat$param$testing.pattern.syph
    syph.annualtest.int <- dat$param$syph.annualtest.int
    syph.6motest.int <- dat$param$syph.6motest.int
    
    tsincelntst.syph <- at - dat$attr$last.neg.test.syph
    tsincelntst.syph[is.na(tsincelntst.syph)] <- at - dat$attr$arrival.time[is.na(tsincelntst.syph)]
    
    ## Process
    selected.ann <- which(tt.traj.syph == 3)
    selected.6mo <- which(tt.traj.syph == 4)
    
    if (testing.pattern.syph == "memoryless") {
        elig.syph.ann <- which(tt.traj.syph == 3 &
                               (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                prepStat == 0)
        rates.syph <- rep(1/syph.annualtest.int, length(elig.syph.ann))
        tst.syph.nprep.ann <- elig.syph.ann[rbinom(length(elig.syph.ann), 1, rates.syph) == 1]
        
        
        elig.syph.6mo <- which(tt.traj.syph == 4 &
                               (diag.status.syph == 0 | is.na(diag.status.syph)) &
                               prepStat == 0)
        rates.syph <- rep(1/syph.6motest.int, length(elig.syph.6mo))
        tst.syph.nprep.6mo <- elig.syph.6mo[rbinom(length(elig.syph.6mo), 1, rates.syph) == 1]
        tst.syph.nprep <- c(tst.syph.nprep.ann, tst.syph.nprep.6mo)
    }
    
    if (testing.pattern.syph == "interval" ) {
        tst.syph.annual.interval <- which(tt.traj.syph == 3 &
                                          (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                          tsincelntst.syph >= 2*(syph.annualtest.int) & 
                                          prepStat == 0)

    
        tst.syph.6mo.interval <- which(tt.traj.syph == 4 &
                                       (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                       tsincelntst.syph >= 2*(syph.6motest.int) &
                                       prepStat == 0)
        tst.syph.nprep <- c(tst.syph.annual.interval, tst.syph.6mo.interval)
    }
    
    
    # PrEP testing
    tst.syph.prep <- which((diag.status.syph == 0 | is.na(diag.status.syph)) &
                           prepStat == 1 &
                           tsincelntst.syph >= prep.tst.int)
    
    tst.all <- c(tst.syph.nprep, tst.syph.prep)
    
    tst.pos <- tst.all[syphstatus[tst.all] == 1 & stage.syph[tst.all] %in% c(2, 3, 4, 5, 6, 7)]
    tst.neg <- setdiff(tst.all, tst.pos)
    
    # Attributes
    dat$attr$last.neg.test.syph[tst.neg] <- at
    dat$attr$diag.status.syph[tst.pos] <- 1
    dat$attr$diag.time.syph[tst.pos] <- at
    
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

