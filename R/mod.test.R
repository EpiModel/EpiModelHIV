
#' @title HIV Testing Module
#'
#' @description Module function for HIV diagnostic testing of infected persons.
#'
#' @inheritParams aging_camplc
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
  active <- dat$attr$active
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  tt.traj <- dat$attr$tt.traj
  status <- dat$attr$status
  inf.time <- dat$attr$inf.time
  asmm<- dat$attr$asmm

  prepStat <- dat$attr$prepStat
  prep.tst.int <- dat$param$prep.tst.int

  # Parameters
  testing.pattern <- dat$param$testing.pattern
  mean.test.B.int <- dat$param$mean.test.B.int
  mean.test.W.int <- dat$param$mean.test.W.int
  mean.test.B.int.asmm <- dat$param$mean.test.B.int.asmm
  mean.test.W.int.asmm <- dat$param$mean.test.W.int.asmm
  twind.int <- dat$param$test.window.int

  ## Process

  if (testing.pattern == "memoryless") {
    elig.B.msm <- which(active == 1 & race == "B" & asmm == 0 & tt.traj != "NN" &
                      (diag.status == 0 | is.na(diag.status)))
    rates.B.msm <- rep(1/mean.test.B.int, length(elig.B.msm))
    rates.B.msm[which(prepStat[elig.B.msm] == 1)] <- 1/prep.tst.int
    tst.B.msm <- elig.B.msm[rbinom(length(elig.B.msm), 1, rates.B.msm) == 1]

    elig.W.msm <- which(active == 1 & race == "W" & asmm == 0 & tt.traj != "NN" &
                      (diag.status == 0 | is.na(diag.status)))
    rates.W.msm <- rep(1/mean.test.W.int, length(elig.W.msm))
    rates.W.msm[which(prepStat[elig.W.msm] == 1)] <- 1/prep.tst.int
    tst.W.msm <- elig.W.msm[rbinom(length(elig.W.msm), 1, rates.W.msm) == 1]
    
    
    elig.B.asmm <- which(active == 1 & race == "B" & asmm == 1 & tt.traj != "NN" &
                      (diag.status == 0 | is.na(diag.status)))
    rates.B.asmm <- rep(1/mean.test.B.int.asmm, length(elig.B.asmm))
    rates.B.asmm[which(prepStat[elig.B.asmm] == 1)] <- 1/prep.tst.int
    tst.B.asmm <- elig.B.asmm[rbinom(length(elig.B.asmm), 1, rates.B.asmm) == 1]
    
    elig.W.asmm <- which(active == 1 & race == "W" & asmm == 0 & tt.traj != "NN" &
                      (diag.status == 0 | is.na(diag.status)))
    rates.W.asmm <- rep(1/mean.test.W.int.asmm, length(elig.W.asmm))
    rates.W.asmm[which(prepStat[elig.W.asmm] == 1)] <- 1/prep.tst.int
    tst.W.asmm <- elig.W.asmm[rbinom(length(elig.W.asmm), 1, rates.W.asmm) == 1]
  }

  if (testing.pattern == "interval") {
    # Time since last neg test MSM
    tsincelntst <- at - dat$attr$last.neg.test
    tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]

    tst.B.nprep.msm <- which(active == 1 & race == "B" & asmm == 0 & tt.traj != "NN" &
                           (diag.status == 0 | is.na(diag.status)) &
                           tsincelntst >= (mean.test.B.int))
    tst.B.prep.msm <- which(active == 1 & race == "B" & asmm == 0 & tt.traj != "NN" &
                          (diag.status == 0 | is.na(diag.status)) &
                          prepStat == 1 & tsincelntst >= prep.tst.int)
    tst.B.msm <- c(tst.B.nprep.msm, tst.B.prep.msm)

    tst.W.nprep.msm <- which(active == 1 & race == "W" & asmm == 0 & tt.traj != "NN" &
                           (diag.status == 0 | is.na(diag.status)) &
                           tsincelntst >= (mean.test.W.int))
    tst.W.prep.msm <- which(active == 1 & race == "W" & asmm == 0 & tt.traj != "NN" &
                          (diag.status == 0 | is.na(diag.status)) &
                          prepStat == 1 & tsincelntst >= prep.tst.int)
    tst.W.msm <- c(tst.W.nprep.msm, tst.W.prep.msm)
    
    # Time since last neg test ASMM
    tsincelntst <- at - dat$attr$last.neg.test
    tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]
    
    tst.B.nprep.asmm <- which(active == 1 & race == "B" & asmm == 1 & tt.traj != "NN" &
                           (diag.status == 0 | is.na(diag.status)) &
                           tsincelntst >= (mean.test.B.int.asmm))
    tst.B.prep.asmm <- which(active == 1 & race == "B" & asmm == 1 & tt.traj != "NN" &
                          (diag.status == 0 | is.na(diag.status)) &
                          prepStat == 1 & tsincelntst >= prep.tst.int)
    tst.B.asmm <- c(tst.B.nprep.asmm, tst.B.prep.asmm)
    
    tst.W.nprep.asmm <- which(active == 1 & race == "W" & asmm == 1 & tt.traj != "NN" &
                           (diag.status == 0 | is.na(diag.status)) &
                           tsincelntst >= (mean.test.W.int.asmm))
    tst.W.prep.asmm <- which(active == 1 & race == "W" & asmm == 1 & tt.traj != "NN" &
                          (diag.status == 0 | is.na(diag.status)) &
                          prepStat == 1 & tsincelntst >= prep.tst.int)
    tst.W.asmm <- c(tst.W.nprep.asmm, tst.W.prep.asmm)
  }

  tst.pos.B.msm <- tst.B.msm[status[tst.B.msm] == 1 & inf.time[tst.B.msm] <= at - twind.int]
  tst.neg.B.msm <- setdiff(tst.B.msm, tst.pos.B.msm)

  tst.pos.W.msm <- tst.W.msm[status[tst.W.msm] == 1 & inf.time[tst.W.msm] <= at - twind.int]
  tst.neg.W.msm <- setdiff(tst.W.msm, tst.pos.W.msm)
  
  tst.pos.B.asmm <- tst.B.asmm[status[tst.B.asmm] == 1 & inf.time[tst.B.asmm] <= at - twind.int]
  tst.neg.B.asmm <- setdiff(tst.B.asmm, tst.pos.B.asmm)
  
  tst.pos.W.asmm <- tst.W.asmm[status[tst.W.asmm] == 1 & inf.time[tst.W.asmm] <= at - twind.int]
  tst.neg.W.asmm <- setdiff(tst.W.asmm, tst.pos.W.asmm)

  tst.pos <- c(tst.pos.B.msm, tst.pos.W.msm, tst.pos.B.asmm, tst.pos.W.asmm)
  tst.neg <- c(tst.neg.B.msm, tst.neg.W.msm, tst.neg.B.asmm, tst.neg.W.asmm)


  ## Output

  # Attributes
  dat$attr$last.neg.test[tst.neg] <- at
  dat$attr$diag.status[tst.pos] <- 1
  dat$attr$diag.time[tst.pos] <- at

  ## Summary statistics
  dat$epi$tst.W.inc.msm[at] <- length(tst.W.msm)
  dat$epi$tst.B.inc.msm[at] <- length(tst.B.msm)
  dat$epi$tst.W.inc.asmm[at] <- length(tst.W.asmm)
  dat$epi$tst.B.inc.asmm[at] <- length(tst.B.asmm)


  return(dat)
}

