
#' @title HIV Testing Module
#'
#' @description Module function for HIV diagnostic testing of infected persons.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This testing module supports memoryless HIV testing for stochastic and
#' geometrically-distributed waiting times to test (constant hazard).
#'
#' @return
#' This function returns the \code{dat} object with updated \code{last.neg.test},
#' \code{diag.status} and \code{diag.time} attributes.
#'
#' @keywords module msm
#'
#' @export
#'
hivtest_msm <- function(dat, at) {

  ## Variables

  # Attributes
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  status <- dat$attr$status
  inf.time <- dat$attr$inf.time
  stage <- dat$attr$stage
  late.tester <- dat$attr$late.tester

  prepStat <- dat$attr$prepStat
  prep.tst.int <- dat$param$prep.tst.int

  # Parameters
  hiv.test.rate <- dat$param$hiv.test.rate
  aids.test.int <- dat$param$vl.aids.int/2
  twind.int <- dat$param$test.window.int

  tsincelntst <- at - dat$attr$last.neg.test
  tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]

  # General interval testing
  elig <- which((diag.status == 0 | is.na(diag.status)) &
                prepStat == 0 & late.tester == 0)

  # Interval testing rates by race
  rates <- hiv.test.rate[race[elig]]
  idsTstGen <- elig[rbinom(length(elig), 1, rates) == 1]

  # Late testing (Neg, then AIDS)
  eligNeg <- which((diag.status == 0 | is.na(diag.status)) &
                   prepStat == 0 & status == 0 & late.tester == 1)
  ratesNeg <- 1/(12*52)
  idsTstLate <- eligNeg[rbinom(length(eligNeg), 1, ratesNeg) == 1]

  eligAIDS <- which((diag.status == 0 | is.na(diag.status)) &
                   prepStat == 0 & stage == 4 & late.tester == 1)
  ratesAIDS <- 1/aids.test.int
  idsTstAIDS <- eligAIDS[rbinom(length(eligAIDS), 1, ratesAIDS) == 1]

  # PrEP testing
  idsTstPrEP <- which((diag.status == 0 | is.na(diag.status)) &
                      prepStat == 1 &
                      tsincelntst >= prep.tst.int)

  tstAll <- c(idsTstGen, idsTstLate, idsTstAIDS, idsTstPrEP)

  tstPos <- tstAll[status[tstAll] == 1 & inf.time[tstAll] <= at - twind.int]
  tstNeg <- setdiff(tstAll, tstPos)

  # Attributes
  dat$attr$last.neg.test[tstNeg] <- at
  dat$attr$diag.status[tstPos] <- 1
  dat$attr$diag.time[tstPos] <- at
  dat$attr$diag.stage[tstPos] <- stage[tstPos]

  # Summary stats
  if (at >= 52*65) {
    dat$attr$num.neg.tests[tstNeg] <- dat$attr$num.neg.tests[tstNeg] + 1
  }
  dat$epi$tot.tests[at] <- length(tstAll)
  dat$epi$tot.tests.B[at] <- length(intersect(tstAll, which(race == 1)))
  dat$epi$tot.tests.H[at] <- length(intersect(tstAll, which(race == 2)))
  dat$epi$tot.tests.W[at] <- length(intersect(tstAll, which(race == 3)))
  dat$epi$tot.tests.nprep[at] <- length(c(idsTstGen, idsTstLate, idsTstAIDS))

  dat$epi$tot.neg.tests[at] <- length(tstNeg)

  return(dat)
}


#' @export
#' @rdname hivtest_msm
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

