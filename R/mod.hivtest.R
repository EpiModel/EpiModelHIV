
#' @title HIV Testing Module
#'
#' @description Module function for HIV diagnostic testing of infected persons.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This testing module supports memoryless HIV testing for stochastic and geometrically-distributed
#' waiting times to test (constant hazard).
#'
#' @return
#' This function returns the `dat` object with updated `last.neg.test`, `diag.status` and
#' `diag.time` attributes.
#'
#' @export
#'
hivtest_msm <- function(dat, at) {

  ## Inputs
  # Attributes
  diag.status   <- get_attr(dat, "diag.status")
  race          <- get_attr(dat, "race")
  status        <- get_attr(dat, "status")
  inf.time      <- get_attr(dat, "inf.time")
  stage         <- get_attr(dat, "stage")
  late.tester   <- get_attr(dat, "late.tester")
  part.ident    <- get_attr(dat, "part.ident")
  entrTime      <- get_attr(dat, "entrTime")
  prepStat      <- get_attr(dat, "prepStat")
  last.neg.test <- get_attr(dat, "last.neg.test")
  diag.time     <- get_attr(dat, "diag.time")
  diag.stage    <- get_attr(dat, "diag.stage")

  # EDP Attributes
  edp.prepStat <- get_attr(dat, "edp.prepStat")

  # Parameters
  hiv.test.rate      <- get_param(dat, "hiv.test.rate")
  part.hiv.test.rate <- get_param(dat, "part.hiv.test.rate")
  vl.aids.int        <- get_param(dat, "vl.aids.int")
  test.window.int    <- get_param(dat, "test.window.int")
  prep.tst.int       <- get_param(dat, "prep.tst.int")

  aids.test.int      <- vl.aids.int / 2

  ## Process
  # Time since last negative test
  tsincelntst <- at - last.neg.test
  tsincelntst[is.na(tsincelntst)] <- at - entrTime[is.na(tsincelntst)]

  # Identified Partner Testing
  tstNeg.part <- numeric()
  tstPos.part <- numeric()
  eligPart <- which(part.ident == at & (diag.status == 0 | is.na(diag.status)))

  if (length(eligPart) > 0) {
    ## Testing: If any partners identified above, test randomly based on testing rate
    # Race of individuals identified for testing
    prob.screen <- part.hiv.test.rate[race[eligPart]]
    # Test screening
    screened <- eligPart[runif(length(eligPart)) < prob.screen]
    dat <- set_attr(dat, "part.scrnd", at, posit_ids = screened)
    # Testing results
    tstPos.part <- eligPart[status[screened] == 1]
    tstNeg.part <- eligPart[status[screened] == 0]
  }

  # General interval testing
  elig <- which((diag.status == 0 | is.na(diag.status)) &
                  (prepStat == 0 | edp.prepStat == 0) & late.tester == 0)

  elig <- setdiff(elig, eligPart)

  # Interval testing rates by race
  rates <- hiv.test.rate[race[elig]]
  idsTstGen <- elig[runif(length(elig)) < rates]

  # Late testing (Neg, then AIDS)
  eligNeg <- which((diag.status == 0 | is.na(diag.status)) &
                     (prepStat == 0 | edp.prepStat == 0) & status == 0 & late.tester == 1)

  eligNeg <- setdiff(eligNeg, eligPart)
  ratesNeg <- 1 / (12 * 52)
  idsTstLate <- eligNeg[runif(length(eligNeg)) < ratesNeg]

  eligAIDS <- which((diag.status == 0 | is.na(diag.status)) &
                      (prepStat == 0 | edp.prepStat == 0) & stage == 4 & late.tester == 1)
  ratesAIDS <- 1 / aids.test.int
  idsTstAIDS <- eligAIDS[runif(length(eligAIDS)) < ratesAIDS]

  # PrEP testing
  idsTstPrEP <- which((diag.status == 0 | is.na(diag.status)) &
                        (prepStat == 1 | edp.prepStat == 1) &
                        tsincelntst >= prep.tst.int)
  idsTstPrEP <- setdiff(idsTstPrEP, eligPart)

  tstAll <- c(idsTstGen, idsTstLate, idsTstAIDS, idsTstPrEP)

  tstPos <- tstAll[status[tstAll] == 1 & inf.time[tstAll] <= at - test.window.int]
  tstNeg <- setdiff(tstAll, tstPos)

  tstPos <- c(tstPos, tstPos.part)
  tstNeg <- c(tstNeg, tstNeg.part)
  tstAll <- c(tstAll, tstPos.part, tstNeg.part)

  # Update attributes
  last.neg.test[tstNeg] <- at
  diag.status[tstPos] <- 1
  diag.time[tstPos] <- at
  diag.stage[tstPos] <- stage[tstPos]

  ## Output
  # Set Attributes
  dat <- set_attr(dat, "last.neg.test", last.neg.test)
  dat <- set_attr(dat, "diag.status", diag.status)
  dat <- set_attr(dat, "diag.time", diag.time)
  dat <- set_attr(dat, "diag.stage", diag.stage)

  return(dat)
}
