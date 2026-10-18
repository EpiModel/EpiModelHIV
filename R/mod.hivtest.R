
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
hivtest_msm   <- function(dat, at) {

  # Attributes
  diag.status             <- get_attr(dat, "diag.status")
  race                    <- get_attr(dat, "race")
  status                  <- get_attr(dat, "status")
  inf.time                <- get_attr(dat, "inf.time")
  stage                   <- get_attr(dat, "stage")
  late.tester             <- get_attr(dat, "late.tester")
  part.ident              <- get_attr(dat, "part.ident")
  entrTime                <- get_attr(dat, "entrTime")
  prepStat                <- get_attr(dat, "prepStat")
  last.neg.test           <- get_attr(dat, "last.neg.test")
  last.pos.test           <- get_attr(dat, "last.pos.test")
  diag.time               <- get_attr(dat, "diag.time")
  diag.stage              <- get_attr(dat, "diag.stage")
  tx.status               <- get_attr(dat, "tx.status")
  cuml.time.on.tx         <- get_attr(dat, "cuml.time.on.tx")

  # Parameters
  hiv.test.rate        <- get_param(dat, "hiv.test.rate")
  part.hiv.test.rate   <- get_param(dat, "part.hiv.test.rate")
  vl.aids.int          <- get_param(dat, "vl.aids.int")
  test.window.int      <- get_param(dat, "test.window.int")
  prep.tst.int         <- get_param(dat, "prep.tst.int")

  aids.test.int        <- vl.aids.int / 2

  ## Process
  # Time since last negative test
  tsincelntst   <- at - last.neg.test
  tsincelntst[is.na(tsincelntst)]   <- at - entrTime[is.na(tsincelntst)]


#  P A R T N E R - B A S E D    H I V    T E S T I N G      (P B T)
#-------------------------------------------------------------------------------

  # Tests both 1st and 2nd generation of partners to ND/PP indexes
  tstNeg.part   <- numeric()
  tstPos.part   <- numeric()

  dat   <- set_epi(dat, "tot.part.ident", at, length(which(part.ident == at)))

  eligPart   <- which(part.ident == at & (diag.status == 0 | is.na(diag.status)))
  dat   <- set_epi(dat, "elig.for.scrn", at, length(which(part.ident == at &
                                                        (diag.status == 0
                                                         | is.na(diag.status)))))

  if (length(eligPart) > 0) {
    ## Testing: If any partners identified above, test randomly based on testing rate
    # Race of individuals identified for testing
    prob.screen   <- part.hiv.test.rate[race[eligPart]]

    # Test screening
    screened   <- unique(eligPart[runif(length(eligPart)) < prob.screen])
    dat   <- set_epi(dat, "part.scrnd.tst", at, length(screened))
    dat   <- set_attr(dat, "part.scrnd", at, posit_ids = screened)

    # Testing results
    tstPos.part   <- screened[which(status[screened] == 1)]
    dat   <- set_epi(dat, "positive.part", at, length(tstPos.part))

    tstNeg.part   <- setdiff(screened, tstPos.part)
  }



#  I N D I V I D U A L - B A S E D    H I V    T E S T I N G    (I B T)
#-------------------------------------------------------------------------------

  #1.General testing
  #-----------------
  #HIV negative/ unknown persons
    elig   <- which((diag.status == 0 | is.na(diag.status)) &
                    prepStat == 0 & late.tester == 0)
    elig   <- setdiff(elig, eligPart)
    rates   <- hiv.test.rate[race[elig]]
    idsTstGen   <- elig[runif(length(elig)) < rates]


  #HIV positive persons (retests)
   if (at < get_param(dat, "prevpos.retest.start")) {
     eligPrevPos0   <- which(status == 1 & diag.status == 1)
     eligPrevPos    <- which(diag.status == 1 & tx.status == 0 & cuml.time.on.tx >= 0)

     idsTstGen.PrevPos  <- NULL
   }

   if (at >= get_param(dat, "prevpos.retest.start")) {
     eligPrevPos0   <- which(status == 1 & diag.status == 1)
     eligPrevPos   <- which(diag.status == 1 & tx.status == 0 & cuml.time.on.tx >= 0)

     #Retest rates among PP
     diag.time.prevpos   <- diag.time[eligPrevPos]
     time.since.diag   <- (at - diag.time.prevpos) / 52
     ratesPrevPos   <- exp(1.185 + (-0.024 * (time.since.diag))) / (2 * 52)

     #retests
     idsTstGen.PrevPos   <- eligPrevPos[runif(length(eligPrevPos)) < ratesPrevPos]
   }



  # 2. Late testing (Neg, then AIDS)
  #---------------------------------
    eligNeg   <- which((diag.status == 0 | is.na(diag.status)) & prepStat == 0 & status == 0 & late.tester == 1)
    dat   <- set_epi(dat, "elig.lateAIDS.pre", at, length(eligNeg))

    eligNeg   <- setdiff(eligNeg, eligPart)
    dat   <- set_epi(dat, "elig.lateAIDS.post", at, length(eligNeg))

    ratesNeg   <- 1 / (12 * 52)
    idsTstLate   <- eligNeg[runif(length(eligNeg)) < ratesNeg]

    eligAIDS   <- which((diag.status == 0 | is.na(diag.status)) & prepStat == 0 & stage == 4 & late.tester == 1)
    ratesAIDS   <- 1 / aids.test.int
    idsTstAIDS   <- eligAIDS[runif(length(eligAIDS)) < ratesAIDS]



  #3. PrEP testing
  #---------------
    idsTstPrEP   <- which((diag.status == 0 | is.na(diag.status)) & prepStat == 1 & tsincelntst >= prep.tst.int)
    idsTstPrEP   <- setdiff(idsTstPrEP, eligPart)



  # Summarize screening and test results
  tstAll   <- c(idsTstGen, idsTstLate, idsTstAIDS, idsTstPrEP)

  tstPos   <- tstAll[status[tstAll] == 1 & inf.time[tstAll] <= at - test.window.int]
  tstNeg   <- setdiff(tstAll, tstPos)
  dat   <- set_epi(dat, "tot.pos.ibtnotPP", at, length(tstPos))

  tstPos   <- c(tstPos, tstPos.part)
  tstNeg   <- c(tstNeg, tstNeg.part)
  tstAll   <- c(tstAll, tstPos.part, tstNeg.part, idsTstGen.PrevPos)

  # Update attributes
  last.neg.test[tstNeg]   <- at
  diag.status[tstPos]   <- 1
  diag.time[tstPos]   <- at
  diag.stage[tstPos]   <- stage[tstPos]
  last.pos.test[idsTstGen.PrevPos]   <- at

  # Set Attributes
  dat   <- set_attr(dat, "last.neg.test", last.neg.test)
  dat   <- set_attr(dat, "diag.status", diag.status)
  dat   <- set_attr(dat, "diag.time", diag.time)
  dat   <- set_attr(dat, "diag.stage", diag.stage)
  dat   <- set_attr(dat, "last.pos.test", last.pos.test)

  # Epi Trackers
  dat   <- set_epi(dat, "tot.tests", at, length(tstAll))
  dat   <- set_epi(dat, "tot.tests.B", at, length(intersect(tstAll, which(race == 1))))
  dat   <- set_epi(dat, "tot.tests.H", at, length(intersect(tstAll, which(race == 2))))
  dat   <- set_epi(dat, "tot.tests.W", at, length(intersect(tstAll, which(race == 3))))
  dat   <- set_epi(dat, "tot.neg.tests", at, length(tstNeg))
  dat   <- set_epi(dat, "tot.tests.ibt", at, length(c(idsTstGen, idsTstLate, idsTstAIDS,
                                                      idsTstPrEP, idsTstGen.PrevPos)))
  dat   <- set_epi(dat, "tot.tests.ibtNegunk", at, length(idsTstGen))
  dat   <- set_epi(dat, "tot.tests.ibtLateAIDs", at, length(c(idsTstLate, idsTstAIDS)))
  dat   <- set_epi(dat, "tot.tests.ibtPrEP", at, length(idsTstPrEP))
  dat   <- set_epi(dat, "tot.tests.ibtnprepnPP", at, length(c(idsTstGen, idsTstLate, idsTstAIDS)))
  dat   <- set_epi(dat, "tot.tests.ibtnotPP", at, length(c(idsTstGen, idsTstLate, idsTstAIDS, idsTstPrEP)))
  dat   <- set_epi(dat, "tot.tests.ibtPP", at, length(idsTstGen.PrevPos))

  dat   <- set_epi(dat, "allPrevPos", at, length(eligPrevPos0))
  dat   <- set_epi(dat, "eligPP.for.retest", at, length(eligPrevPos))
  dat   <- set_epi(dat, "eligPPforRetest.rxnaive", at,
                 length(intersect(eligPrevPos, which(tx.status == 0 & cuml.time.on.tx == 0))))

  return(dat)
}
