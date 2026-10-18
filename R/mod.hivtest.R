
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

  sexStigmaClass<- get_attr(dat, "sexStigmaClass")
  mdd.tx.status     <- get_attr(dat, "mdd.tx.status")
  mdd.txintv.status <- get_attr(dat, "mdd.txintv.status")
  mdd.txintvinit.time<- get_attr(dat, "mdd.txintvinit.time")
  mdd.txintv.last.endtime<- get_attr(dat, "mdd.txintv.last.endtime")


  # Parameters
  hiv.test.rate      <- get_param(dat, "hiv.test.rate")
  part.hiv.test.rate <- get_param(dat, "part.hiv.test.rate")
  vl.aids.int        <- get_param(dat, "vl.aids.int")
  test.window.int    <- get_param(dat, "test.window.int")
  prep.tst.int       <- get_param(dat, "prep.tst.int")

  stigma.hivtest.mult     <- get_param(dat, "stigma.hivtest.mult")
  minadeq.txdur           <- get_param(dat, "minadeq.txdur")
  stigcounse.efxdur       <- get_param(dat, "stigcounse.efxdur")
  mh.stigtxefx.hivtst     <- get_param(dat, "mh.stigtxefx.hivtst")


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
    prob.screen <- part.hiv.test.rate[race[eligPart]]

        #add sexual stigma effects on HIV testing
        if(at >= get_param(dat, "mhefx.start")) {
          txintv.multelig        <- rep(0, length(eligPart))

          #get stigma class, intv tx (mdd) status, tottime spent in mdd therapy or
          #time since therapy ended (if the indiv participated in one recently) to
          #identify those for whom testing impacts of stigma will be mitigated. (Mitigating
          #stigma impacts on testing for those with anticipated healthcare predominant stigma
          #[stigmaClass==2] only, since those with family/gen soc stigma [stigmaClass==3]
          #already had positive HIV testing uptake impacts)
          stigmaClass            <- sexStigmaClass[eligPart]

          txintv.status          <- mdd.txintv.status[eligPart]
          txintv.tottime         <- at - mdd.txintvinit.time[eligPart]
          txintv.lastendtime     <- at - mdd.txintv.last.endtime[eligPart]

          #id those eligible for multiplier adjustment because they are on stigma-impacting mdd tx
          #here this would be individual's with anticipated healthcare predominant stigma
          #who have received minimally adequate tx or recently stopped tx after remission
          txintv.multelig[which(txintv.status == 1 & txintv.tottime >= minadeq.txdur |
                                   txintv.lastendtime <= stigcounse.efxdur)] <- 1

          #get multiplier (pulled from param) based on individual's stigmaClass
          hivtst.mult <- stigma.hivtest.mult[stigmaClass]
          #increase screening prob for those with anticipated healthcare stigma to 1 if they have
          #received minimally adequate stig interv] or
          #have recently completed mdd intvtx [within the past 6 mos]
          #set to increase odds of testing by 50%
          hivtst.mult [txintv.multelig == 1] <- mh.stigtxefx.hivtst

          #apply final multiplier to test rate
          prob.screen <- plogis(qlogis(prob.screen) + log(hivtst.mult))
        }

    screened <- unique(eligPart[runif(length(eligPart)) < prob.screen])
    dat   <- set_epi(dat, "tot.tests.pbt", at, length(screened))

    dat <- set_attr(dat, "part.scrnd", at, posit_ids = screened)

    # Testing results
    tstPos.part <- screened[which(status[screened] == 1)]
    tstNeg.part <- setdiff(screened, tstPos.part)
  }



  # General interval testing
  elig <- which((diag.status == 0 | is.na(diag.status)) &
                  prepStat == 0 & late.tester == 0)
  elig <- setdiff(elig, eligPart)

  rates <- hiv.test.rate[race[elig]]
  #add sexual stigma effects on HIV testing
  if(at >= get_param(dat, "mhefx.start")) {
    txintv.multelig        <- rep(0, length(elig))

    #get stigma class, intv tx (mdd) status, tottime spent in mdd therapy or
    #time since therapy ended (if the indiv participated in one recently) to
    #identify those for whom testing impacts of stigma will be mitigated. (Mitigating
    #stigma impacts on testing for those with anticipated healthcare predominant stigma
    #[stigmaClass==2] only, since those with family/gen soc stigma [stigmaClass==3]
    #already had positive HIV testing uptake impacts)
    stigmaClass            <- sexStigmaClass[elig]

    txintv.status          <- mdd.txintv.status[elig]
    txintv.tottime         <- at - mdd.txintvinit.time[elig]
    txintv.lastendtime     <- at - mdd.txintv.last.endtime[elig]

    #id those eligible for multiplier adjustment because they are on stigma-impacting mdd tx
    #here this would be individual's with anticipated healthcare predominant stigma
    #who have received minimally adequate tx or recently stopped tx after remission
    txintv.multelig[which(txintv.status == 1 & txintv.tottime >= minadeq.txdur |
                            txintv.lastendtime <= stigcounse.efxdur)] <- 1

    #get multiplier (pulled from param) based on individual's stigmaClass
    hivtst.mult <- stigma.hivtest.mult[stigmaClass]
    #adjust prepinit prob for those eligible for multiplier adjustment (here adjusting by
    #set to increase odds of HIV testing by 50%
    hivtst.mult [txintv.multelig == 1] <- mh.stigtxefx.hivtst

    #apply adjusted multiplier to test rate
    rates <- exp(log(rates) + log(hivtst.mult))
  }

  idsTstGen <- elig[runif(length(elig)) < rates]



  # Late testing (Neg, then AIDS)
  eligNeg <- which((diag.status == 0 | is.na(diag.status)) &
                     prepStat == 0 & status == 0 & late.tester == 1)
  eligNeg <- setdiff(eligNeg, eligPart)

  ratesNeg <- 1 / (12 * 52)
  idsTstLate <- eligNeg[runif(length(eligNeg)) < ratesNeg]


  eligAIDS <- which((diag.status == 0 | is.na(diag.status)) &
                      prepStat == 0 & stage == 4 & late.tester == 1)
  ratesAIDS <- 1 / aids.test.int
  idsTstAIDS <- eligAIDS[runif(length(eligAIDS)) < ratesAIDS]



  # PrEP testing
  idsTstPrEP <- which((diag.status == 0 | is.na(diag.status)) &
                        prepStat == 1 &
                        tsincelntst >= prep.tst.int)
  idsTstPrEP <- setdiff(idsTstPrEP, eligPart)



  #all testing
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

  # Epi Trackers
  dat   <- set_epi(dat, "tot.tests", at, length(tstAll))
  dat   <- set_epi(dat, "tot.tests.ibt", at, length(c(idsTstGen, idsTstLate, idsTstAIDS,
                                                      idsTstPrEP)))


  return(dat)
}
