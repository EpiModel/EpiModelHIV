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
  diag.status <- get_attr(dat, "diag.status")
  race <- get_attr(dat, "race")
  status <- get_attr(dat, "status")
  inf.time <- get_attr(dat, "inf.time")
  stage <- get_attr(dat, "stage")
  late.tester <- get_attr(dat, "late.tester")
  prep <- get_attr(dat, "prep")
  prep.otc <- get_attr(dat, "prep.otc")
  prep.any <- get_attr(dat, "prep.any")
  last.neg.test <- get_attr(dat, "last.neg.test")
  prep.start.last <- get_attr(dat, "prep.start.last")
  prep.otc.start.last <- get_attr(dat, "prep.otc.start.last")

  # Parameters
  hiv.test.rate <- get_param(dat, "hiv.test.rate")
  vl.aids.int <- get_param(dat, "vl.aids.int")
  test.window.int <- get_param(dat, "test.window.int")
  prep.tst.int <- get_param(dat, "prep.tst.int")

  prep.otc.tst.int <- get_param(dat, "prep.otc.tst.int")
  prep.otc.always.hiv.tst <- get_param(dat, "prep.otc.always.hiv.tst")

  aids.test.int <- vl.aids.int / 2

  ## Process
  # Time since last negative test
  tsincelntst <- at - last.neg.test

  undiag_base <- diag.status == 0

  # General interval testing
  gen_elig_ids <- which(
    undiag_base &
      prep.any == 0 &
      late.tester == 0
  )

  # Interval testing rates by race
  rates <- hiv.test.rate[race[gen_elig_ids]]
  idsTstGen <- gen_elig_ids[runif(length(gen_elig_ids)) < rates]

  # Late testing (Neg, then AIDS)
  eligNeg <- which(
    undiag_base &
      prep.any == 0 &
      status == 0 &
      late.tester == 1
  )

  ratesNeg <- 1 / (12 * 52)
  idsTstLate <- eligNeg[runif(length(eligNeg)) < ratesNeg]

  eligAIDS <- which(
    undiag_base &
      prep.any == 0 &
      stage == 4 &
      late.tester == 1
  )
  ratesAIDS <- 1 / aids.test.int
  idsTstAIDS <- eligAIDS[runif(length(eligAIDS)) < ratesAIDS]

  # PrEP testing
  idsTstPrEP <- which(
    undiag_base &
      prep == 1 &
      # PrEP just started or test interval
      (prep.start.last == at | tsincelntst >= prep.tst.int)
  )

  if (prep.otc.always.hiv.tst == 1) {
    # PrEP OTC testing - interval
    prep_otc_test_ids <- which(
      undiag_base &
        prep.otc == 1 &
        # PrEP just started or test interval
        (prep.otc.start.last == at | tsincelntst >= prep.otc.tst.int)
    )
  } else {
    # PrEP OTC testing - rate based
    rate <- 1 / prep.otc.tst.int
    elig_ids <- which(undiag_base & prep.otc == 1)
    prep_otc_test_ids <- elig_ids[runif(length(elig_ids)) < rate]
    if (prep.otc.always.hiv.tst > 0) {
      # Some OTC users test on startup
      rate <- prep.otc.always.hiv.tst
      elig_ids <- which(prep.otc.start.last == at)
      always_test_ids <- elig_ids[runif(length(elig_ids)) < rate]
      prep_otc_test_ids <- union(prep_otc_test_ids, always_test_ids)
    }
  }
  tstAll <- c(idsTstGen, idsTstLate, idsTstAIDS, idsTstPrEP, prep_otc_test_ids)

  tstPos <- tstAll[
    status[tstAll] == 1 & inf.time[tstAll] <= at - test.window.int
  ]
  tstNeg <- setdiff(tstAll, tstPos)

  ## Output
  dat <- set_pos_hiv_test_attrs(dat, tstPos)
  dat <- set_neg_hiv_test_attrs(dat, tstNeg)

  if (get_control(dat, "debug")) {
    dat <- set_epi(dat, "dbg_hiv_tst_general_elig", at, length(gen_elig_ids))
    dat <- set_epi(dat, "dbg_hiv_tst_general", at, length(idsTstGen))
    dat <- set_epi(dat, "dbg_hiv_tst_late_elig", at, length(eligNeg))
    dat <- set_epi(dat, "dbg_hiv_tst_late", at, length(idsTstLate))
    dat <- set_epi(dat, "dbg_hiv_tst_aids_elig", at, length(eligAIDS))
    dat <- set_epi(dat, "dbg_hiv_tst_aids", at, length(idsTstAIDS))
    dat <- set_epi(dat, "dbg_hiv_tst_prep_std", at, length(idsTstPrEP))
    dat <- set_epi(dat, "dbg_hiv_tst_prep_otc", at, length(prep_otc_test_ids))
  }

  return(dat)
}

set_pos_hiv_test_attrs <- function(dat, ids) {
  at <- get_current_timestep(dat)
  dat <- set_attr(dat, "diag.status", 1, posit_ids = ids)
  dat <- set_attr(dat, "diag.time", at, posit_ids = ids)
  new_diag_stage <- get_attr(dat, "stage", posit_ids = ids)
  dat <- set_attr(dat, "diag.stage", new_diag_stage, posit_ids = ids)

  # Stop PrEP if diagnosed with HIV
  prep_ids <- ids[get_attr(dat, "prep", posit_ids = ids) == 1]
  dat <- set_stop_prep_attrs(dat, prep_ids)

  # Stop PrEP OTC if diagnosed with HIV
  prep_otc_ids <- ids[get_attr(dat, "prep.otc", posit_ids = ids) == 1]
  dat <- set_stop_prep_otc_attrs(dat, prep_otc_ids)

  return(dat)
}

set_neg_hiv_test_attrs <- function(dat, ids) {
  at <- get_current_timestep(dat)
  dat <- set_attr(dat, "last.neg.test", at, posit_ids = ids)
  return(dat)
}
