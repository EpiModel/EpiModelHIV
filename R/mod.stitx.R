#' @title STI Treatment Module
#'
#' @description Stochastically simulates GC/CT and syphilis diagnosis and
#'              treatment.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
stitx_msm <- function(dat, at) {

  ## Inputs
  # Attributes
  recentpartners    <- get_attr(dat, "recentpartners")
  race              <- get_attr(dat, "race")

  # Infection status
  status            <- get_attr(dat, "status")
  role.class        <- get_attr(dat, "role.class")
  rGC               <- get_attr(dat, "rGC")
  uGC               <- get_attr(dat, "uGC")
  rCT               <- get_attr(dat, "rCT")
  uCT               <- get_attr(dat, "uCT")
  # Infection time
  rGC.infTime       <- get_attr(dat, "rGC.infTime")
  uGC.infTime       <- get_attr(dat, "uGC.infTime")
  rCT.infTime       <- get_attr(dat, "rCT.infTime")
  uCT.infTime       <- get_attr(dat, "uCT.infTime")
  # Symptoms
  rGC.sympt         <- get_attr(dat, "rGC.sympt")
  uGC.sympt         <- get_attr(dat, "uGC.sympt")
  rCT.sympt         <- get_attr(dat, "rCT.sympt")
  uCT.sympt         <- get_attr(dat, "uCT.sympt")
  # Tx
  uGC.tx            <- get_attr(dat, "uGC.tx")
  uGC.tx.prep       <- get_attr(dat, "uGC.tx.prep")
  uGC.tx.ept        <- get_attr(dat, "uGC.tx.ept")
  rGC.tx            <- get_attr(dat, "rGC.tx")
  rGC.tx.prep       <- get_attr(dat, "rGC.tx.prep")
  rGC.tx.ept        <- get_attr(dat, "rGC.tx.ept")
  uCT.tx            <- get_attr(dat, "uCT.tx")
  uCT.tx.prep       <- get_attr(dat, "uCT.tx.prep")
  uCT.tx.ept        <- get_attr(dat, "uCT.tx.ept")
  rCT.tx            <- get_attr(dat, "rCT.tx")
  rCT.tx.prep       <- get_attr(dat, "rCT.tx.prep")
  rCT.tx.ept        <- get_attr(dat, "rCT.tx.ept")
  prepStartTime     <- get_attr(dat, "prepStartTime")
  prepLastStiScreen <- get_attr(dat, "prepLastStiScreen")
  # tt.traj
  tt.traj.gc.hivpos <- get_attr(dat, "tt.traj.gc.hivpos")
  tt.traj.gc.hivneg <- get_attr(dat, "tt.traj.gc.hivneg")
  tt.traj.ct.hivpos <- get_attr(dat, "tt.traj.ct.hivpos")
  tt.traj.ct.hivneg <- get_attr(dat, "tt.traj.ct.hivneg")
  # Diagnosis/ Attributes from testing
  diag.status.gc    <- get_attr(dat, "diag.status.gc")
  diag.status.ct    <- get_attr(dat, "diag.status.ct")
  tsinceltst.rgc    <- get_attr(dat, "time.since.last.test.rgc")
  tsinceltst.ugc    <- get_attr(dat, "time.since.last.test.ugc")
  tsinceltst.rct    <- get_attr(dat, "time.since.last.test.rct")
  tsinceltst.uct    <- get_attr(dat, "time.since.last.test.uct")
  # EPT
  eptpartEligTx_GC  <- get_attr(dat, "eptpartEligTx_GC")
  eptpartEligTx_CT  <- get_attr(dat, "eptpartEligTx_CT")
  eptpartEligTxdate <- get_attr(dat, "eptpartEligTxdate")
  eptindexEligdate  <- get_attr(dat, "eptindexEligdate")
  eptindexElig      <- get_attr(dat, "eptindexElig")

  # Parameters
  gc.sympt.prob.tx     <- get_param(dat, "gc.sympt.prob.tx")
  ct.sympt.prob.tx     <- get_param(dat, "ct.sympt.prob.tx")
  gc.asympt.prob.tx    <- get_param(dat, "gc.asympt.prob.tx")
  ct.asympt.prob.tx    <- get_param(dat, "ct.asympt.prob.tx")
  prep.sti.screen.int  <- get_param(dat, "prep.sti.screen.int")
  prep.sti.prob.tx     <- get_param(dat, "prep.sti.prob.tx")
  ept.gc.success       <- get_param(dat, "ept.gc.success")
  ept.ct.success       <- get_param(dat, "ept.ct.success")
  ept.coverage         <- get_param(dat, "ept.coverage")
  ept.cov.rate         <- get_param(dat, "ept.cov.rate")
  sti.correlation.time <- get_param(dat, "sti.correlation.time")

  # Gonorrhea ---------------------------------------------------------------

  # symptomatic gc treatment
  idsRGC_tx_sympt <- which(
    rGC == 1 &
    rGC.infTime < at &
    rGC.sympt == 1 &
    (is.na(rGC.tx) | rGC.tx == 0)
  )

  idsUGC_tx_sympt <- which(
    uGC == 1 &
    uGC.infTime < at &
    uGC.sympt == 1 &
    (is.na(uGC.tx) | uGC.tx == 0)
  )

  idsGC_tx_sympt <- c(idsRGC_tx_sympt, idsUGC_tx_sympt)

  probs <- gc.sympt.prob.tx[race[idsGC_tx_sympt]]
  txGC_sympt <- idsGC_tx_sympt[runif(length(idsGC_tx_sympt)) < probs]
  txRGC_sympt <- intersect(idsRGC_tx_sympt, txGC_sympt)
  txUGC_sympt <- intersect(idsUGC_tx_sympt, txGC_sympt)

  # asymptomatic gc treatment
  idsRGC_tx_asympt <- which(
    rGC == 1 &
    rGC.infTime < at &
    rGC.sympt == 0 &
    diag.status.gc == 1 &
    (is.na(rGC.tx) | rGC.tx == 0)
  )

  idsUGC_tx_asympt <- which(
    uGC == 1 &
    uGC.infTime < at &
    uGC.sympt == 0 &
    diag.status.gc == 1 &
    (is.na(uGC.tx) | uGC.tx == 0)
  )

  idsGC_tx_asympt <- c(idsRGC_tx_asympt, idsUGC_tx_asympt)

  probs <- gc.asympt.prob.tx[race[idsGC_tx_asympt]]
  txGC_asympt <- idsGC_tx_asympt[runif(length(idsGC_tx_asympt)) < probs]

  txRGC_asympt <- intersect(idsRGC_tx_asympt, txGC_asympt)
  txUGC_asympt <- intersect(idsUGC_tx_asympt, txGC_asympt)

  # all treated GC
  txRGC <- union(txRGC_sympt, txRGC_asympt)
  txUGC <- union(txUGC_sympt, txUGC_asympt)

  idsRGC_tx <- union(idsRGC_tx_sympt, idsRGC_tx_asympt)
  idsUGC_tx <- union(idsUGC_tx_sympt, idsUGC_tx_asympt)


  # Chlamydia ---------------------------------------------------------------

  # symptomatic ct treatment
  idsRCT_tx_sympt <- which(
    rCT == 1 &
    rCT.infTime < at &
    rCT.sympt == 1 &
    (is.na(rCT.tx) | rCT.tx == 0)
  )

  idsUCT_tx_sympt <- which(
    uCT == 1 &
    uCT.infTime < at &
    uCT.sympt == 1 &
    (is.na(uCT.tx) | uCT.tx == 0)
  )

  idsCT_tx_sympt <- c(idsRCT_tx_sympt, idsUCT_tx_sympt)

  probs <- ct.sympt.prob.tx[race[idsCT_tx_sympt]]
  txCT_sympt <- idsCT_tx_sympt[runif(length(idsCT_tx_sympt)) < probs]
  txRCT_sympt <- intersect(idsRCT_tx_sympt, txCT_sympt)
  txUCT_sympt <- intersect(idsUCT_tx_sympt, txCT_sympt)

  # asymptomatic ct treatment
  idsRCT_tx_asympt <- which(
    rCT == 1 &
    rCT.infTime < at &
    rCT.sympt == 0 &
    diag.status.ct == 1 &
    (is.na(rCT.tx) | rCT.tx == 0)
  )

  idsUCT_tx_asympt <- which(
    uCT == 1 &
    uCT.infTime < at &
    uCT.sympt == 0 &
    diag.status.ct == 1 &
    (is.na(uCT.tx) | uCT.tx == 0)
  )

  idsCT_tx_asympt <- c(idsRCT_tx_asympt, idsUCT_tx_asympt)

  probs <- ct.asympt.prob.tx[race[idsCT_tx_asympt]]
  txCT_asympt <- idsCT_tx_asympt[runif(length(idsCT_tx_asympt)) < probs]

  txRCT_asympt <- intersect(idsRCT_tx_asympt, txCT_asympt)
  txUCT_asympt <- intersect(idsUCT_tx_asympt, txCT_asympt)

  # all treated CT
  txRCT <- union(txRCT_sympt, txRCT_asympt)
  txUCT <- union(txUCT_sympt, txUCT_asympt)

  idsRCT_tx <- union(idsRCT_tx_sympt, idsRCT_tx_asympt)
  idsUCT_tx <- union(idsUCT_tx_sympt, idsUCT_tx_asympt)


  # PrEP-related treatment for all STIs --------------------------------------

  # Interval-based treatment for MSM on PrEP
  idsSTI_screen <- which(
    prepStartTime == at |
    (at - prepLastStiScreen >= prep.sti.screen.int)
  )

  dat <- set_epi(dat, "sti.screening.ep.prep", at, length(idsSTI_screen))

  prepLastStiScreen[idsSTI_screen] <- at

  idsRGC_prep_tx <- intersect(
    idsSTI_screen,
    which(
      rGC == 1 &
      rGC.infTime < at &
      (is.na(rGC.tx.prep) | rGC.tx.prep == 0)
    )
  )

  idsUGC_prep_tx <- intersect(
    idsSTI_screen,
    which(
      uGC == 1 &
      uGC.infTime < at &
      (is.na(uGC.tx.prep) | uGC.tx.prep == 0)
    )
  )

  idsRCT_prep_tx <- intersect(
    idsSTI_screen,
    which(
      rCT == 1 &
      rCT.infTime < at &
      (is.na(rCT.tx.prep) | rCT.tx.prep == 0)
    )
  )

  idsUCT_prep_tx <- intersect(
    idsSTI_screen,
    which(
      uCT == 1 &
      uCT.infTime < at &
      (is.na(uCT.tx.prep) | uCT.tx.prep == 0)
    )
  )

  txRGC_prep <- idsRGC_prep_tx[runif(length(idsRGC_prep_tx)) < prep.sti.prob.tx]
  txUGC_prep <- idsUGC_prep_tx[runif(length(idsUGC_prep_tx)) < prep.sti.prob.tx]
  txRCT_prep <- idsRCT_prep_tx[runif(length(idsRCT_prep_tx)) < prep.sti.prob.tx]
  txUCT_prep <- idsUCT_prep_tx[runif(length(idsUCT_prep_tx)) < prep.sti.prob.tx]

  # Summarize all treated for each STI - EPT-treated not eligible to provide EPT to their partners
  txRGC_all <- c(txRGC, txRGC_prep)
  txUGC_all <- c(txUGC, txUGC_prep)
  txRCT_all <- c(txRCT, txRCT_prep)
  txUCT_all <- c(txUCT, txUCT_prep)

  # Subset all treated for GC/CT to treated with partners (for EPT)
  ept_txRGC_all <- txRGC_all[recentpartners[txRGC_all] > 0]
  ept_txUGC_all <- txUGC_all[recentpartners[txUGC_all] > 0]
  ept_txRCT_all <- txRCT_all[recentpartners[txRCT_all] > 0]
  ept_txUCT_all <- txUCT_all[recentpartners[txUCT_all] > 0]
  ept_txGC_all <- c(ept_txRGC_all, ept_txUGC_all)
  ept_txCT_all <- c(ept_txRCT_all, ept_txUCT_all)
  ept_tx_all <- unique(c(
    ept_txRGC_all,
    ept_txUGC_all,
    ept_txRCT_all,
    ept_txUCT_all
  ))

  # Update EPT index status and eligibility for GC/CT treated with partners
  if (length(ept_tx_all > 0)) {
    dat <- set_attr(dat, "eptindexElig", 1, posit_ids = ept_tx_all)
    dat <- set_attr(dat, "eptindexStat", 0, posit_ids = ept_tx_all)
    dat <- set_attr(dat, "eptindexEligdate", at, posit_ids = ept_tx_all)
  }

  # EPT Treatment for Non-index (no test is done) ------------------------------

  # Have prevalent infection, are eligible for tx through EPT, are untreated,
  # are not previously assigned for EPT tx, and were provided/uptake EPT last
  # time step
  idsRGC_tx_ept <- which(
    rGC == 1 &
    rGC.infTime < at &
    eptpartEligTx_GC == 1 &
    eptpartEligTxdate == (at - 1) &
    (is.na(rGC.tx) | rGC.tx == 0) &
    (is.na(rGC.tx.prep) | rGC.tx.prep == 0) &
    is.na(rGC.tx.ept)
  )

  idsUGC_tx_ept <- which(
    uGC == 1 &
    uGC.infTime < at &
    eptpartEligTx_GC == 1 &
    eptpartEligTxdate == (at - 1) &
    (is.na(uGC.tx) | uGC.tx == 0) &
    (is.na(uGC.tx.prep) | uGC.tx.prep == 0) &
    is.na(uGC.tx.ept)
  )

  idsGC_tx_ept <- c(idsRGC_tx_ept, idsUGC_tx_ept)
  txGC_ept <- idsGC_tx_ept[runif(length(idsGC_tx_ept)) < ept.gc.success]
  txRGC_ept <- intersect(idsRGC_tx_ept, txGC_ept)
  txUGC_ept <- intersect(idsUGC_tx_ept, txGC_ept)

  idsRCT_tx_ept <- which(
    rCT == 1 &
    rCT.infTime < at &
    eptpartEligTx_CT == 1 &
    eptpartEligTxdate == (at - 1) &
    (is.na(rCT.tx)  | rCT.tx == 0) &
    (is.na(rCT.tx.prep) | rCT.tx.prep == 0) &
    is.na(rCT.tx.ept)
  )

  idsUCT_tx_ept <- which(
    uCT == 1 &
    uCT.infTime < at &
    eptpartEligTx_CT == 1 &
    eptpartEligTxdate == (at - 1) &
    (is.na(uCT.tx)  | uCT.tx == 0) &
    (is.na(uCT.tx.prep) | uCT.tx.prep == 0) &
    is.na(uCT.tx.ept)
  )

  idsCT_tx_ept <- c(idsRCT_tx_ept, idsUCT_tx_ept)

  txCT_ept <- idsCT_tx_ept[runif(length(idsCT_tx_ept)) < ept.ct.success]
  txRCT_ept <- intersect(idsRCT_tx_ept, txCT_ept)
  txUCT_ept <- intersect(idsUCT_tx_ept, txCT_ept)

  # All EPT-treated index ids
  allidsept <- unique(c(idsCT_tx_ept, idsGC_tx_ept))

  # Summarize all successfully treated for each STI, now including EPT
  alltxRGC <- c(txRGC, txRGC_prep, txRGC_ept)
  alltxUGC <- c(txUGC, txUGC_prep, txUGC_ept)
  alltxRCT <- c(txRCT, txRCT_prep, txRCT_ept)
  alltxUCT <- c(txUCT, txUCT_prep, txUCT_ept)
  alltxEPT <- c(txRGC_ept, txUGC_ept, txRCT_ept, txUCT_ept)

  # EPT Initiation for Index Partner -------------------------------------------

  # Eligibility only lasts one time step - so coverage is 0 for current eligibles
  eptCov <- 0
  idsEligSt <- which(
    eptindexElig == 1 &
    eptindexEligdate == at
  )
  nEligSt <- length(idsEligSt)

  nStart <- max(0, min(nEligSt, round((ept.coverage - eptCov) * nEligSt)))
  ept_idsStart <- numeric()
  if (nStart > 0) {
    if (ept.cov.rate >= 1) {
      ept_idsStart <- ssample(idsEligSt, nStart)
    } else {
      ept_idsStart <- idsEligSt[runif(nStart) <  ept.cov.rate]
    }
  }
  eptCov <- (length(ept_idsStart)) / nEligSt

  # Update EPT index status for those selected to receive EPT for their partners
  dat <- set_attr(dat, "eptindexStat", 1, posit_ids = ept_idsStart)

  index_gc <- intersect(ept_txGC_all, ept_idsStart)
  index_ct <- intersect(ept_txCT_all, ept_idsStart)

  # Correlated testing for other STIs if symptomatic for one -------------------

  # All treated for other site of STIs, minus those getting treated for STI through EPT
  tst.rgc <- setdiff(txRCT_sympt, txRGC_ept)
  tst.ugc <- setdiff(txUCT_sympt, txUGC_ept)
  tst.rct <- setdiff(txRGC_sympt, txRCT_ept)
  tst.uct <- setdiff(txUGC_sympt, txUCT_ept)

  # Remove those just treated for STI (either sympt/asympt) from testing
  tst.rgc <- setdiff(tst.rgc, txRGC)
  tst.ugc <- setdiff(tst.ugc, txUGC)
  tst.rct <- setdiff(tst.rct, txRCT)
  tst.uct <- setdiff(tst.uct, txUCT)

  # Subset to those not tested for particular STI recently (in last 12 weeks)
  tst.rgc <- tst.rgc[which(
    tsinceltst.rgc[tst.rgc] > sti.correlation.time &
      (
        is.na(diag.status.gc[tst.rgc]) |
        diag.status.gc[tst.rgc]
      ) &
      role.class[tst.rgc] %in% c("R", "V")
    )]

  tst.ugc <- tst.ugc[which(
    tsinceltst.ugc[tst.ugc] > sti.correlation.time &
    (
      is.na(diag.status.gc[tst.ugc]) |
      diag.status.gc[tst.ugc]
    ) &
    role.class[tst.ugc] %in% c("I", "V")
  )]

  tst.rct <- tst.rct[which(
    tsinceltst.rct[tst.rct] > sti.correlation.time &
    (
      is.na(diag.status.ct[tst.rct]) |
      diag.status.ct[tst.rct]
    ) &
    role.class[tst.rct] %in% c("R", "V")
  )]

  tst.uct <- tst.uct[which(
    tsinceltst.uct[tst.uct] > sti.correlation.time &
    (
      is.na(diag.status.ct[tst.uct]) |
      diag.status.ct[tst.uct]
    ) &
    role.class[tst.uct] %in% c("I", "V")
  )]

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

  # Output ---------------------------------------------------------------------

  # PrEP
  dat$attr$prepLastStiScreen <- prepLastStiScreen

  # Gonorrhea
  dat <- set_attr(dat, "rGC.tx", 0, posit_ids = idsRGC_tx)
  dat <- set_attr(dat, "rGC.tx", 1, posit_ids = txRGC)
  dat <- set_attr(dat, "uGC.tx", 0, posit_ids = idsUGC_tx)
  dat <- set_attr(dat, "uGC.tx", 1, posit_ids = txUGC)
  dat <- set_attr(dat, "rGC.tx.prep", 0, posit_ids = idsRGC_prep_tx)
  dat <- set_attr(dat, "rGC.tx.prep", 1, posit_ids = txRGC_prep)
  dat <- set_attr(dat, "rGC.tx.ept", 0, posit_ids = idsRGC_tx_ept)
  dat <- set_attr(dat, "rGC.tx.ept", 1, posit_ids = txRGC_ept)
  dat <- set_attr(dat, "last.tx.time.rgc", at, posit_ids = txRGC_all)
  dat <- set_attr(dat, "last.tx.time.ugc", at, posit_ids = txUGC_all)
  dat <- set_attr(dat, "last.tx.time.rgc.prep", at, posit_ids = txRGC_prep)
  dat <- set_attr(dat, "last.tx.time.ugc.prep", at, posit_ids = txUGC_prep)
  dat <- set_attr(dat, "uGC.tx.prep", 0, posit_ids = idsUGC_prep_tx)
  dat <- set_attr(dat, "uGC.tx.prep", 1, posit_ids = txUGC_prep)
  dat <- set_attr(dat, "uGC.tx.ept", 0, posit_ids = idsUGC_tx_ept)
  dat <- set_attr(dat, "uGC.tx.ept", 1, posit_ids = txUGC_ept)

  uGC.tx <- get_attr(dat, "uGC.tx")
  rGC.tx <- get_attr(dat, "rGC.tx")
  uGC.tx.prep <- get_attr(dat, "uGC.tx.prep")
  rGC.tx.prep <- get_attr(dat, "rGC.tx.prep")
  uGC.tx.ept <- get_attr(dat, "uGC.tx.ept")
  rGC.tx.ept <- get_attr(dat, "rGC.tx.ept")
  rGC <- get_attr(dat, "rGC")
  uGC <- get_attr(dat, "uGC")

  dat <- set_attr(dat, "rGC.tx", 1, posit_ids = (uGC.tx == 1 | uGC.tx.prep == 1 | uGC.tx.ept == 1) & rGC == 1)
  dat <- set_attr(dat, "uGC.tx", 1, posit_ids = (rGC.tx == 1 | rGC.tx.prep == 1 | rGC.tx.ept == 1) & uGC == 1)

  dat <- set_attr(dat, "last.diag.time.gc", at, posit_ids = txGC_sympt)
  dat <- set_attr(dat, "last.neg.test.rgc", NA, posit_ids = txRGC_sympt)
  dat <- set_attr(dat, "last.neg.test.ugc", NA, posit_ids = txUGC_sympt)
  dat <- set_attr(dat, "time.since.last.test.rgc", 0, posit_ids = txRGC_sympt)
  dat <- set_attr(dat, "time.since.last.test.ugc", 0, posit_ids = txUGC_sympt)
  dat <- set_attr(dat, "diag.status.gc", 0, posit_ids = idsGC_tx_sympt)
  dat <- set_attr(dat, "diag.status.gc", 1, posit_ids = txGC_sympt)

  ## Correlated Gonorrhea testing
  dat <- set_attr(dat, "last.neg.test.rgc", at, posit_ids = tst.rgc.neg)
  dat <- set_attr(dat, "last.neg.test.ugc", at, posit_ids = tst.ugc.neg)
  dat <- set_attr(dat, "last.neg.test.rgc", NA, posit_ids = tst.rgc.pos)
  dat <- set_attr(dat, "last.neg.test.ugc", NA, posit_ids = tst.ugc.pos)
  dat <- set_attr(dat, "diag.status.gc", 1, posit_ids = tst.gc.pos)
  dat <- set_attr(dat, "last.diag.time.gc", at, posit_ids = tst.gc.pos)
  dat <- set_attr(dat, "time.since.last.test.rgc", 0, posit_ids = tst.rgc)
  dat <- set_attr(dat, "time.since.last.test.ugc", 0, posit_ids = tst.ugc)

  # Chlamydia
  dat <- set_attr(dat, "rCT.tx", 0, posit_ids = idsRCT_tx)
  dat <- set_attr(dat, "rCT.tx", 1, posit_ids = txRCT)
  dat <- set_attr(dat, "uCT.tx", 0, posit_ids = idsUCT_tx)
  dat <- set_attr(dat, "uCT.tx", 1, posit_ids = txUCT)
  dat <- set_attr(dat, "rCT.tx.prep", 0, posit_ids = idsRCT_prep_tx)
  dat <- set_attr(dat, "rCT.tx.prep", 1, posit_ids = txRCT_prep)
  dat <- set_attr(dat, "rCT.tx.ept", 0, posit_ids = idsRCT_tx_ept)
  dat <- set_attr(dat, "rCT.tx.ept", 1, posit_ids = txRCT_ept)
  dat <- set_attr(dat, "last.tx.time.rct", at, posit_ids = txRCT_all)
  dat <- set_attr(dat, "last.tx.time.uct", at, posit_ids = txUCT_all)
  dat <- set_attr(dat, "last.tx.time.rct.prep", at, posit_ids = txRCT_prep)
  dat <- set_attr(dat, "last.tx.time.uct.prep", at, posit_ids = txUCT_prep)
  dat <- set_attr(dat, "uCT.tx.prep", 0, posit_ids = idsUCT_prep_tx)
  dat <- set_attr(dat, "uCT.tx.prep", 1, posit_ids = txUCT_prep)
  dat <- set_attr(dat, "uCT.tx.ept", 0, posit_ids = idsUCT_tx_ept)
  dat <- set_attr(dat, "uCT.tx.ept", 1, posit_ids = txUCT_ept)
  dat <- set_attr(dat, "last.diag.time.ct", at, posit_ids = txCT_sympt)
  dat <- set_attr(dat, "last.neg.test.rct", NA, posit_ids = txRCT_sympt)
  dat <- set_attr(dat, "last.neg.test.uct", NA, posit_ids = txUCT_sympt)
  dat <- set_attr(dat, "time.since.last.test.rct", 0, posit_ids = txRCT_sympt)
  dat <- set_attr(dat, "time.since.last.test.uct", 0, posit_ids = txUCT_sympt)
  dat <- set_attr(dat, "diag.status.ct", 0, posit_ids = idsCT_tx_sympt)
  dat <- set_attr(dat, "diag.status.ct", 1, posit_ids = txCT_sympt)

  uCT.tx <- get_attr(dat, "uCT.tx")
  rCT.tx <- get_attr(dat, "rCT.tx")
  uCT.tx.prep <- get_attr(dat, "uCT.tx.prep")
  rCT.tx.prep <- get_attr(dat, "rCT.tx.prep")
  uCT.tx.ept <- get_attr(dat, "uCT.tx.ept")
  rCT.tx.ept <- get_attr(dat, "rCT.tx.ept")
  rCT <- get_attr(dat, "rCT")
  uCT <- get_attr(dat, "uCT")

  dat <- set_attr(dat, "rCT.tx", 1, posit_ids = (uCT.tx == 1 | uCT.tx.prep == 1 | uCT.tx.ept == 1) & rCT == 1)
  dat <- set_attr(dat, "uCT.tx", 1, posit_ids = (rCT.tx == 1 | rCT.tx.prep == 1 | rCT.tx.ept == 1) & uCT == 1)

  ## Correlated Chlamydia testing
  dat <- set_attr(dat, "last.neg.test.rct", at, posit_ids = tst.rct.neg)
  dat <- set_attr(dat, "last.neg.test.uct", at, posit_ids = tst.uct.neg)
  dat <- set_attr(dat, "last.neg.test.rct", NA, posit_ids = tst.rct.pos)
  dat <- set_attr(dat, "last.neg.test.uct", NA, posit_ids = tst.uct.pos)
  dat <- set_attr(dat, "diag.status.ct", 1, posit_ids = tst.ct.pos)
  dat <- set_attr(dat, "last.diag.time.ct", at, posit_ids = tst.ct.pos)
  dat <- set_attr(dat, "time.since.last.test.rct", 0, posit_ids = tst.rct)
  dat <- set_attr(dat, "time.since.last.test.uct", 0, posit_ids = tst.uct)

  ## Proportion of infections treated in past year
  #dat$epi$tx.gc.prop[at] <- ifelse(length(which(at - dat$attr$last.rGC.infTime <= 52)) +
  #                                   length(which(at - dat$attr$last.uGC.infTime <= 52)) == 0 |
  #                                   is.na(length(which(at - dat$attr$last.rGC.infTime <= 52)) +
  #                                           length(which(at - dat$attr$last.uGC.infTime <= 52))) |
  #                                   is.nan(length(which(at - dat$attr$last.rGC.infTime <= 52)) +
  #                                            length(which(at - dat$attr$last.uGC.infTime <= 52))) |
  #                                   is.null(length(which(at - dat$attr$last.rGC.infTime <= 52)) +
  #                                             length(which(at - dat$attr$last.uGC.infTime <= 52))),
  #                                 0,
  #                                 (length(which(at - dat$attr$last.tx.time.rgc <= 52 &
  #                                                 at - dat$attr$last.rGC.infTime <= 52)) +
  #                                    length(which(at - dat$attr$last.tx.time.ugc <= 52 &
  #                                                   at - dat$attr$last.uGC.infTime <= 52))) /
  #                                   (length(which(at - dat$attr$last.rGC.infTime <= 52)) +
  #                                      length(which(at - dat$attr$last.uGC.infTime <= 52))))

  #dat$epi$tx.ct.prop[at] <- ifelse(length(which(at - dat$attr$last.rCT.infTime <= 52)) +
  #                                   length(which(at - dat$attr$last.uCT.infTime <= 52)) == 0 |
  #                                   is.na(length(which(at - dat$attr$last.rCT.infTime <= 52)) +
  #                                           length(which(at - dat$attr$last.uCT.infTime <= 52))) |
  #                                   is.nan(length(which(at - dat$attr$last.rCT.infTime <= 52)) +
  #                                            length(which(at - dat$attr$last.uCT.infTime <= 52))) |
  #                                   is.null(length(which(at - dat$attr$last.rCT.infTime <= 52)) +
  #                                             length(which(at - dat$attr$last.uCT.infTime <= 52))),
  #                                 0,
  #                                 (length(which(at - dat$attr$last.tx.time.rct <= 52 &
  #                                                 at - dat$attr$last.rCT.infTime <= 52)) +
  #                                    length(which(at - dat$attr$last.tx.time.uct <= 52 &
  #                                                   at - dat$attr$last.uCT.infTime <= 52))) /
  #                                   (length(which(at - dat$attr$last.rCT.infTime <= 52)) +
  #                                      length(which(at - dat$attr$last.uCTC.infTime <= 52))))

  #dat$epi$tx.gcct.prop[at] <- ifelse((length(which(at - dat$attr$last.rGC.infTime <= 52)) +
  #                                      length(which(at - dat$attr$last.uGC.infTime <= 52)) +
  #                                      length(which(at - dat$attr$last.rCT.infTime <= 52)) +
  #                                      length(which(at - dat$attr$last.uCT.infTime <= 52))) == 0 |
  #                                     is.na(length(which(at - dat$attr$last.rGC.infTime <= 52)) +
  #                                             length(which(at - dat$attr$last.uGC.infTime <= 52)) +
  #                                             length(which(at - dat$attr$last.rCT.infTime <= 52)) +
  #                                             length(which(at - dat$attr$last.uCT.infTime <= 52))) |
  #                                     is.nan(length(which(at - dat$attr$last.rGC.infTime <= 52)) +
  #                                              length(which(at - dat$attr$last.uGC.infTime <= 52)) +
  #                                              length(which(at - dat$attr$last.rCT.infTime <= 52)) +
  #                                              length(which(at - dat$attr$last.uCT.infTime <= 52))) |
  #                                     is.null(length(which(at - dat$attr$last.rGC.infTime <= 52)) +
  #                                               length(which(at - dat$attr$last.uGC.infTime <= 52)) +
  #                                               length(which(at - dat$attr$last.rCT.infTime <= 52)) +
  #                                               length(which(at - dat$attr$last.uCT.infTime <= 52))), 0,
  #                                   (length(which(at - dat$attr$last.tx.time.rgc <= 52 &
  #                                                   at - dat$attr$last.rGC.infTime <= 52)) +
  #                                      length(which(at - dat$attr$last.tx.time.ugc <= 52 &
  #                                                     at - dat$attr$last.uGC.infTime <= 52)) +
  #                                      length(which(at - dat$attr$last.tx.time.rct <= 52 &
  #                                                     at - dat$attr$last.rCT.infTime <= 52)) +
  #                                      length(which(at - dat$attr$last.tx.time.uct <= 52 &
  #                                                     at - dat$attr$last.uCT.infTime <= 52))) /
  #                                     (length(which(at - dat$attr$last.rGC.infTime <= 52)) +
  #                                        length(which(at - dat$attr$last.uGC.infTime <= 52)) +
  #                                        length(which(at - dat$attr$last.rCT.infTime <= 52)) +
  #                                        length(which(at - dat$attr$last.uCT.infTime <= 52))))

  # Non-index EPT-treated
  dat <- set_attr(dat, "eptpartEligTx_GC", NA, posit_ids = txGC_ept)
  dat <- set_attr(dat, "eptpartEligTx_CT", NA, posit_ids = txCT_ept)
  dat <- set_attr(dat, "eptpartEligTxdate", NA, posit_ids = alltxEPT)
  dat <- set_attr(dat, "eptpartTx", 0, posit_ids = allidsept)
  dat <- set_attr(dat, "eptpartTx", 1, posit_ids = alltxEPT)

  # Epi Trackers - Tx count
  hivpos.ids <- which(status == 1)
  hivneg.ids <- which(status == 0)

  dat <- set_epi(dat, "gc.tx.sympt.hivpos", at, length(intersect(txGC_sympt, hivpos.ids)))
  dat <- set_epi(dat, "gc.tx.sympt.hivneg", at, length(intersect(txGC_sympt, hivneg.ids)))
  dat <- set_epi(dat, "gc.tx.asympt.hivpos", at, length(intersect(txGC_asympt, hivpos.ids)))
  dat <- set_epi(dat, "gc.tx.asympt.hivneg", at, length(intersect(txGC_asympt, hivneg.ids)))
  dat <- set_epi(dat, "gc.tx.prep", at, length(union(txRGC_prep, txUGC_prep)))

  dat <- set_epi(dat, "ct.tx.sympt.hivpos", at, length(intersect(txCT_sympt, hivpos.ids)))
  dat <- set_epi(dat, "ct.tx.sympt.hivneg", at, length(intersect(txCT_sympt, hivneg.ids)))
  dat <- set_epi(dat, "ct.tx.asympt.hivpos", at, length(intersect(txCT_asympt, hivpos.ids)))
  dat <- set_epi(dat, "ct.tx.asympt.hivneg", at, length(intersect(txCT_asympt, hivneg.ids)))
  dat <- set_epi(dat, "ct.tx.prep", at, length(union(txRCT_prep, txUCT_prep)))

  ## summary statistics
  #if (is.null(dat$epi$num.asympt.tx)) {
  #  dat$epi$num.asympt.tx <- rep(NA, length(dat$control$nsteps))
  #  dat$epi$num.asympt.cases <- rep(NA, length(dat$control$nsteps))
  #  dat$epi$num.asympt.tx.prep <- rep(NA, length(dat$control$nsteps))
  #  dat$epi$num.asympt.cases.prep <- rep(NA, length(dat$control$nsteps))
  #  dat$epi$num.rect.tx <- rep(NA, length(dat$epi$control$nsteps))
  #  dat$epi$num.rect.cases <- rep(NA, length(dat$epi$control$nsteps))
  #  dat$epi$num.rect.tx.prep <- rep(NA, length(dat$epi$control$nsteps))
  #  dat$epi$num.rect.cases.prep <- rep(NA, length(dat$epi$control$nsteps))
  #}

  ## Update number tested due to urethral or rectal symptoms
  #dat$epi$rGC_symptstidxtime[at] <- length(tst.rgc)
  #dat$epi$uGC_symptstidxtime[at] <- length(tst.ugc)
  #dat$epi$rCT_symptstidxtime[at] <- length(tst.rct)
  #dat$epi$uCT_symptstidxtime[at] <- length(tst.uct)

  ## Update number testing positive for other STI due to symptoms
  #dat$epi$rGC_pos_symptstidxtime[at] <- length(tst.rgc.pos)
  #dat$epi$uGC_pos_symptstidxtime[at] <- length(tst.ugc.pos)
  #dat$epi$rCT_pos_symptstidxtime[at] <- length(tst.rct.pos)
  #dat$epi$uCT_pos_symptstidxtime[at] <- length(tst.uct.pos)

  ## Update symptomatic tests to now include symptomatic STI-correlated testing
  #dat$epi$rGCsympttests[at] <- sum(length(unique(txRGC_sympt)),
  #                                 dat$epi$rGC_symptstidxtime[at],
  #                                 na.rm = TRUE)
  #dat$epi$uGCsympttests[at] <- sum(length(unique(txUGC_sympt)),
  #                                 dat$epi$uGC_symptstidxtime[at],
  #                                 na.rm = TRUE)
  #dat$epi$GCsympttests[at] <- sum(dat$epi$rGCsympttests[at],
  #                                dat$epi$uGCsympttests[at],
  #                                na.rm = TRUE)

  #dat$epi$rGCsympttests.pos[at] <- sum(length(unique(txRGC_sympt)),
  #                                     dat$epi$rGC_pos_symptstidxtime[at],
  #                                     na.rm = TRUE)
  #dat$epi$uGCsympttests.pos[at] <- sum(length(unique(txRGC_sympt)),
  #                                     dat$epi$uGC_pos_symptstidxtime[at],
  #                                     na.rm = TRUE)
  #dat$epi$GCsympttests.pos[at] <- sum(dat$epi$rGCsympttests.pos[at],
  #                                    dat$epi$uGCsympttests.pos[at],
  #                                    na.rm = TRUE)

  #dat$epi$rCTsympttests[at] <- sum(length(unique(txRCT_sympt)),
  #                                 dat$epi$RCT_symptstidxtime[at],
  #                                 na.rm = TRUE)
  #dat$epi$uCTsympttests[at] <- sum(length(unique(txUCT_sympt)),
  #                                 dat$epi$uCT_symptstidxtime[at],
  #                                 na.rm = TRUE)
  #dat$epi$CTsympttests[at] <- sum(dat$epi$rCTsympttests[at],
  #                                dat$epi$rCT_symptstidxtime[at],
  #                                dat$epi$uCTsympttests[at],
  #                                dat$epi$uCT_symptstidxtime[at],
  #                                na.rm = TRUE)

  #dat$epi$rCTsympttests.pos[at] <-  sum(length(unique(txRGC_sympt)),
  #                                      dat$epi$rCT_pos_symptstidxtime[at],
  #                                      na.rm = TRUE)
  #dat$epi$uCTsympttests.pos[at] <- sum(length(unique(txRGC_sympt)),
  #                                     dat$epi$uCT_pos_symptstidxtime[at],
  #                                     na.rm = TRUE)
  #dat$epi$CTsympttests.pos[at] <- sum(dat$epi$rCTsympttests.pos[at],
  #                                    dat$epi$uCTsympttests.pos[at],
  #                                    na.rm = TRUE)

  #dat$epi$stisympttests[at] <- sum(dat$epi$GCsympttests[at],
  #                                 dat$epi$rCTsympttests[at],
  #                                 na.rm = TRUE)

  #dat$epi$stisympttests.pos[at] <- sum(dat$epi$GCsympttests.pos[at],
  #                                     dat$epi$CTsympttests.pos[at],
  #                                     na.rm = TRUE)

  ## Asymptomatic treated
  #asympt.tx <- c(intersect(txRGC_all, which(dat$attr$rGC.sympt == 0)),
  #               intersect(txUGC_all, which(dat$attr$uGC.sympt == 0)),
  #               intersect(txRCT_all, which(dat$attr$rCT.sympt == 0)),
  #               intersect(txUCT_all, which(dat$attr$uCT.sympt == 0)))
  #dat$epi$num.asympt.tx[at] <- length(unique(asympt.tx))

  #asympt.cases <- c(idsRGC_tx_asympt, intersect(idsRGC_prep_tx, which(dat$attr$rGC.sympt == 0)),
  #                  idsUGC_tx_asympt, intersect(idsUGC_prep_tx, which(dat$attr$uGC.sympt == 0)),
  #                  idsRCT_tx_asympt, intersect(idsRCT_prep_tx, which(dat$attr$rCT.sympt == 0)),
  #                  idsUCT_tx_asympt, intersect(idsUCT_prep_tx, which(dat$attr$uCT.sympt == 0)))
  #dat$epi$num.asympt.cases[at] <- length(unique(asympt.cases))

  #asympt.tx.prep <- c(intersect(txRGC_prep, which(dat$attr$rGC.sympt == 0)),
  #                    intersect(txUGC_prep, which(dat$attr$uGC.sympt == 0)),
  #                    intersect(txRCT_prep, which(dat$attr$rCT.sympt == 0)),
  #                    intersect(txUCT_prep, which(dat$attr$uCT.sympt == 0)))
  #dat$epi$num.asympt.tx.prep[at] <- length(unique(asympt.tx.prep))

  #asympt.cases.prep <- c(intersect(idsRGC_prep_tx, which(dat$attr$rGC.sympt == 0)),
  #                       intersect(idsUGC_prep_tx, which(dat$attr$uGC.sympt == 0)),
  #                       intersect(idsRCT_prep_tx, which(dat$attr$rCT.sympt == 0)),
  #                       intersect(idsUCT_prep_tx, which(dat$attr$uCT.sympt == 0)))
  #dat$epi$num.asympt.cases.prep[at] <- length(unique(asympt.cases.prep))

  #rect.tx <- c(txRGC_all, txRCT_all)
  #dat$epi$num.rect.tx[at] <- length(unique(rect.tx))
  #rect.cases <- c(idsRGC_tx, idsRGC_prep_tx, idsRCT_tx, idsRCT_prep_tx)
  #dat$epi$num.rect.cases[at] <- length(unique(rect.cases))

  #rect.tx.prep <- c(txRGC_prep, txRCT_prep)
  #dat$epi$num.rect.tx.prep[at] <- length(unique(rect.tx.prep))
  #rect.cases.prep <- c(idsRGC_prep_tx, idsRCT_prep_tx)
  #dat$epi$num.rect.cases.prep[at] <- length(unique(rect.cases.prep))

  ## Track total number treated
  #dat$epi$txGC[at] <- length(unique(txRGC)) + length(unique(txUGC))
  #dat$epi$txGC_asympt[at] <- length(unique(txGC_asympt))
  #dat$epi$txCT[at] <- length(unique(txRCT)) + length(unique(txUCT))
  #dat$epi$txCT_asympt[at] <- length(unique(txCT_asympt))
  #dat$epi$txSTI_asympt[at] <- dat$epi$txGC_asympt[at] + dat$epi$txCT_asympt[at]
  #dat$epi$txSTI[at] <- dat$epi$txGC[at] + dat$epi$txCT[at]

  ##number of treatment episodes - cost-effectiveness analysis
  ##all
  #dat$epi$stitxepisodes[at] <- length(unique(c(txRGC, txUGC, txRCT, txUCT)))
  ##asymptomatic
  #dat$epi$stitxepisodes_asympt[at] <- length(unique(c(txRGC_asympt, txUGC_asympt, txRCT_asympt, txUCT_asympt)))

  ## Risk group-specific treatment and test counters
  #dat$epi$txGC.tttraj1[at] <- length(which(tt.traj.gc.hivneg[unique(c(txRGC, txUGC))] == 1 |
  #                                           tt.traj.gc.hivpos[unique(c(txRGC, txUGC))] == 1))
  #dat$epi$txGC_asympt.tttraj1[at] <- length(which(tt.traj.gc.hivneg[unique(txGC_asympt)] == 1 |
  #                                                  tt.traj.gc.hivpos[unique(txGC_asympt)] == 1))
  #dat$epi$txGC.tttraj2[at] <- length(which(tt.traj.gc.hivneg[unique(c(txRGC, txUGC))] == 2 |
  #                                           tt.traj.gc.hivpos[unique(c(txRGC, txUGC))] == 2))
  #dat$epi$txGC_asympt.tttraj2[at] <- length(which(tt.traj.gc.hivneg[unique(txGC_asympt)] == 2 |
  #                                                  tt.traj.gc.hivpos[unique(txGC_asympt)] == 2))

  #dat$epi$txCT.tttraj1[at] <- length(which(tt.traj.ct.hivpos[unique(c(txRCT, txUCT))] == 1 |
  #                                           tt.traj.ct.hivneg[unique(c(txRCT, txUCT))] == 1))
  #dat$epi$txCT_asympt.tttraj1[at] <- length(which(tt.traj.ct.hivpos[unique(txCT_asympt)] == 1 |
  #                                                  tt.traj.ct.hivneg[unique(txCT_asympt)] == 1))
  #dat$epi$txCT.tttraj2[at] <- length(which(tt.traj.ct.hivpos[unique(c(txRCT, txUCT))] == 2 |
  #                                           tt.traj.ct.hivneg[unique(c(txRCT, txUCT))] == 2))
  #dat$epi$txCT_asympt.tttraj2[at] <- length(which(tt.traj.ct.hivpos[unique(txCT_asympt)] == 2 |
  #                                                  tt.traj.ct.hivneg[unique(txCT_asympt)] == 2))

  #dat$epi$txSTI_asympt.tttraj1[at] <- dat$epi$txGC_asympt.tttraj1[at] + dat$epi$txCT_asympt.tttraj1[at]

  #dat$epi$txSTI_asympt.tttraj2[at] <- dat$epi$txGC_asympt.tttraj2[at] + dat$epi$txCT_asympt.tttraj2[at]

  #dat$epi$txSTI.tttraj1[at] <- dat$epi$txGC.tttraj1[at] + dat$epi$txCT.tttraj1[at]

  #dat$epi$txSTI.tttraj2[at] <- dat$epi$txGC.tttraj2[at] + dat$epi$txCT.tttraj2[at]

  #dat$epi$rGCsympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 1)) +
  #  length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 1))
  #dat$epi$uGCsympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 1)) +
  #  length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 1))
  #dat$epi$GCsympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 1)) +
  #  length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 1)) +
  #  length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 1)) +
  #  length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 1))

  #dat$epi$rGCsympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 2)) +
  #  length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 2))
  #dat$epi$uGCsympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 2)) +
  #  length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 2))
  #dat$epi$GCsympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 2)) +
  #  length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 2)) +
  #  length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 2)) +
  #  length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 2))

  #dat$epi$rCTsympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 1)) +
  #  length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 1))
  #dat$epi$uCTsympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 1)) +
  #  length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 1))
  #dat$epi$CTsympttests.tttraj1[at] <- length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 1)) +
  #  length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 1)) +
  #  length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 1)) +
  #  length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 1))

  #dat$epi$rCTsympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 2)) +
  #  length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 2))
  #dat$epi$uCTsympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 2)) +
  #  length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 2))
  #dat$epi$CTsympttests.tttraj2[at] <- length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 2)) +
  #  length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 2)) +
  #  length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 2)) +
  #  length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 2))

  #dat$epi$stisympttests.tttraj1[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 1)) +
  #  length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 1)) +
  #  length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 1)) +
  #  length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 1)) +
  #  length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 1)) +
  #  length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 1)) +
  #  length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 1)) +
  #  length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 1))

  #dat$epi$stisympttests.tttraj2[at] <- length(which(tt.traj.gc.hivpos[unique(txRGC_sympt)] == 2)) +
  #  length(which(tt.traj.gc.hivneg[unique(txRGC_sympt)] == 2)) +
  #  length(which(tt.traj.gc.hivpos[unique(txUGC_sympt)] == 2)) +
  #  length(which(tt.traj.gc.hivneg[unique(txUGC_sympt)] == 2)) +
  #  length(which(tt.traj.ct.hivpos[unique(txRCT_sympt)] == 2)) +
  #  length(which(tt.traj.ct.hivneg[unique(txRCT_sympt)] == 2)) +
  #  length(which(tt.traj.ct.hivpos[unique(txUCT_sympt)] == 2)) +
  #  length(which(tt.traj.ct.hivneg[unique(txUCT_sympt)] == 2))




  ## EPT
  ## Number of index provided with EPT
  #dat$epi$eptindexprovided_gc[at] <- length(index_gc)
  #dat$epi$eptindexprovided_ct[at] <- length(index_ct)

  ## Proportion of treated GC/CT index who have current partners - e.g. eligibility for EPT
  #dat$epi$propindexeptElig[at] <- ifelse(length(unique(c(txRGC_all, txUGC_all, txRCT_all, txUCT_all))) > 0,
  #                                       length(unique(ept_tx_all)) /
  #                                         length(unique(c(txRGC_all, txUGC_all, txRCT_all, txUCT_all))),
  #                                       NA)

  ## Proportion of eligible index who will receive EPT - varies by sim scenario
  #dat$epi$eptCov[at] <- eptCov

  ## Number of non-index treated due to EPT - mix of provision, uptake, tx success
  #dat$epi$eptTx[at] <- length(unique(alltxEPT))

  ## Proportion of all non-index eligible to be treated who had treatment success
  ## dat$epi$eptprop_tx[at] <- ifelse(length(unique(allidsept)) > 0,
  ##                                  length(unique(alltxEPT)) / length(unique(allidsept)),
  ##                                  0)

  return(dat)
}
