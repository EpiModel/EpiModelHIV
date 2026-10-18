
#' @title STI Treatment Module
#'
#' @description Stochastically simulates GC/CT diagnosis and treatment.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
stitx_msm <- function(dat, at) {

  ## Input
  # Attributes
  race <- get_attr(dat, "race")
  prepStat <- get_attr(dat, "prepStat")
  prepStartTime <- get_attr(dat, "prepStartTime")
  prepLastStiScreen <- get_attr(dat, "prepLastStiScreen")

  rGC <- get_attr(dat, "rGC")
  rGC.infTime <- get_attr(dat, "rGC.infTime")
  rGC.sympt <- get_attr(dat, "rGC.sympt")
  rGC.tx <- get_attr(dat, "rGC.tx")
  rGC.tx.prep <- get_attr(dat, "rGC.tx.prep")

  uGC <- get_attr(dat, "uGC")
  uGC.infTime <- get_attr(dat, "uGC.infTime")
  uGC.sympt <- get_attr(dat, "uGC.sympt")
  uGC.tx <- get_attr(dat, "uGC.tx")
  uGC.tx.prep <- get_attr(dat, "uGC.tx.prep")

  rCT <- get_attr(dat, "rCT")
  rCT.infTime <- get_attr(dat, "rCT.infTime")
  rCT.sympt <- get_attr(dat, "rCT.sympt")
  rCT.tx <- get_attr(dat, "rCT.tx")
  rCT.tx.prep <- get_attr(dat, "rCT.tx.prep")

  uCT <- get_attr(dat, "uCT")
  uCT.infTime <- get_attr(dat, "uCT.infTime")
  uCT.sympt <- get_attr(dat, "uCT.sympt")
  uCT.tx <- get_attr(dat, "uCT.tx")
  uCT.tx.prep <- get_attr(dat, "uCT.tx.prep")

  # Parameters
  gc.sympt.prob.tx <- get_param(dat, "gc.sympt.prob.tx")
  ct.sympt.prob.tx <- get_param(dat, "ct.sympt.prob.tx")
  gc.asympt.prob.tx <- get_param(dat, "gc.asympt.prob.tx")
  ct.asympt.prob.tx <- get_param(dat, "ct.asympt.prob.tx")

  prep.sti.screen.int <- get_param(dat, "prep.sti.screen.int")
  prep.sti.prob.tx <- get_param(dat, "prep.sti.prob.tx")


  ## Symptomatic GC Treatment ##
  idsRGC_tx_sympt <- which(
    rGC == 1 &
    rGC.infTime < at &
    rGC.sympt == 1 &
    is.na(rGC.tx)
  )
  idsUGC_tx_sympt <- which(
    uGC == 1 &
    uGC.infTime < at &
    uGC.sympt == 1 &
    is.na(uGC.tx)
  )

  # Subset by race
  idsGC_tx_sympt <- union(idsRGC_tx_sympt, idsUGC_tx_sympt)
  races <- sort(unique(race[idsGC_tx_sympt]))
  txGC_sympt <- rep(NA, length(idsGC_tx_sympt))
  for (i in races) {
    ids.race <- which(race[idsGC_tx_sympt] == i)
    txGC_sympt[ids.race] <- runif(length(ids.race)) < gc.sympt.prob.tx[i]
  }
  ids_txGC_sympt <- idsGC_tx_sympt[which(txGC_sympt == 1)]

  # Subset by site
  txRGC_sympt <- intersect(idsRGC_tx_sympt, ids_txGC_sympt)
  txUGC_sympt <- intersect(idsUGC_tx_sympt, ids_txGC_sympt)

  ## Asymptomatic GC Treatment ##
  idsRGC_tx_asympt <- which(
    rGC == 1 &
    rGC.infTime < at &
    rGC.sympt == 0 &
    is.na(rGC.tx) &
    prepStat == 0
  )
  idsUGC_tx_asympt <- which(
    uGC == 1 &
    uGC.infTime < at &
    uGC.sympt == 0 &
    is.na(uGC.tx) &
    prepStat == 0
  )

  # Subset by race
  idsGC_tx_asympt <- union(idsRGC_tx_asympt, idsUGC_tx_asympt)
  races <- sort(unique(race[idsGC_tx_asympt]))
  txGC_asympt <- rep(NA, length(idsGC_tx_asympt))
  for (i in races) {
    ids.race <- which(race[idsGC_tx_asympt] == i)
    txGC_asympt[ids.race] <- runif(length(ids.race)) < gc.asympt.prob.tx[i]
  }
  ids_txGC_asympt <- idsGC_tx_asympt[which(txGC_asympt == 1)]

  # Subset by site
  txRGC_asympt <- intersect(idsRGC_tx_asympt, ids_txGC_asympt)
  txUGC_asympt <- intersect(idsUGC_tx_asympt, ids_txGC_asympt)

  ## All Treated GC ##

  # IDs of men sucessfully treated
  txRGC <- union(txRGC_sympt, txRGC_asympt)
  txUGC <- union(txUGC_sympt, txUGC_asympt)

  # IDs of men eligible for treatment
  idsRGC_tx <- union(idsRGC_tx_sympt, idsRGC_tx_asympt)
  idsUGC_tx <- union(idsUGC_tx_sympt, idsUGC_tx_asympt)


  ## Symptomatic CT Treatment ##
  idsRCT_tx_sympt <- which(
    rCT == 1 &
    rCT.infTime < at &
    rCT.sympt == 1 &
    is.na(rCT.tx)
  )
  idsUCT_tx_sympt <- which(
    uCT == 1 &
    uCT.infTime < at &
    uCT.sympt == 1 &
    is.na(uCT.tx)
  )

  # Subset by race
  idsCT_tx_sympt <- union(idsRCT_tx_sympt, idsUCT_tx_sympt)
  races <- sort(unique(race[idsCT_tx_sympt]))
  txCT_sympt <- rep(NA, length(idsCT_tx_sympt))
  for (i in races) {
    ids.race <- which(race[idsCT_tx_sympt] == i)
    txCT_sympt[ids.race] <- runif(length(ids.race)) < ct.sympt.prob.tx[i]
  }
  ids_txCT_sympt <- idsCT_tx_sympt[which(txCT_sympt == 1)]

  # Subset by site
  txRCT_sympt <- intersect(idsRCT_tx_sympt, ids_txCT_sympt)
  txUCT_sympt <- intersect(idsUCT_tx_sympt, ids_txCT_sympt)


  ## Asymptomatic CT Treatment ##
  idsRCT_tx_asympt <- which(
    rCT == 1 &
    rCT.infTime < at &
    rCT.sympt == 0 &
    is.na(rCT.tx) &
    prepStat == 0
  )
  idsUCT_tx_asympt <- which(
    uCT == 1 &
    uCT.infTime < at &
    uCT.sympt == 0 &
    is.na(uCT.tx) &
    prepStat == 0
  )

  # Subset by race
  idsCT_tx_asympt <- union(idsRCT_tx_asympt, idsUCT_tx_asympt)
  races <- sort(unique(race[idsCT_tx_asympt]))
  txCT_asympt <- rep(NA, length(idsCT_tx_asympt))
  for (i in races) {
    ids.race <- which(race[idsCT_tx_asympt] == i)
    txCT_asympt[ids.race] <- runif(length(ids.race)) < ct.asympt.prob.tx[i]
  }
  ids_txCT_asympt <- idsCT_tx_asympt[which(txCT_asympt == 1)]

  # Subset by site
  txRCT_asympt <- intersect(idsRCT_tx_asympt, ids_txCT_asympt)
  txUCT_asympt <- intersect(idsUCT_tx_asympt, ids_txCT_asympt)

  ## All Treated CT ##
  txRCT <- union(txRCT_sympt, txRCT_asympt)
  txUCT <- union(txUCT_sympt, txUCT_asympt)

  idsRCT_tx <- union(idsRCT_tx_sympt, idsRCT_tx_asympt)
  idsUCT_tx <- union(idsUCT_tx_sympt, idsUCT_tx_asympt)


  ## Interval-based treatment for MSM on PrEP ##
  idsSTI_screen <- which(prepStartTime == at |
                           (at - prepLastStiScreen >= prep.sti.screen.int))

  dat <- set_attr(dat, "prepLastStiScreen", at, posit_ids = idsSTI_screen)


  idsRGC_prep_tx <- intersect(
    idsSTI_screen,
    which(
      rGC == 1 &
      rGC.infTime < at &
      is.na(rGC.tx.prep)
    )
  )
  idsUGC_prep_tx <- intersect(
    idsSTI_screen,
    which(
      uGC == 1 &
      uGC.infTime < at &
      is.na(uGC.tx.prep)
    )
  )
  idsRCT_prep_tx <- intersect(
    idsSTI_screen,
    which(
      rCT == 1 &
      rCT.infTime < at &
      is.na(rCT.tx.prep)
    )
  )
  idsUCT_prep_tx <- intersect(
    idsSTI_screen,
    which(
      uCT == 1 &
      uCT.infTime < at &
      is.na(uCT.tx.prep)
    )
  )

  txRGC_prep <- idsRGC_prep_tx[runif(length(idsRGC_prep_tx)) < prep.sti.prob.tx]
  txUGC_prep <- idsUGC_prep_tx[runif(length(idsUGC_prep_tx)) < prep.sti.prob.tx]
  txRCT_prep <- idsRCT_prep_tx[runif(length(idsRCT_prep_tx)) < prep.sti.prob.tx]
  txUCT_prep <- idsUCT_prep_tx[runif(length(idsUCT_prep_tx)) < prep.sti.prob.tx]


  ## Update Attributes ##
  rGC.tx[idsRGC_tx] <- 0
  rGC.tx[txRGC] <- 1

  uGC.tx[idsUGC_tx] <- 0
  uGC.tx[txUGC] <- 1

  rCT.tx[idsRCT_tx] <- 0
  rCT.tx[txRCT] <- 1

  uCT.tx[idsUCT_tx] <- 0
  uCT.tx[txUCT] <- 1

  rGC.tx.prep[idsRGC_prep_tx] <- 0
  rGC.tx.prep[txRGC_prep] <- 1

  uGC.tx.prep[idsUGC_prep_tx] <- 0
  uGC.tx.prep[txUGC_prep] <- 1

  rCT.tx.prep[idsRCT_prep_tx] <- 0
  rCT.tx.prep[txRCT_prep] <- 1

  uCT.tx.prep[idsUCT_prep_tx] <- 0
  uCT.tx.prep[txUCT_prep] <- 1

  ## Add tx at other anatomical site ##
  rGC.tx[(uGC.tx == 1 | uGC.tx.prep == 1) & rGC == 1] <- 1
  uGC.tx[(rGC.tx == 1 | rGC.tx.prep == 1) & uGC == 1] <- 1

  rCT.tx[(uCT.tx == 1 | uCT.tx.prep == 1) & rCT == 1] <- 1
  uCT.tx[(rCT.tx == 1 | rCT.tx.prep == 1) & uCT == 1] <- 1

  dat <- set_attr(dat, "rGC.tx", rGC.tx)
  dat <- set_attr(dat, "uGC.tx", uGC.tx)
  dat <- set_attr(dat, "rCT.tx", rCT.tx)
  dat <- set_attr(dat, "uCT.tx", uCT.tx)
  dat <- set_attr(dat, "rGC.tx.prep", rGC.tx.prep)
  dat <- set_attr(dat, "uGC.tx.prep", uGC.tx.prep)
  dat <- set_attr(dat, "rCT.tx.prep", rCT.tx.prep)
  dat <- set_attr(dat, "uCT.tx.prep", uCT.tx.prep)

  return(dat)
}
