
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

  # Parameters
  gc.sympt.prob.tx <- dat$param$gc.sympt.prob.tx
  ct.sympt.prob.tx <- dat$param$ct.sympt.prob.tx
  gc.asympt.prob.tx <- dat$param$gc.asympt.prob.tx
  ct.asympt.prob.tx <- dat$param$ct.asympt.prob.tx

  prep.sti.screen.int <- dat$param$prep.sti.screen.int
  prep.sti.prob.tx <- dat$param$prep.sti.prob.tx

  # Attributes
  race <- dat$attr$race

  ## Symptomatic GC Treatment ##
  idsRGC_tx_sympt <- which(dat$attr$rGC == 1 &
                           dat$attr$rGC.infTime < at &
                           dat$attr$rGC.sympt == 1 &
                           is.na(dat$attr$rGC.tx))
  idsUGC_tx_sympt <- which(dat$attr$uGC == 1 &
                           dat$attr$uGC.infTime < at &
                           dat$attr$uGC.sympt == 1 &
                           is.na(dat$attr$uGC.tx))

  # Subset by race
  idsGC_tx_sympt <- union(idsRGC_tx_sympt, idsUGC_tx_sympt)
  races <- sort(unique(race[idsGC_tx_sympt]))
  txGC_sympt <- rep(NA, length(idsGC_tx_sympt))
  for (i in races) {
    ids.race <- which(race[idsGC_tx_sympt] == i)
    txGC_sympt[ids.race] <- rbinom(length(ids.race), 1, gc.sympt.prob.tx[i])
  }
  ids_txGC_sympt <- idsGC_tx_sympt[which(txGC_sympt == 1)]

  # Subset by site
  txRGC_sympt <- intersect(idsRGC_tx_sympt, ids_txGC_sympt)
  txUGC_sympt <- intersect(idsUGC_tx_sympt, ids_txGC_sympt)

  ## Asymptomatic GC Treatment ##
  idsRGC_tx_asympt <- which(dat$attr$rGC == 1 &
                            dat$attr$rGC.infTime < at &
                            dat$attr$rGC.sympt == 0 &
                            is.na(dat$attr$rGC.tx) &
                            dat$attr$prepStat == 0)
  idsUGC_tx_asympt <- which(dat$attr$uGC == 1 &
                            dat$attr$uGC.infTime < at &
                            dat$attr$uGC.sympt == 0 &
                            is.na(dat$attr$uGC.tx) &
                            dat$attr$prepStat == 0)

  # Subset by race
  idsGC_tx_asympt <- union(idsRGC_tx_asympt, idsUGC_tx_asympt)
  races <- sort(unique(race[idsGC_tx_asympt]))
  txGC_asympt <- rep(NA, length(idsGC_tx_asympt))
  for (i in races) {
    ids.race <- which(race[idsGC_tx_asympt] == i)
    txGC_asympt[ids.race] <- rbinom(length(ids.race), 1, gc.asympt.prob.tx[i])
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
  idsRCT_tx_sympt <- which(dat$attr$rCT == 1 &
                           dat$attr$rCT.infTime < at &
                           dat$attr$rCT.sympt == 1 &
                           is.na(dat$attr$rCT.tx))
  idsUCT_tx_sympt <- which(dat$attr$uCT == 1 &
                           dat$attr$uCT.infTime < at &
                           dat$attr$uCT.sympt == 1 &
                           is.na(dat$attr$uCT.tx))

  # Subset by race
  idsCT_tx_sympt <- union(idsRCT_tx_sympt, idsUCT_tx_sympt)
  races <- sort(unique(race[idsCT_tx_sympt]))
  txCT_sympt <- rep(NA, length(idsCT_tx_sympt))
  for (i in races) {
    ids.race <- which(race[idsCT_tx_sympt] == i)
    txCT_sympt[ids.race] <- rbinom(length(ids.race), 1, ct.sympt.prob.tx[i])
  }
  ids_txCT_sympt <- idsCT_tx_sympt[which(txCT_sympt == 1)]

  # Subset by site
  txRCT_sympt <- intersect(idsRCT_tx_sympt, ids_txCT_sympt)
  txUCT_sympt <- intersect(idsUCT_tx_sympt, ids_txCT_sympt)


  ## Asymptomatic CT Treatment ##
  idsRCT_tx_asympt <- which(dat$attr$rCT == 1 &
                            dat$attr$rCT.infTime < at &
                            dat$attr$rCT.sympt == 0 &
                            is.na(dat$attr$rCT.tx) &
                            dat$attr$prepStat == 0)
  idsUCT_tx_asympt <- which(dat$attr$uCT == 1 &
                            dat$attr$uCT.infTime < at &
                            dat$attr$uCT.sympt == 0 &
                            is.na(dat$attr$uCT.tx) &
                            dat$attr$prepStat == 0)

  # Subset by race
  idsCT_tx_asympt <- union(idsRCT_tx_asympt, idsUCT_tx_asympt)
  races <- sort(unique(race[idsCT_tx_asympt]))
  txCT_asympt <- rep(NA, length(idsCT_tx_asympt))
  for (i in races) {
    ids.race <- which(race[idsCT_tx_asympt] == i)
    txCT_asympt[ids.race] <- rbinom(length(ids.race), 1, ct.asympt.prob.tx[i])
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
  idsSTI_screen <- which(dat$attr$prepStartTime == at |
                           (at - dat$attr$prepLastStiScreen >= prep.sti.screen.int))

  dat$attr$prepLastStiScreen[idsSTI_screen] <- at


  idsRGC_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$rGC == 1 &
                                    dat$attr$rGC.infTime < at &
                                    is.na(dat$attr$rGC.tx.prep)))
  idsUGC_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$uGC == 1 &
                                    dat$attr$uGC.infTime < at &
                                    is.na(dat$attr$uGC.tx.prep)))
  idsRCT_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$rCT == 1 &
                                    dat$attr$rCT.infTime < at &
                                    is.na(dat$attr$rCT.tx.prep)))
  idsUCT_prep_tx <- intersect(idsSTI_screen,
                              which(dat$attr$uCT == 1 &
                                    dat$attr$uCT.infTime < at &
                                    is.na(dat$attr$uCT.tx.prep)))

  txRGC_prep <- idsRGC_prep_tx[which(rbinom(length(idsRGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txUGC_prep <- idsUGC_prep_tx[which(rbinom(length(idsUGC_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txRCT_prep <- idsRCT_prep_tx[which(rbinom(length(idsRCT_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]
  txUCT_prep <- idsUCT_prep_tx[which(rbinom(length(idsUCT_prep_tx), 1,
                                            prep.sti.prob.tx) == 1)]


  ## Update Attributes ##
  dat$attr$rGC.tx[idsRGC_tx] <- 0
  dat$attr$rGC.tx[txRGC] <- 1

  dat$attr$uGC.tx[idsUGC_tx] <- 0
  dat$attr$uGC.tx[txUGC] <- 1

  dat$attr$rCT.tx[idsRCT_tx] <- 0
  dat$attr$rCT.tx[txRCT] <- 1

  dat$attr$uCT.tx[idsUCT_tx] <- 0
  dat$attr$uCT.tx[txUCT] <- 1

  dat$attr$rGC.tx.prep[idsRGC_prep_tx] <- 0
  dat$attr$rGC.tx.prep[txRGC_prep] <- 1

  dat$attr$uGC.tx.prep[idsUGC_prep_tx] <- 0
  dat$attr$uGC.tx.prep[txUGC_prep] <- 1

  dat$attr$rCT.tx.prep[idsRCT_prep_tx] <- 0
  dat$attr$rCT.tx.prep[txRCT_prep] <- 1

  dat$attr$uCT.tx.prep[idsUCT_prep_tx] <- 0
  dat$attr$uCT.tx.prep[txUCT_prep] <- 1


  ## Add tx at other anatomical site ##
  dat$attr$rGC.tx[which((dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1) &
                          dat$attr$rGC == 1)] <- 1
  dat$attr$uGC.tx[which((dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1) &
                          dat$attr$uGC == 1)] <- 1

  dat$attr$rCT.tx[which((dat$attr$uCT.tx == 1 | dat$attr$uCT.tx.prep == 1) &
                          dat$attr$rCT == 1)] <- 1
  dat$attr$uCT.tx[which((dat$attr$rCT.tx == 1 | dat$attr$rCT.tx.prep == 1) &
                          dat$attr$uCT == 1)] <- 1

  return(dat)
}
