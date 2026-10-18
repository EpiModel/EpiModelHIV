
#' @title STI Treatment Module
#'
#' @description Stochastically simulates GC/CT treatment.
#'
#' @inheritParams aging_msm
#'
#' @export
#'
stitx_msm <- function(dat, at) {


  # Setup ---------------------

  ## Input
  # Attributes
  rGC <- get_attr(dat, "rGC")
  rGC.infTime <- get_attr(dat, "rGC.infTime")
  rGC.sympt <- get_attr(dat, "rGC.sympt")
  uGC <- get_attr(dat, "uGC")
  uGC.infTime <- get_attr(dat, "uGC.infTime")
  uGC.sympt <- get_attr(dat, "uGC.sympt")
  GC.tx <- get_attr(dat, "GC.tx")

  rCT <- get_attr(dat, "rCT")
  rCT.infTime <- get_attr(dat, "rCT.infTime")
  rCT.sympt <- get_attr(dat, "rCT.sympt")
  uCT <- get_attr(dat, "uCT")
  uCT.infTime <- get_attr(dat, "uCT.infTime")
  uCT.sympt <- get_attr(dat, "uCT.sympt")
  CT.tx <- get_attr(dat, "CT.tx")

  diag.status.gc <- get_attr(dat, "diag.status.gc")
  diag.status.gc.prep <- get_attr(dat, "diag.status.gc.prep")
  gc.last.screen.time <- get_attr(dat, "gc.last.screen.time")
  gc.last.dxpos.time <- get_attr(dat, "gc.last.dxpos.time")

  diag.status.ct <- get_attr(dat, "diag.status.ct")
  diag.status.ct.prep <- get_attr(dat, "diag.status.ct.prep")
  ct.last.screen.time <- get_attr(dat, "ct.last.screen.time")
  ct.last.dxpos.time <- get_attr(dat, "ct.last.dxpos.time")

  prepStat <- get_attr(dat, "prepStat")

  # Parameters
  gc.sympt.tx.prob  <- get_param(dat, "gc.sympt.tx.prob")[1]
  ct.sympt.tx.prob  <- get_param(dat, "ct.sympt.tx.prob")[1]
  gc.asympt.tx.prob <- get_param(dat, "gc.asympt.tx.prob")[1]
  ct.asympt.tx.prob <- get_param(dat, "ct.asympt.tx.prob")[1]

  sti.prep.tx.prob <- get_param(dat, "sti.prep.tx.prob")



  # GC Treatment ---------------------

  ## Symptomatic treatment eligible
  idsRGC_txElig_sympt <- which(
    rGC == 1 &
    rGC.infTime == at - 1 &
    rGC.sympt == 1 &
    is.na(GC.tx)
  )
  idsUGC_txElig_sympt <- which(
    uGC == 1 &
    uGC.infTime == at - 1 &
    uGC.sympt == 1 &
    is.na(GC.tx)
  )
  idsGC_txElig_sympt <- union(idsRGC_txElig_sympt, idsUGC_txElig_sympt)
  idsGC_tx_sympt <- idsGC_txElig_sympt[runif(length(idsGC_txElig_sympt)) < gc.sympt.tx.prob]


  ## Asymptomatic treatment
  idsGC_txElig_asympt <- which(
    diag.status.gc == 1 &
    gc.last.dxpos.time == at &
    prepStat == 0
  )
  idsGC_tx_asympt <- idsGC_txElig_asympt[runif(length(idsGC_txElig_asympt)) < gc.asympt.tx.prob]


  ## PrEP treatment
  idsGC_txElig_prep <- which(
    diag.status.gc.prep == 1 &
    gc.last.dxpos.time == at &
    prepStat == 1
  )
  idsGC_tx_prep <- idsGC_txElig_prep[runif(length(idsGC_txElig_prep)) < sti.prep.tx.prob]


  ## Combining treatment
  idsGC_txElig <- unique(c(idsGC_txElig_sympt, idsGC_txElig_asympt, idsGC_txElig_prep))
  idsGC_tx <- unique(c(idsGC_tx_sympt, idsGC_tx_asympt, idsGC_tx_prep))



  # CT treatment ---------------------

  ## Symptomatic treatment eligible
  idsRCT_txElig_sympt <- which(
    rCT == 1 &
    rCT.infTime == at - 1 &
    rCT.sympt == 1 &
    is.na(CT.tx)
  )
  idsUCT_txElig_sympt <- which(
    uCT == 1 &
    uCT.infTime == at - 1 &
    uCT.sympt == 1 &
    is.na(CT.tx)
  )
  idsCT_txElig_sympt <- union(idsRCT_txElig_sympt, idsUCT_txElig_sympt)
  idsCT_tx_sympt <- idsCT_txElig_sympt[runif(length(idsCT_txElig_sympt)) < ct.sympt.tx.prob]


  ## Asymptomatic treatment
  idsCT_txElig_asympt <- which(
    diag.status.ct == 1 &
    ct.last.dxpos.time == at &
    prepStat == 0
  )
  idsCT_tx_asympt <- idsCT_txElig_asympt[runif(length(idsCT_txElig_asympt)) < ct.asympt.tx.prob]


  ## PrEP treatment
  idsCT_txElig_prep <- which(
    diag.status.ct.prep == 1 &
    ct.last.dxpos.time == at &
    prepStat == 1
  )
  idsCT_tx_prep <- idsCT_txElig_prep[runif(length(idsCT_txElig_prep)) < sti.prep.tx.prob]


  ## Combining treatment
  idsCT_txElig <- unique(c(idsCT_txElig_sympt, idsCT_txElig_asympt, idsCT_txElig_prep))
  idsCT_tx <- unique(c(idsCT_tx_sympt, idsCT_tx_asympt, idsCT_tx_prep))



  # Update attributes ---------------------
  GC.tx[idsGC_txElig] <- 0
  GC.tx[idsGC_tx] <- 1

  CT.tx[idsCT_txElig] <- 0
  CT.tx[idsCT_tx] <- 1

  # Set last screen/dxpos times for symptomatic Dx/Tx
  gc.last.screen.time[idsGC_tx_sympt] <- at
  gc.last.dxpos.time[idsGC_tx_sympt] <- at
  ct.last.screen.time[idsGC_tx_sympt] <- at
  ct.last.dxpos.time[idsCT_tx_sympt] <- at

  dat <- set_attr(dat, "GC.tx", GC.tx)
  dat <- set_attr(dat, "CT.tx", CT.tx)
  dat <- set_attr(dat, "gc.last.screen.time", gc.last.screen.time)
  dat <- set_attr(dat, "gc.last.dxpos.time", gc.last.dxpos.time)
  dat <- set_attr(dat, "ct.last.screen.time", ct.last.screen.time)
  dat <- set_attr(dat, "ct.last.dxpos.time", ct.last.dxpos.time)

  return(dat)
}
