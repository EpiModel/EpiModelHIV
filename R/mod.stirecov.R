
#' @title STI Recovery Module
#'
#' @description Stochastically simulates GC/CT recovery.
#'
#' @inheritParams aging_msm
#'
#' @export
#'
stirecov_msm <- function(dat, at) {

  ## Input
  # Attributes
  rGC <- get_attr(dat, "rGC")
  rGC.infTime <- get_attr(dat, "rGC.infTime")
  uGC <- get_attr(dat, "uGC")
  uGC.infTime <- get_attr(dat, "uGC.infTime")
  GC.tx <- get_attr(dat, "GC.tx")

  rCT <- get_attr(dat, "rCT")
  rCT.infTime <- get_attr(dat, "rCT.infTime")
  uCT <- get_attr(dat, "uCT")
  uCT.infTime <- get_attr(dat, "uCT.infTime")
  CT.tx <- get_attr(dat, "CT.tx")


  # Parameters
  gc.tx.rate   <- 1 / get_param(dat, "gc.tx.int")
  gc.ntx.rate <- 1 / get_param(dat, "gc.ntx.int")
  ct.tx.rate   <- 1 / get_param(dat, "ct.tx.int")
  ct.ntx.rate <- 1 / get_param(dat, "ct.ntx.int")


  # GC Recovery ---------------------------------------------------------

  # Untreated (asymptomatic and symptomatic)
  idsRGC_recovElig_ntx <- which(
    rGC == 1 &
    rGC.infTime < at &
    (is.na(GC.tx) | GC.tx == 0)
  )
  idsUGC_recovElig_ntx <- which(
    uGC == 1 &
    uGC.infTime < at &
    (is.na(GC.tx) | GC.tx == 0)
  )
  idsGC_recovElig_ntx <- union(idsRGC_recovElig_ntx, idsUGC_recovElig_ntx)
  idsGC_recov_ntx <- idsGC_recovElig_ntx[runif(length(idsGC_recovElig_ntx)) < gc.ntx.rate]


  # Treated (asymptomatic and symptomatic)
  idsGC_recovElig_tx <- which(GC.tx == 1)
  idsGC_recov_tx <- idsGC_recovElig_tx[runif(length(idsGC_recovElig_tx)) < gc.tx.rate]
  idsGC_recov <- union(idsGC_recov_ntx, idsGC_recov_tx)


  # Update attr
  dat <- set_attr(dat, "rGC", 0, idsGC_recov)
  dat <- set_attr(dat, "rGC.sympt", NA, idsGC_recov)
  dat <- set_attr(dat, "rGC.infTime", NA, idsGC_recov)
  dat <- set_attr(dat, "uGC", 0, idsGC_recov)
  dat <- set_attr(dat, "uGC.sympt", NA, idsGC_recov)
  dat <- set_attr(dat, "uGC.infTime", NA, idsGC_recov)
  dat <- set_attr(dat, "GC.tx", NA, idsGC_recov)
  dat <- set_attr(dat, "diag.status.gc", NA, idsGC_recov)
  dat <- set_attr(dat, "diag.status.gc.prep", NA, idsGC_recov)


  # CT Recovery ---------------------------------------------------------

  # Untreated (asymptomatic and symptomatic)
  idsRCT_recovElig_ntx <- which(
    rCT == 1 &
    rCT.infTime < at &
    (is.na(CT.tx) | CT.tx == 0)
  )
  idsUCT_recovElig_ntx <- which(
    uCT == 1 &
    uCT.infTime < at &
    (is.na(CT.tx) | CT.tx == 0)
  )
  idsCT_recovElig_ntx <- union(idsRCT_recovElig_ntx, idsUCT_recovElig_ntx)
  idsCT_recov_ntx <- idsCT_recovElig_ntx[runif(length(idsCT_recovElig_ntx)) < ct.ntx.rate]


  # Treated (asymptomatic and symptomatic)
  idsCT_recovElig_tx <- which(CT.tx == 1)
  idsCT_recov_tx <- idsCT_recovElig_tx[runif(length(idsCT_recovElig_tx)) < ct.tx.rate]
  idsCT_recov <- union(idsCT_recov_ntx, idsCT_recov_tx)


  # Update attr
  dat <- set_attr(dat, "rCT", 0, idsCT_recov)
  dat <- set_attr(dat, "rCT.sympt", NA, idsCT_recov)
  dat <- set_attr(dat, "rCT.infTime", NA, idsCT_recov)
  dat <- set_attr(dat, "uCT", 0, idsCT_recov)
  dat <- set_attr(dat, "uCT.sympt", NA, idsCT_recov)
  dat <- set_attr(dat, "uCT.infTime", NA, idsCT_recov)
  dat <- set_attr(dat, "CT.tx", NA, idsCT_recov)
  dat <- set_attr(dat, "diag.status.ct", NA, idsCT_recov)
  dat <- set_attr(dat, "diag.status.ct.prep", NA, idsCT_recov)

  return(dat)
}
