
#' @title STI Recovery Module
#'
#' @description Stochastically simulates GC/CT recovery.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
stirecov_msm <- function(dat, at) {

  # Parameters ----------------------------------------------------------

  rgc.ntx.int <- dat$param$rgc.ntx.int
  ugc.ntx.int <- dat$param$ugc.ntx.int
  gc.tx.int <- dat$param$gc.tx.int

  rct.ntx.int <- dat$param$rct.ntx.int
  uct.ntx.int <- dat$param$uct.ntx.int
  ct.tx.int <- dat$param$ct.tx.int

  # GC Recovery ---------------------------------------------------------

  # Untreated (asymptomatic and symptomatic)
  idsRGC_ntx <- which(dat$attr$rGC == 1 &
                        dat$attr$rGC.infTime < at &
                        (is.na(dat$attr$rGC.tx) | dat$attr$rGC.tx == 0) &
                        (is.na(dat$attr$rGC.tx.prep) | dat$attr$rGC.tx.prep == 0))
  idsUGC_ntx <- which(dat$attr$uGC == 1 &
                        dat$attr$uGC.infTime < at &
                        (is.na(dat$attr$uGC.tx) | dat$attr$uGC.tx == 0) &
                        (is.na(dat$attr$uGC.tx.prep) | dat$attr$uGC.tx.prep == 0))

  # recovRGC_ntx <- idsRGC_ntx[which(rbinom(length(idsRGC_ntx), 1,
  #                                         1/rgc.ntx.int) == 1)]
  # recovUGC_ntx <- idsUGC_ntx[which(rbinom(length(idsUGC_ntx), 1,
  #                                         1/ugc.ntx.int) == 1)]
  recovRGC_ntx <- idsRGC_ntx[at - dat$attr$rGC.infTime[idsRGC_ntx] >= rgc.ntx.int]
  recovUGC_ntx <- idsUGC_ntx[at - dat$attr$uGC.infTime[idsUGC_ntx] >= ugc.ntx.int]


  # Treated (asymptomatic and symptomatic)
  idsRGC_tx <- which(dat$attr$rGC == 1 &
                       dat$attr$rGC.infTime < at &
                       (dat$attr$rGC.tx == 1 | dat$attr$rGC.tx.prep == 1))
  idsUGC_tx <- which(dat$attr$uGC == 1 &
                       dat$attr$uGC.infTime < at &
                       (dat$attr$uGC.tx == 1 | dat$attr$uGC.tx.prep == 1))

  # recovRGC_tx <- idsRGC_tx[which(rbinom(length(idsRGC_tx), 1,
  #                                       1/gc.tx.int) == 1)]
  # recovUGC_tx <- idsUGC_tx[which(rbinom(length(idsUGC_tx), 1,
  #                                       1/gc.tx.int) == 1)]
  recovRGC_tx <- idsRGC_tx[at - dat$attr$rGC.infTime[idsRGC_tx] >= gc.tx.int]
  recovUGC_tx <- idsUGC_tx[at - dat$attr$uGC.infTime[idsUGC_tx] >= gc.tx.int]

  recovRGC <- c(recovRGC_ntx, recovRGC_tx)
  recovUGC <- c(recovUGC_ntx, recovUGC_tx)

  dat$attr$rGC[recovRGC] <- 0
  dat$attr$rGC.sympt[recovRGC] <- NA
  dat$attr$rGC.infTime[recovRGC] <- NA
  dat$attr$rGC.tx[recovRGC] <- NA
  dat$attr$rGC.tx.prep[recovRGC] <- NA

  dat$attr$uGC[recovUGC] <- 0
  dat$attr$uGC.sympt[recovUGC] <- NA
  dat$attr$uGC.infTime[recovUGC] <- NA
  dat$attr$uGC.tx[recovUGC] <- NA
  dat$attr$uGC.tx.prep[recovUGC] <- NA



  # CT Recovery ---------------------------------------------------------

  # Untreated (asymptomatic and symptomatic)
  idsRCT_ntx <- which(dat$attr$rCT == 1 &
                        dat$attr$rCT.infTime < at &
                        (is.na(dat$attr$rCT.tx) | dat$attr$rCT.tx == 0) &
                        (is.na(dat$attr$rCT.tx.prep) | dat$attr$rCT.tx.prep == 0))
  idsUCT_ntx <- which(dat$attr$uCT == 1 &
                        dat$attr$uCT.infTime < at &
                        (is.na(dat$attr$uCT.tx) | dat$attr$uCT.tx == 0) &
                        (is.na(dat$attr$uCT.tx.prep) | dat$attr$uCT.tx.prep == 0))

  # recovRCT_ntx <- idsRCT_ntx[which(rbinom(length(idsRCT_ntx),
  #                                         1, 1/rct.ntx.int) == 1)]
  # recovUCT_ntx <- idsUCT_ntx[which(rbinom(length(idsUCT_ntx),
  #                                         1, 1/uct.ntx.int) == 1)]
  recovRCT_ntx <- idsRCT_ntx[at - dat$attr$rCT.infTime[idsRCT_ntx] >= rct.ntx.int]
  recovUCT_ntx <- idsUCT_ntx[at - dat$attr$uCT.infTime[idsUCT_ntx] >= uct.ntx.int]

  # Treated (asymptomatic and symptomatic)
  idsRCT_tx <- which(dat$attr$rCT == 1 &
                       dat$attr$rCT.infTime < at &
                       (dat$attr$rCT.tx == 1 | dat$attr$rCT.tx.prep == 1))
  idsUCT_tx <- which(dat$attr$uCT == 1 &
                       dat$attr$uCT.infTime < at &
                       (dat$attr$uCT.tx == 1 | dat$attr$uCT.tx.prep == 1))

  # recovRCT_tx <- idsRCT_tx[which(rbinom(length(idsRCT_tx),
  #                                       1, 1/ct.tx.int) == 1)]
  # recovUCT_tx <- idsUCT_tx[which(rbinom(length(idsUCT_tx),
  #                                       1, 1/ct.tx.int) == 1)]
  recovRCT_tx <- idsRCT_tx[at - dat$attr$rCT.infTime[idsRCT_tx] >= ct.tx.int]
  recovUCT_tx <- idsUCT_tx[at - dat$attr$uCT.infTime[idsUCT_tx] >= ct.tx.int]

  recovRCT <- c(recovRCT_ntx, recovRCT_tx)
  recovUCT <- c(recovUCT_ntx, recovUCT_tx)


  # Output ------------------------------------------------------------------

  dat$attr$rCT[recovRCT] <- 0
  dat$attr$rCT.sympt[recovRCT] <- NA
  dat$attr$rCT.infTime[recovRCT] <- NA
  dat$attr$rCT.tx[recovRCT] <- NA
  dat$attr$rCT.tx.prep[recovRCT] <- NA

  dat$attr$uCT[recovUCT] <- 0
  dat$attr$uCT.sympt[recovUCT] <- NA
  dat$attr$uCT.infTime[recovUCT] <- NA
  dat$attr$uCT.tx[recovUCT] <- NA
  dat$attr$uCT.tx.prep[recovUCT] <- NA

  # dat$epi$recov.rgc[at] <- length(recovRGC)
  # dat$epi$recov.ugc[at] <- length(recovUGC)
  # dat$epi$recov.rct[at] <- length(recovRCT)
  # dat$epi$recov.uct[at] <- length(recovUCT)

  return(dat)
}
