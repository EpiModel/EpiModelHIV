
#' @title PrEP Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
prep_msm <- function(dat, at) {

  if (at < dat$param$prep.start) {

      # Update # of PrEP asymptomatic STI tests to 0
      dat$epi$rGCasympttests.prep[at] <- 0
      dat$epi$uGCasympttests.prep[at] <- 0
      dat$epi$GCasympttests.prep[at] <- 0
      dat$epi$rGCasympttests.pos.prep[at] <- 0
      dat$epi$uGCasympttests.pos.prep[at] <- 0
      dat$epi$GCasympttests.pos.prep[at] <- 0

      dat$epi$rCTasympttests.prep[at] <- 0
      dat$epi$uCTasympttests.prep[at] <- 0
      dat$epi$CTasympttests.prep[at] <- 0
      dat$epi$rCTasympttests.pos.prep[at] <- 0
      dat$epi$uCTasympttests.pos.prep[at] <- 0
      dat$epi$CTasympttests.pos.prep[at] <- 0

      dat$epi$syphasympttests.prep[at] <- 0
      dat$epi$syphasympttests.pos.prep[at] <- 0

    return(dat)
  }

  ## Variables

  # Attributes
  active <- dat$attr$active
  race <- dat$attr$race
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  diag.status.syph <- dat$attr$diag.status.syph
  diag.status.gc <- dat$attr$diag.status.gc
  diag.status.ct <- dat$attr$diag.status.ct
  tst.rect.sti.rr <- dat$param$tst.rect.sti.rr
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT
  syphilis <- dat$attr$syphilis
  stage.syph <- dat$attr$stage.syph
  lnt <- dat$attr$last.neg.test
  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  prepLastRisk <- dat$attr$prepLastRisk
  prepStartTime <- dat$attr$prepStartTime
  prepLastStiScreen <- dat$attr$prepLastStiScreen
  prep.tst.int <- dat$param$prep.tst.int
  time.on.prep <- dat$attr$time.on.prep


  # Parameters

  prep.coverage <- dat$param$prep.coverage
  prep.cov.rate <- dat$param$prep.cov.rate
  prep.class.prob <- dat$param$prep.class.prob

  if (at == dat$param$prep.start) {
      dat$attr$time.hivneg[status == 0] <- 0
      dat$attr$time.off.prep[status == 0] <- 0
      dat$attr$stage.time[status == 1] <- 0
      dat$attr$stage.time.ar.ndx[status == 1] <- 0
      dat$attr$stage.time.ar.dx[status == 1] <- 0
      dat$attr$stage.time.af.ndx[status == 1] <- 0
      dat$attr$stage.time.af.dx[status == 1] <- 0
      dat$epi$stage.time.early.chronic.ndx[status == 1] <- 0
      dat$epi$stage.time.early.chronic.dx.yrone[status == 1] <- 0
      dat$epi$stage.time.early.chronic.dx.yrstwotolate[status == 1] <- 0
      dat$epi$stage.time.early.chronic.art[status == 1] <- 0
      dat$epi$stage.time.late.chronic.ndx[status == 1] <- 0
      dat$epi$stage.time.late.chronic.dx[status == 1] <- 0
      dat$epi$stage.time.late.chronic.art[status == 1] <- 0
      dat$attr$stage.time.aids.ndx[status == 1] <- 0
      dat$attr$stage.time.aids.dx[status == 1] <- 0
      dat$attr$stage.time.aids.art[status == 1] <- 0

  }

  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(status == 0 & prepStat == 0 & lnt == at)

  # Core eligiblity
  ind1 <- dat$attr$prep.ind.uai.mono
  ind2 <- dat$attr$prep.ind.uai.nmain
  ind3 <- dat$attr$prep.ind.ai.sd
  ind4 <- dat$attr$prep.ind.sti

  twind <- at - dat$param$prep.risk.int
  idsEligStart <- intersect(which(ind1 >= twind | ind2 >= twind |
                                    ind3 >= twind | ind4 >= twind),
                            idsEligStart)

  prepElig[idsEligStart] <- 1


  ## Stoppage -----------------------------------------------------------------

  # No indications
  idsRiskAssess <- which(prepStat == 1 & lnt == at & (at - prepLastRisk) >= 52)
  prepLastRisk[idsRiskAssess] <- at

  idsEligStop <- intersect(which(ind1 < twind & ind2 < twind &
                                   ind3 < twind & ind4 < twind),
                           idsRiskAssess)

  prepElig[idsEligStop] <- 0

  # Diagnosis
  idsStpDx <- which(prepStat == 1 & diag.status == 1)

  # Reset PrEP status
  idsStp <- c(idsStpDx, idsEligStop)
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  prepLastStiScreen[idsStp] <- NA

  # Update time on PrEP after people stop
  time.on.prep[prepStat == 1] <- time.on.prep[prepStat == 1] + 1


  ## Initiation ---------------------------------------------------------------

  prepCov <- sum(prepStat == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  prepCov <- ifelse(is.nan(prepCov), 0, prepCov)

  idsEligSt <- which(prepElig == 1)
  nEligSt <- length(idsEligSt)

  nStart <- max(0, min(nEligSt, round((prep.coverage - prepCov) *
                                        sum(prepElig == 1, na.rm = TRUE))))
  idsStart <- NULL
  if (nStart > 0) {
    if (prep.cov.rate >= 1) {
      idsStart <- ssample(idsEligSt, nStart)
    } else {
      idsStart <- idsEligSt[rbinom(nStart, 1, prep.cov.rate) == 1]
    }
  }

  # Attributes
  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1
    prepStartTime[idsStart] <- at
    prepLastRisk[idsStart] <- at

    # PrEP class
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
  }

  # Update time off PrEP for those not starting PrEP (housed in PrEP module starting at PrEP start time)
  if (at > dat$param$prep.start) {
      dat$attr$time.off.prep[prepStat == 0] <-
          dat$attr$time.off.prep[prepStat == 0] + 1
  }

  ## STI Testing on PrEP ------------------------------------------------------

  if (is.null(dat$epi$num.asympt.tx)) {
      dat$epi$rGCasympttests.prep <- rep(0, length(dat$control$nsteps))
      dat$epi$uGCasympttests.prep <- rep(0, length(dat$control$nsteps))
      dat$epi$GCasympttests.prep <- rep(0, length(dat$control$nsteps))
      dat$epi$rCTasympttests.prep <- rep(0, length(dat$control$nsteps))
      dat$epi$uCTasympttests.prep <- rep(0, length(dat$control$nsteps))
      dat$epi$CTasympttests.prep <- rep(0, length(dat$control$nsteps))
      dat$epi$syphasympttests.prep <- rep(0, length(dat$control$nsteps))
  }

    ## Testing
  tsincelntst.syph <- at - dat$attr$last.neg.test.syph
  tsincelntst.syph[is.na(tsincelntst.syph)] <- at - dat$attr$arrival.time[is.na(tsincelntst.syph)]

  tsincelntst.rgc <- at - dat$attr$last.neg.test.rgc
  tsincelntst.ugc <- at - dat$attr$last.neg.test.ugc
  tsincelntst.rgc[is.na(tsincelntst.rgc)] <- at - dat$attr$arrival.time[is.na(tsincelntst.rgc)]
  tsincelntst.ugc[is.na(tsincelntst.ugc)] <- at - dat$attr$arrival.time[is.na(tsincelntst.ugc)]
  tsincelntst.gc <- min(tsincelntst.rgc, tsincelntst.ugc)

  tsincelntst.rct <- at - dat$attr$last.neg.test.rct
  tsincelntst.uct <- at - dat$attr$last.neg.test.uct
  tsincelntst.rct[is.na(tsincelntst.rct)] <- at - dat$attr$arrival.time[is.na(tsincelntst.rct)]
  tsincelntst.uct[is.na(tsincelntst.uct)] <- at - dat$attr$arrival.time[is.na(tsincelntst.uct)]
  tsincelntst.ct <- min(tsincelntst.rct, tsincelntst.uct)

  # PrEP STI testing
  tst.syph.prep <- which((diag.status.syph == 0 | is.na(diag.status.syph)) &
                             prepStat == 1 &
                             tsincelntst.syph >= prep.tst.int)
  tst.gc.prep <- which((diag.status.gc == 0 | is.na(diag.status.gc)) &
                           prepStat == 1 &
                           tsincelntst.gc >= prep.tst.int)
  tst.ct.prep <- which((diag.status.ct == 0 | is.na(diag.status.ct)) &
                           prepStat == 1 &
                           tsincelntst.ct >= prep.tst.int)

  # Syphilis PrEP testing
  tst.syph.pos <- tst.syph.prep[syphilis[tst.syph.prep] == 1 & stage.syph[tst.syph.prep] %in% c(2, 3, 4, 5, 6, 7)]
  tst.syph.neg <- setdiff(tst.syph.prep, tst.syph.pos)

  # GC PrEP testing
  tst.rgc <- tst.gc.prep[dat$attr$role.class %in% c("R", "V")]
  tst.rgc <- sample(tst.rgc, tst.rect.sti.rr * length(tst.rgc))
  tst.ugc <- tst.gc.prep[dat$attr$role.class %in% c("I", "V")]
  tst.rgc.pos <- tst.rgc[rGC == 1]
  tst.ugc.pos <- tst.ugc[uGC == 1]
  tst.rgc.neg <- setdiff(tst.rgc, tst.rgc.pos)
  tst.ugc.neg <- setdiff(tst.ugc, tst.ugc.pos)
  tst.gc.pos <- unique(c(tst.rgc.pos, tst.ugc.pos))

  # CT PrEP testing
  tst.rct <- tst.ct.prep[dat$attr$role.class %in% c("R", "V")]
  tst.rct <- sample(tst.rct, tst.rect.sti.rr * length(tst.rct))
  tst.uct <- tst.ct.prep[dat$attr$role.class %in% c("I", "V")]
  tst.rct.pos <- tst.rct[rCT == 1]
  tst.uct.pos <- tst.uct[uCT == 1]
  tst.rct.neg <- setdiff(tst.rct, tst.rct.pos)
  tst.uct.neg <- setdiff(tst.uct, tst.uct.pos)
  tst.ct.pos <- unique(c(tst.rct.pos, tst.uct.pos))

  # Syphilis Attributes
  dat$attr$last.neg.test.syph[tst.syph.neg] <- at
  dat$attr$last.neg.test.syph[tst.syph.pos] <- NA
  dat$attr$diag.status.syph[tst.syph.pos] <- 1
  dat$attr$lastdiag.time.syph[tst.syph.pos] <- at

  # GC Attributes
  dat$attr$last.neg.test.rgc[tst.rgc.neg] <- at
  dat$attr$last.neg.test.ugc[tst.ugc.neg] <- at
  dat$attr$last.neg.test.rgc[tst.rgc.pos] <- NA
  dat$attr$last.neg.test.ugc[tst.ugc.pos] <- NA
  dat$attr$diag.status.gc[tst.gc.pos] <- 1
  dat$attr$lastdiag.time.gc[tst.gc.pos] <- at

  # CT Attributes
  dat$attr$last.neg.test.rct[tst.rct.neg] <- at
  dat$attr$last.neg.test.uct[tst.uct.neg] <- at
  dat$attr$last.neg.test.rct[tst.rct.pos] <- NA
  dat$attr$last.neg.test.uct[tst.uct.pos] <- NA
  dat$attr$diag.status.ct[tst.ct.pos] <- 1
  dat$attr$lastdiag.time.ct[tst.ct.pos] <- at

  # Count number of tests due to PrEP
  dat$epi$rGCasympttests.prep[at] <- length(tst.rgc)
  dat$epi$uGCasympttests.prep[at] <- length(tst.ugc)
  dat$epi$GCasympttests.prep[at] <- length(c(tst.rgc, tst.ugc))

  dat$epi$rGCasympttests.pos.prep[at] <- length(tst.rgc)
  dat$epi$uGCasympttests.pos.prep[at] <- length(tst.ugc)
  dat$epi$GCasympttests.pos.prep[at] <- length(c(tst.rgc.pos, tst.ugc.pos))

  dat$epi$rCTasympttests.prep[at] <- length(tst.rct)
  dat$epi$uCTasympttests.prep[at] <- length(tst.uct)
  dat$epi$CTasympttests.prep[at] <- length(c(tst.rct, tst.uct))

  dat$epi$rCTasympttests.pos.prep[at] <- length(tst.rct)
  dat$epi$uCTasympttests.pos.prep[at] <- length(tst.uct)
  dat$epi$CTasympttests.pos.prep[at] <- length(c(tst.rct.pos, tst.uct.pos))

  dat$epi$syphasympttests.prep[at] <- length(c(tst.syph.prep))
  dat$epi$syphasympttests.pos.prep[at] <- length(c(tst.syph.pos))

  ## Output -------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepClass <- prepClass
  dat$attr$prepLastRisk <- prepLastRisk
  dat$attr$prepLastStiScreen <- prepLastStiScreen
  dat$attr$time.on.prep <- time.on.prep

  # Summary Statistics
  dat$epi$prepCov[at] <- prepCov
  dat$epi$prepStart[at] <- length(idsStart)




  return(dat)
}
