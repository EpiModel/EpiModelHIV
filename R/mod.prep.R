
#' @export
prep.mard <- function(dat, at) {

  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  male <- dat$attr$male
  status <- dat$attr$status

  prep.start <- dat$param$prep.start
  prep.rec.rate <- dat$param$prep.rec.rate
  prep.pop <- dat$param$prep.pop

  if (is.null(dat$attr$prepElig)) dat$attr$prepElig <- rep(NA, length(active))
  prepElig <- dat$attr$prepElig

  if (is.null(dat$attr$prepEligTime)) dat$attr$prepEligTime <- rep(NA, length(active))
  prepEligTime <- dat$attr$prepEligTime


  # Recruitment Models ------------------------------------------------------

  vecEligSt <- NULL

  if (prep.pop == "gwomen" & at >= prep.start) {
    idsEligSt <- which(active == 1 & status == "s" & male == 0 & is.na(prepElig))
    if (length(idsEligSt) > 0) {
      vecEligSt <- which(rbinom(length(idsEligSt), 1, prep.rec.rate) == 1)
    }
  }


  # Update attributes -------------------------------------------------------

  if (length(vecEligSt) > 0) {
    idsStart <- idsEligSt[vecEligSt]
    prepElig[idsStart] <- 1
    prepEligTime[idsStart] <- at
  }

  dat$attr$prepElig <- prepElig
  dat$attr$prepEligTime <- prepEligTime

  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  n <- length(active)

  prepElig <- dat$attr$prepElig

  prep.coverage <- dat$param$prep.coverage
  prep.adhere.full <- dat$param$prep.adhere.full
  prep.adhere.part <- dat$param$prep.adhere.part


  # Process -----------------------------------------------------------------

  ## Start PreP

  # Calculate prep coverage
  prepCov <- sum(!is.na(prepStartTime[which(prepElig == 1)]))/
    sum(prepElig == 1, na.rm = TRUE)
  prepCov <- ifelse(is.nan(prepCov), 0, prepCov)

  idsEligSt <- which(active == 1 & prepElig == 1 & is.na(prepStartTime))
  nEligSt <- length(idsEligSt)
  idsStart <- NULL

  ## Start PreP
  nStart <- max(0, min(nEligSt, round((prep.coverage - prepCov) *
                                        sum(prepElig == 1, na.rm = TRUE))))
  if (nStart > 0) {
    idsStart <- ssample(idsEligSt, nStart)
  }

  ## Prep adherence type assignment
  if (length(idsStart) > 0) {
    needprepType <- which(is.na(prepType[idsStart]))
    if (length(needprepType) > 0) {
      prepType[idsStart[needprepType]] <- rbinom(length(needprepType), 1,
                                                 prep.adhere.full)
    }
    if (prep.adhere.part == 0) {
      idsStart <- intersect(idsStart, which(prepType == 1))
    }
  }

  ## Update starting attributes
  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1
    prepStartTime[idsStart] <- at
    prepStops[idsStart] <- 0
    prepTimeOn[idsStart] <- 0
    prepTimeOff[idsStart] <- 0
  }




  # Output ------------------------------------------------------------------

  idsOnPrep <- which(prepStat == 1)
  idsOffPrep <- which(prepStat == 0 & !is.na(prepStartTime))
  prepTimeOn[idsOnPrep] <- prepTimeOn[idsOnPrep] + 1
  prepTimeOff[idsOffPrep] <- prepTimeOff[idsOffPrep] + 1

  dat$attr$prepStat <- prepStat
  dat$attr$prepType <- prepType
  dat$attr$prepElig <- prepElig

  dat$epi$prepCov[at] <- prepCov
  dat$epi$prepStart[at] <- length(idsStart)
  dat$epi$prepAdhereNum[at] <- sum(prepStat == 1, na.rm = TRUE)


  return(dat)
}
