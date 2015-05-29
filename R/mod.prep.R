
#' @export
prep <- function(dat, at) {

  dat <- prep.elig(dat, at)
  dat <- prep.update(dat, at)

  return(dat)
}


#' @export
prep.elig <- function(dat, at) {
  # if (at >= 5) browser()
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
  if (prep.pop == "gsdcoup" & at >= prep.start) {
    del <- as.data.frame(discord_edgelist.hiv(dat, at))
    dxStatInf <- dat$attr$dxStat[del$inf]
    idsEligSt <- del$sus[which(dxStatInf == 1)]
    if (at > 2) {
      # Removes ids who are already on PreP
      idsEligSt <- idsEligSt[-which(dat$attr$prepStat[idsEligSt] == 1)]
      # Removes ids who have cycled off Prep (non-adherent)
      idsEligSt <- idsEligSt[-which(!is.na(dat$attr$prepStartTime[idsEligSt]))]
    }
    if (length(idsEligSt) > 0) {
      vecEligSt <- which(rbinom(length(idsEligSt), 1, prep.rec.rate) == 1)
    }
  }

  if (prep.pop == "gwomen" & at >= prep.start) {
    idsEligSt <- which(active == 1 & status == "s" & male == 0 & is.na(prepElig))
    if (length(idsEligSt) > 0) {
      vecEligSt <- which(rbinom(length(idsEligSt), 1, prep.rec.rate) == 1)
    }
  }

  if (prep.pop == "rctsdcoup" & at >= prep.start) {
    rctStat <- dat$param$rctStat
    del <- as.data.frame(discord_edgelist.hiv(dat, at))
    idsEligSt <- which()
  }

  if (prep.pop == "rctwomen" & at >= prep.start) {
    rctStat <- dat$attr$rctStat
    rctTx <- dat$attr$rctTx
    idsEligSt <- which(active == 1 & rctStat == 1 & rctTx == 1 & is.na(prepElig))
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

  return(dat)
}


#' @export
prep.update <- function(dat, at) {
  # if (at >= 5) browser()
  # Variables ---------------------------------------------------------------
  active <- dat$attr$active
  n <- length(active)

  prepElig <- dat$attr$prepElig

  ## Pull or set attr
  if (is.null(dat$attr$prepStat)) dat$attr$prepStat <- rep(NA, n)
  prepStat <- dat$attr$prepStat

  if (is.null(dat$attr$prepStartTime)) dat$attr$prepStartTime <- rep(NA, n)
  prepStartTime <- dat$attr$prepStartTime

  if (is.null(dat$attr$prepStopTime)) dat$attr$prepStopTime <- rep(NA, n)
  prepStopTime <- dat$attr$prepStopTime

  if (is.null(dat$attr$prepStopReas)) dat$attr$prepStopReas <- rep(NA, n)
  prepStopReas <- dat$attr$prepStopReas

  if (is.null(dat$attr$prepType)) dat$attr$prepType <- rep(NA, n)
  prepType <- dat$attr$prepType

  if (is.null(dat$attr$prepStops)) dat$attr$prepStops <- rep(NA, n)
  prepStops <- dat$attr$prepStops

  if (is.null(dat$attr$prepTimeOn)) dat$attr$prepTimeOn <- rep(NA, n)
  prepTimeOn <- dat$attr$prepTimeOn

  if (is.null(dat$attr$prepTimeOff)) dat$attr$prepTimeOff <- rep(NA, n)
  prepTimeOff <- dat$attr$prepTimeOff

  prep.coverage <- dat$param$prep.coverage
  prep.adhere.full <- dat$param$prep.adhere.full
  prep.adhere.part <- dat$param$prep.adhere.part


  # Process -----------------------------------------------------------------

  ## Start PreP

  # Calculate prep coverage
  prepCov <- sum(!is.na(prepStartTime[which(prepElig == 1)]))/sum(prepElig == 1, na.rm = TRUE)
  prepCov <- ifelse(is.nan(prepCov), 0, prepCov)

  idsEligSt <- which(active == 1 & prepElig == 1 & is.na(prepStartTime))
  nEligSt <- length(idsEligSt)
  idsStart <- NULL

  ## Start PreP
  nStart <- max(0, min(nEligSt, round((prep.coverage - prepCov) * sum(prepElig == 1, na.rm = TRUE))))
  if (nStart > 0) {
    idsStart <- ssample(idsEligSt, nStart)
  }

  ## Prep adherence type assignment
  if (length(idsStart) > 0) {
    needprepType <- which(is.na(prepType[idsStart]))
    if (length(needprepType) > 0) {
      prepType[idsStart[needprepType]] <- rbinom(length(needprepType), 1, prep.adhere.full)
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

  ## Stop Prep (sus becomes infected)
  idsStopConv <- which(active == 1 & dat$attr$prepStat == 1 & dat$attr$status == "i")
  if (length(idsStopConv) > 0) {
    prepStat[idsStopConv] <- 0
    prepStopTime[idsStopConv] <- at
    prepStopReas[idsStopConv] <- "sc"
    prepElig[idsStopConv] <- 0
  }

  ## Stop Prep (serodis couple dissolved)
  idsStopDiss <- NULL
  if (dat$param$prep.pop %in% c("gsdcoup", "rctsdcoup")) {
    del <- as.data.frame(discord_edgelist.hiv(dat, at))
    idsOn <- which(dat$attr$prepStat == 1)
    idsStopDiss <- setdiff(idsOn, del$sus)
    if (length(idsStopDiss) > 0) {
      prepStat[idsStopDiss] <- 0
      prepStopTime[idsStopDiss] <- at
      prepStopReas[idsStopDiss] <- "ds"
      prepElig[idsStopDiss] <- 0
      prepStartTime[idsStopDiss] <- NA
    }
  }

  ## Stop Prep (non-adherence)
  idsStopNad <- NULL
  idsEligStopNad <- which(active == 1 & dat$attr$prepStat == 1 & prepType == 0)
  nEligStopNad <- length(idsEligStopNad)
  if (nEligStopNad > 0) {
    vecStopNad <- which(rbinom(nEligStopNad, 1, (1 - prep.adhere.part)) == 1)
    if (length(vecStopNad) > 0) {
      idsStopNad <- idsEligStopNad[vecStopNad]
      prepStat[idsStopNad] <- 0
      prepStops[idsStopNad] <- prepStops[idsStopNad] + 1
      prepStopTime[idsStopNad] <- at
      prepStopReas[idsStopNad] <- "nad"
    }
  }


  ## Restart Prep (resume adherence)
  idsRest <- NULL
  idsEligRest <- which(active == 1 & dat$attr$prepStat == 0 & prepStopReas == "nad")
  nEligRest <- length(idsEligRest)
  if (nEligRest > 0) {
    vecRest <- which(rbinom(nEligRest, 1, prep.adhere.part) == 1)
    if (length(vecRest) > 0) {
      idsRest <- idsEligRest[vecRest]
      prepStat[idsRest] <- 1
    }
  }


  # Output ------------------------------------------------------------------

  idsOnPrep <- which(prepStat == 1)
  idsOffPrep <- which(prepStat == 0 & !is.na(prepStartTime))
  prepTimeOn[idsOnPrep] <- prepTimeOn[idsOnPrep] + 1
  prepTimeOff[idsOffPrep] <- prepTimeOff[idsOffPrep] + 1

  dat$attr$prepStat <- prepStat
  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepStopTime <- prepStopTime
  dat$attr$prepStops <- prepStops
  dat$attr$prepTimeOn <- prepTimeOn
  dat$attr$prepTimeOff <- prepTimeOff
  dat$attr$prepType <- prepType
  dat$attr$prepElig <- prepElig
  dat$attr$prepStopReas <- prepStopReas

  if (at == 2) {
    dat$epi$prepCov <- c(NA, prepCov)
    dat$epi$prepStart <- c(0, length(idsStart))
    dat$epi$prepStopNad <- c(0, length(idsStopNad))
    dat$epi$prepStopConv <- c(0, length(idsStopConv))
    dat$epi$prepStopDiss <- c(0, length(idsStopDiss))
    dat$epi$prepAdhereNum <- c(0, sum(prepStat == 1, na.rm = TRUE))
    dat$epi$prepAdherePct <- c(0, sum(active == 1 & prepStat == 1, na.rm = TRUE) /
                                 sum(active == 1 & !is.na(prepStartTime) &
                                       (is.na(prepStopReas) | prepStopReas == "nad"),
                                     na.rm = TRUE))
  } else {
    dat$epi$prepCov[at] <- prepCov
    dat$epi$prepStart[at] <- length(idsStart)
    dat$epi$prepStopNad[at] <- length(idsStopNad)
    dat$epi$prepStopConv[at] <- length(idsStopConv)
    dat$epi$prepStopDiss[at] <- length(idsStopDiss)
    dat$epi$prepAdhereNum[at] <- sum(prepStat == 1, na.rm = TRUE)
    dat$epi$prepAdherePct[at] <- sum(active == 1 & prepStat == 1, na.rm = TRUE) /
      sum(active == 1 & !is.na(prepStartTime) &
            (is.na(prepStopReas) | prepStopReas == "nad"),
          na.rm = TRUE)
  }

  return(dat)
}
