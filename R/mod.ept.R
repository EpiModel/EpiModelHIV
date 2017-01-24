
#' @title EPT Module
#'
#' @description Module function for implementation and uptake of expedited
#'              partner therapy (EPT) to prevent STI infection.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
ept_msm <- function(dat, at) {
    
    if (at < dat$param$ept.start) {
        return(dat)
    }
    
    ## Variables
    
    # Attributes
    rGC <- dat$attr$rGC
    uGC <- dat$attr$uGC
    rCT <- dat$attr$rCT
    uCT <- dat$attr$uCT
    
    # Fix this - these tx vars are NA after recovery - can only be eligible to provide medication at one time step?
    # Or, change this variable to a day and have a x-week interval where EPT could be provided?
    rGC.tx <- dat$attr$rGC.tx
    uGC.tx <- dat$attr$uGC.tx
    rCT.tx <- dat$attr$rCT.tx
    uCT.tx <- dat$attr$uCT.tx
    
    active <- dat$attr$active
    sexactive <- dat$attr$sexactive
    status <- dat$attr$status
    diag.status <- dat$attr$diag.status
    diag.status.syph <- dat$attr$diag.status.syph
    diag.status.gc <- dat$attr$diag.status.gc
    diag.status.ct <- dat$attr$diag.status.ct
    
    lnt <- dat$attr$last.neg.test
    lnt.rgc <- dat$attr$last.neg.test.rgc
    lnt.ugc <- dat$attr$last.neg.test.ugc
    lnt.rct <- dat$attr$last.neg.test.rct
    lnt.uct <- dat$attr$last.neg.test.uct
    lnt.syph <- dat$attr$last.neg.test.syph
    
    eptElig <- dat$attr$eptElig
    eptStat <- dat$attr$eptStat
    eptEligdate <- dat$attr$eptEligdate
    #eptLastRisk <- dat$attr$eptLastRisk
    #eptStartTime <- dat$attr$eptStartTime

    # Parameters
    ept.coverage <- dat$param$ept.coverage
    ept.cov.rate <- dat$param$ept.cov.rate
    
    
    ## Eligibility ---------------------------------------------------------------
    
    # Base eligibility (index partners who were just treated)
    idsEligStart <- which(active == 1 & eptEligdate >= at & eptStat == 0)
    
    # Interval - Last sexually active date: Need to be sexually active in last X days to have partners to give meds to
    ind1 <- dat$attr$ept.ind1
    
    # Draw from edge list?

    twind <- at - dat$param$ept.risk.int
    idsEligStart <- intersect(which(ind1 >= twind), idsEligStart)
    
    eptElig[idsEligStart] <- 1

    # ## EPT
    # idsept <- which((at - sexactive) <= ept.risk.int)
    # dat$attr$ept.ind1[dat$param$ept.risk.int] <- at
    
    
    ## Stoppage ------------------------------------------------------------------
    
    # Change EPT eligibility back to NA?
    
    # No indications
    #idsRiskAssess <- which(active == 1 & eptStat == 1) #& lnt == at & (at - eptLastRisk) >= 52)
    #eptLastRisk[idsRiskAssess] <- at
    # 
    # idsEligStop <- intersect(which(ind1 < twind & ind2 < twind &
    #                                   ind3 < twind & ind4 < twind),
    #                         idsRiskAssess)
    # 
    # eptElig[idsEligStop] <- NA
    # 
    # Diagnosis - should be partner's status?
    # idsStpDx <- which(active == 1 & eptStat == 1 & )
    # 
    # # Death
    #idsStpDth <- which(active == 0 & eptStat == 1)
    # 
    # # Reset EPT status
    # idsStp <- c(idsStpDx, idsStpDth, idsEligStop)
    # eptStat[idsStp] <- NA

    ## Initiation ----------------------------------------------------------------
    
    eptCov <- sum(eptStat == 1, na.rm = TRUE)/sum(eptElig == 1, na.rm = TRUE)
    eptCov <- ifelse(is.nan(eptCov), 0, eptCov)
    
    idsEligSt <- which(eptElig == 1)
    nEligSt <- length(idsEligSt)
    
    nStart <- max(0, min(nEligSt, round((ept.coverage - eptCov) *
                                            sum(eptElig == 1, na.rm = TRUE))))
     idsStart <- NULL
     if (nStart > 0) {
        if (ept.cov.rate >= 1) {
            idsStart <- ssample(idsEligSt, nStart)
        } else {
            idsStart <- idsEligSt[rbinom(nStart, 1, ept.cov.rate) == 1]
        }
     }
    
    # Attributes
     if (length(idsStart) > 0) {
         eptStat[idsStart] <- 1
         #eptStartTime[idsStart] <- at
         #eptLastRisk[idsStart] <- at
     }
    
    
    ## Output --------------------------------------------------------------------
    
    # Attributes
    dat$attr$eptElig <- eptElig
    dat$attr$eptStat <- eptStat
    #dat$attr$eptStartTime <- eptStartTime
    #dat$attr$eptLastRisk <- eptLastRisk
    
    # Summary Statistics
    dat$epi$eptCov[at] <- eptCov
    dat$epi$eptStart[at] <- length(idsStart)
    
    #intvars <- grep(names(p), pattern = ".int", fixed = TRUE)
    #p[intvars] <- lapply(p[intvars], FUN = function(x) round(x / p$time.unit))
    
    #ratevars <- grep(names(p), pattern = ".rate", fixed = TRUE)
    #p[ratevars] <- lapply(p[ratevars], FUN = function(x) x * p$time.unit)
    
    return(dat)
}
