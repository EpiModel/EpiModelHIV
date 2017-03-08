
#' @title EPT Module
#'
#' @description Module function for implementation of expedited
#'              partner therapy (EPT) in index partner to prevent 
#'              STI infection. Eligbility is handled in the risk history
#'              module.
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
    # Add a last tx time attribute and set to at
    
    active <- dat$attr$active
    
    eptElig <- dat$attr$eptElig
    eptStat <- dat$attr$eptStat
    eptEligdate <- dat$attr$eptEligdate
    eptTx <- dat$attr$eptTx
    eptStartTime <- dat$attr$eptStartTime

    # Parameters
    ept.risk.int <- dat$param$ept.risk.int
    ept.coverage <- dat$param$ept.coverage
    ept.cov.rate <- dat$param$ept.cov.rate    

    part.list <- dat$temp$part.list
    
    ## Stoppage (Index) ---------------------------------------------------------------
    
    # Index no longer eligible( > 60 days since treatment time)
    idseptExpired <- which(at - eptEligdate > ept.risk.int)
    
    # Death
    idsStpDth <- which(active == 0 & eptStat == 1)
    
    # Reset EPT status
    idsStp <- c(idseptExpired, idsStpDth)
    eptStat[idsStp] <- NA
    
    ## Initiation (index) -------------------------------------------------------------
    
    eptCov <- sum(eptStat == 1, na.rm = TRUE) / sum(eptElig == 1, na.rm = TRUE)
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
    
    # Update attributes of index
    if (length(idsStart) > 0) {
        eptStat[idsStart] <- 1
        eptStartTime[idsStart] <- at
    }
    
    ## Output --------------------------------------------------------------------
    
    # Attributes
    dat$attr$eptElig <- eptElig
    dat$attr$eptStat <- eptStat
    dat$attr$eptTX <- eptTx
    dat$attr$eptStartTime <- eptStartTime
    
    # Summary Statistics
    dat$epi$eptCov[at] <- eptCov
    dat$epi$eptStart[at] <- length(idsStart)

    return(dat)
}
