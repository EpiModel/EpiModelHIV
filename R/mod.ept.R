
#' @title EPT Module
#'
#' @description Module function for implementation of expedited
#'              partner therapy (EPT) in index partner to prevent 
#'              STI infection. Eligibility is handled in the risk history
#'              module.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
sti_ept_msm <- function(dat, at) {
    
    if (at < dat$param$ept.start) {
        return(dat)
    }
    
    ## Variables
  
#Need to convert uid to regular id?
    
    # Attributes
    rGC <- dat$attr$rGC
    uGC <- dat$attr$uGC
    rCT <- dat$attr$rCT
    uCT <- dat$attr$uCT
    
    # Fix this - these tx vars are NA after recovery - can only be eligible to 
    # provide medication at one time step?
    # Or, change this variable to a day and have a x-week interval where EPT 
    # could be provided?
    # Add a last tx time attribute and set to at
    
    eptindexElig <- dat$attr$eptindexElig
    eptindexStat <- dat$attr$eptindexStat
    eptindexEligdate <- dat$attr$eptindexEligdate
    eptpartEligTx <- dat$attr$eptpartEligTx
    eptpartTx <- dat$attr$eptTx
    eptindexStartTime <- dat$attr$eptindexStartTime

    # Parameters
    ept.risk.int <- dat$param$ept.risk.int
    ept.coverage <- dat$param$ept.coverage
    ept.cov.rate <- dat$param$ept.cov.rate    

    part.list <- dat$temp$part.list
    
    ## Stoppage (Index) -------------------------------------------------------
    
    # Index no longer eligible(> 60 days since treatment time)
    idseptExpired <- which(((at - eptindexEligdate) > ept.risk.int) & eptindexStat == 1)
    
    # Reset EPT status
    idsStp <- c(idseptExpired)
    eptindexStat[idsStp] <- NA
    eptindexElig[idsStp] <- NA
    
    ## Initiation (index) -----------------------------------------------------
    
    eptCov <- sum(eptindexStat == 1, na.rm = TRUE) / sum(eptindexElig == 1, na.rm = TRUE)
    eptCov <- ifelse(is.nan(eptCov), 0, eptCov)
    
    idsEligSt <- which(eptindexElig == 1)
    nEligSt <- length(idsEligSt)
    
    nStart <- max(0, min(nEligSt, round((ept.coverage - eptCov) * sum(eptindexElig == 1, na.rm = TRUE))))
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
            eptindexStat[idsStart] <- 1
            eptindexStartTime[idsStart] <- at
    }
    
    ## Output -----------------------------------------------------------------
    
    # Attributes
    dat$attr$eptindexElig <- eptindexElig
    dat$attr$eptindexStat <- eptindexStat
    dat$attr$eptpartEligTX <- eptpartEligTx
    dat$attr$eptpartTX <- eptpartTx
    dat$attr$eptindexStartTime <- eptindexStartTime
    
    # Summary Statistics
    dat$epi$eptCov[at] <- eptCov
    dat$epi$eptindexStart[at] <- length(idsStart)

    return(dat)
}
