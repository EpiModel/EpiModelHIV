
#' @title EPT Module
#'
#' @description Module function for provision  of expedited
#'              partner therapy (EPT) from index partner to non-index partner
#'              and uptake by non-index partner to prevent STI infection.
#'              Eligibility for index partner is handled in the STI treatment
#'              module and eligibility for non-index handled in risk history
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

    # Partnership list
    part.list <- dat$temp$part.list

    ## Stoppage for Index) -------------------------------------------------------

    # Index no longer eligible(> 60 days since treatment time)
    idseptExpired <- which(((at - eptindexEligdate) > ept.risk.int) & eptindexStat == 1)

    # Reset EPT status
    idsStp <- c(idseptExpired)
    eptindexStat[idsStp] <- NA
    eptindexElig[idsStp] <- NA


    ## Output -----------------------------------------------------------------

    # Attributes
    dat$attr$eptindexElig <- eptindexElig
    dat$attr$eptindexStat <- eptindexStat
    dat$attr$eptpartEligTX <- eptpartEligTx
    dat$attr$eptpartTX <- eptpartTx
    dat$attr$eptindexStartTime <- eptindexStartTime

    return(dat)
}
