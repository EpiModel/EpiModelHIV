
#' @title EPT Module
#'
#' @description Module function for eligibility of non-index partner, provision
#'              of expedited partner therapy (EPT) from index partner to
#'              non-index partner, and uptake by non-index partner to prevent
#'              STI infection. Eligibility for index partner is handled in the
#'              STI treatment module.
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

    ## Variables ---------------------------------------------------------------

    # Attributes
    rGC <- dat$attr$rGC
    uGC <- dat$attr$uGC
    rCT <- dat$attr$rCT
    uCT <- dat$attr$uCT

    ## Attributes
    #uid <- dat$attr$uid
    rGC.tx <- dat$attr$rGC.tx
    uGC.tx <- dat$attr$uGC.tx
    rCT.tx <- dat$attr$rCT.tx
    uCT.tx <- dat$attr$uCT.tx
    rGC.tx.prep <- dat$attr$rGC.tx.prep
    uGC.tx.prep <- dat$attr$uGC.tx.prep
    rCT.tx.prep <- dat$attr$rCT.tx.prep
    uCT.tx.prep <- dat$attr$uCT.tx.prep
    rGC.tx.ept <- dat$attr$rGC.tx.ept
    uGC.tx.ept <- dat$attr$uGC.tx.ept
    rCT.tx.ept <- dat$attr$rCT.tx.ept
    uCT.tx.ept <- dat$attr$uCT.tx.ept
    eptindexElig <- dat$attr$eptindexElig
    eptindexStat <- dat$attr$eptindexStat
    eptindexEligdate <- dat$attr$eptindexEligdate

    ## Parameters
    ept.risk.int <- dat$param$ept.risk.int
    ept.provision.main.ong <- dat$param$ept.provision.partner.main.ong
    ept.provision.pers.ong <- dat$param$ept.provision.partner.pers.ong
    ept.provision.main.end <- dat$param$ept.provision.partner.main.end
    ept.provision.pers.end <- dat$param$ept.provision.partner.pers.end
    ept.provision.inst <- dat$param$ept.provision.partner.inst
    ept.uptake.main <- dat$param$ept.uptake.partner.main
    ept.uptake.pers <- dat$param$ept.uptake.partner.pers
    ept.uptake.inst <- dat$param$ept.uptake.partner.inst

    # Partnership list
    part.list <- dat$temp$part.list

    ## Stoppage for Index ------------------------------------------------------

    # Index no longer eligible(> 1 time step since treatment time)
    idseptExpired <- which((at - eptindexEligdate) > 1)

    # Reset EPT status
    idsStp <- c(idseptExpired)
    eptindexStat[idsStp] <- NA
    eptindexElig[idsStp] <- NA


    ## Indications for non-index-------------------------------------------------

    ## Eligibility of partners
    part.list <- dat$temp$part.list

    # Subset partner list to partnerships active within an EPT interval - last active date within risk interval
    part.list <- part.list[which((at - (part.list[, "last.active.time"]) <= ept.risk.int)), , drop = FALSE]

    # Subset partner list to where both partners are alive (a dead index can't provide EPT to alive non-index)
    part.list <- part.list[which(part.list[, "uid1"] %in% dat$attr$uid & part.list[, "uid2"] %in% dat$attr$uid), , drop = FALSE]

    # Different partnership subsets
    part.listept.main.ong <- part.list[which((part.list[, "ptype"] == 1) &
                                                (part.list[, "last.active.time"] == at)), , drop = FALSE]
    part.listept.pers.ong <- part.list[which((part.list[, "ptype"] == 2) &
                                               (part.list[, "last.active.time"] == at)), , drop = FALSE]
    part.listept.main.end <- part.list[which((part.list[, "ptype"] == 1) &
                                                (part.list[, "last.active.time"] < at)), , drop = FALSE]
    part.listept.pers.end <- part.list[which((part.list[, "ptype"] == 2) &
                                                (part.list[, "last.active.time"] < at)), , drop = FALSE]
    part.listept.inst <- part.list[which((part.list[, "ptype"] == 3)), , drop = FALSE]


    ### Partner 1 has been given EPT, so partner 2 eligible
    ## Main, ongoing
    # List Partner 1 IDs
    idspartlist.col1.main.ong <- which(dat$attr$uid %in% part.listept.main.ong[, "uid1"])

    # Return ID for partner 1 who has been given EPT
    idspartlist.col1.ept.main.ong <- idspartlist.col1.main.ong[which(eptindexStat[idspartlist.col1.main.ong] == 1)]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.main.ong <- part.listept.main.ong[which(part.listept.main.ong[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.main.ong]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.main.ong <- which(dat$attr$uid %in% partlist.col1.ept.main.ong[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.main.ong <- idspartlistsept1.main.ong[which(rGC.tx[idspartlistsept1.main.ong] %in% c(0, NA) |
                                                                   uGC.tx[idspartlistsept1.main.ong] %in% c(0, NA) |
                                                                   rCT.tx[idspartlistsept1.main.ong] %in% c(0, NA) |
                                                                   uCT.tx[idspartlistsept1.main.ong] %in% c(0, NA) |
                                                                   rGC.tx.prep[idspartlistsept1.main.ong] %in% c(0, NA) |
                                                                   uGC.tx.prep[idspartlistsept1.main.ong] %in% c(0, NA) |
                                                                   rCT.tx.prep[idspartlistsept1.main.ong] %in% c(0, NA) |
                                                                   uCT.tx.prep[idspartlistsept1.main.ong] %in% c(0, NA) |
                                                                   is.na(rGC.tx.ept[idspartlistsept1.main.ong]) |
                                                                   is.na(uGC.tx.ept[idspartlistsept1.main.ong]) |
                                                                   is.na(rCT.tx.ept[idspartlistsept1.main.ong]) |
                                                                   is.na(uCT.tx.ept[idspartlistsept1.main.ong]))]
    ## Casual, ongoing
    # List Partner 1 IDs
    idspartlist.col1.pers.ong <- which(dat$attr$uid %in% part.listept.pers.ong[, "uid1"])

    # Return ID for partner 1 who has been given EPT
    idspartlist.col1.ept.pers.ong <- idspartlist.col1.pers.ong[which(eptindexStat[idspartlist.col1.pers.ong] == 1)]
    eptindexStat[idspartlist.col1.ept.pers.ong]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.pers.ong <- part.listept.pers.ong[which(part.listept.pers.ong[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.pers.ong]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.pers.ong <- which(dat$attr$uid %in% partlist.col1.ept.pers.ong[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.pers.ong <- idspartlistsept1.pers.ong[which(rGC.tx[idspartlistsept1.pers.ong] %in% c(0, NA) |
                                                                   uGC.tx[idspartlistsept1.pers.ong] %in% c(0, NA) |
                                                                   rCT.tx[idspartlistsept1.pers.ong] %in% c(0, NA) |
                                                                   uCT.tx[idspartlistsept1.pers.ong] %in% c(0, NA) |
                                                                   rGC.tx.prep[idspartlistsept1.pers.ong] %in% c(0, NA) |
                                                                   uGC.tx.prep[idspartlistsept1.pers.ong] %in% c(0, NA) |
                                                                   rCT.tx.prep[idspartlistsept1.pers.ong] %in% c(0, NA) |
                                                                   uCT.tx.prep[idspartlistsept1.pers.ong] %in% c(0, NA) |
                                                                   is.na(rGC.tx.ept[idspartlistsept1.pers.ong]) |
                                                                   is.na(uGC.tx.ept[idspartlistsept1.pers.ong]) |
                                                                   is.na(rCT.tx.ept[idspartlistsept1.pers.ong]) |
                                                                   is.na(uCT.tx.ept[idspartlistsept1.pers.ong]))]

    ## Main, ended
    # List Partner 1 IDs
    idspartlist.col1.main.end <- which(dat$attr$uid %in% part.listept.main.end[, "uid1"])

    # Return ID for partner 1 who has been given EPT
    idspartlist.col1.ept.main.end <- idspartlist.col1.main.end[which(eptindexStat[idspartlist.col1.main.end] == 1)]
    eptindexStat[idspartlist.col1.ept.main.end]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.main.end <- part.listept.main.end[which(part.listept.main.end[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.main.end]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.main.end <- which(dat$attr$uid %in% partlist.col1.ept.main.end[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.main.end <- idspartlistsept1.main.end[which(rGC.tx[idspartlistsept1.main.end] %in% c(0, NA) |
                                                                   uGC.tx[idspartlistsept1.main.end] %in% c(0, NA) |
                                                                   rCT.tx[idspartlistsept1.main.end] %in% c(0, NA) |
                                                                   uCT.tx[idspartlistsept1.main.end] %in% c(0, NA) |
                                                                   rGC.tx.prep[idspartlistsept1.main.end] %in% c(0, NA) |
                                                                   uGC.tx.prep[idspartlistsept1.main.end] %in% c(0, NA) |
                                                                   rCT.tx.prep[idspartlistsept1.main.end] %in% c(0, NA) |
                                                                   uCT.tx.prep[idspartlistsept1.main.end] %in% c(0, NA) |
                                                                   is.na(rGC.tx.ept[idspartlistsept1.main.end]) |
                                                                   is.na(uGC.tx.ept[idspartlistsept1.main.end]) |
                                                                   is.na(rCT.tx.ept[idspartlistsept1.main.end]) |
                                                                   is.na(uCT.tx.ept[idspartlistsept1.main.end]))]
    ## Casual, ended
    # List Partner 1 IDs
    idspartlist.col1.pers.end <- which(dat$attr$uid %in% part.listept.pers.end[, "uid1"])

    # Return ID for partner 1 who has been given EPT
    idspartlist.col1.ept.pers.end <- idspartlist.col1.pers.end[which(eptindexStat[idspartlist.col1.pers.end] == 1)]
    eptindexStat[idspartlist.col1.ept.pers.end]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.pers.end <- part.listept.pers.end[which(part.listept.pers.end[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.pers.end]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.pers.end <- which(dat$attr$uid %in% partlist.col1.ept.pers.end[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.pers.end <- idspartlistsept1.pers.end[which(rGC.tx[idspartlistsept1.pers.end] %in% c(0, NA) |
                                                                   uGC.tx[idspartlistsept1.pers.end] %in% c(0, NA) |
                                                                   rCT.tx[idspartlistsept1.pers.end] %in% c(0, NA) |
                                                                   uCT.tx[idspartlistsept1.pers.end] %in% c(0, NA) |
                                                                   rGC.tx.prep[idspartlistsept1.pers.end] %in% c(0, NA) |
                                                                   uGC.tx.prep[idspartlistsept1.pers.end] %in% c(0, NA) |
                                                                   rCT.tx.prep[idspartlistsept1.pers.end] %in% c(0, NA) |
                                                                   uCT.tx.prep[idspartlistsept1.pers.end] %in% c(0, NA) |
                                                                   is.na(rGC.tx.ept[idspartlistsept1.pers.end]) |
                                                                   is.na(uGC.tx.ept[idspartlistsept1.pers.end]) |
                                                                   is.na(rCT.tx.ept[idspartlistsept1.pers.end]) |
                                                                   is.na(uCT.tx.ept[idspartlistsept1.pers.end]))]

    ## Instantaneous
    # List Partner 1 IDs
    idspartlist.col1.inst <- which(dat$attr$uid %in% part.listept.inst[, "uid1"])

    # Return ID for partner 1 who has been given EPT
    idspartlist.col1.ept.inst <- idspartlist.col1.inst[which(eptindexStat[idspartlist.col1.inst] == 1)]
    eptindexStat[idspartlist.col1.ept.inst]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.inst <- part.listept.inst[which(part.listept.inst[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.inst]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.inst <- which(dat$attr$uid %in% partlist.col1.ept.inst[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.inst <- idspartlistsept1.inst[which(rGC.tx[idspartlistsept1.inst] %in% c(0, NA) |
                                                                   uGC.tx[idspartlistsept1.inst] %in% c(0, NA) |
                                                                   rCT.tx[idspartlistsept1.inst] %in% c(0, NA) |
                                                                   uCT.tx[idspartlistsept1.inst] %in% c(0, NA) |
                                                                   rGC.tx.prep[idspartlistsept1.inst] %in% c(0, NA) |
                                                                   uGC.tx.prep[idspartlistsept1.inst] %in% c(0, NA) |
                                                                   rCT.tx.prep[idspartlistsept1.inst] %in% c(0, NA) |
                                                                   uCT.tx.prep[idspartlistsept1.inst] %in% c(0, NA) |
                                                                   is.na(rGC.tx.ept[idspartlistsept1.inst]) |
                                                                   is.na(uGC.tx.ept[idspartlistsept1.inst]) |
                                                                   is.na(rCT.tx.ept[idspartlistsept1.inst]) |
                                                                   is.na(uCT.tx.ept[idspartlistsept1.inst]))]


    ### Partner 2 has been given EPT, so partner 1 eligible
    ## Main, ongoing
    idspartlist.col2.main.ong <- which(dat$attr$uid %in% part.listept.main.ong[, "uid2"])

    # Return ID for partner 2 who has been given EPT
    idspartlist.col2.ept.main.ong <- idspartlist.col2.main.ong[which(eptindexStat[idspartlist.col2.main.ong] == 1)]
    eptindexStat[idspartlist.col2.ept.main.ong]

    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.main.ong <- part.listept.main.ong[which(part.listept.main.ong[, "uid2"] %in% dat$attr$uid[idspartlist.col1.ept.main.ong]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.main.ong <- which(dat$attr$uid %in% partlist.col2.ept.main.ong[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.main.ong <- idspartlistsept2.main.ong[which(rGC.tx[idspartlistsept2.main.ong] %in% c(0, NA) |
                                                                   uGC.tx[idspartlistsept2.main.ong] %in% c(0, NA) |
                                                                   rCT.tx[idspartlistsept2.main.ong] %in% c(0, NA) |
                                                                   uCT.tx[idspartlistsept2.main.ong] %in% c(0, NA) |
                                                                   rGC.tx.prep[idspartlistsept2.main.ong] %in% c(0, NA) |
                                                                   uGC.tx.prep[idspartlistsept2.main.ong] %in% c(0, NA) |
                                                                   rCT.tx.prep[idspartlistsept2.main.ong] %in% c(0, NA) |
                                                                   uCT.tx.prep[idspartlistsept2.main.ong] %in% c(0, NA) |
                                                                   is.na(rGC.tx.ept[idspartlistsept2.main.ong]) |
                                                                   is.na(uGC.tx.ept[idspartlistsept2.main.ong]) |
                                                                   is.na(rCT.tx.ept[idspartlistsept2.main.ong]) |
                                                                   is.na(uCT.tx.ept[idspartlistsept2.main.ong]))]
    ## Casual, ongoing
    # List Partner 2 IDs
    idspartlist.col2.pers.ong <- which(dat$attr$uid %in% part.listept.pers.ong[, "uid2"])

    # Return ID for partner 2 who has been given EPT
    idspartlist.col2.ept.pers.ong <- idspartlist.col2.pers.ong[which(eptindexStat[idspartlist.col2.pers.ong] == 1)]
    eptindexStat[idspartlist.col2.ept.pers.ong]

    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.pers.ong <- part.listept.pers.ong[which(part.listept.pers.ong[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.pers.ong]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.pers.ong <- which(dat$attr$uid %in% partlist.col2.ept.pers.ong[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.pers.ong <- idspartlistsept2.pers.ong[which(rGC.tx[idspartlistsept2.pers.ong] %in% c(0, NA) |
                                                                   uGC.tx[idspartlistsept2.pers.ong] %in% c(0, NA) |
                                                                   rCT.tx[idspartlistsept2.pers.ong] %in% c(0, NA) |
                                                                   uCT.tx[idspartlistsept2.pers.ong] %in% c(0, NA) |
                                                                   rGC.tx.prep[idspartlistsept2.pers.ong] %in% c(0, NA) |
                                                                   uGC.tx.prep[idspartlistsept2.pers.ong] %in% c(0, NA) |
                                                                   rCT.tx.prep[idspartlistsept2.pers.ong] %in% c(0, NA) |
                                                                   uCT.tx.prep[idspartlistsept2.pers.ong] %in% c(0, NA) |
                                                                   is.na(rGC.tx.ept[idspartlistsept2.pers.ong]) |
                                                                   is.na(uGC.tx.ept[idspartlistsept2.pers.ong]) |
                                                                   is.na(rCT.tx.ept[idspartlistsept2.pers.ong]) |
                                                                   is.na(uCT.tx.ept[idspartlistsept2.pers.ong]))]

    ## Main, ended
    # List Partner 2 IDs
    idspartlist.col2.main.end <- which(dat$attr$uid %in% part.listept.main.end[, "uid2"])

    # Return ID for partner 2 who has been given EPT
    idspartlist.col2.ept.main.end <- idspartlist.col2.main.end[which(eptindexStat[idspartlist.col2.main.end] == 1)]
    eptindexStat[idspartlist.col2.ept.main.end]

    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.main.end <- part.listept.main.end[which(part.listept.main.end[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.main.end]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.main.end <- which(dat$attr$uid %in% partlist.col2.ept.main.end[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.main.end <- idspartlistsept1.main.end[which(rGC.tx[idspartlistsept2.main.end] %in% c(0, NA) |
                                                                   uGC.tx[idspartlistsept2.main.end] %in% c(0, NA) |
                                                                   rCT.tx[idspartlistsept2.main.end] %in% c(0, NA) |
                                                                   uCT.tx[idspartlistsept2.main.end] %in% c(0, NA) |
                                                                   rGC.tx.prep[idspartlistsept2.main.end] %in% c(0, NA) |
                                                                   uGC.tx.prep[idspartlistsept2.main.end] %in% c(0, NA) |
                                                                   rCT.tx.prep[idspartlistsept2.main.end] %in% c(0, NA) |
                                                                   uCT.tx.prep[idspartlistsept2.main.end] %in% c(0, NA) |
                                                                   is.na(rGC.tx.ept[idspartlistsept2.main.end]) |
                                                                   is.na(uGC.tx.ept[idspartlistsept2.main.end]) |
                                                                   is.na(rCT.tx.ept[idspartlistsept2.main.end]) |
                                                                   is.na(uCT.tx.ept[idspartlistsept2.main.end]))]
    ## Casual, ended
    # List Partner 2 IDs
    idspartlist.col2.pers.end <- which(dat$attr$uid %in% part.listept.pers.end[, "uid2"])

    # Return ID for partner 2 who has been given EPT
    idspartlist.col2.ept.pers.end <- idspartlist.col2.pers.end[which(eptindexStat[idspartlist.col2.pers.end] == 1)]
    eptindexStat[idspartlist.col2.ept.pers.end]

    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.pers.end <- part.listept.pers.end[which(part.listept.pers.end[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.pers.end]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.pers.end <- which(dat$attr$uid %in% partlist.col2.ept.pers.end[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.pers.end <- idspartlistsept2.pers.end[which(rGC.tx[idspartlistsept2.pers.end] %in% c(0, NA) |
                                                                   uGC.tx[idspartlistsept2.pers.end] %in% c(0, NA) |
                                                                   rCT.tx[idspartlistsept2.pers.end] %in% c(0, NA) |
                                                                   uCT.tx[idspartlistsept2.pers.end] %in% c(0, NA) |
                                                                   rGC.tx.prep[idspartlistsept2.pers.end] %in% c(0, NA) |
                                                                   uGC.tx.prep[idspartlistsept2.pers.end] %in% c(0, NA) |
                                                                   rCT.tx.prep[idspartlistsept2.pers.end] %in% c(0, NA) |
                                                                   uCT.tx.prep[idspartlistsept2.pers.end] %in% c(0, NA) |
                                                                   is.na(rGC.tx.ept[idspartlistsept2.pers.end]) |
                                                                   is.na(uGC.tx.ept[idspartlistsept2.pers.end]) |
                                                                   is.na(rCT.tx.ept[idspartlistsept2.pers.end]) |
                                                                   is.na(uCT.tx.ept[idspartlistsept2.pers.end]))]

    ## Instantaneous
    # List Partner 2 IDs
    idspartlist.col2.inst <- which(dat$attr$uid %in% part.listept.inst[, "uid2"])

    # Return ID for partner 2 who has been given EPT
    idspartlist.col2.ept.inst <- idspartlist.col1.inst[which(eptindexStat[idspartlist.col2.inst] == 1)]
    eptindexStat[idspartlist.col2.ept.inst]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col2.ept.inst <- part.listept.inst[which(part.listept.inst[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.inst]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept2.inst <- which(dat$attr$uid %in% partlist.col2.ept.inst[, "uid1"])

    # Check STI Tx status of partner 2
    idspartlistsept2.inst <- idspartlistsept2.inst[which(rGC.tx[idspartlistsept2.inst] %in% c(0, NA) |
                                                           uGC.tx[idspartlistsept2.inst] %in% c(0, NA) |
                                                           rCT.tx[idspartlistsept2.inst] %in% c(0, NA) |
                                                           uCT.tx[idspartlistsept2.inst] %in% c(0, NA) |
                                                           rGC.tx.prep[idspartlistsept2.inst] %in% c(0, NA) |
                                                           uGC.tx.prep[idspartlistsept2.inst] %in% c(0, NA) |
                                                           rCT.tx.prep[idspartlistsept2.inst] %in% c(0, NA) |
                                                           uCT.tx.prep[idspartlistsept2.inst] %in% c(0, NA) |
                                                           is.na(rGC.tx.ept[idspartlistsept2.inst]) |
                                                           is.na(uGC.tx.ept[idspartlistsept2.inst]) |
                                                           is.na(rCT.tx.ept[idspartlistsept2.inst]) |
                                                           is.na(uCT.tx.ept[idspartlistsept2.inst]))]

    # All EPT-tx eligible IDs (partners of index)
    idsept <- unique(c(idspartlistsept1.main.ong, idspartlistsept2.main.ong,
                       idspartlistsept1.pers.ong, idspartlistsept2.pers.ong,
                       idspartlistsept1.main.end, idspartlistsept2.main.end,
                       idspartlistsept1.pers.end, idspartlistsept2.pers.end,
                       idspartlistsept1.inst, idspartlistsept2.inst))

    idsept.main.ong <- unique(c(idspartlistsept1.main.ong,
                                idspartlistsept2.main.ong))
    idsept.pers.ong <- unique(c(idspartlistsept1.pers.ong,
                                idspartlistsept2.pers.ong))
    idsept.main.end <- unique(c(idspartlistsept1.main.end,
                                idspartlistsept2.main.end))
    idsept.pers.end <- unique(c(idspartlistsept1.pers.end,
                                idspartlistsept2.pers.end))
    idsept.inst <- unique(c(idspartlistsept1.inst,
                            idspartlistsept2.inst))

    ## Provision to non-index partners -----------------------------------------
    ##(to be treated at next time step)
    idsprovided.main.ong <- idsept.main.ong[which(rbinom(length(idsept.main.ong), 1,
                                                         ept.provision.main.ong) == 1)]
    idsprovided.pers.ong <- idsept.pers.ong[which(rbinom(length(idsept.pers.ong), 1,
                                                         ept.provision.pers.ong) == 1)]
    idsprovided.main.end <- idsept.main.end[which(rbinom(length(idsept.main.end), 1,
                                                         ept.provision.main.end) == 1)]
    idsprovided.pers.end <- idsept.pers.end[which(rbinom(length(idsept.pers.end), 1,
                                                         ept.provision.pers.end) == 1)]
    idsprovided.inst <- idsept.inst[which(rbinom(length(idsept.inst), 1,
                                                 ept.provision.inst) == 1)]

    idsprovided_ept <- unique(c(idsprovided.main.ong, idsprovided.pers.ong,
                                idsprovided.main.end, idsprovided.pers.end,
                                idsprovided.inst))

    idsprovided.main_ept <- unique(c(idsprovided.main.ong, idsprovided.main.end))

    idsprovided.pers_ept <- unique(c(idsprovided.pers.ong, idsprovided.pers.end))

    idsprovided.inst_ept <- unique(c(idsprovided.inst))

    # Uptake by non-index ------------------------------------------------------
    # Uptake occurs in same time step (or before next step) but non-index is
    # actually treated at next time step

    idsept_tx.main <- idsprovided.main_ept[which(rbinom(length(idsprovided.main_ept), 1,
                                                        ept.uptake.main) == 1)]
    idsept_tx.pers <- idsprovided.pers_ept[which(rbinom(length(idsprovided.pers_ept), 1,
                                                        ept.uptake.pers) == 1)]
    idsept_tx.inst <- idsprovided.inst_ept[which(rbinom(length(idsprovided.inst_ept), 1,
                                                        ept.uptake.inst) == 1)]
    idsuptake_ept <- unique(c(idsept_tx.main, idsept_tx.pers, idsept_tx.inst))

    ## Output -----------------------------------------------------------------

    # Index attributes
    dat$attr$eptindexElig <- eptindexElig
    dat$attr$eptindexStat <- eptindexStat
    dat$attr$eptindexEligdate <- eptindexEligdate

    # Non-index attributes
    dat$attr$eptpartEligReceive[idsept] <- 1
    dat$attr$eptpartEligTx[idsprovided_ept] <- 0
    dat$attr$eptpartEligTx[idsuptake_ept] <- 1
    dat$attr$eptpartEligTxdate[idsprovided_ept] <- at

    # Update Epi
    dat$epi$eptpartelig[at] <- length(idsept)
    dat$epi$eptpartprovided[at] <- length(idsprovided_ept)
    dat$epi$eptpartuptake[at] <- length(idsuptake_ept)
    dat$epi$eptprop_provided[at] <- dat$epi$eptpartprovided[at] / dat$epi$eptpartelig[at]
    dat$epi$eptuninfectedprovided[at] <- sum(rGC[idsprovided_ept] == 0 &
                                             uGC[idsprovided_ept] == 0 &
                                             rCT[idsprovided_ept] == 0 &
                                             uCT[idsprovided_ept] == 0) /
                                         length(idsprovided_ept)
    dat$epi$eptuninfecteduptake[at] <- sum(rGC[idsuptake_ept] == 0 &
                                           uGC[idsuptake_ept] == 0 &
                                           rCT[idsuptake_ept] == 0 &
                                           uCT[idsuptake_ept] == 0) /
                                       length(idsuptake_ept)

    return(dat)
}
