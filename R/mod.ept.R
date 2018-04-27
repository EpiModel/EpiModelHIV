
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
    eptpartEligReceive <- dat$attr$eptpartEligReceive

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


    ## Indications for non-index--------------------------------------
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


    ## Gonorrhea--------------------------------------

    ### Partner 1 has been given EPT, so partner 2 eligible
    ## Main, ongoing
    # List Partner 1 IDs
    idspartlist.col1.main.ong.gc <- which(dat$attr$uid %in% part.listept.main.ong[, "uid1"])

    # Return ID for partner 1 who has been given EPT and is currently being treated for GC
    idspartlist.col1.ept.main.ong.gc <- idspartlist.col1.main.ong.gc[which(eptindexStat[idspartlist.col1.main.ong.gc] == 1 &
                                                                 (rGC.tx[idspartlist.col1.main.ong.gc] == 1 | uGC.tx[idspartlist.col1.main.ong.gc] == 1))]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.main.ong.gc <- part.listept.main.ong[which(part.listept.main.ong[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.main.ong.gc]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.main.ong.gc <- which(dat$attr$uid %in% partlist.col1.ept.main.ong.gc[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.main.ong.gc <- idspartlistsept1.main.ong.gc[which(rGC.tx[idspartlistsept1.main.ong.gc] %in% c(0, NA) &
                                                                   uGC.tx[idspartlistsept1.main.ong.gc] %in% c(0, NA) &
                                                                   rGC.tx.prep[idspartlistsept1.main.ong.gc] %in% c(0, NA) &
                                                                   uGC.tx.prep[idspartlistsept1.main.ong.gc] %in% c(0, NA) &
                                                                   is.na(rGC.tx.ept[idspartlistsept1.main.ong.gc]) &
                                                                   is.na(uGC.tx.ept[idspartlistsept1.main.ong.gc]) &
                                                                   is.na(eptpartEligReceive[idspartlistsept1.main.ong.gc]))]
    ## Casual, ongoing
    # List Partner 1 IDs
    idspartlist.col1.pers.ong.gc <- which(dat$attr$uid %in% part.listept.pers.ong[, "uid1"])

    # Return ID for partner 1 who has been given EPT and is currently being treated for GC
    idspartlist.col1.ept.pers.ong.gc <- idspartlist.col1.pers.ong.gc[which(eptindexStat[idspartlist.col1.pers.ong.gc] == 1 &
                                                                       (rGC.tx[idspartlist.col1.pers.ong.gc] == 1 | uGC.tx[idspartlist.col1.pers.ong.gc] == 1))]
    eptindexStat[idspartlist.col1.ept.pers.ong.gc]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.pers.ong.gc <- part.listept.pers.ong[which(part.listept.pers.ong[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.pers.ong.gc]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.pers.ong.gc <- which(dat$attr$uid %in% partlist.col1.ept.pers.ong.gc[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.pers.ong.gc <- idspartlistsept1.pers.ong.gc[which(rGC.tx[idspartlistsept1.pers.ong.gc] %in% c(0, NA) &
                                                                   uGC.tx[idspartlistsept1.pers.ong.gc] %in% c(0, NA) &
                                                                   rGC.tx.prep[idspartlistsept1.pers.ong.gc] %in% c(0, NA) &
                                                                   uGC.tx.prep[idspartlistsept1.pers.ong.gc] %in% c(0, NA) &
                                                                   is.na(rGC.tx.ept[idspartlistsept1.pers.ong.gc]) &
                                                                   is.na(uGC.tx.ept[idspartlistsept1.pers.ong.gc]) &
                                                                   is.na(eptpartEligReceive[idspartlistsept1.pers.ong.gc]))]

    ## Main, ended
    # List Partner 1 IDs
    idspartlist.col1.main.end.gc <- which(dat$attr$uid %in% part.listept.main.end[, "uid1"])

    # Return ID for partner 1 who has been given EPT and is currently being treated for GC
    idspartlist.col1.ept.main.end.gc <- idspartlist.col1.main.end.gc[which(eptindexStat[idspartlist.col1.main.end.gc] == 1 &
                                                                                 (rGC.tx[idspartlist.col1.main.end.gc] == 1 | uGC.tx[idspartlist.col1.main.end.gc] == 1))]
    eptindexStat[idspartlist.col1.ept.main.end.gc]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.main.end.gc <- part.listept.main.end[which(part.listept.main.end[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.main.end.gc]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.main.end.gc <- which(dat$attr$uid %in% partlist.col1.ept.main.end.gc[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.main.end.gc <- idspartlistsept1.main.end.gc[which(rGC.tx[idspartlistsept1.main.end.gc] %in% c(0, NA) &
                                                                   uGC.tx[idspartlistsept1.main.end.gc] %in% c(0, NA) &
                                                                   rGC.tx.prep[idspartlistsept1.main.end.gc] %in% c(0, NA) &
                                                                   uGC.tx.prep[idspartlistsept1.main.end.gc] %in% c(0, NA) &
                                                                   is.na(rGC.tx.ept[idspartlistsept1.main.end.gc]) &
                                                                   is.na(uGC.tx.ept[idspartlistsept1.main.end.gc]) &
                                                                   is.na(eptpartEligReceive[idspartlistsept1.main.end.gc]))]
    ## Casual, ended
    # List Partner 1 IDs
    idspartlist.col1.pers.end.gc <- which(dat$attr$uid %in% part.listept.pers.end[, "uid1"])

    # Return ID for partner 1 who has been given EPT and is currently being treated for GC
    idspartlist.col1.ept.pers.end.gc <- idspartlist.col1.pers.end.gc[which(eptindexStat[idspartlist.col1.pers.end.gc] == 1 &
                                                                             (rGC.tx[idspartlist.col1.pers.end.gc] == 1 | uGC.tx[idspartlist.col1.pers.end.gc] == 1))]
    eptindexStat[idspartlist.col1.ept.pers.end.gc]


    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.pers.end.gc <- part.listept.pers.end[which(part.listept.pers.end[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.pers.end.gc]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.pers.end.gc <- which(dat$attr$uid %in% partlist.col1.ept.pers.end.gc[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.pers.end.gc <- idspartlistsept1.pers.end.gc[which(rGC.tx[idspartlistsept1.pers.end.gc] %in% c(0, NA) &
                                                                   uGC.tx[idspartlistsept1.pers.end.gc] %in% c(0, NA) &
                                                                   rGC.tx.prep[idspartlistsept1.pers.end.gc] %in% c(0, NA) &
                                                                   uGC.tx.prep[idspartlistsept1.pers.end.gc] %in% c(0, NA) &
                                                                   is.na(rGC.tx.ept[idspartlistsept1.pers.end.gc]) &
                                                                   is.na(uGC.tx.ept[idspartlistsept1.pers.end.gc]) &
                                                                   is.na(eptpartEligReceive[idspartlistsept1.pers.end.gc]))]

    ## Instantaneous
    # List Partner 1 IDs
    idspartlist.col1.inst.gc <- which(dat$attr$uid %in% part.listept.inst[, "uid1"])

    # Return ID for partner 1 who has been given EPT and is currently being treated for GC
    idspartlist.col1.ept.inst.gc <- idspartlist.col1.inst.gc[which(eptindexStat[idspartlist.col1.inst.gc] == 1 &
                                                                             (rGC.tx[idspartlist.col1.inst.gc] == 1 | uGC.tx[idspartlist.col1.inst.gc] == 1))]
    eptindexStat[idspartlist.col1.ept.inst.gc]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.inst.gc <- part.listept.inst[which(part.listept.inst[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.inst.gc]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.inst.gc <- which(dat$attr$uid %in% partlist.col1.ept.inst.gc[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.inst.gc <- idspartlistsept1.inst.gc[which(rGC.tx[idspartlistsept1.inst.gc] %in% c(0, NA) &
                                                                   uGC.tx[idspartlistsept1.inst.gc] %in% c(0, NA) &
                                                                   rGC.tx.prep[idspartlistsept1.inst.gc] %in% c(0, NA) &
                                                                   uGC.tx.prep[idspartlistsept1.inst.gc] %in% c(0, NA) &
                                                                   is.na(rGC.tx.ept[idspartlistsept1.inst.gc]) &
                                                                   is.na(uGC.tx.ept[idspartlistsept1.inst.gc]) &
                                                                   is.na(eptpartEligReceive[idspartlistsept1.inst.gc]))]


    ### Partner 2 has been given EPT, so partner 1 eligible
    ## Main, ongoing
    idspartlist.col2.main.ong.gc <- which(dat$attr$uid %in% part.listept.main.ong[, "uid2"])

    # Return ID for partner 2 who has been given EPT and is currently being treated for GC
    idspartlist.col2.ept.main.ong.gc <- idspartlist.col2.main.ong.gc[which(eptindexStat[idspartlist.col2.main.ong.gc] == 1 &
                                                                    (rGC.tx[idspartlist.col2.main.ong.gc] == 1 | uGC.tx[idspartlist.col2.main.ong.gc] == 1))]
    eptindexStat[idspartlist.col2.ept.main.ong.gc]

    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.main.ong.gc <- part.listept.main.ong[which(part.listept.main.ong[, "uid2"] %in% dat$attr$uid[idspartlist.col1.ept.main.ong.gc]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.main.ong.gc <- which(dat$attr$uid %in% partlist.col2.ept.main.ong.gc[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.main.ong.gc <- idspartlistsept2.main.ong.gc[which(rGC.tx[idspartlistsept2.main.ong.gc] %in% c(0, NA) &
                                                                   uGC.tx[idspartlistsept2.main.ong.gc] %in% c(0, NA) &
                                                                   rGC.tx.prep[idspartlistsept2.main.ong.gc] %in% c(0, NA) &
                                                                   uGC.tx.prep[idspartlistsept2.main.ong.gc] %in% c(0, NA) &
                                                                   is.na(rGC.tx.ept[idspartlistsept2.main.ong.gc]) &
                                                                   is.na(uGC.tx.ept[idspartlistsept2.main.ong.gc]) &
                                                                   is.na(eptpartEligReceive[idspartlistsept2.main.ong.gc]))]
    ## Casual, ongoing
    # List Partner 2 IDs
    idspartlist.col2.pers.ong.gc <- which(dat$attr$uid %in% part.listept.pers.ong[, "uid2"])

    # Return ID for partner 2 who has been given EPT and is currently being treated for GC
    idspartlist.col2.ept.pers.ong.gc <- idspartlist.col2.pers.ong.gc[which(eptindexStat[idspartlist.col2.pers.ong.gc] == 1 &
                                                                             (rGC.tx[idspartlist.col2.pers.ong.gc] == 1 | uGC.tx[idspartlist.col2.pers.ong.gc] == 1))]
    eptindexStat[idspartlist.col2.ept.pers.ong.gc]

    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.pers.ong.gc <- part.listept.pers.ong[which(part.listept.pers.ong[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.pers.ong.gc]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.pers.ong.gc <- which(dat$attr$uid %in% partlist.col2.ept.pers.ong.gc[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.pers.ong.gc <- idspartlistsept2.pers.ong.gc[which(rGC.tx[idspartlistsept2.pers.ong.gc] %in% c(0, NA) &
                                                                   uGC.tx[idspartlistsept2.pers.ong.gc] %in% c(0, NA) &
                                                                   rGC.tx.prep[idspartlistsept2.pers.ong.gc] %in% c(0, NA) &
                                                                   uGC.tx.prep[idspartlistsept2.pers.ong.gc] %in% c(0, NA) &
                                                                   is.na(rGC.tx.ept[idspartlistsept2.pers.ong.gc]) &
                                                                   is.na(uGC.tx.ept[idspartlistsept2.pers.ong.gc]) &
                                                                   is.na(eptpartEligReceive[idspartlistsept2.pers.ong.gc]))]

    ## Main, ended
    # List Partner 2 IDs
    idspartlist.col2.main.end.gc <- which(dat$attr$uid %in% part.listept.main.end[, "uid2"])

    # Return ID for partner 2 who has been given EPT and is currently being treated for GC
    idspartlist.col2.ept.main.end.gc <- idspartlist.col2.main.end.gc[which(eptindexStat[idspartlist.col2.main.end.gc] == 1 &
                                                                             (rGC.tx[idspartlist.col2.main.end.gc] == 1 | uGC.tx[idspartlist.col2.main.end.gc] == 1))]
    eptindexStat[idspartlist.col2.ept.main.end.gc]


    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.main.end.gc <- part.listept.main.end[which(part.listept.main.end[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.main.end.gc]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.main.end.gc <- which(dat$attr$uid %in% partlist.col2.ept.main.end.gc[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.main.end.gc <- idspartlistsept2.main.end.gc[which(rGC.tx[idspartlistsept2.main.end.gc] %in% c(0, NA) &
                                                                   uGC.tx[idspartlistsept2.main.end.gc] %in% c(0, NA) &
                                                                   rGC.tx.prep[idspartlistsept2.main.end.gc] %in% c(0, NA) &
                                                                   uGC.tx.prep[idspartlistsept2.main.end.gc] %in% c(0, NA) &
                                                                   is.na(rGC.tx.ept[idspartlistsept2.main.end.gc]) &
                                                                   is.na(uGC.tx.ept[idspartlistsept2.main.end.gc]) &
                                                                   is.na(eptpartEligReceive[idspartlistsept2.main.end.gc]))]
    ## Casual, ended
    # List Partner 2 IDs
    idspartlist.col2.pers.end.gc <- which(dat$attr$uid %in% part.listept.pers.end[, "uid2"])

    # Return ID for partner 2 who has been given EPT and is currently being treated for GC
    idspartlist.col2.ept.pers.end.gc <- idspartlist.col2.pers.end.gc[which(eptindexStat[idspartlist.col2.pers.end.gc] == 1 &
                                                                             (rGC.tx[idspartlist.col2.pers.end.gc] == 1 | uGC.tx[idspartlist.col2.pers.end.gc] == 1))]
    eptindexStat[idspartlist.col2.ept.pers.end.gc]

    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.pers.end.gc <- part.listept.pers.end[which(part.listept.pers.end[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.pers.end.gc]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.pers.end.gc <- which(dat$attr$uid %in% partlist.col2.ept.pers.end.gc[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.pers.end.gc <- idspartlistsept2.pers.end.gc[which(rGC.tx[idspartlistsept2.pers.end.gc] %in% c(0, NA) &
                                                                   uGC.tx[idspartlistsept2.pers.end.gc] %in% c(0, NA) &
                                                                   rGC.tx.prep[idspartlistsept2.pers.end.gc] %in% c(0, NA) &
                                                                   uGC.tx.prep[idspartlistsept2.pers.end.gc] %in% c(0, NA) &
                                                                   is.na(rGC.tx.ept[idspartlistsept2.pers.end.gc]) &
                                                                   is.na(uGC.tx.ept[idspartlistsept2.pers.end.gc]) &
                                                                   is.na(eptpartEligReceive[idspartlistsept2.pers.end.gc]))]

    ## Instantaneous
    # List Partner 2 IDs
    idspartlist.col2.inst.gc <- which(dat$attr$uid %in% part.listept.inst[, "uid2"])

    # Return ID for partner 2 who has been given EPT and is currently being treated for GC
    idspartlist.col2.ept.inst.gc <- idspartlist.col2.inst.gc[which(eptindexStat[idspartlist.col2.inst.gc] == 1 &
                                                                             (rGC.tx[idspartlist.col2.inst.gc] == 1 | uGC.tx[idspartlist.col2.inst.gc] == 1))]
    eptindexStat[idspartlist.col2.ept.inst.gc]


    # Return rows in each subset where partner 1 has been given EPT
    partlist.col2.ept.inst.gc <- part.listept.inst[which(part.listept.inst[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.inst.gc]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept2.inst.gc <- which(dat$attr$uid %in% partlist.col2.ept.inst.gc[, "uid1"])

    # Check STI Tx status of partner 2
    idspartlistsept2.inst.gc <- idspartlistsept2.inst.gc[which(rGC.tx[idspartlistsept2.inst.gc] %in% c(0, NA) &
                                                           uGC.tx[idspartlistsept2.inst.gc] %in% c(0, NA) &
                                                           rGC.tx.prep[idspartlistsept2.inst.gc] %in% c(0, NA) &
                                                           uGC.tx.prep[idspartlistsept2.inst.gc] %in% c(0, NA) &
                                                           is.na(rGC.tx.ept[idspartlistsept2.inst.gc]) &
                                                           is.na(uGC.tx.ept[idspartlistsept2.inst.gc]) &
                                                           is.na(eptpartEligReceive[idspartlistsept2.inst.gc]))]


    ## Chlamydia--------------------------------------

    ### Partner 1 has been given EPT, so partner 2 eligible
    ## Main, ongoing
    # List Partner 1 IDs
    idspartlist.col1.main.ong.ct <- which(dat$attr$uid %in% part.listept.main.ong[, "uid1"])

    # Return ID for partner 1 who has been given EPT and is being treated for CT
    idspartlist.col1.ept.main.ong.ct <- idspartlist.col1.main.ong.ct[which(eptindexStat[idspartlist.col1.main.ong.ct] == 1 &
                                                                       (rCT.tx[idspartlist.col1.main.ong.ct] == 1 | uCT.tx[idspartlist.col1.main.ong.ct] == 1))]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.main.ong.ct <- part.listept.main.ong[which(part.listept.main.ong[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.main.ong.ct]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.main.ong.ct <- which(dat$attr$uid %in% partlist.col1.ept.main.ong.ct[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.main.ong.ct <- idspartlistsept1.main.ong.ct[which(rCT.tx[idspartlistsept1.main.ong.ct] %in% c(0, NA) &
                                                                         uCT.tx[idspartlistsept1.main.ong.ct] %in% c(0, NA) &
                                                                         rCT.tx.prep[idspartlistsept1.main.ong.ct] %in% c(0, NA) &
                                                                         uCT.tx.prep[idspartlistsept1.main.ong.ct] %in% c(0, NA) &
                                                                         is.na(rCT.tx.ept[idspartlistsept1.main.ong.ct]) &
                                                                         is.na(uCT.tx.ept[idspartlistsept1.main.ong.ct]) &
                                                                         is.na(eptpartEligReceive[idspartlistsept1.main.ong.ct]))]
    ## Casual, ongoing
    # List Partner 1 IDs
    idspartlist.col1.pers.ong.ct <- which(dat$attr$uid %in% part.listept.pers.ong[, "uid1"])

    # Return ID for partner 1 who has been given EPT and is being treated for CT
    idspartlist.col1.ept.pers.ong.ct <- idspartlist.col1.pers.ong.ct[which(eptindexStat[idspartlist.col1.pers.ong.ct] == 1 &
                                                                             (rCT.tx[idspartlist.col1.pers.ong.ct] == 1 | uCT.tx[idspartlist.col1.pers.ong.ct] == 1))]
    eptindexStat[idspartlist.col1.ept.pers.ong.ct]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.pers.ong.ct <- part.listept.pers.ong[which(part.listept.pers.ong[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.pers.ong.ct]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.pers.ong.ct <- which(dat$attr$uid %in% partlist.col1.ept.pers.ong.ct[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.pers.ong.ct <- idspartlistsept1.pers.ong.ct[which(rCT.tx[idspartlistsept1.pers.ong.ct] %in% c(0, NA) &
                                                                         uCT.tx[idspartlistsept1.pers.ong.ct] %in% c(0, NA) &
                                                                         rCT.tx.prep[idspartlistsept1.pers.ong.ct] %in% c(0, NA) &
                                                                         uCT.tx.prep[idspartlistsept1.pers.ong.ct] %in% c(0, NA) &
                                                                         is.na(rCT.tx.ept[idspartlistsept1.pers.ong.ct]) &
                                                                         is.na(uCT.tx.ept[idspartlistsept1.pers.ong.ct]) &
                                                                         is.na(eptpartEligReceive[idspartlistsept1.pers.ong.ct]))]

    ## Main, ended
    # List Partner 1 IDs
    idspartlist.col1.main.end.ct <- which(dat$attr$uid %in% part.listept.main.end[, "uid1"])

    # Return ID for partner 1 who has been given EPT and is being treated for CT
    idspartlist.col1.ept.main.end.ct <- idspartlist.col1.main.end.ct[which(eptindexStat[idspartlist.col1.main.end.ct] == 1 &
                                                                             (rCT.tx[idspartlist.col1.main.end.ct] == 1 | uCT.tx[idspartlist.col1.main.end.ct] == 1))]
    eptindexStat[idspartlist.col1.ept.main.end.ct]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.main.end.ct <- part.listept.main.end[which(part.listept.main.end[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.main.end.ct]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.main.end.ct <- which(dat$attr$uid %in% partlist.col1.ept.main.end.ct[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.main.end.ct <- idspartlistsept1.main.end.ct[which(rCT.tx[idspartlistsept1.main.end.ct] %in% c(0, NA) &
                                                                         uCT.tx[idspartlistsept1.main.end.ct] %in% c(0, NA) &
                                                                         rCT.tx.prep[idspartlistsept1.main.end.ct] %in% c(0, NA) &
                                                                         uCT.tx.prep[idspartlistsept1.main.end.ct] %in% c(0, NA) &
                                                                         is.na(rCT.tx.ept[idspartlistsept1.main.end.ct]) &
                                                                         is.na(uCT.tx.ept[idspartlistsept1.main.end.ct]) &
                                                                         is.na(eptpartEligReceive[idspartlistsept1.main.end.ct]))]
    ## Casual, ended
    # List Partner 1 IDs
    idspartlist.col1.pers.end.ct <- which(dat$attr$uid %in% part.listept.pers.end[, "uid1"])

    # Return ID for partner 1 who has been given EPT and is being treated for CT
    idspartlist.col1.ept.pers.end.ct <- idspartlist.col1.pers.end.ct[which(eptindexStat[idspartlist.col1.pers.end.ct] == 1 &
                                                                             (rCT.tx[idspartlist.col1.pers.end.ct] == 1 | uCT.tx[idspartlist.col1.pers.end.ct] == 1))]
    eptindexStat[idspartlist.col1.ept.pers.end.ct]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.pers.end.ct <- part.listept.pers.end[which(part.listept.pers.end[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.pers.end.ct]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.pers.end.ct <- which(dat$attr$uid %in% partlist.col1.ept.pers.end.ct[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.pers.end.ct <- idspartlistsept1.pers.end.ct[which(rCT.tx[idspartlistsept1.pers.end.ct] %in% c(0, NA) &
                                                                         uCT.tx[idspartlistsept1.pers.end.ct] %in% c(0, NA) &
                                                                         rCT.tx.prep[idspartlistsept1.pers.end.ct] %in% c(0, NA) &
                                                                         uCT.tx.prep[idspartlistsept1.pers.end.ct] %in% c(0, NA) &
                                                                         is.na(rCT.tx.ept[idspartlistsept1.pers.end.ct]) &
                                                                         is.na(uCT.tx.ept[idspartlistsept1.pers.end.ct]) &
                                                                         is.na(eptpartEligReceive[idspartlistsept1.pers.end.ct]))]

    ## Instantaneous
    # List Partner 1 IDs
    idspartlist.col1.inst.ct <- which(dat$attr$uid %in% part.listept.inst[, "uid1"])

    # Return ID for partner 1 who has been given EPT and is being treated for CT
    idspartlist.col1.ept.inst.ct <- idspartlist.col1.inst.ct[which(eptindexStat[idspartlist.col1.inst.ct] == 1 &
                                                                             (rCT.tx[idspartlist.col1.inst.ct] == 1 | uCT.tx[idspartlist.col1.inst.ct] == 1))]
    eptindexStat[idspartlist.col1.ept.inst.ct]

    # Return rows in each subset where partner 1 has been given EPT
    partlist.col1.ept.inst.ct <- part.listept.inst[which(part.listept.inst[, "uid1"] %in% dat$attr$uid[idspartlist.col1.ept.inst.ct]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept1.inst.ct <- which(dat$attr$uid %in% partlist.col1.ept.inst.ct[, "uid2"])

    # Check STI Tx status of partner 2
    idspartlistsept1.inst.ct <- idspartlistsept1.inst.ct[which(rCT.tx[idspartlistsept1.inst.ct] %in% c(0, NA) &
                                                                 uCT.tx[idspartlistsept1.inst.ct] %in% c(0, NA) &
                                                                 rCT.tx.prep[idspartlistsept1.inst.ct] %in% c(0, NA) &
                                                                 uCT.tx.prep[idspartlistsept1.inst.ct] %in% c(0, NA) &
                                                                 is.na(rCT.tx.ept[idspartlistsept1.inst.ct]) &
                                                                 is.na(uCT.tx.ept[idspartlistsept1.inst.ct]) &
                                                                 is.na(eptpartEligReceive[idspartlistsept1.inst.ct]))]


    ### Partner 2 has been given EPT, so partner 1 eligible
    ## Main, ongoing
    idspartlist.col2.main.ong.ct <- which(dat$attr$uid %in% part.listept.main.ong[, "uid2"])

    # Return ID for partner 2 who has been given EPT and is being treated for CT
    idspartlist.col2.ept.main.ong.ct <- idspartlist.col2.main.ong.ct[which(eptindexStat[idspartlist.col2.main.ong.ct] == 1 &
                                                                     (rCT.tx[idspartlist.col2.main.ong.ct] == 1 | uCT.tx[idspartlist.col2.main.ong.ct] == 1))]
    eptindexStat[idspartlist.col2.ept.main.ong.ct]

    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.main.ong.ct <- part.listept.main.ong[which(part.listept.main.ong[, "uid2"] %in% dat$attr$uid[idspartlist.col1.ept.main.ong.ct]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.main.ong.ct <- which(dat$attr$uid %in% partlist.col2.ept.main.ong.ct[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.main.ong.ct <- idspartlistsept2.main.ong.ct[which(rCT.tx[idspartlistsept2.main.ong.ct] %in% c(0, NA) &
                                                                      uCT.tx[idspartlistsept2.main.ong.ct] %in% c(0, NA) &
                                                                      rCT.tx.prep[idspartlistsept2.main.ong.ct] %in% c(0, NA) &
                                                                      uCT.tx.prep[idspartlistsept2.main.ong.ct] %in% c(0, NA) &
                                                                      is.na(rCT.tx.ept[idspartlistsept2.main.ong.ct]) &
                                                                      is.na(uCT.tx.ept[idspartlistsept2.main.ong.ct]) &
                                                                      is.na(eptpartEligReceive[idspartlistsept2.main.ong.ct]))]
    ## Casual, ongoing
    # List Partner 2 IDs
    idspartlist.col2.pers.ong.ct <- which(dat$attr$uid %in% part.listept.pers.ong[, "uid2"])

    # Return ID for partner 2 who has been given EPT and is being treated for CT
    idspartlist.col2.ept.pers.ong.ct <- idspartlist.col2.pers.ong.ct[which(eptindexStat[idspartlist.col2.pers.ong.ct] == 1 &
                                                                             (rCT.tx[idspartlist.col2.pers.ong.ct] == 1 | uCT.tx[idspartlist.col2.pers.ong.ct] == 1))]
    eptindexStat[idspartlist.col2.ept.pers.ong.ct]

    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.pers.ong.ct <- part.listept.pers.ong[which(part.listept.pers.ong[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.pers.ong.ct]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.pers.ong.ct <- which(dat$attr$uid %in% partlist.col2.ept.pers.ong.ct[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.pers.ong.ct <- idspartlistsept2.pers.ong.ct[which(rCT.tx[idspartlistsept2.pers.ong.ct] %in% c(0, NA) &
                                                                         uCT.tx[idspartlistsept2.pers.ong.ct] %in% c(0, NA) &
                                                                         rCT.tx.prep[idspartlistsept2.pers.ong.ct] %in% c(0, NA) &
                                                                         uCT.tx.prep[idspartlistsept2.pers.ong.ct] %in% c(0, NA) &
                                                                         is.na(rCT.tx.ept[idspartlistsept2.pers.ong.ct]) &
                                                                         is.na(uCT.tx.ept[idspartlistsept2.pers.ong.ct]) &
                                                                         is.na(eptpartEligReceive[idspartlistsept2.pers.ong.ct]))]

    ## Main, ended
    # List Partner 2 IDs
    idspartlist.col2.main.end.ct <- which(dat$attr$uid %in% part.listept.main.end[, "uid2"])

    # Return ID for partner 2 who has been given EPT and is being treated for CT
    idspartlist.col2.ept.main.end.ct <- idspartlist.col2.main.end.ct[which(eptindexStat[idspartlist.col2.main.end.ct] == 1 &
                                                                             (rCT.tx[idspartlist.col2.main.end.ct] == 1 | uCT.tx[idspartlist.col2.main.end.ct] == 1))]
    eptindexStat[idspartlist.col2.ept.main.end.ct]

    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.main.end.ct <- part.listept.main.end[which(part.listept.main.end[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.main.end.ct]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.main.end.ct <- which(dat$attr$uid %in% partlist.col2.ept.main.end.ct[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.main.end.ct <- idspartlistsept2.main.end.ct[which(rCT.tx[idspartlistsept2.main.end.ct] %in% c(0, NA) &
                                                                         uCT.tx[idspartlistsept2.main.end.ct] %in% c(0, NA) &
                                                                         rCT.tx.prep[idspartlistsept2.main.end.ct] %in% c(0, NA) &
                                                                         uCT.tx.prep[idspartlistsept2.main.end.ct] %in% c(0, NA) &
                                                                         is.na(rCT.tx.ept[idspartlistsept2.main.end.ct]) &
                                                                         is.na(uCT.tx.ept[idspartlistsept2.main.end.ct]) &
                                                                         is.na(eptpartEligReceive[idspartlistsept2.main.end.ct]))]
    ## Casual, ended
    # List Partner 2 IDs
    idspartlist.col2.pers.end.ct <- which(dat$attr$uid %in% part.listept.pers.end[, "uid2"])

    # Return ID for partner 2 who has been given EPT and is being treated for CT
    idspartlist.col2.ept.pers.end.ct <- idspartlist.col2.pers.end.ct[which(eptindexStat[idspartlist.col2.pers.end.ct] == 1 &
                                                                             (rCT.tx[idspartlist.col2.pers.end.ct] == 1 | uCT.tx[idspartlist.col2.pers.end.ct] == 1))]
    eptindexStat[idspartlist.col2.ept.pers.end.ct]

    # Return rows in each subset where partner 2 has been given EPT
    partlist.col2.ept.pers.end.ct <- part.listept.pers.end[which(part.listept.pers.end[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.pers.end.ct]), , drop = FALSE]

    # Select IDs of partner 1
    idspartlistsept2.pers.end.ct <- which(dat$attr$uid %in% partlist.col2.ept.pers.end.ct[, "uid1"])

    # Check STI Tx status of partner 1
    idspartlistsept2.pers.end.ct <- idspartlistsept2.pers.end.ct[which(rCT.tx[idspartlistsept2.pers.end.ct] %in% c(0, NA) &
                                                                         uCT.tx[idspartlistsept2.pers.end.ct] %in% c(0, NA) &
                                                                         rCT.tx.prep[idspartlistsept2.pers.end.ct] %in% c(0, NA) &
                                                                         uCT.tx.prep[idspartlistsept2.pers.end.ct] %in% c(0, NA) &
                                                                         is.na(rCT.tx.ept[idspartlistsept2.pers.end.ct]) &
                                                                         is.na(uCT.tx.ept[idspartlistsept2.pers.end.ct]) &
                                                                         is.na(eptpartEligReceive[idspartlistsept2.pers.end.ct]))]

    ## Instantaneous
    # List Partner 2 IDs
    idspartlist.col2.inst.ct <- which(dat$attr$uid %in% part.listept.inst[, "uid2"])

    # Return ID for partner 2 who has been given EPT and is being treated for CT
    idspartlist.col2.ept.inst.ct <- idspartlist.col2.inst.ct[which(eptindexStat[idspartlist.col2.inst.ct] == 1 &
                                                                             (rCT.tx[idspartlist.col2.inst.ct] == 1 | uCT.tx[idspartlist.col2.inst.ct] == 1))]
    eptindexStat[idspartlist.col2.ept.inst.ct]


    # Return rows in each subset where partner 1 has been given EPT
    partlist.col2.ept.inst.ct <- part.listept.inst[which(part.listept.inst[, "uid2"] %in% dat$attr$uid[idspartlist.col2.ept.inst.ct]), , drop = FALSE]

    # Select IDs of partner 2
    idspartlistsept2.inst.ct <- which(dat$attr$uid %in% partlist.col2.ept.inst.ct[, "uid1"])

    # Check STI Tx status of partner 2
    idspartlistsept2.inst.ct <- idspartlistsept2.inst.ct[which(rCT.tx[idspartlistsept2.inst.ct] %in% c(0, NA) &
                                                              uCT.tx[idspartlistsept2.inst.ct] %in% c(0, NA) &
                                                              rCT.tx.prep[idspartlistsept2.inst.ct] %in% c(0, NA) &
                                                              uCT.tx.prep[idspartlistsept2.inst.ct] %in% c(0, NA) &
                                                              is.na(rCT.tx.ept[idspartlistsept2.inst.ct]) &
                                                              is.na(uCT.tx.ept[idspartlistsept2.inst.ct]) &
                                                              is.na(eptpartEligReceive[idspartlistsept2.inst.ct]))]

    # All EPT-tx eligible IDs (partners of index)
    idsept <- unique(c(idspartlistsept1.main.ong.gc, idspartlistsept2.main.ong.gc,
                       idspartlistsept1.pers.ong.gc, idspartlistsept2.pers.ong.gc,
                       idspartlistsept1.main.end.gc, idspartlistsept2.main.end.gc,
                       idspartlistsept1.pers.end.gc, idspartlistsept2.pers.end.gc,
                       idspartlistsept1.inst.gc, idspartlistsept2.inst.gc,
                       idspartlistsept1.main.ong.ct, idspartlistsept2.main.ong.ct,
                       idspartlistsept1.pers.ong.ct, idspartlistsept2.pers.ong.ct,
                       idspartlistsept1.main.end.ct, idspartlistsept2.main.end.ct,
                       idspartlistsept1.pers.end.ct, idspartlistsept2.pers.end.ct,
                       idspartlistsept1.inst.ct, idspartlistsept2.inst.ct))

    idsept.main.ong <- unique(c(idspartlistsept1.main.ong.gc,
                                idspartlistsept2.main.ong.gc,
                                idspartlistsept1.main.ong.ct,
                                idspartlistsept2.main.ong.ct))
    idsept.pers.ong <- unique(c(idspartlistsept1.pers.ong.gc,
                                idspartlistsept2.pers.ong.gc,
                                idspartlistsept1.pers.ong.ct,
                                idspartlistsept2.pers.ong.ct))
    idsept.main.end <- unique(c(idspartlistsept1.main.end.gc,
                                idspartlistsept2.main.end.gc,
                                idspartlistsept1.main.end.ct,
                                idspartlistsept2.main.end.ct))
    idsept.pers.end <- unique(c(idspartlistsept1.pers.end.gc,
                                idspartlistsept2.pers.end.gc,
                                idspartlistsept1.pers.end.ct,
                                idspartlistsept2.pers.end.ct))
    idsept.inst <- unique(c(idspartlistsept1.inst.gc,
                            idspartlistsept2.inst.gc,
                            idspartlistsept1.inst.ct,
                            idspartlistsept2.inst.ct))
    ids.ept.gc <- unique(c(idspartlistsept1.main.ong.gc, idspartlistsept2.main.ong.gc,
                                         idspartlistsept1.pers.ong.gc, idspartlistsept2.pers.ong.gc,
                                         idspartlistsept1.main.end.gc, idspartlistsept2.main.end.gc,
                                         idspartlistsept1.pers.end.gc, idspartlistsept2.pers.end.gc,
                                         idspartlistsept1.inst.gc, idspartlistsept2.inst.gc))
    ids.ept.ct <- unique(c(idspartlistsept1.main.ong.ct, idspartlistsept2.main.ong.ct,
                           idspartlistsept1.pers.ong.ct, idspartlistsept2.pers.ong.ct,
                           idspartlistsept1.main.end.ct, idspartlistsept2.main.end.ct,
                           idspartlistsept1.pers.end.ct, idspartlistsept2.pers.end.ct,
                           idspartlistsept1.inst.ct, idspartlistsept2.inst.ct))

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

    idsprovided_ept <- c(idsprovided.main.ong, idsprovided.pers.ong,
                                idsprovided.main.end, idsprovided.pers.end,
                                idsprovided.inst)

    idsprovided.main_ept <- c(idsprovided.main.ong, idsprovided.main.end)

    idsprovided.pers_ept <- c(idsprovided.pers.ong, idsprovided.pers.end)

    idsprovided.inst_ept <- c(idsprovided.inst)

    # Need to further refine
    idsprovided_gc <- intersect(idsprovided_ept, ids.ept.gc)
    idsprovided_ct <- intersect(idsprovided_ept, ids.ept.ct)

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

    idsept_tx.gc <- intersect(idsuptake_ept, ids.ept.gc)
    idsept_tx.ct <- intersect(idsuptake_ept, ids.ept.ct)


    ## Output -----------------------------------------------------------------

    # Update with new trackers
    if (is.null(dat$epi$eptpartprovided_gc)) {
      dat$epi$eptpartprovided_gc <- rep(NA, length(dat$control$nsteps))
      dat$epi$eptpartprovided_ct <- rep(NA, length(dat$control$nsteps))
      dat$epi$eptpartprovided_main <- rep(NA, length(dat$control$nsteps))
      dat$epi$eptpartprovided_pers <- rep(NA, length(dat$control$nsteps))
      dat$epi$eptpartprovided_inst <- rep(NA, length(dat$control$nsteps))
      dat$epi$eptpartuptake_main <- rep(NA, length(dat$control$nsteps))
      dat$epi$eptpartuptake_pers <- rep(NA, length(dat$control$nsteps))
      dat$epi$eptpartuptake_inst <- rep(NA, length(dat$control$nsteps))
      dat$epi$eptpartuptake_gc <- rep(NA, length(dat$control$nsteps))
      dat$epi$eptgcinfectundiaghiv <- rep(NA, length(dat$control$nsteps))
      dat$epi$eptctinfectundiaghiv <- rep(NA, length(dat$control$nsteps))
      dat$epi$eptgcctinfectundiaghiv <- rep(NA, length(dat$control$nsteps))
    }

    # Index attributes
    dat$attr$eptindexElig <- eptindexElig
    dat$attr$eptindexStat <- eptindexStat
    dat$attr$eptindexEligdate <- eptindexEligdate

    # Non-index attributes
    dat$attr$eptpartEligReceive[idsept] <- 1
    dat$attr$eptpartEligTx_GC[idsprovided_ept] <- 0
    dat$attr$eptpartEligTx_CT[idsprovided_ept] <- 0
    dat$attr$eptpartEligTx_GC[idsept_tx.gc] <- 1
    dat$attr$eptpartEligTx_CT[idsept_tx.ct] <- 1
    dat$attr$eptpartEligTxdate[idsprovided_ept] <- at

    # Update Epi
    dat$epi$eptpartelig[at] <- length(idsept)
    dat$epi$eptpartprovided[at] <- length(idsprovided_ept)
    dat$epi$eptpartprovided_gc[at] <- length(idsprovided_gc)
    dat$epi$eptpartprovided_ct[at] <- length(idsprovided_ct)
    dat$epi$eptpartprovided_main[at] <- length(idsprovided.main_ept)
    dat$epi$eptpartprovided_pers[at] <- length(idsprovided.pers_ept)
    dat$epi$eptpartprovided_inst[at] <- length(idsprovided.inst_ept)
    dat$epi$eptpartuptake[at] <- length(idsuptake_ept)
    dat$epi$eptpartuptake_main[at] <- length(idsept_tx.main)
    dat$epi$eptpartuptake_pers[at] <- length(idsept_tx.pers)
    dat$epi$eptpartuptake_inst[at] <- length(idsept_tx.inst)
    dat$epi$eptpartuptake_gc[at] <- length(idsept_tx.gc)
    dat$epi$eptpartuptake_ct[at] <- length(idsept_tx.ct)

    # Wasted EPT
    dat$epi$eptuninfectedprovided[at] <- length(which(rGC[idsprovided_ept] == 0 &
                                             uGC[idsprovided_ept] == 0 &
                                             rCT[idsprovided_ept] == 0 &
                                             uCT[idsprovided_ept] == 0))
    dat$epi$eptuninfecteduptake[at] <- length(which(rGC[idsuptake_ept] == 0 &
                                                  uGC[idsuptake_ept] == 0 &
                                                  rCT[idsuptake_ept] == 0 &
                                                  uCT[idsuptake_ept] == 0))

    # Missed opportunities EPT
    dat$epi$eptgcinfectsti[at] <-  length(idsept_tx.gc[which(dat$attr$diag.status.gc[idsept_tx.gc] == 1 |
                                                           dat$attr$diag.status.ct[idsept_tx.gc] == 1 |
                                                           dat$attr$diag.status.syph[idsept_tx.gc] == 1 |
                                                           dat$attr$diag.status[idsept_tx.gc] == 1)])
    dat$epi$eptctinfectsti[at] <- length(idsept_tx.ct[which(dat$attr$diag.status.gc[idsept_tx.ct] == 1 |
                                                          dat$attr$diag.status.ct[idsept_tx.ct] == 1 |
                                                          dat$attr$diag.status.syph[idsept_tx.ct] == 1 |
                                                          dat$attr$diag.status[idsept_tx.ct] == 1)])
    dat$epi$eptgcinfectundiaghiv[at] <- length(idsept_tx.gc[which(dat$attr$status[idsept_tx.gc] == 1 &
                                                                     dat$attr$diag.status[idsept_tx.gc] == 0)])
    dat$epi$eptctinfectundiaghiv[at] <- length(idsept_tx.ct[which(dat$attr$status[idsept_tx.ct] == 1 &
                                                               dat$attr$diag.status[idsept_tx.ct] == 0)])
    dat$epi$eptgcctinfectundiaghiv[at] <- length(idsuptake_ept[which(dat$attr$status[idsuptake_ept] == 1 &
                                                                           dat$attr$diag.status[idsuptake_ept] == 0)])


    return(dat)
}
