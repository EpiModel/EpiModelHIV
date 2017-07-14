
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
# browser()
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
    # Return rows in each subset where partner 1 has been given EPT
    part.listept1.main.ong <- part.listept.main.ong[which(eptindexStat[which(dat$attr$uid %in% part.listept.main.ong[, "uid1"])] == 1), , drop = FALSE]
    part.listept1.pers.ong <- part.listept.pers.ong[which(eptindexStat[which(dat$attr$uid %in% part.listept.pers.ong[, "uid1"])] == 1), , drop = FALSE]
    part.listept1.main.end <- part.listept.main.end[which(eptindexStat[which(dat$attr$uid %in% part.listept.main.end[, "uid1"])] == 1), , drop = FALSE]
    part.listept1.pers.end <- part.listept.pers.end[which(eptindexStat[which(dat$attr$uid %in% part.listept.pers.end[, "uid1"])] == 1), , drop = FALSE]
    part.listept1.inst <- part.listept.inst[which(eptindexStat[which(dat$attr$uid %in% part.listept.inst[, "uid1"])] == 1), , drop = FALSE]

    # Choose other partner in those rows
    part2.ept.main.ong <- which(dat$attr$uid %in% part.listept1.main.ong[, "uid2"])
    part2.ept.pers.ong <- which(dat$attr$uid %in% part.listept1.pers.ong[, "uid2"])
    part2.ept.main.end <- which(dat$attr$uid %in% part.listept1.main.end[, "uid2"])
    part2.ept.pers.end <- which(dat$attr$uid %in% part.listept1.pers.end[, "uid2"])
    part2.ept.inst <- which(dat$attr$uid %in% part.listept1.inst[, "uid2"])

    # Check current Tx status for those partners - remove those currently being treated
    # idspartlistsept2.main.ong <- part2.ept.main.ong[]
    # idspartlistsept2.pers.ong <- part2.ept.pers.ong[]
    # idspartlistsept2.main.end <- part2.ept.main.end[]
    # idspartlistsept2.pers.end <- part2.ept.pers.end[]
    # idspartlistsept2.inst <- part2.ept.inst[]

    ### Partner 2 has been given EPT, so partner 1 eligible
    # Return rows where partner 2 has been given EPT
    part.listept2.main.ong <- part.listept.main.ong[which(eptindexStat[which(dat$attr$uid %in% part.listept.main.ong[, "uid2"])] == 1), , drop = FALSE]
    part.listept2.pers.ong <- part.listept.pers.ong[which(eptindexStat[which(dat$attr$uid %in% part.listept.pers.ong[, "uid2"])] == 1), , drop = FALSE]
    part.listept2.main.end <- part.listept.main.end[which(eptindexStat[which(dat$attr$uid %in% part.listept.main.end[, "uid2"])] == 1), , drop = FALSE]
    part.listept2.pers.end <- part.listept.pers.end[which(eptindexStat[which(dat$attr$uid %in% part.listept.pers.end[, "uid2"])] == 1), , drop = FALSE]
    part.listept2.inst <- part.listept.inst[which(eptindexStat[which(dat$attr$uid %in% part.listept.inst[, "uid2"])] == 1), , drop = FALSE]

    # Choose other partner in those rows
    part1.ept.main.ong <- which(dat$attr$uid %in% part.listept2.main.ong[, "uid1"])
    part1.ept.pers.ong <- which(dat$attr$uid %in% part.listept2.pers.ong[, "uid1"])
    part1.ept.main.end <- which(dat$attr$uid %in% part.listept2.main.end[, "uid1"])
    part1.ept.pers.end <- which(dat$attr$uid %in% part.listept2.pers.end[, "uid1"])
    part1.ept.inst <- which(dat$attr$uid %in% part.listept2.inst[, "uid1"])

    # Check current Tx status for those partners  - remove those currently being treated
    # idspartlistsept1.main.ong <- part1.ept.main.ong[]
    # idspartlistsept1.pers.ong <- part1.ept.pers.ong[]
    # idspartlistsept1.main.end <- part1.ept.main.end[]
    # idspartlistsept1.pers.end <- part1.ept.pers.end[]
    # idspartlistsept1.inst <- part1.ept.inst[]

    # part.listept1.main.ong <- part.list[which(((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
    #                                                (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
    #                                                (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
    #                                                (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
    #                                                (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
    #                                                (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
    #                                                (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
    #                                                (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
    #                                                is.na(rGC.tx.ept[part.list[, "uid2"]]) |
    #                                                is.na(uGC.tx.ept[part.list[, "uid2"]]) |
    #                                                is.na(rCT.tx.ept[part.list[, "uid2"]]) |
    #                                                is.na(uCT.tx.ept[part.list[, "uid2"]])) &
    #                                             eptindexStat[part.list[, "uid1"]] == 1), , drop = FALSE]

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
    # Has to be a time component here - can only uptake at next time step
    # Create partEligTxdate?
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

    # Update Epi
    dat$epi$eptpartelig[at] <- length(idsept)
    dat$epi$eptpartprovided[at] <- length(idsprovided_ept)
    dat$epi$eptpartuptake[at] <- length(idsuptake_ept)
    dat$epi$eptprop_provided[at] <- dat$epi$eptpartprovided[at] / dat$epi$eptpartelig[at]
    #dat$epi$eptuninfectedprovided[at] <- / length(idsprovided_ept)
    #dat$epi$eptuninfecteduptake[at] <- / length(idsprovided_ept)## ADD IN UNINFECTED MEN UPTAKE TX


    return(dat)
}
