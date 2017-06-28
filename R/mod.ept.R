
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
    eptindexStartTime <- dat$attr$eptindexStartTime # DBNU

    ## Parameters
    ept.risk.int <- dat$param$ept.risk.int
    ept.provision.main.ong <- dat$param$ept.provision.partner.main.ong
    ept.provision.casl.ong <- dat$param$ept.provision.partner.casl.ong
    ept.provision.main.end <- dat$param$ept.provision.partner.main.end
    ept.provision.casl.end <- dat$param$ept.provision.partner.casl.end
    ept.provision.inst <- dat$param$ept.provision.partner.inst
    ept.uptake.main <- dat$param$ept.uptake.partner.main
    ept.uptake.casl <- dat$param$ept.uptake.partner.casl
    ept.uptake.inst <- dat$param$ept.uptake.partner.inst

    # Partnership list
    part.list <- dat$temp$part.list

    ## Stoppage for Index ------------------------------------------------------

    # Index no longer eligible(> 60 days since treatment time)
    idseptExpired <- which(((at - eptindexEligdate) > ept.risk.int) & eptindexStat == 1)

    # Reset EPT status
    idsStp <- c(idseptExpired)
    eptindexStat[idsStp] <- NA
    eptindexElig[idsStp] <- NA


    # Indications for non-index-------------------------------------------------

    ## Eligibility of partners
    part.list <- dat$temp$part.list

    # Subset partner list to partnerships active within an EPT interval - last active date within risk interval
    part.list <- part.list[which((at - (part.list[, "last.active.time"]) <= ept.risk.int)), , drop = FALSE]

    # Convert uid to regular ids in partnership list
    #idspartlist <- which(uid %in% part.list[, c("uid1", "uid2")])

    #### Partner 1  recently treated, so partner 2 eligible for EPT
    #### Currently only eligible for EPT once

    # Criteria: eligible within EPT risk interval, currently untreated, partner
    # received EPT, main partnership, and last active today
    part.listept1.main.ong <- part.list[which((at - eptindexEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   is.na(rGC.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(uGC.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(rCT.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(uCT.tx.ept[part.list[, "uid2"]])) &
                                                eptindexStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 1 &
                                                (part.list[, "last.active.time"]) == at), , drop = FALSE]

    # Criteria: eligible within EPT risk interval, currently untreated, partner
    # received EPT, casual partnership, and last active today
    part.listept1.casl.ong <- part.list[which((at - eptindexEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   is.na(rGC.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(uGC.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(rCT.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(uCT.tx.ept[part.list[, "uid2"]])) &
                                                eptindexStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 2 &
                                                (part.list[, "last.active.time"]) == at), , drop = FALSE]

    # Criteria: eligible within EPT risk interval, currently untreated, partner
    # received EPT, one-off partnership, and not active today
    part.listept1.main.end <- part.list[which((at - eptindexEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   is.na(rGC.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(uGC.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(rCT.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(uCT.tx.ept[part.list[, "uid2"]])) &
                                                eptindexStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 1 &
                                                (at - part.list[, "last.active.time"]) > 0), , drop = FALSE]

    # Criteria: eligible within EPT risk interval, currently untreated, index partner
    # received EPT, main partnership, and not active today
    part.listept1.casl.end <- part.list[which((at - eptindexEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                   is.na(rGC.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(uGC.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(rCT.tx.ept[part.list[, "uid2"]]) |
                                                   is.na(uCT.tx.ept[part.list[, "uid2"]])) &
                                                eptindexStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 2 &
                                                (at - part.list[, "last.active.time"]) > 0), , drop = FALSE]

    # Criteria: eligible within EPT risk interval, currently untreated, index partner
    # received EPT, one-off
    part.listept1.inst <- part.list[which((at - eptindexEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                            ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                               (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                               (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                               (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                               (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                               (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                               (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                               (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                               is.na(rGC.tx.ept[part.list[, "uid2"]]) |
                                               is.na(uGC.tx.ept[part.list[, "uid2"]]) |
                                               is.na(rCT.tx.ept[part.list[, "uid2"]]) |
                                               is.na(uCT.tx.ept[part.list[, "uid2"]])) &
                                            eptindexStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 3), , drop = FALSE]


    idspartlistsept1.main.ong <- part.listept1.main.ong[, "uid2"]
    idspartlistsept1.casl.ong <- part.listept1.casl.ong[, "uid2"]
    idspartlistsept1.main.end <- part.listept1.main.end[, "uid2"]
    idspartlistsept1.casl.end <- part.listept1.casl.end[, "uid2"]
    idspartlistsept1.inst <- part.listept1.inst[, "uid2"]


    ### Partner 1  recently treated, so partner 2 eligible for EPT

    # Criteria: eligible within EPT risk interval, currently untreated, partner
    # received EPT, main partnership, and last active today
    part.listept2.main.ong <- part.list[which((at - eptindexEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   is.na(rGC.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(uGC.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(rCT.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(uCT.tx.ept[part.list[, "uid1"]])) &
                                                eptindexStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 1 &
                                                (part.list[, "last.active.time"]) == at), , drop = FALSE]

    # Criteria: eligible within EPT risk interval, currently untreated, partner
    # received EPT, casual partnership, and last active today
    part.listept2.casl.ong <- part.list[which((at - eptindexEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   is.na(rGC.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(uGC.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(rCT.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(uCT.tx.ept[part.list[, "uid1"]])) &
                                                eptindexStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 2 &
                                                (part.list[, "last.active.time"]) == at), , drop = FALSE]

    # Criteria: eligible within EPT risk interval, currently untreated, partner
    # received EPT, one-off partnership, and not active today
    part.listept2.main.end <- part.list[which((at - eptindexEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   is.na(rGC.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(uGC.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(rCT.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(uCT.tx.ept[part.list[, "uid1"]])) &
                                                eptindexStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 1 &
                                                (at - part.list[, "last.active.time"]) > 0), , drop = FALSE]

    # Criteria: eligible within EPT risk interval, currently untreated, index partner
    # received EPT, main partnership, and not active today
    part.listept2.casl.end <- part.list[which((at - eptindexEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                   is.na(rGC.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(uGC.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(rCT.tx.ept[part.list[, "uid1"]]) |
                                                   is.na(uCT.tx.ept[part.list[, "uid1"]])) &
                                                eptindexStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 2 &
                                                (at - part.list[, "last.active.time"]) > 0), , drop = FALSE]

    # Criteria: eligible within EPT risk interval, currently untreated, index partner
    # received EPT, one-off
    part.listept2.inst <- part.list[which((at - eptindexEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                            ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                               (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                               (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                               (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                               (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                               (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                               (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                               (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                               is.na(rGC.tx.ept[part.list[, "uid1"]]) |
                                               is.na(uGC.tx.ept[part.list[, "uid1"]]) |
                                               is.na(rCT.tx.ept[part.list[, "uid1"]]) |
                                               is.na(uCT.tx.ept[part.list[, "uid1"]])) &
                                            eptindexStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 3), , drop = FALSE]


    idspartlistsept2.main.ong <- part.listept2.main.ong[, "uid1"]
    idspartlistsept2.casl.ong <- part.listept2.casl.ong[, "uid1"]
    idspartlistsept2.main.end <- part.listept2.main.end[, "uid1"]
    idspartlistsept2.casl.end <- part.listept2.casl.end[, "uid1"]
    idspartlistsept2.inst <- part.listept2.inst[, "uid1"]

    # All EPT-tx eligible IDs (partners of index)
    idsept <- unique(c(idspartlistsept1.main.ong, idspartlistsept2.main.ong,
                       idspartlistsept1.casl.ong, idspartlistsept2.casl.ong,
                       idspartlistsept1.main.end, idspartlistsept2.main.end,
                       idspartlistsept1.casl.end, idspartlistsept2.casl.end,
                       idspartlistsept1.inst, idspartlistsept2.inst))

    idsept.main.ong <- unique(c(idspartlistsept1.main.ong,
                                idspartlistsept2.main.ong))
    idsept.casl.ong <- unique(c(idspartlistsept1.casl.ong,
                                idspartlistsept2.casl.ong))
    idsept.main.end <- unique(c(idspartlistsept1.main.end,
                                idspartlistsept2.main.end))
    idsept.casl.end <- unique(c(idspartlistsept1.casl.end,
                                idspartlistsept2.casl.end))
    idsept.inst <- unique(c(idspartlistsept1.inst,
                            idspartlistsept2.inst))

    ## Provision to and uptake of partners (to be treated at next time step)
    idsprovided.main.ong <- idsept.main.ong[which(rbinom(length(idsept.main.ong), 1,
                                                         ept.provision.main.ong) == 1)]
    idsprovided.casl.ong <- idsept.casl.ong[which(rbinom(length(idsept.casl.ong), 1,
                                                         ept.provision.casl.ong) == 1)]
    idsprovided.main.end <- idsept.main.end[which(rbinom(length(idsept.main.end), 1,
                                                         ept.provision.main.end) == 1)]
    idsprovided.casl.end <- idsept.casl.end[which(rbinom(length(idsept.casl.end), 1,
                                                         ept.provision.casl.end) == 1)]
    idsprovided.inst <- idsept.inst[which(rbinom(length(idsept.inst), 1,
                                                 ept.provision.inst) == 1)]

    idsprovided_ept <- unique(c(idsprovided.main.ong, idsprovided.casl.ong,
                                idsprovided.main.end, idsprovided.casl.end,
                                idsprovided.inst))

    idsprovided.main_ept <- unique(c(idsprovided.main.ong, idsprovided.main.end))

    idsprovided.casl_ept <- unique(c(idsprovided.casl.ong, idsprovided.casl.end))

    idsprovided.inst_ept <- unique(c(idsprovided.inst))

    # Uptake
    idsept_tx.main <- idsprovided.main_ept[which(rbinom(length(idsprovided.main_ept), 1,
                                                        ept.uptake.main) == 1)]
    idsept_tx.casl <- idsprovided.casl_ept[which(rbinom(length(idsprovided.casl_ept), 1,
                                                        ept.uptake.casl) == 1)]
    idsept_tx.inst <- idsprovided.inst_ept[which(rbinom(length(idsprovided.inst_ept), 1,
                                                        ept.uptake.inst) == 1)]
    idsuptake_ept <- unique(c(idsept_tx.main, idsept_tx.casl, idsept_tx.inst))

    ## Output -----------------------------------------------------------------

    # Index attributes
    dat$attr$eptindexElig <- eptindexElig
    dat$attr$eptindexStat <- eptindexStat
    dat$attr$eptindexStartTime <- eptindexStartTime

    # Non-index attributes
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
