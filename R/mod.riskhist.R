
#' @title Risk History for PrEP Module
#'
#' @description Module function to track the risk history of uninfected persons
#'              for purpose of PrEP prevention intervention targeting.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_prep_msm <- function(dat, at) {

  if (at < dat$param$riskh.prep.start) {
    return(dat)
  }

  ## Attributes
  uid <- dat$attr$uid
  dx <- dat$attr$diag.status
  since.test <- at - dat$attr$last.neg.test
  rGC.tx <- dat$attr$rGC.tx
  uGC.tx <- dat$attr$uGC.tx
  rCT.tx <- dat$attr$rCT.tx
  uCT.tx <- dat$attr$uCT.tx
  syph.tx <- dat$attr$syph.tx

  ## Parameters
  time.unit <- dat$param$time.unit

  ## Edgelist, adds uai summation per partnership from act list
  pid <- NULL # For R CMD Check
  al <- as.data.frame(dat$temp$al)
  by_pid <- group_by(al, pid)
  uai <- summarise(by_pid, uai = sum(uai))[, 2]
  el <- as.data.frame(cbind(dat$temp$el, uai))

  # Remove concordant positive edges
  el2 <- el[el$st2 == 0, ]

  # Initialize attributes
  if (is.null(dat$attr$prep.ind.uai.mono)) {
    dat$attr$prep.ind.uai.mono <- rep(NA, length(uid))
    dat$attr$prep.ind.uai.nmain <- rep(NA, length(uid))
    dat$attr$prep.ind.ai.sd <- rep(NA, length(uid))
    dat$attr$prep.ind.sti <- rep(NA, length(uid))
  }

  ## Degree ##
  main.deg <- get_degree(dat$el[[1]])
  casl.deg <- get_degree(dat$el[[2]])
  inst.deg <- get_degree(dat$el[[3]])


  # Indications -------------------------------------------------------------

  # Any UAI
  uai.any <- unique(c(el2$p1[el2$uai > 0],
                      el2$p2[el2$uai > 0]))

  # Monogamous partnerships: 1-sided
  tot.deg <- main.deg + casl.deg + inst.deg
  uai.mono1 <- intersect(which(tot.deg == 1), uai.any)

  # "Negative" partnerships
  tneg <- unique(c(el2$p1[el2$st1 == 0], el2$p2[el2$st1 == 0]))
  fneg <- unique(c(el2$p1[which(dx[el2$p1] == 0)], el2$p2[which(dx[el2$p1] == 0)]))
  all.neg <- c(tneg, fneg)

  ## Condition 1b: UAI in 1-sided "monogamous" "negative" partnership,
  ##               partner not tested in past 6 months
  uai.mono1.neg <- intersect(uai.mono1, all.neg)
  part.id1 <- c(el2[el2$p1 %in% uai.mono1.neg, 2], el2[el2$p2 %in% uai.mono1.neg, 1])
  not.tested.6mo <- since.test[part.id1] > (180/time.unit)
  part.not.tested.6mo <- uai.mono1.neg[which(not.tested.6mo == TRUE)]
  dat$attr$prep.ind.uai.mono[part.not.tested.6mo] <- at

  ## Condition 2b: UAI in non-main partnerships
  uai.nmain <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$ptype %in% 2:3],
                        el2$p2[el2$uai > 0 & el2$ptype %in% 2:3]))
  dat$attr$prep.ind.uai.nmain[uai.nmain] <- at

  ## Condition 3a: AI within known serodiscordant partnerships
  el2.cond3 <- el2[el2$st1 == 1 & el2$ptype %in% 1:2, ]

  # Disclosure
  discl.list <- dat$temp$discl.list
  disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
  delt.cdl <- uid[el2.cond3[, 1]] * 1e7 + uid[el2.cond3[, 2]]
  discl <- (delt.cdl %in% disclose.cdl)
  ai.sd <- el2.cond3$p2[discl == TRUE]
  dat$attr$prep.ind.ai.sd[ai.sd] <- at

  ## Condition 4, any current STI diagnosis (before recovery)
  idsDx <- which(rGC.tx == 1 | uGC.tx == 1 |
                   rCT.tx == 1 | uCT.tx == 1 | syph.tx == 1)
  dat$attr$prep.ind.sti[idsDx] <- at


  return(dat)
}


#' @title Risk History for STI Testing Module
#'
#' @description Module function to track the risk history of uninfected persons
#'              for purpose of STI testing prevention intervention targeting.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_stitest_msm <- function(dat, at) {

  if (at < dat$param$riskh.stitest.start) {
    return(dat)
  }

  ## Parameters
  partnercutoff <- dat$param$partnercutoff

  ## Attributes
  uid <- dat$attr$uid
  race <- dat$attr$race

  # Indications -------------------------------------------------------------
  part.list <- dat$temp$part.list

  ### Lower risk - existing in partner list (6 month look back)
  idspartlist <- which(uid %in% part.list[, c("uid1", "uid2")])
  idsnotpartlist <- setdiff(which(race %in% c("B","W")), idspartlist)
  # these are relative ids of nodes in partner list

  ### High-risk: Have more than one sex partner in last x months
  # Reset # of partners at each time step- length of "recent" interval is drawn from interval of partner list lookback
  dat$attr$recentpartners <- rep(0, length(which(race %in% c("B","W"))))

  # For those who had partners, calculate # of occurrences in partner list
  part.count <- as.data.frame(table(part.list[, c("uid1", "uid2")]))

  # Calculate # of recent partners: 0 for those not in part list, update numbers for only actives in part list
  dat$attr$recentpartners[idspartlist] <- part.count[which(part.count[, "Var1"] %in% uid), 2]

  # Choose those who have had more than X partners in last x months
  idsrecentpartners <- which(dat$attr$recentpartners > partnercutoff)
  idsnotrecentpartners <- setdiff(which(race %in% c("B","W")), idsrecentpartners)

  ### Update STI indication attributes
  dat$attr$stitest.ind.active[idspartlist] <- 1
  dat$attr$stitest.ind.active[idsnotpartlist] <- 0
  dat$attr$stitest.ind.recentpartners[idsrecentpartners] <- 1
  dat$attr$stitest.ind.recentpartners[idsnotrecentpartners] <- 0

  return(dat)
}



#' @title Risk History for EPT Module
#'
#' @description Module function to track the risk history of uninfected persons
#'              for purpose of EPT prevention intervention targeting.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_ept_msm <- function(dat, at) {

  if (at < dat$param$riskh.ept.start) {
      return(dat)
  }

  ## Attributes
  rGC.tx <- dat$attr$rGC.tx
  uGC.tx <- dat$attr$uGC.tx
  rCT.tx <- dat$attr$rCT.tx
  uCT.tx <- dat$attr$uCT.tx
  rGC.tx.prep <- dat$attr$rGC.tx.prep
  uGC.tx.prep <- dat$attr$uGC.tx.prep
  rCT.tx.prep <- dat$attr$rCT.tx.prep
  uCT.tx.prep <- dat$attr$uCT.tx.prep
  eptStat <- dat$attr$eptStat
  eptEligdate <- dat$attr$eptEligdate
  eptEligTx <- dat$attr$eptEligTx       # DBNU
  eptStartTime <- dat$attr$eptStartTime # DBNU
  eptTx <- dat$attr$eptTx

  ## Parameters
  ept.risk.int <- dat$param$ept.risk.int
  ept.provision.main <- dat$param$ept.provision.partner.main
  ept.provision.casl <- dat$param$ept.provision.partner.casl
  ept.provision.inst <- dat$param$ept.provision.partner.inst
  ept.uptake.main <- dat$param$ept.uptake.partner.main
  ept.uptake.casl <- dat$param$ept.uptake.partner.casl
  ept.uptake.inst <- dat$param$ept.uptake.partner.inst
  ept.provision.short.main.rr <- dat$param$ept.provision.short.main.rr
  ept.provision.short.casl.rr <- dat$param$ept.provision.short.casl.rr
  ept.provision.short.inst.rr <- dat$param$ept.provision.short.inst.rr
  ept.provision.med.main.rr <- dat$param$ept.provision.med.main.rr
  ept.provision.med.casl.rr <- dat$param$ept.provision.med.casl.rr
  ept.provision.med.inst.rr <- dat$param$ept.provision.med.inst.rr
  ept.provision.long.main.rr <- dat$param$ept.provision.long.main.rr
  ept.provision.long.casl.rr <- dat$param$ept.provision.long.casl.rr
  ept.provision.long.inst.rr <- dat$param$ept.provision.long.inst.rr

  ## Edgelist, adds uai summation per partnership from act list
  pid <- NULL # For R CMD Check
  al <- as.data.frame(dat$temp$al)
  by_pid <- group_by(al, pid)
  uai <- summarise(by_pid, uai = sum(uai))[, 2]


  # Indications -------------------------------------------------------------

  ## Eligibility of partners
  part.list <- dat$temp$part.list

  # Subset partner list to those active within an EPT interval - last active date within 60 days
  part.list <- part.list[which((at - (part.list[, "last.active.time"]) <= ept.risk.int)), , drop = FALSE]

  # Subset partner list to alive ---> may need to remove!
  # Could still be indicated even if partner is dead
  # part.list <- part.list[which(active[part.list[, "uid1"]] == 1 & active[part.list[, "uid2"]] == 1), , drop = FALSE]

  # Partner 1 was recently treated, so partner 2 would be eligible to be treated
  # through EPT if not currently being treated for anything
  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, main partnership, and last active within last 20 days
  part.listept1.main.short <- part.list[which((at - eptEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                  ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA)) &
                                                  eptStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 1 &
                                                  (at - part.list[, "last.active.time"]) <= 2), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, casual partnership, and last active within last 20 days
  part.listept1.casl.short <- part.list[which((at - eptEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                  ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA)) &
                                                  eptStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 2 &
                                                  (at - part.list[, "last.active.time"]) <= 2), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, one-off partnership, and last active within last 20 days
  part.listept1.inst.short <- part.list[which((at - eptEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                  ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                       (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA)) &
                                                  eptStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 3 &
                                                  (at - part.list[, "last.active.time"]) <= 2), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, main partnership, and last active between 21 and 41 days ago
  part.listept1.main.med <- part.list[which((at - eptEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA)) &
                                                eptStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 1 &
                                                (at - part.list[, "last.active.time"]) >= 3 &
                                                (at - part.list[, "last.active.time"]) < 6), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, casual partnership, and last active between 21 and 41 days ago
  part.listept1.casl.med <- part.list[which((at - eptEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA)) &
                                                eptStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 2 &
                                                (at - part.list[, "last.active.time"]) >= 3 &
                                                (at - part.list[, "last.active.time"]) < 6), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, one-off partnership, and last active between 21 and 41 days ago
  part.listept1.inst.med <- part.list[which((at - eptEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                     (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA)) &
                                                eptStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 3 &
                                                (at - part.list[, "last.active.time"]) >= 3 &
                                                (at - part.list[, "last.active.time"]) < 6), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, main partnership, and last active between 42 and 60 days ago
  part.listept1.main.long <- part.list[which((at - eptEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                 ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA)) &
                                                 eptStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 1 &
                                                 (at - part.list[, "last.active.time"]) >= 6 &
                                                 (at - part.list[, "last.active.time"]) <= 9), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, casual partnership, and last active between 42 and 60 days ago
  part.listept1.casl.long <- part.list[which((at - eptEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                 ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA)) &
                                                 eptStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 2 &
                                                 (at - part.list[, "last.active.time"]) >= 6 &
                                                 (at - part.list[, "last.active.time"]) <= 9), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, one-off partnership, and last active between 42 and 60 days ago
  part.listept1.inst.long <- part.list[which((at - eptEligdate[part.list[, "uid1"]]) <= ept.risk.int &
                                                 ((rGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uGC.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (rCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uCT.tx[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (rGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uGC.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (rCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA) |
                                                      (uCT.tx.prep[part.list[, "uid2"]]) %in% c(0, NA)) &
                                                 eptStat[part.list[, "uid1"]] == 1 & part.list[, "ptype"] == 3 &
                                                 (at - part.list[, "last.active.time"]) >= 6 &
                                                 (at - part.list[, "last.active.time"]) <= 9), , drop = FALSE]

  idspartlistsept1.main.short <- part.listept1.main.short[, "uid2"]
  idspartlistsept1.casl.short <- part.listept1.casl.short[, "uid2"]
  idspartlistsept1.inst.short <- part.listept1.inst.short[, "uid2"]

  idspartlistsept1.main.med <- part.listept1.main.med[, "uid2"]
  idspartlistsept1.casl.med <- part.listept1.casl.med[, "uid2"]
  idspartlistsept1.inst.med <- part.listept1.inst.med[, "uid2"]

  idspartlistsept1.main.long <- part.listept1.main.long[, "uid2"]
  idspartlistsept1.casl.long <- part.listept1.casl.long[, "uid2"]
  idspartlistsept1.inst.long <- part.listept1.inst.long[, "uid2"]


  # Partner 2 was recently treated, so partner 1 would be eligible to be
  # treated through EPT if not currently being treated for anything
  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, main partnership, and last active within last 20 days
  part.listept2.main.short <- part.list[which((at - eptEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                  ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA)) &
                                                  eptStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 1 &
                                                  (at - part.list[, "last.active.time"]) <= 2), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, casual partnership, and last active within last 20 days
  part.listept2.casl.short <- part.list[which((at - eptEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                  ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA)) &
                                                  eptStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 2 &
                                                  (at - part.list[, "last.active.time"]) <= 2), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, one-off partnership, and last active within last 20 days
  part.listept2.inst.short <- part.list[which((at - eptEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                  ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                       (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA)) &
                                                  eptStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 3 &
                                                  (at - part.list[, "last.active.time"]) <= 2), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, main partnership, and last active between 21 and 41 days ago
  part.listept2.main.med <- part.list[which((at - eptEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA)) &
                                                eptStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 1 &
                                                (at - part.list[, "last.active.time"]) >= 3 &
                                                (at - part.list[, "last.active.time"]) < 6), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, casual partnership, and last active between 21 and 41 days ago
  part.listept2.casl.med <- part.list[which((at - eptEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA)) &
                                                eptStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 2 &
                                                (at - part.list[, "last.active.time"]) >= 3 &
                                                (at - part.list[, "last.active.time"]) < 6), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, one-off partnership, and last active between 21 and 41 days ago
  part.listept2.inst.med <- part.list[which((at - eptEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                     (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA)) &
                                                eptStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 3 &
                                                (at - part.list[, "last.active.time"]) >= 3 &
                                                (at - part.list[, "last.active.time"]) < 6), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, main partnership, and last active between 42 and 63 days ago
  part.listept2.main.long <- part.list[which((at - eptEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                 ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA)) &
                                                 eptStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 1 &
                                                 (at - part.list[, "last.active.time"]) >= 6 &
                                                 (at - part.list[, "last.active.time"]) <= 9), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, casual partnership, and last active between 42 and 63 days ago
  part.listept2.casl.long <- part.list[which((at - eptEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                 ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA)) &
                                                 eptStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 2 &
                                                 (at - part.list[, "last.active.time"]) >= 6 &
                                                 (at - part.list[, "last.active.time"]) <= 9), , drop = FALSE]

  # Criteria: eligible within EPT risk interval, currently untreated, partner
  # received EPT, one-off partnership, and last active between 42 and 63 days ago
  part.listept2.inst.long <- part.list[which((at - eptEligdate[part.list[, "uid2"]]) <= ept.risk.int &
                                                 ((rGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uGC.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (rCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uCT.tx[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (rGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uGC.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (rCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA) |
                                                      (uCT.tx.prep[part.list[, "uid1"]]) %in% c(0, NA)) &
                                                 eptStat[part.list[, "uid2"]] == 1 & part.list[, "ptype"] == 3 &
                                                 (at - part.list[, "last.active.time"]) >= 6 &
                                                 (at - part.list[, "last.active.time"]) <= 9), , drop = FALSE]

  idspartlistsept2.main.short <- part.listept2.main.short[, "uid1"]
  idspartlistsept2.casl.short <- part.listept2.casl.short[, "uid1"]
  idspartlistsept2.inst.short <- part.listept2.inst.short[, "uid1"]

  idspartlistsept2.main.med <- part.listept2.main.med[, "uid1"]
  idspartlistsept2.casl.med <- part.listept2.casl.med[, "uid1"]
  idspartlistsept2.inst.med <- part.listept2.inst.med[, "uid1"]

  idspartlistsept2.main.long <- part.listept2.main.long[, "uid1"]
  idspartlistsept2.casl.long <- part.listept2.casl.long[, "uid1"]
  idspartlistsept2.inst.long <- part.listept2.inst.long[, "uid1"]

  # All EPT eligible IDs (partners of index)
  idsept <- unique(c(idspartlistsept1.main.short, idspartlistsept1.casl.short,
                     idspartlistsept1.inst.short, idspartlistsept2.main.short,
                     idspartlistsept2.casl.short, idspartlistsept2.inst.short,
                     idspartlistsept1.main.med, idspartlistsept1.casl.med,
                     idspartlistsept1.inst.med, idspartlistsept2.main.med,
                     idspartlistsept2.casl.med, idspartlistsept2.inst.med,
                     idspartlistsept1.main.long, idspartlistsept1.casl.long,
                     idspartlistsept1.inst.long, idspartlistsept2.main.long,
                     idspartlistsept2.casl.long, idspartlistsept2.inst.long))

  idsept.main.short <- unique(c(idspartlistsept1.main.short,
                                idspartlistsept2.main.short))
  idsept.casl.short <- unique(c(idspartlistsept1.casl.short,
                                idspartlistsept2.casl.short))
  idsept.inst.short <- unique(c(idspartlistsept1.inst.short,
                                idspartlistsept2.inst.short))

  idsept.main.med <- unique(c(idspartlistsept1.main.med,
                              idspartlistsept2.main.med))
  idsept.casl.med <- unique(c(idspartlistsept1.casl.med,
                              idspartlistsept2.casl.med))
  idsept.inst.med <- unique(c(idspartlistsept1.inst.med,
                              idspartlistsept2.inst.med))

  idsept.main.long <- unique(c(idspartlistsept1.main.long,
                               idspartlistsept2.main.long))
  idsept.casl.long <- unique(c(idspartlistsept1.casl.long,
                               idspartlistsept2.casl.long))
  idsept.inst.long <- unique(c(idspartlistsept1.inst.long,
                               idspartlistsept2.inst.long))

  ## Provision to and uptake of partners (to be treated at next time step)
  idsprovided.main.short <- idsept.main.short[which(rbinom(length(idsept.main.short), 1,
                                                           ept.provision.main * ept.provision.short.main.rr) == 1)]
  idsprovided.casl.short <- idsept.casl.short[which(rbinom(length(idsept.casl.short), 1,
                                                           ept.provision.casl * ept.provision.short.casl.rr) == 1)]
  idsprovided.inst.short <- idsept.inst.short[which(rbinom(length(idsept.inst.short), 1,
                                                           ept.provision.inst * ept.provision.short.inst.rr) == 1)]
  idsprovided.main.med <- idsept.main.med[which(rbinom(length(idsept.main.med), 1,
                                                       ept.provision.main * ept.provision.med.main.rr) == 1)]
  idsprovided.casl.med <- idsept.casl.med[which(rbinom(length(idsept.casl.med), 1,
                                                       ept.provision.casl * ept.provision.med.casl.rr) == 1)]
  idsprovided.inst.med <- idsept.inst.med[which(rbinom(length(idsept.inst.med), 1,
                                                       ept.provision.inst * ept.provision.med.inst.rr) == 1)]
  idsprovided.main.long <- idsept.main.long[which(rbinom(length(idsept.main.long), 1,
                                                         ept.provision.main * ept.provision.long.main.rr) == 1)]
  idsprovided.casl.long <- idsept.casl.long[which(rbinom(length(idsept.casl.long), 1,
                                                         ept.provision.casl * ept.provision.long.casl.rr) == 1)]
  idsprovided.inst.long <- idsept.inst.long[which(rbinom(length(idsept.inst.long), 1,
                                                         ept.provision.inst * ept.provision.long.inst.rr) == 1)]

  idsprovided_ept <- unique(c(idsprovided.main.short, idsprovided.casl.short,
                              idsprovided.inst.short, idsprovided.main.med,
                              idsprovided.casl.med, idsprovided.inst.med,
                              idsprovided.main.long, idsprovided.casl.long,
                              idsprovided.inst.long))

  idsprovided.main_ept <- unique(c(idsprovided.main.short, idsprovided.main.med,
                                   idsprovided.main.long))

  idsprovided.casl_ept <- unique(c(idsprovided.casl.short, idsprovided.casl.med,
                                   idsprovided.casl.long))

  idsprovided.inst_ept <- unique(c(idsprovided.inst.short, idsprovided.inst.med,
                                   idsprovided.inst.long))

  # Uptake
  idsept_tx.main <- idsprovided.main_ept[which(rbinom(length(idsprovided.main_ept), 1,
                                                      ept.uptake.main) == 1)]
  idsept_tx.casl <- idsprovided.casl_ept[which(rbinom(length(idsprovided.casl_ept), 1,
                                                      ept.uptake.casl) == 1)]
  idsept_tx.inst <- idsprovided.inst_ept[which(rbinom(length(idsprovided.inst_ept), 1,
                                                      ept.uptake.inst) == 1)]
  idsuptake_ept <- unique(c(idsept_tx.main, idsept_tx.casl, idsept_tx.inst))

  # Update attributes for partner who is assigned to uptake EPT
  eptTx[idsuptake_ept] <- 1

  # Update Epi
  if (at >= dat$param$riskh.ept.start) {
    dat$epi$eptpartelig[at] <- length(idsept)
    dat$epi$eptprovided[at] <- length(idsprovided_ept)
    dat$epi$eptTx[at] <- length(idsuptake_ept)
    dat$epi$eptprop_provided[at] <- dat$epi$eptprovided[at] / dat$epi$eptpartelig[at]
    dat$epi$eptprop_tx[at] <- dat$epi$eptTx[at] / dat$epi$eptprovided[at]
  }

  return(dat)
}
