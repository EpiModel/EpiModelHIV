
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
  syph.incub.tx <- dat$attr$syph.incub.tx
  syph.prim.tx <- dat$attr$syph.prim.tx
  syph.seco.tx <- dat$attr$syph.seco.tx
  syph.earlat.tx <- dat$attr$syph.earlat.tx
  syph.latelat.tx <- dat$attr$syph.latelat.tx
  syph.tert.tx <- dat$attr$syph.tert.tx

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
  pers.deg <- get_degree(dat$el[[2]])
  inst.deg <- get_degree(dat$el[[3]])


  # Indications -------------------------------------------------------------

  # Any UAI
  uai.any <- unique(c(el2$p1[el2$uai > 0],
                      el2$p2[el2$uai > 0]))

  # Monogamous partnerships: 1-sided
  tot.deg <- main.deg + pers.deg + inst.deg
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
                   rCT.tx == 1 | uCT.tx == 1 | syph.incub.tx == 1 |
                   syph.prim.tx == 1 | syph.seco.tx == 1 |
                   syph.earlat.tx == 1 | syph.latelat.tx == 1 |
                   syph.tert.tx == 1)
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

  ## Parameters
  partnercutoff <- dat$param$partnercutoff
  stitest.active.int <- dat$param$stitest.active.int

  ## Attributes
  uid <- dat$attr$uid
  race <- dat$attr$race
  stitestind1 <- dat$attr$time.last.sex

  # Indications -------------------------------------------------------------

  part.list <- dat$temp$part.list

  ### Lower risk - sexually active in the last year
  idsactive <- which((at - stitestind1) <= stitest.active.int)
  idsnotactive <- setdiff(which(race %in% c("B","W")), idsactive)
  # these are relative ids of nodes in partner list

  ### High-risk: Have more than one sex partner in last x months
  # Reset # of partners at each time step- length of "recent" interval is drawn from interval of partner list lookback
  dat$attr$recentpartners <- rep(0, length(dat$attr$race))

  if (at >= dat$param$stitest.start | at >= dat$param$ept.start) {

  # For those who had partners, calculate # of occurrences in partner list
  part.count <- as.data.frame(table(part.list[, c("uid1", "uid2")]))

  if (nrow(part.count) > 1) {

    # Calculate # of recent partners: 0 for those not in part list, update numbers for only actives in part list
    dat$attr$recentpartners[which(part.count[, "Var1"] %in% uid)] <- part.count[which(part.count[, "Var1"] %in% uid), 2]

    # Choose those who have had more than X partners in last x months
    idsrecentpartners <- which(dat$attr$recentpartners > partnercutoff)
    idsnotrecentpartners <- setdiff(which(race %in% c("B","W")), idsrecentpartners)

    dat$attr$stitest.ind.recentpartners[idsrecentpartners] <- 1
    dat$attr$stitest.ind.recentpartners[idsnotrecentpartners] <- 0
    dat$epi$recentpartners[at] <- length(idsrecentpartners)

  }

  }
  ### Update STI indication attributes
  dat$attr$stitest.ind.active[idsactive] <- 1
  dat$attr$stitest.ind.active[idsnotactive] <- 0

  dat$epi$stiactiveind[at] <- length(idsactive)


  return(dat)
}
