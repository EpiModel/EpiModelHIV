
#' @title Risk History Module
#'
#' @description Module function to track the risk history of uninfected persons
#'              for purpose of intervention targeting.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_msm <- function(dat, at) {

  # if (at < dat$param$riskh.stitest.start) {
  #   return(dat)
  # }
    
  if (at < min(dat$param$riskh.start, dat$param$riskh.stitest.start, dat$param$riskh.ept.start)) {
        return(dat)
  }

  ## Attributes
  uid <- dat$attr$uid
  active <- dat$attr$active
  dx <- dat$attr$diag.status
  since.test <- at - dat$attr$last.neg.test
  rGC.tx <- dat$attr$rGC.tx
  uGC.tx <- dat$attr$uGC.tx
  rCT.tx <- dat$attr$rCT.tx
  uCT.tx <- dat$attr$uCT.tx
  syph.tx <- dat$attr$syph.tx
  rGC.tx.prep <- dat$attr$rGC.tx.prep
  uGC.tx.prep <- dat$attr$uGC.tx.prep
  rCT.tx.prep <- dat$attr$rCT.tx.prep
  uCT.tx.prep <- dat$attr$uCT.tx.prep
  syph.tx.prep <- dat$attr$syph.tx.prep
  sexactive <- dat$attr$sexactive
  sexnewedge <- dat$attr$sexnewedge
  tt.traj.ct <- dat$attr$tt.traj.ct
  tt.traj.gc <- dat$attr$tt.traj.gc
  tt.traj.syph <- dat$attr$tt.traj.syph
  eptStat <- dat$attr$eptStat
  eptEligdate <- dat$attr$eptEligdate
  eptEligTx <- dat$attr$eptEligTx
  eptStartTime <- dat$attr$eptStartTime
  eptTx <- dat$attr$eptTx

  ## Parameters
  time.unit <- dat$param$time.unit
  stitest.active.int <- dat$param$stitest.active.int
  sti.highrisktest.int <- dat$param$sti.highrisktest.int
  ept.risk.int <- dat$param$ept.risk.int
  ept.provision.main <- dat$param$ept.provision.partner.main
  ept.provision.casl <- dat$param$ept.provision.partner.casl
  ept.provision.inst <- dat$param$ept.provision.partner.inst
  ept.uptake.main <- dat$param$ept.uptake.partner.main
  ept.uptake.casl <- dat$param$ept.uptake.partner.casl
  ept.uptake.inst <- dat$param$ept.uptake.partner.inst

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
    dat$attr$stitest.ind.active <- rep(NA, length(uid))
    dat$attr$stitest.ind.recentpartners <- rep(NA, length(uid))
    dat$attr$stitest.ind.sti <- rep(NA, length(uid))
    dat$attr$stitest.ind.newpartners <- rep(NA, length(uid))
    dat$attr$stitest.ind.concurrpartner <- rep(NA, length(uid))
    dat$attr$stitest.ind.partnersti <- rep(NA, length(uid))
    dat$attr$stitest.ind.uai.nmain <- rep(NA, length(uid))
    dat$attr$stitest.ind.uai.any <- rep(NA, length(uid))
  }

  ## Degree ##
  main.deg <- get_degree(dat$el[[1]])
  casl.deg <- get_degree(dat$el[[2]])
  inst.deg <- get_degree(dat$el[[3]])


  ## Preconditions ##

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
  
  
  ############################
  ## STI Testing Conditions ##
  ############################
  
  part.list <- dat$temp$part.list
 
  # Sexually active - annual testing for syphilis, CT, GC
  idsactive <- which((at - sexactive) <= stitest.active.int)
  dat$attr$stitest.ind.active[idsactive] <- at
  
  # High-risk: CDC definition of increased risk for women - Add these indications to top of this module
  
  ### Have a new sex partner in last x months
  idsnewpartners <- which((at - sexnewedge) <= sti.highrisktest.int)
  dat$attr$stitest.ind.newpartners[idsnewpartners] <- at
  
  ### Multiple sex partners
  #	Have more than one sex partner in last x months
  part.count <- as.data.frame(table(part.list[, 1:2]))
  idspartlist <- which(uid %in% part.list[, 1:2])
  idsnotpartlist <- which(!(uid %in% part.list[, 1:2]))
  dat$attr$recentpartners[idspartlist] <- part.count[idspartlist, 2]
  idsrecentpartners <- which(dat$attr$recentpartners > 1)
  dat$attr$stitest.ind.recentpartners[idsrecentpartners] <- at

  # Reset # of partners to 0 for those no longer in part.list
  dat$attr$recentpartners[idsnotpartlist] <- 0
  
  ### Had/have a sex partner with concurrent partners
  # Partner 1 has multiple partners, Partner 2 indicated
  part.listmult1 <- part.list[which((dat$attr$recentpartners[part.list[, 1]]) > 1) , , drop = FALSE]
  # Partner 2 indicated
  idspartlistmult1 <- part.listmult1[, 2]
  
  # Partner 2 has multiple partners, so partner 1 is indicated
  part.listmult2 <- part.list[which((dat$attr$recentpartners[part.list[, 2]]) > 1), , drop = FALSE]
  # Partner 1 indicated
  idspartlistmult2 <- part.listmult2[, 1]
  
  # Combine into one list for indication
  idspartmult <- unique(c(idspartlistmult1, idspartlistmult2))
  idspartmult <-  uid[which(uid %in% idspartmult)]
  dat$attr$stitest.ind.concurrpartner[which(uid %in% idspartmult)] <- at

  
  ### Have a sex partner who has a treated sexually transmitted infection in the last interval
  # Partner 1 has a STI, Partner 2 indicated
  part.liststi1 <- part.list[which((at - dat$attr$last.tx.time.rct[part.list[, 1]]) <= sti.highrisktest.int |
                                       (at - dat$attr$last.tx.time.uct[part.list[, 1]]) <= sti.highrisktest.int |
                                       (at - dat$attr$last.tx.time.rgc[part.list[, 1]]) <= sti.highrisktest.int |
                                       (at - dat$attr$last.tx.time.ugc[part.list[, 1]]) <= sti.highrisktest.int |
                                       (at - dat$attr$last.tx.time.syph[part.list[, 1]]) <= sti.highrisktest.int), , drop = FALSE]
  # Partner 2 indicated
  idspartliststi1 <- part.liststi1[, 2]
  
  # Partner 2 has a STI, so partner 1 is indicated
  part.liststi2 <- part.list[which((at - dat$attr$last.tx.time.rct[part.list[, 2]]) <= sti.highrisktest.int |
                                       (at - dat$attr$last.tx.time.uct[part.list[, 2]]) <= sti.highrisktest.int |
                                       (at - dat$attr$last.tx.time.rgc[part.list[, 2]]) <= sti.highrisktest.int |
                                       (at - dat$attr$last.tx.time.ugc[part.list[, 2]]) <= sti.highrisktest.int |
                                       (at - dat$attr$last.tx.time.syph[part.list[, 2]]) <= sti.highrisktest.int), , drop = FALSE]
  # Partner 1 indicated
  idspartliststi2 <- part.liststi2[, 1]
  
  # Combine into one list for indication
  idspartsti <- unique(c(idspartliststi1, idspartliststi2))
  idspartsti <-  uid[which(uid %in% idspartsti)]
  dat$attr$stitest.ind.partnersti[which(uid %in% idspartsti)] <- at
  
  ### Inconsistent condom use among persons who are not in mutually monogamous relationships - includes concordant HIV
  # Using PrEP logic for UAI in non-main relationships: Could there be a closer approximation?
  uai.nmain <- unique(c(el$p1[el$uai > 0 & el$ptype %in% 2:3],
                        el$p2[el$uai > 0 & el$ptype %in% 2:3]))
  dat$attr$stitest.ind.uai.nmain[uai.nmain] <- at
  
  ## Any CAI
  uai.any <- unique(c(el$p1[el$uai > 0], el$p2[el$uai > 0]))
  dat$attr$stitest.ind.uai.any[uai.any] <- at
  
  
  ### Previous or coexisting STIs (treated) in a time interval
  idsSTI <- which((at - dat$attr$last.tx.time.syph) <= sti.highrisktest.int | (at - dat$attr$last.tx.time.rgc) <= sti.highrisktest.int |
                      (at - dat$attr$last.tx.time.ugc) <= sti.highrisktest.int | (at - dat$attr$last.tx.time.rct) <= sti.highrisktest.int | 
                      (at - dat$attr$last.tx.time.uct) <= sti.highrisktest.int)
  dat$attr$stitest.ind.sti[idsSTI] <- at
  
  # Update attributes
  dat$attr$tt.traj.ct <- tt.traj.ct
  dat$attr$tt.traj.gc <- tt.traj.gc
  dat$attr$tt.traj.syph <- tt.traj.syph
  
  # Update epi to show prevalence of STI testing indications
  if (at >= dat$param$riskh.stitest.start) {
      dat$epi$stiactiveind[at] <- length(idsactive) / dat$epi$num[at]
      dat$epi$newpartner[at] <- length(idsnewpartners) / dat$epi$num[at]
      dat$epi$recentpartners[at] <- length(idsrecentpartners) / dat$epi$num[at]
      dat$epi$concurrpart[at] <- length(which(uid %in% idspartmult)) / dat$epi$num[at]
      dat$epi$partnersti[at] <- length(which(uid %in% idspartsti)) / dat$epi$num[at]
      dat$epi$uai.nmain[at] <- length(uai.nmain) / dat$epi$num[at]
      dat$epi$uai.any[at] <- length(uai.any) / dat$epi$num[at]
      dat$epi$recentSTI[at] <- length(idsSTI) / dat$epi$num[at]
  }
  
  
  #####################
  ## EPT  Conditions ##
  #####################
  
  ## Eligibility of partners---------------------------------------------------------------
  
  # Subset partner list to those active within an EPT interval - last active date within 60 days
  part.list <- part.list[which((at - (part.list[, 5]) <= ept.risk.int)), , drop = FALSE]
  
  # Subset partner list to alive
  part.list <- part.list[which(active[part.list[, 1]] == 1 & active[part.list[, 2]] == 1), , drop = FALSE]
  
  # Partner 1 was recently treated, so partner 2 would be eligible to be treated through EPT if not currently being treated for anything
  part.listept1.main <- part.list[which((at - eptEligdate[part.list[, 1]]) <= ept.risk.int & 
                                            ((rGC.tx[part.list[, 2]]) %in% c(0, NA) | (uGC.tx[part.list[, 2]]) %in% c(0, NA) |
                                            (rCT.tx[part.list[, 2]]) %in% c(0, NA) | (uCT.tx[part.list[, 2]]) %in% c(0, NA) |
                                            (rGC.tx.prep[part.list[, 2]]) %in% c(0, NA) | (uGC.tx.prep[part.list[, 2]]) %in% c(0, NA) |
                                            (rCT.tx.prep[part.list[, 2]]) %in% c(0, NA) | (uCT.tx.prep[part.list[, 2]]) %in% c(0, NA)) &
                                            eptStat[part.list[, 1]] == 1 & part.list[, 3] == 1), , drop = FALSE]
  
  part.listept1.casl <- part.list[which((at - eptEligdate[part.list[, 1]]) <= ept.risk.int & 
                                            ((rGC.tx[part.list[, 2]]) %in% c(0, NA) | (uGC.tx[part.list[, 2]]) %in% c(0, NA) |
                                            (rCT.tx[part.list[, 2]]) %in% c(0, NA) | (uCT.tx[part.list[, 2]]) %in% c(0, NA) |
                                            (rGC.tx.prep[part.list[, 2]]) %in% c(0, NA) | (uGC.tx.prep[part.list[, 2]]) %in% c(0, NA) |
                                            (rCT.tx.prep[part.list[, 2]]) %in% c(0, NA) | (uCT.tx.prep[part.list[, 2]]) %in% c(0, NA)) &
                                            eptStat[part.list[, 1]] == 1 & part.list[, 3] == 2), , drop = FALSE]
  
  part.listept1.inst <- part.list[which((at - eptEligdate[part.list[, 1]]) <= ept.risk.int & 
                                            ((rGC.tx[part.list[, 2]]) %in% c(0, NA) | (uGC.tx[part.list[, 2]]) %in% c(0, NA) |
                                            (rCT.tx[part.list[, 2]]) %in% c(0, NA) | (uCT.tx[part.list[, 2]]) %in% c(0, NA) |
                                            (rGC.tx.prep[part.list[, 2]]) %in% c(0, NA) | (uGC.tx.prep[part.list[, 2]]) %in% c(0, NA) |
                                            (rCT.tx.prep[part.list[, 2]]) %in% c(0, NA) | (uCT.tx.prep[part.list[, 2]]) %in% c(0, NA)) &
                                            eptStat[part.list[, 1]] == 1 & part.list[, 3] == 3), , drop = FALSE]
  
  rGC.tx.prep <- dat$attr$rGC.tx.prep
  uGC.tx.prep <- dat$attr$uGC.tx.prep
  rCT.tx.prep <- dat$attr$rCT.tx.prep
  uCT.tx.prep <- dat$attr$uCT.tx.prep
  syph.tx.prep <- dat$attr$syph.tx.prep
  
  idspartlistsept1.main <- part.listept1.main[, 2]
  idspartlistsept1.casl <- part.listept1.casl[, 2]
  idspartlistsept1.inst <- part.listept1.inst[, 2]
  
  
  # Partner 2 was recently treated, so partner 1 would be eligible to be treated through EPT if not currently being treated for anything
  part.listept2.main <- part.list[which((at - eptEligdate[part.list[, 2]]) <= ept.risk.int & 
                                            ((rGC.tx[part.list[, 1]]) %in% c(0, NA) | (uGC.tx[part.list[, 1]]) %in% c(0, NA) |
                                            (rCT.tx[part.list[, 1]]) %in% c(0, NA) | (uCT.tx[part.list[, 1]]) %in% c(0, NA) |
                                            (rGC.tx.prep[part.list[, 1]]) %in% c(0, NA) | (uGC.tx.prep[part.list[, 1]]) %in% c(0, NA) |
                                            (rCT.tx.prep[part.list[, 1]]) %in% c(0, NA) | (uCT.tx.prep[part.list[, 1]]) %in% c(0, NA)) &
                                            eptStat[part.list[, 2]] == 1 & part.list[, 3] == 1), , drop = FALSE]
  
  part.listept2.casl <- part.list[which((at - eptEligdate[part.list[, 2]]) <= ept.risk.int & 
                                            ((rGC.tx[part.list[, 1]]) %in% c(0, NA) | (uGC.tx[part.list[, 1]]) %in% c(0, NA) |
                                            (rCT.tx[part.list[, 1]]) %in% c(0, NA) | (uCT.tx[part.list[, 1]]) %in% c(0, NA)|
                                            (rGC.tx.prep[part.list[, 1]]) %in% c(0, NA) | (uGC.tx.prep[part.list[, 1]]) %in% c(0, NA) |
                                            (rCT.tx.prep[part.list[, 1]]) %in% c(0, NA) | (uCT.tx.prep[part.list[, 1]]) %in% c(0, NA)) &
                                            eptStat[part.list[, 2]] == 1 & part.list[, 3] == 2), , drop = FALSE]
  
  part.listept2.inst <- part.list[which((at - eptEligdate[part.list[, 2]]) <= ept.risk.int & 
                                            ((rGC.tx[part.list[, 1]]) %in% c(0, NA) | (uGC.tx[part.list[, 1]]) %in% c(0, NA) |
                                            (rCT.tx[part.list[, 1]]) %in% c(0, NA) | (uCT.tx[part.list[, 1]]) %in% c(0, NA) |
                                            (rGC.tx.prep[part.list[, 1]]) %in% c(0, NA) | (uGC.tx.prep[part.list[, 1]]) %in% c(0, NA) |
                                            (rCT.tx.prep[part.list[, 1]]) %in% c(0, NA) | (uCT.tx.prep[part.list[, 1]]) %in% c(0, NA)) &
                                            eptStat[part.list[, 2]] == 1 & part.list[, 3] == 3), , drop = FALSE]
  
  idspartlistsept2.main <- part.listept2.main[, 1]
  idspartlistsept2.casl <- part.listept2.casl[, 1]
  idspartlistsept2.inst <- part.listept2.inst[, 1]
  
  # All EPT eligible IDs (partners of index)
  idsept <- unique(c(idspartlistsept1.main, idspartlistsept1.casl, idspartlistsept1.inst, 
                     idspartlistsept2.main, idspartlistsept2.casl, idspartlistsept2.inst))
  idsept.main <- unique(c(idspartlistsept1.main, idspartlistsept2.main))
  idsept.casl <- unique(c(idspartlistsept1.casl, idspartlistsept2.casl))
  idsept.inst <- unique(c(idspartlistsept1.inst, idspartlistsept2.inst))
  
  ## Provision to and uptake of partners (to be treated at next time step)
  idsprovided.main <- idsept.main[which(rbinom(length(idsept.main), 1,
                                               ept.provision.main) == 1)]
  idsprovided.casl <- idsept.casl[which(rbinom(length(idsept.casl), 1,
                                               ept.provision.casl) == 1)]
  idsprovided.inst <- idsept.inst[which(rbinom(length(idsept.inst), 1,
                                               ept.provision.inst) == 1)]
  idsprovided_ept <- unique(c(idsprovided.main, idsprovided.casl, idsprovided.inst))
  
  # Uptake
  idsept_tx.main <- idsprovided.main[which(rbinom(length(idsprovided.main), 1,
                                                  ept.uptake.main) == 1)]
  idsept_tx.casl <- idsprovided.casl[which(rbinom(length(idsprovided.casl), 1,
                                                  ept.uptake.casl) == 1)]
  idsept_tx.inst <- idsprovided.inst[which(rbinom(length(idsprovided.inst), 1,
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
