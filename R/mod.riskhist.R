
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
    
  if (at < min(dat$param$riskh.start, dat$param$riskh.stitest.start)) {
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
  sexactive <- dat$attr$sexactive
  sexnewedge <- dat$attr$sexnewedge
  tt.traj.ct <- dat$attr$tt.traj.ct
  tt.traj.gc <- dat$attr$tt.traj.gc
  tt.traj.syph <- dat$attr$tt.traj.syph

  ## Parameters
  time.unit <- dat$param$time.unit
  stitest.active.int <- dat$param$stitest.active.int
  sti.highrisktest.int <- dat$param$sti.highrisktest.int
  ept.risk.int <- dat$param$ept.risk.int

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

  ## Condition 4, any STI diagnosis
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
    

  # Reset attributes
  dat$attr$tt.traj.ct <- tt.traj.ct
  dat$attr$tt.traj.gc <- tt.traj.gc
  dat$attr$tt.traj.syph <- tt.traj.syph
  
  #####################
  ## EPT  Conditions ##
  #####################
  
  
  return(dat)
}
