
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

  if (at < dat$param$riskh.start) {
    return(dat)
  }

  ## Attributes
  uid <- dat$attr$uid

  ## Parameters

  ## Edgelist, adds uai summation per partnership from act list
  pid <- NULL # For R CMD Check
  al <- as.data.frame(dat$temp$al)
  by_pid <- group_by(al, pid)
  uai <- summarise(by_pid, uai = sum(uai))[, 2]
  el <- as.data.frame(cbind(dat$temp$el, uai))

  # Remove concordant positive edges
  el2 <- el[el$st2 == 0, ]

  if (is.null(dat$attr$uai.mono2.nt.3mo)) {
    dat$attr$uai.mono2.nt.3mo <- rep(NA, length(uid))
    dat$attr$uai.mono1.nt.3mo <- rep(NA, length(uid))
    dat$attr$uai.mono2.nt.6mo <- rep(NA, length(uid))
    dat$attr$uai.mono1.nt.6mo <- rep(NA, length(uid))
    dat$attr$uai.nonmonog <- rep(NA, length(uid))
    dat$attr$uai.nmain <- rep(NA, length(uid))
    dat$attr$ai.sd.mc <- rep(NA, length(uid))
    dat$attr$uai.sd.mc <- rep(NA, length(uid))
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

  # Monogamous partnerships: 2-sided
  mono.2s <- tot.deg[el2$p1] == 1 & tot.deg[el2$p2] == 1
  ai.mono2 <- sort(unname(do.call("c", c(el2[mono.2s, 1:2]))))
  uai.mono2 <- intersect(ai.mono2, uai.any)

  # "Negative" partnerships
  tneg <- unique(c(el2$p1[el2$st1 == 0], el2$p2[el2$st1 == 0]))
  dx <- dat$attr$diag.status
  fneg <- unique(c(el2$p1[which(dx[el2$p1] == 0)], el2$p2[which(dx[el2$p1] == 0)]))
  all.neg <- c(tneg, fneg)
  since.test <- at - dat$attr$last.neg.test


  ## Condition 1a: UAI in 2-sided monogamous "negative" partnership,
  ##               partner not tested in past 3, 6 months
  uai.mono2.neg <- intersect(uai.mono2, all.neg)
  part.id2 <- c(el2[el2$p1 %in% uai.mono2.neg, 2], el2[el2$p2 %in% uai.mono2.neg, 1])
  not.tested.3mo <- since.test[part.id2] > (90/dat$param$time.unit)
  part.not.tested.3mo <- uai.mono2.neg[which(not.tested.3mo == TRUE)]
  dat$attr$uai.mono2.nt.3mo[part.not.tested.3mo] <- at

  not.tested.6mo <- since.test[part.id2] > (180/dat$param$time.unit)
  part.not.tested.6mo <- uai.mono2.neg[which(not.tested.6mo == TRUE)]
  dat$attr$uai.mono2.nt.6mo[part.not.tested.6mo] <- at


  ## Condition 1b: UAI in 1-sided "monogamous" "negative" partnership,
  ##               partner not tested in past 3, 6 months
  uai.mono1.neg <- intersect(uai.mono1, all.neg)
  part.id1 <- c(el2[el2$p1 %in% uai.mono1.neg, 2], el2[el2$p2 %in% uai.mono1.neg, 1])
  not.tested.3mo <- since.test[part.id1] > (90/dat$param$time.unit)
  part.not.tested.3mo <- uai.mono1.neg[which(not.tested.3mo == TRUE)]
  dat$attr$uai.mono1.nt.3mo[part.not.tested.3mo] <- at

  not.tested.6mo <- since.test[part.id1] > (180/dat$param$time.unit)
  part.not.tested.6mo <- uai.mono1.neg[which(not.tested.6mo == TRUE)]
  dat$attr$uai.mono1.nt.6mo[part.not.tested.6mo] <- at


  ## Condition 2a: UAI in non-monogamous partnerships
  el2.uai <- el2[el2$uai > 0, ]
  vec <- c(el2.uai[, 1], el2.uai[, 2])
  uai.nonmonog <- unique(vec[duplicated(vec)])
  dat$attr$uai.nonmonog[uai.nonmonog] <- at


  ## Condition 2b: UAI in non-main partnerships
  uai.nmain <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$ptype %in% 2:3],
                        el2$p2[el2$uai > 0 & el2$ptype %in% 2:3]))
  dat$attr$uai.nmain[uai.nmain] <- at


  ## Condition 3a: AI within known serodiscordant partnerships
  el2.cond3 <- el2[el2$st1 == 1 & el2$ptype %in% 1:2, ]

  # Disclosure
  discl.list <- dat$temp$discl.list
  disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
  delt.cdl <- uid[el2.cond3[, 1]] * 1e7 + uid[el2.cond3[, 2]]
  discl <- (delt.cdl %in% disclose.cdl)

  ai.sd.mc <- el2.cond3$p2[discl == TRUE]
  dat$attr$ai.sd.mc[ai.sd.mc] <- at


  ## Condition 3b: UAI within known serodiscordant partnerships
  uai.sd.mc <- el2.cond3$p2[discl == TRUE & el2.cond3$uai > 0]
  dat$attr$uai.sd.mc[uai.sd.mc] <- at

  return(dat)
}
