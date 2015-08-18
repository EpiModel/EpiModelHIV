
#' @title Risk History Module
#'
#' @description Module function to track the risk history of uninfected persons
#'              for purpose of intervention targeting.
#'
#' @inheritParams aging.mard
#'
#' @keywords module
#' @export
#'
riskhist.mard <- function(dat, at) {

  ## Attributes
  uid <- dat$attr$uid

  ## Parameters
  prep.risk.int <- ceiling(dat$param$prep.risk.int)


  ## Edgelist
  el <- dat$temp$el

  # Remove concordant positive edges
  el2 <- el[el$st2 == 0, ]


  ## Any UAI
  uai.any <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0], el2$p2[el2$uai > 0]))
  uai.main <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$type == "main"],
                       el2$p2[el2$uai > 0 & el2$type == "main"]))
  uai.nmain <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$type %in% c("pers", "inst")],
                        el2$p2[el2$uai > 0 & el2$type %in% c("pers", "inst")]))


  ## Degree > 1
  main.deg <- as.numeric(summary(dat$nw$m ~ sociality(base = 0), at = at))
  casl.deg <- as.numeric(summary(dat$nw$p ~ sociality(base = 0), at = at))
  inst.deg <- as.numeric(summary(dat$nw$i ~ sociality(base = 0), at = at))


  # UAI within known serodiscordant partnerships
  dlist <- dat$temp$discl.list
  discl <- sapply(1:nrow(el2), function(x) {
    length(intersect(which(uid[el2$p1[x]] == dlist$pos),
                     which(uid[el2$p2[x]] == dlist$neg))) != 0
  })
  uai.sd <- el2$p2[discl == TRUE & el2$uai > 0]
  uai.sd.main <- el2$p2[discl == TRUE & el2$uai > 0 & el2$type == "main"]
  uai.sd.nmain <- el2$p2[discl == TRUE & el2$uai > 0 & el2$type %in% c("pers", "inst")]


  ## Truncate riskh matrices
  for (i in 1:length(dat$riskh)) {
    dat$riskh[[i]] <- dat$riskh[[i]][, -1]
    dat$riskh[[i]] <- cbind(dat$riskh[[i]], rep(NA, nrow(dat$riskh[[i]])))
  }


  ## Update attributes
  dat$riskh$mdeg[, prep.risk.int] <- main.deg
  dat$riskh$pdeg[, prep.risk.int] <- casl.deg
  dat$riskh$ideg[, prep.risk.int] <- inst.deg

  dat$riskh$uai[, prep.risk.int] <- 0
  dat$riskh$uai[uai.any, prep.risk.int] <- 1

  dat$riskh$uai.main[, prep.risk.int] <- 0
  dat$riskh$uai.main[uai.main, prep.risk.int] <- 1

  dat$riskh$uai.nmain[, prep.risk.int] <- 0
  dat$riskh$uai.nmain[uai.nmain, prep.risk.int] <- 1

  dat$riskh$uai.sd[, prep.risk.int] <- 0
  dat$riskh$uai.sd[uai.sd, prep.risk.int] <- 1

  dat$riskh$uai.sd.main[, prep.risk.int] <- 0
  dat$riskh$uai.sd.main[uai.sd.main, prep.risk.int] <- 1

  dat$riskh$uai.sd.nmain[, prep.risk.int] <- 0
  dat$riskh$uai.sd.nmain[uai.sd.nmain, prep.risk.int] <- 1

  ## Summary risk criteria for PrEP


  return(dat)
}
