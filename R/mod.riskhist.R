
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

  ## Edgelist
  el <- dat$temp$el

  # Remove concordant positive edges
  el2 <- el[el$st2 == 0, ]


  ## Any UAI
  uai.any <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0], el2$p2[el2$uai > 0]))
  uai.main <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$type == "main"],
                       el2$p2[el2$uai > 0 & el2$type == "main"]))
  uai.casl <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$type == "pers"],
                       el2$p2[el2$uai > 0 & el2$type == "pers"]))
  uai.inst <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$type == "inst"],
                       el2$p2[el2$uai > 0 & el2$type == "inst"]))
  uai.nmain <- union(uai.casl, uai.inst)


  # UAI within known serodiscordant partnerships
  dlist <- dat$temp$discl.list
  discl <- sapply(1:nrow(el2), function(x) {
    length(intersect(which(uid[el2$p1[x]] == dlist$pos),
                     which(uid[el2$p2[x]] == dlist$neg))) != 0
  })

  uai.serodis <- el2$p2[discl == TRUE & el2$uai > 0]
  uai.serodis.main <- el2$p2[discl == TRUE & el2$uai > 0 & el2$type == "main"]
  uai.serodis.casl <- el2$p2[discl == TRUE & el2$uai > 0 & el2$type == "pers"]
  uai.serodis.inst <- el2$p2[discl == TRUE & el2$uai > 0 & el2$type == "inst"]
  uai.serodis.nmain <- union(uai.serodis.casl, uai.serodis.inst)


  ## Summary statistics
  dat$attr$uai.last[uai.any] <- at
  dat$attr$uai.main[uai.main] <- at
  dat$attr$uai.casl.last[uai.casl] <- at
  dat$attr$uai.inst.last[uai.inst] <- at
  dat$attr$uai.nmain.last[uai.nmain] <- at
  dat$attr$uai.sd[uai.serodis] <- at
  dat$attr$uai.sd.main[uai.serodis.main] <- at
  dat$attr$uai.sd.casl[uai.serodis.casl] <- at
  dat$attr$uai.sd.inst[uai.serodis.inst] <- at
  dat$attr$uai.sd.nmain[uai.serodis.nmain] <- at

  return(dat)
}
