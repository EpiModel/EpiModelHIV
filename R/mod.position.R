
#' @title Position Module
#'
#' @description Module function for establishing sexual role or position in each
#'              act on the discordant edgelist.
#'
#' @inheritParams aging.mard
#'
#' @details
#' The sexual role within each act is determined by each nodes "role identity"
#' as exclusively receptive, exclusively insertive, or versatile. This function
#' determines whether the infected or the susceptible partner is the insertive
#' partner for that act. For the first two role identity types, that is
#' deterministic based on identity. For versatile-versatile pairs, this is
#' determined stochastically for each act.
#'
#' @return
#' This function returns the updated discordant edgelist with a \code{ins}
#' attribute for values of whether the infected node is insertive or the
#' susceptible node is insertive for that act.
#'
#' @keywords module
#' @export
#'
position.mard <- function(dat, at) {

  ## Variables
  el <- dat$temp$el

  role.class <- dat$attr$role.class
  ins.quot <- dat$attr$ins.quot
  race <- dat$attr$race

  vv.iev.BB.prob <- dat$param$vv.iev.BB.prob
  vv.iev.BW.prob <- dat$param$vv.iev.BW.prob
  vv.iev.WW.prob <- dat$param$vv.iev.WW.prob


  ## Process
  el$rc1 <- role.class[el$p1]
  el$rc2 <- role.class[el$p2]

  el$p1ins <- 0
  el$p1ins[el$rc1 == "I"] <- el$ai[el$rc1 == "I"]
  el$p1ins[el$rc1 == "V" & el$rc2 == "R"] <- el$ai[el$rc1 == "V" & el$rc2 == "R"]

  vv <- which(el$rc1 == "V" & el$rc2 == "V" & el$ai > 0)
  vv.race.combo <- paste0(race[el$p1[vv]], race[el$p2[vv]])
  vv.race.combo[vv.race.combo == "WB"] <- "BW"
  vv.iev.prob <- (vv.race.combo == "BB") * vv.iev.BB.prob +
                 (vv.race.combo == "BW") * vv.iev.BW.prob +
                 (vv.race.combo == "WW") * vv.iev.WW.prob

  iev <- rbinom(length(vv), el$ai[vv], vv.iev.prob)
  el$ai[vv[which(iev > 0)]] <- el$ai[vv[which(iev > 0)]] + iev[which(iev > 0)]
  el$p1ins[vv[which(iev > 0)]] <- rbinom(sum(iev > 0), el$ai[vv[which(iev > 0)]], 0.5)

  vv.remaining <- vv[iev == 0]

  p1.prob <- ins.quot[el$p1[vv.remaining]] /
             (ins.quot[el$p1[vv.remaining]] + ins.quot[el$p2[vv.remaining]])
  p1insall <- rbinom(length(vv.remaining), 1, p1.prob)
  el$p1ins[vv.remaining[p1insall == 1]] <- el$ai[vv.remaining[p1insall == 1]]

  ## Output
  dat$temp$el <- el

  return(dat)
}
