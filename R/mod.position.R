
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
  dal <- dat$temp$dal
  if (nrow(dal) == 0) {
    return(dat)
  }

  role.class <- dat$attr$role.class
  ins.quot <- dat$attr$ins.quot
  race <- dat$attr$race

  vv.iev.BB.prob <- dat$param$vv.iev.BB.prob
  vv.iev.BW.prob <- dat$param$vv.iev.BW.prob
  vv.iev.WW.prob <- dat$param$vv.iev.WW.prob


  ## Process
  pos.role.class <- role.class[dal$pos]
  neg.role.class <- role.class[dal$neg]

  dal$ins <- NA
  dal$ins[pos.role.class == "I"] <- "P"
  dal$ins[pos.role.class == "R"] <- "N"
  dal$ins[neg.role.class == "I"] <- "N"
  dal$ins[neg.role.class == "R"] <- "P"

  vv <- which(pos.role.class == "V" & neg.role.class == "V")
  vv.race.combo <- paste0(race[dal$pos[vv]], race[dal$neg[vv]])
  vv.race.combo[vv.race.combo == "WB"] <- "BW"
  vv.iev.prob <- (vv.race.combo == "BB") * vv.iev.BB.prob +
                 (vv.race.combo == "BW") * vv.iev.BW.prob +
                 (vv.race.combo == "WW") * vv.iev.WW.prob

  iev <- rbinom(length(vv), 1, vv.iev.prob)
  dal$ins[vv[iev == 1]] <- "B"
  vv.remaining <- vv[iev == 0]

  inspos.prob <- ins.quot[dal$pos[vv.remaining]] /
                 (ins.quot[dal$pos[vv.remaining]] + ins.quot[dal$neg[vv.remaining]])
  inspos <- rbinom(length(vv.remaining), 1, inspos.prob)
  dal$ins[vv.remaining[inspos == 1]] <- "P"
  dal$ins[vv.remaining[inspos == 0]] <- "N"


  ## Output
  dat$temp$dal <- dal

  return(dat)
}
