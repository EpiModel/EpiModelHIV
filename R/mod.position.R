
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

  # Variables ---------------------------------------------------------------

  if (dat$control$save.dal == TRUE) {
    dal <- dat$temp$dal[[at]]
  } else {
    dal <- dat$temp$dal
  }

  role.class <- dat$attr$role.class
  ins.quot <- dat$attr$ins.quot
  race <- dat$attr$race

  vv.prob.iev.BB <- dat$param$vv.prob.iev.BB
  vv.prob.iev.BW <- dat$param$vv.prob.iev.BW
  vv.prob.iev.WW <- dat$param$vv.prob.iev.WW


  # Process -----------------------------------------------------------------

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
  vv.prob.iev <- (vv.race.combo == "BB") * vv.prob.iev.BB +
                 (vv.race.combo == "BW") * vv.prob.iev.BW +
                 (vv.race.combo == "WW") * vv.prob.iev.WW

  iev <- rbinom(length(vv), 1, vv.prob.iev)
  dal$ins[vv[iev == 1]] <- "B"
  vv.remaining <- vv[iev == 0]

  prob.inspos <- ins.quot[dal$pos[vv.remaining]] /
                 (ins.quot[dal$pos[vv.remaining]] + ins.quot[dal$neg[vv.remaining]])
  inspos <- rbinom(length(vv.remaining), 1, prob.inspos)
  dal$ins[vv.remaining[inspos == 1]] <- "P"
  dal$ins[vv.remaining[inspos == 0]] <- "N"


  # Output ------------------------------------------------------------------
  if (dat$control$save.dal == TRUE) {
    dat$temp$dal[[at]] <- dal
  } else {
    dat$temp$dal <- dal
  }



  return(dat)
}
