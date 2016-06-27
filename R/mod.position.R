
#' @title Position Module
#'
#' @description Module function for establishing sexual role or position in each
#'              act on the discordant edgelist.
#'
#' @inheritParams aging_msm
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
#' @keywords module msm
#' 
#' @export
#'
position_msm <- function(dat, at) {

  ## Variables
  al <- dat$temp$al
  if (nrow(al) == 0) {
    return(dat)
  }

  status <- dat$attr$status
  dal <- al[which(status[al[, 1]] == 1 & status[al[, 2]] == 0), ]
  dat$temp$al <- NULL

  role.class <- dat$attr$role.class
  ins.quot <- dat$attr$ins.quot
  race <- dat$attr$race

  vv.iev.BB.prob <- dat$param$vv.iev.BB.prob
  vv.iev.BW.prob <- dat$param$vv.iev.BW.prob
  vv.iev.WW.prob <- dat$param$vv.iev.WW.prob


  ## Process
  pos.role.class <- role.class[dal[, 1]]
  neg.role.class <- role.class[dal[, 2]]

  ins <- rep(NA, length(pos.role.class))
  ins[which(pos.role.class == "I")] <- 1  # "P"
  ins[which(pos.role.class == "R")] <- 0  # "N"
  ins[which(neg.role.class == "I")] <- 0  # "N"
  ins[which(neg.role.class == "R")] <- 1  # "P"

  vv <- which(pos.role.class == "V" & neg.role.class == "V")
  vv.race.combo <- paste0(race[dal[, 1]][vv], race[dal[, 2]][vv])
  vv.race.combo[vv.race.combo == "WB"] <- "BW"
  vv.iev.prob <- (vv.race.combo == "BB") * vv.iev.BB.prob +
                 (vv.race.combo == "BW") * vv.iev.BW.prob +
                 (vv.race.combo == "WW") * vv.iev.WW.prob

  iev <- rbinom(length(vv), 1, vv.iev.prob)
  ins[vv[iev == 1]] <- 2 # "B"
  vv.remaining <- vv[iev == 0]

  inspos.prob <- ins.quot[dal[, 1][vv.remaining]] /
                 (ins.quot[dal[, 1][vv.remaining]] + ins.quot[dal[, 2][vv.remaining]])
  inspos <- rbinom(length(vv.remaining), 1, inspos.prob)
  ins[vv.remaining[inspos == 1]] <- 1  # "P"
  ins[vv.remaining[inspos == 0]] <- 0  # "N"


  ## Output
  dat$temp$dal <- cbind(dal, ins)

  return(dat)
}
