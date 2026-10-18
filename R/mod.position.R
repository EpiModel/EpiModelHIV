
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

  al <- dat$temp$al
  if (nrow(al) == 0) {
    return(dat)
  }

  # Attributes
  role.class <- dat$attr$role.class
  ins.quot <- dat$attr$ins.quot

  # Parameters

  ## Process
  p1.role.class <- role.class[al[, "p1"]]
  p2.role.class <- role.class[al[, "p2"]]

  ins <- rep(NA, length(p1.role.class))
  ins[which(p1.role.class == 0)] <- 1
  ins[which(p1.role.class == 1)] <- 0
  ins[which(p2.role.class == 0)] <- 0
  ins[which(p2.role.class == 1)] <- 1

  # Versatile MSM
  vv <- which(p1.role.class == 2 & p2.role.class == 2)
  p1.ins.prob <- ins.quot[al[, 1][vv]] /
                 (ins.quot[al[, 1][vv]] + ins.quot[al[, 2][vv]])
  p1.ins <- rbinom(length(vv), 1, p1.ins.prob)
  ins[vv[p1.ins == 1]] <- 1
  ins[vv[p1.ins == 0]] <- 0

  ## Output
  dat$temp$al <- cbind(al, ins)

  return(dat)
}
