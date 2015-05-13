
#' @title Sexual Acts Module
#'
#' @description Module function for setting the number of sexual acts on the
#'              discordant edgelist.
#'
#' @inheritParams aging.mard
#'
#' @details
#' The number of acts at each time step is specified as a function of the race of
#' both members in a pair and the expected values within black-black, black-white,
#' and white-white combinations. For one-off partnerships, this is deterministically
#' set at 1, whereas for main and causal partnerships it is a stochastic draw
#' from a Poisson distribution. The number of total acts may further be modified
#' by the level of HIV viral suppression in an infected person.
#'
#' @return
#' This function returns the \code{dat} object with the updated discordant act
#' list (\code{dal}). Each element of \code{dal} is a data frame with the ids of the
#' discordant pair repeated the number of times they have AI.
#'
#' @keywords module
#' @export
#'
acts.mard <- function(dat, at) {

  for (type in c("main", "pers", "inst")) {
    # Variables ---------------------------------------------------------------

    # Attributes
    active <- dat$attr$active
    uid <- dat$attr$uid
    status <- dat$attr$status
    race <- dat$attr$race
    diag.status <- dat$attr$diag.status
    tx.status <- dat$attr$tx.status
    tt.traj <- dat$attr$tt.traj

    # Parameters
    if (type == "main") {
      base.ai.BB.rate <- dat$param$base.ai.main.BB.rate
      base.ai.BW.rate <- dat$param$base.ai.main.BW.rate
      base.ai.WW.rate <- dat$param$base.ai.main.WW.rate
      redux.ai.diag.rr <- dat$param$redux.ai.diag.main.rr
      redux.ai.discl.rr <- dat$param$redux.ai.discl.main.rr
      incr.ai.full.supp.rr <- dat$param$incr.ai.full.supp.main.rr
      incr.ai.part.supp.rr <- dat$param$incr.ai.part.supp.main.rr
      fixed <- FALSE
      if (dat$control$delete.nodes == TRUE) {
        el <- matrix(as.edgelist(dat$nw$m), ncol = 2)
      } else {
        el <- get.dyads.active(dat$nw$m, at = at)
      }
    }
    if (type == "pers") {
      base.ai.BB.rate <- dat$param$base.ai.pers.BB.rate
      base.ai.BW.rate <- dat$param$base.ai.pers.BW.rate
      base.ai.WW.rate <- dat$param$base.ai.pers.WW.rate
      redux.ai.diag.rr <- dat$param$redux.ai.diag.pers.rr
      redux.ai.discl.rr <- dat$param$redux.ai.discl.pers.rr
      incr.ai.full.supp.rr <- dat$param$incr.ai.full.supp.pers.rr
      incr.ai.part.supp.rr <- dat$param$incr.ai.part.supp.pers.rr
      fixed <- FALSE
      if (dat$control$delete.nodes == TRUE) {
        el <- matrix(as.edgelist(dat$nw$p), ncol = 2)
      } else {
        el <- get.dyads.active(dat$nw$p, at = at)
      }
    }
    if (type == "inst") {
      base.ai.BB.rate <- 1
      base.ai.BW.rate <- 1
      base.ai.WW.rate <- 1
      redux.ai.diag.rr <- 0
      redux.ai.discl.rr <- 0
      incr.ai.full.supp.rr <- 0
      incr.ai.part.supp.rr <- 0
      fixed <- TRUE
      el <- matrix(as.edgelist(dat$nw$i), ncol = 2)
    }


    # Processes ---------------------------------------------------------------

    # Construct discordant edgelist
    disc.el <- el[status[el[, 1]] - status[el[, 2]] == 1, , drop = FALSE]
    disc.el <- rbind(disc.el, el[status[el[, 2]] - status[el[, 1]] == 1, 2:1, drop = FALSE])


    if (nrow(disc.el) > 0) {

      exp.ai <- rep(NA, dim(disc.el)[1])

      race.1 <- race[disc.el[, 1]]
      race.2 <- race[disc.el[, 2]]
      num.B <- (race.1 == "B") + (race.2 == "B")

      exp.ai <- (num.B == 2) * base.ai.BB.rate +
                (num.B == 1) * base.ai.BW.rate +
                (num.B == 0) * base.ai.WW.rate

      pos.diag <- diag.status[disc.el[, 1]]
      pos.tx    <- tx.status[disc.el[, 1]]
      pos.tt.traj <- tt.traj[disc.el[, 1]]

      dlist <- dat$temp$discl.list
      disclosed <- sapply(1:nrow(disc.el), function(x) {
        length(intersect(which(uid[disc.el[x, 1]] == dlist$pos),
                         which(uid[disc.el[x, 2]] == dlist$neg))) != 0
      })

      exp.ai[pos.diag == 1] <- exp.ai[pos.diag == 1] * (1 - redux.ai.diag.rr)
      exp.ai[disclosed == TRUE] <- exp.ai[disclosed == TRUE] *
                                   (1 - redux.ai.discl.rr)
      exp.ai[pos.tx == 1 & pos.tt.traj == "YF"] <- exp.ai[pos.tx == 1 & pos.tt.traj == "YF"] *
                                                   (1 + incr.ai.full.supp.rr)
      exp.ai[pos.tx == 1 & pos.tt.traj == "YP"] <- exp.ai[pos.tx == 1 & pos.tt.traj == "YP"] *
                                                   (1 + incr.ai.part.supp.rr)

      if (fixed == FALSE) {
        ai <- rpois(length(exp.ai), exp.ai)
      } else {
        ai <- round(exp.ai)
      }

      if (sum(ai) > 0) {
        result <- data.frame(pos = rep(disc.el[, 1], ai),
                             neg = rep(disc.el[, 2], ai),
                             type = toupper(substr(type, 1, 1)),
                             uai = NA,
                             ins = NA,
                             stringsAsFactors = FALSE)
      } else {
        result <- data.frame(pos = NULL,
                             neg = NULL,
                             type = NULL,
                             uai = NULL,
                             ins = NULL,
                             stringsAsFactors = FALSE)
      }

      if (type == "main") {
        if (dat$control$save.dal == TRUE) {
          dat$temp$dal[[at]] <- result
        } else {
          dat$temp$dal <- result
        }
      } else {
        if (dat$control$save.dal == TRUE) {
          dat$temp$dal[[at]] <- rbind(dat$temp$dal[[at]], result)
        } else {
          dat$temp$dal <- rbind(dat$temp$dal, result)
        }
      }

    }
  }

  return(dat)
}
