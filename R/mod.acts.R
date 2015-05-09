
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
      base.exp.ai.BB <- dat$param$base.exp.ai.main.BB
      base.exp.ai.BW <- dat$param$base.exp.ai.main.BW
      base.exp.ai.WW <- dat$param$base.exp.ai.main.WW
      redux.exp.ai.diag <- dat$param$redux.exp.ai.diag.main
      redux.exp.ai.discl <- dat$param$redux.exp.ai.discl.main
      incr.exp.ai.full.supp <- dat$param$incr.exp.ai.full.supp.main
      incr.exp.ai.part.supp <- dat$param$incr.exp.ai.part.supp.main
      fixed <- FALSE
      if (dat$control$delete.nodes == TRUE) {
        el <- matrix(as.edgelist(dat$nw$m), ncol = 2)
      } else {
        el <- get.dyads.active(dat$nw$m, at = at)
      }
    }
    if (type == "pers") {
      base.exp.ai.BB <- dat$param$base.exp.ai.pers.BB
      base.exp.ai.BW <- dat$param$base.exp.ai.pers.BW
      base.exp.ai.WW <- dat$param$base.exp.ai.pers.WW
      redux.exp.ai.diag <- dat$param$redux.exp.ai.diag.pers
      redux.exp.ai.discl <- dat$param$redux.exp.ai.discl.pers
      incr.exp.ai.full.supp <- dat$param$incr.exp.ai.full.supp.pers
      incr.exp.ai.part.supp <- dat$param$incr.exp.ai.part.supp.pers
      fixed <- FALSE
      if (dat$control$delete.nodes == TRUE) {
        el <- matrix(as.edgelist(dat$nw$p), ncol = 2)
      } else {
        el <- get.dyads.active(dat$nw$p, at = at)
      }
    }
    if (type == "inst") {
      base.exp.ai.BB <- 1
      base.exp.ai.BW <- 1
      base.exp.ai.WW <- 1
      redux.exp.ai.diag <- 0
      redux.exp.ai.discl <- 0
      incr.exp.ai.full.supp <- 0
      incr.exp.ai.part.supp <- 0
      fixed <- TRUE
      el <- matrix(as.edgelist(dat$nw$i), ncol = 2)
    }


    # Processes ---------------------------------------------------------------

    disc.el <- disc_check(el, status)

    if (nrow(disc.el) > 0) {

      exp.ai <- rep(NA, dim(disc.el)[1])

      race.1 <- race[disc.el[, 1]]
      race.2 <- race[disc.el[, 2]]
      num.B <- (race.1 == "B") + (race.2 == "B")

      exp.ai <- (num.B == 2) * base.exp.ai.BB +
                (num.B == 1) * base.exp.ai.BW +
                (num.B == 0) * base.exp.ai.WW

      pos.diag <- diag.status[disc.el[, 1]]
      pos.tx    <- tx.status[disc.el[, 1]]
      pos.tt.traj <- tt.traj[disc.el[, 1]]

      dlist <- dat$temp$discl.list
      disclosed <- sapply(1:nrow(disc.el), function(x) {
        length(intersect(which(uid[disc.el[x, 1]] == dlist$pos),
                         which(uid[disc.el[x, 2]] == dlist$neg))) != 0
      })

      exp.ai[pos.diag == 1] <- exp.ai[pos.diag == 1] * (1 - redux.exp.ai.diag)
      exp.ai[disclosed == TRUE] <- exp.ai[disclosed == TRUE] *
                                   (1 - redux.exp.ai.discl)
      exp.ai[pos.tx == 1 & pos.tt.traj == "YF"] <- exp.ai[pos.tx == 1 & pos.tt.traj == "YF"] *
                                                   (1 + incr.exp.ai.full.supp)
      exp.ai[pos.tx == 1 & pos.tt.traj == "YP"] <- exp.ai[pos.tx == 1 & pos.tt.traj == "YP"] *
                                                   (1 + incr.exp.ai.part.supp)

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


#' @title Check for Discordant Pairs
#'
#' @description Checks that the discordant edgelist is in a proper form before
#'              returning its values in a matrix.
#'
#' @param el Current edgelist.
#' @param status Disease status vector.
#'
#' @return
#' This function returns an edgelist of only the discordant relations, with the
#' seropositive man in the first column.
#'
#' @keywords submodule
#' @export
#'
disc_check <- function(el, status) {

  if (class(el) != "matrix") {
    stop("argument el must be a matrix")
  }
  if (ncol(el) != 2) {
    stop("argument el must have two columns")
  }
  if (mode(el) != "numeric") {
    stop("argument el must be numeric")
  }
  if (any(el != round(el)) == "TRUE") {
    stop("argument el must be integer values")
  }
  if ((0 %in% el) == "TRUE") {
    stop("argument el cannot contain 0")
  }
  if (any(el > length(status)) == "TRUE") {
    stop("argument el cannot contain values greater than length of status")
  }

  out <- el[status[el[, 1]] - status[el[, 2]] == 1, , drop = FALSE]
  out <- rbind(out, el[status[el[, 2]] - status[el[, 1]] == 1, 2:1, drop = FALSE])

  return(out)
}

