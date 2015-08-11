
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

    ## Variables ##

    # Attributes
    status <- dat$attr$status
    race <- dat$attr$race

    # Parameters
    ai.scale <- dat$param$ai.scale
    if (type == "main") {
      base.ai.BB.rate <- dat$param$base.ai.main.BB.rate
      base.ai.BW.rate <- dat$param$base.ai.main.BW.rate
      base.ai.WW.rate <- dat$param$base.ai.main.WW.rate
      fixed <- FALSE
      el <- get.dyads.active(dat$nw$m, at = at)
    }
    if (type == "pers") {
      base.ai.BB.rate <- dat$param$base.ai.pers.BB.rate
      base.ai.BW.rate <- dat$param$base.ai.pers.BW.rate
      base.ai.WW.rate <- dat$param$base.ai.pers.WW.rate
      fixed <- FALSE
      el <- get.dyads.active(dat$nw$p, at = at)
    }
    if (type == "inst") {
      base.ai.BB.rate <- 1
      base.ai.BW.rate <- 1
      base.ai.WW.rate <- 1
      fixed <- ifelse(ai.scale != 1, FALSE, TRUE)
      el <- matrix(as.edgelist(dat$nw$i), ncol = 2)
    }

    ## Processes ##
    # Construct discordant edgelist
    disc.el <- el[status[el[, 1]] - status[el[, 2]] == 1, , drop = FALSE]
    disc.el <- rbind(disc.el, el[status[el[, 2]] - status[el[, 1]] == 1, 2:1, drop = FALSE])

    if (nrow(disc.el) > 0) {

      ai.rate <- rep(NA, dim(disc.el)[1])

      race.1 <- race[disc.el[, 1]]
      race.2 <- race[disc.el[, 2]]
      num.B <- (race.1 == "B") + (race.2 == "B")

      ai.rate <- (num.B == 2) * base.ai.BB.rate +
                 (num.B == 1) * base.ai.BW.rate +
                 (num.B == 0) * base.ai.WW.rate
      ai.rate <- ai.rate * ai.scale

      if (fixed == FALSE) {
        ai <- rpois(length(ai.rate), ai.rate)
      } else {
        ai <- round(ai.rate)
      }

      if (sum(ai) > 0) {
        result <- data.frame(pos = rep(disc.el[, 1], ai),
                             neg = rep(disc.el[, 2], ai),
                             type = toupper(substr(type, 1, 1)),
                             uai = NA, ins = NA, stringsAsFactors = FALSE)
      } else {
        result <- data.frame(pos = NULL, neg = NULL, type = NULL,
                             uai = NULL, ins = NULL, stringsAsFactors = FALSE)
      }

      if (type == "main") {
        dat$temp$dal <- result
      } else {
        dat$temp$dal <- rbind(dat$temp$dal, result)
      }

    }
  }

  return(dat)
}
