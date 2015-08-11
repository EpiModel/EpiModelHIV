
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
    el <- as.data.frame(el)
    names(el) <- c("p1", "p2")
    st1 <- status[el$p1]
    st2 <- status[el$p2]
    disc <- abs(st1 - st2) == 1
    el[which(disc == 1 & st2 == 1), ] <- el[which(disc == 1 & st2 == 1), 2:1]
    el$st1 <- status[el$p1]
    el$st2 <- status[el$p2]

    if (nrow(el) > 0) {

      # Base AI rates
      ai.rate <- rep(NA, nrow(el))
      race.p1 <- race[el$p1]
      race.p2 <- race[el$p2]
      num.B <- (race.p1 == "B") + (race.p2 == "B")
      ai.rate <- (num.B == 2) * base.ai.BB.rate +
                 (num.B == 1) * base.ai.BW.rate +
                 (num.B == 0) * base.ai.WW.rate
      ai.rate <- ai.rate * ai.scale

      # Final act number
      if (fixed == FALSE) {
        ai <- rpois(length(ai.rate), ai.rate)
      } else {
        ai <- round(ai.rate)
      }

      # Full edge list
      el$ai <- ai
      el$type <- toupper(substr(type, 1, 1))

      # Discordant act list
      dal <- el[(st1 - st2) == 1, c(1, 2, 5)]
      dal <- data.frame(pos = rep(dal$p1, dal$ai), neg = rep(dal$p2, dal$ai),
                        type = toupper(substr(type, 1, 1)), uai = NA, ins = NA)

      if (type == "main") {
        dat$temp$dal <- dal
        dat$temp$el <- el
      } else {
        dat$temp$dal <- rbind(dat$temp$dal, dal)
        dat$temp$el <- rbind(dat$temp$el, el)
      }
    }

  } # loop over type end

  return(dat)
}
