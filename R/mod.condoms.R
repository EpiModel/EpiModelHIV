
#' @title Condom Use Module
#'
#' @description Module function stochastically simulates potential condom use
#'              for each act on the discordant edgelist.
#'
#' @inheritParams aging.mard
#'
#' @details
#' For each act on the discordant edgelist, condom use is stochastically simulated
#' based on the partnership type and racial combination of the dyad. Other
#' modifiers for the probability of condom use in that pair are diagnosis of
#' disease, disclosure of status, and full or partial HIV viral suppression
#' given HIV anti-retroviral therapy.
#'
#' @return
#' Updates the discordant edgelist with a \code{uai} variable indicating whether
#' condoms were used in that act.
#'
#' @keywords module
#' @export
#'
condoms.mard <- function(dat, at) {

  for (type in c("main", "pers", "inst")) {

    # Variables ---------------------------------------------------------------

    # Attributes
    uid <- dat$attr$uid
    status <- dat$attr$status
    diag.status <- dat$attr$diag.status
    race <- dat$attr$race
    tx.status <- dat$attr$tx.status
    tt.traj <- dat$attr$tt.traj

    # Parameters and data
    if (dat$control$save.dal == TRUE) {
      dal <- dat$temp$dal[[at]]
    } else {
      dal <- dat$temp$dal
    }

    if (type == "main") {
      c.BB.prob <- dat$param$c.main.BB.prob
      c.BW.prob <- dat$param$c.main.BW.prob
      c.WW.prob <- dat$param$c.main.WW.prob
      diag.beta <- dat$param$cond.diag.main.beta
      discl.beta <- dat$param$cond.discl.main.beta
      fsupp.beta <- dat$param$cond.fsupp.main.beta
      psupp.beta <- dat$param$cond.psupp.main.beta
      dal <- dal[dal$type == "M", ]
    }
    if (type == "pers") {
      c.BB.prob <- dat$param$c.pers.BB.prob
      c.BW.prob <- dat$param$c.pers.BW.prob
      c.WW.prob <- dat$param$c.pers.WW.prob
      diag.beta <- dat$param$cond.diag.pers.beta
      discl.beta <- dat$param$cond.discl.pers.beta
      fsupp.beta <- dat$param$cond.fsupp.pers.beta
      psupp.beta <- dat$param$cond.psupp.pers.beta
      dal <- dal[dal$type == "P", ]
    }
    if (type == "inst") {
      c.BB.prob <- dat$param$c.inst.BB.prob
      c.BW.prob <- dat$param$c.inst.BW.prob
      c.WW.prob <- dat$param$c.inst.WW.prob
      diag.beta <- dat$param$cond.diag.inst.beta
      discl.beta <- dat$param$cond.discl.inst.beta
      fsupp.beta <- dat$param$cond.fsupp.inst.beta
      psupp.beta <- dat$param$cond.psupp.inst.beta
      dal <- dal[dal$type == "I", ]
    }



    # Processes ---------------------------------------------------------------

    if (nrow(dal) > 0) {
      cond.prob <- rep(NA, dim(dal)[1])
    }

    race.1 <- race[dal[, 1]]
    race.2 <- race[dal[, 2]]
    num.B <- (race.1 == "B") + (race.2 == "B")

    cond.prob <- (num.B == 2) * c.BB.prob +
                 (num.B == 1) * c.BW.prob +
                 (num.B == 0) * c.WW.prob

    logodds.cond <- log(cond.prob / (1 - cond.prob))

    pos.diag <- diag.status[dal[, 1]]

    dlist <- dat$temp$discl.list
    discl <- sapply(1:nrow(dal), function(x) {
      sum(dlist$pos == uid[dal[x, 1]] &
            dlist$neg == uid[dal[x, 2]]) > 0
    })

    pos.tx <- tx.status[dal[, 1]]
    pos.tt.traj <- tt.traj[dal[, 1]]

    # Odds, Diagnosed
    isDx <- which(pos.diag == 1)
    logodds.cond[isDx] <- logodds.cond[isDx] * (1 - diag.beta)

    # Odds, Disclosed
    isDisc <- which(discl == 1)
    logodds.cond[isDisc] <- logodds.cond[isDisc] * (1 - discl.beta)

    # Odds, Tx Full Suppress Type
    isFS <- which(pos.tx == 1 & pos.tt.traj == "YF")
    logodds.cond[isFS] <- logodds.cond[isFS] * (1 + fsupp.beta)

    # Odds, Tx Part Supress Type
    isPS <- which(pos.tx == 1 & pos.tt.traj == "YP")
    logodds.cond[isPS] <- logodds.cond[isPS] * (1 + psupp.beta)

    old.cond.prob <- cond.prob
    cond.prob <- exp(logodds.cond) / (1 + exp(logodds.cond))

    cond.prob[is.na(cond.prob) & old.cond.prob == 0] <- 0
    cond.prob[is.na(cond.prob) & old.cond.prob == 1] <- 1

    uai <- rbinom(length(cond.prob), 1, 1 - cond.prob)


    # Output ------------------------------------------------------------------

    if (dat$control$save.dal == TRUE) {
      if (type == "main") {
        dat$temp$dal[[at]]$uai[dat$temp$dal[[at]]$type == "M"] <- uai
      }
      if (type == "pers") {
        dat$temp$dal[[at]]$uai[dat$temp$dal[[at]]$type == "P"] <- uai
      }
      if (type == "inst") {
        dat$temp$dal[[at]]$uai[dat$temp$dal[[at]]$type == "I"] <- uai
      }
    } else {
      if (type == "main") {
        dat$temp$dal$uai[dat$temp$dal$type == "M"] <- uai
      }
      if (type == "pers") {
        dat$temp$dal$uai[dat$temp$dal$type == "P"] <- uai
      }
      if (type == "inst") {
        dat$temp$dal$uai[dat$temp$dal$type == "I"] <- uai
      }
    }


  }

  return(dat)
}

