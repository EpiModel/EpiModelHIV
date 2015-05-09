
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
      cprob.BB <- dat$param$cprob.main.BB
      cprob.BW <- dat$param$cprob.main.BW
      cprob.WW <- dat$param$cprob.main.WW
      beta.diag <- dat$param$beta.cond.diag.main
      beta.discl <- dat$param$beta.cond.discl.main
      beta.fsupp <- dat$param$beta.cond.fsupp.main
      beta.psupp <- dat$param$beta.cond.psupp.main
      dal <- dal[dal$type == "M", ]
    }
    if (type == "pers") {
      cprob.BB <- dat$param$cprob.pers.BB
      cprob.BW <- dat$param$cprob.pers.BW
      cprob.WW <- dat$param$cprob.pers.WW
      beta.diag <- dat$param$beta.cond.diag.pers
      beta.discl <- dat$param$beta.cond.discl.pers
      beta.fsupp <- dat$param$beta.cond.fsupp.pers
      beta.psupp <- dat$param$beta.cond.psupp.pers
      dal <- dal[dal$type == "P", ]
    }
    if (type == "inst") {
      cprob.BB <- dat$param$cprob.inst.BB
      cprob.BW <- dat$param$cprob.inst.BW
      cprob.WW <- dat$param$cprob.inst.WW
      beta.diag <- dat$param$beta.cond.diag.inst
      beta.discl <- dat$param$beta.cond.discl.inst
      beta.fsupp <- dat$param$beta.cond.fsupp.inst
      beta.psupp <- dat$param$beta.cond.psupp.inst
      dal <- dal[dal$type == "I", ]
    }



    # Processes ---------------------------------------------------------------

    if (nrow(dal) > 0) {
      prob.cond <- rep(NA, dim(dal)[1])
    }

    race.1 <- race[dal[, 1]]
    race.2 <- race[dal[, 2]]
    num.B <- (race.1 == "B") + (race.2 == "B")

    prob.cond <- (num.B == 2) * cprob.BB +
                 (num.B == 1) * cprob.BW +
                 (num.B == 0) * cprob.WW

    logodds.cond <- log(prob.cond / (1 - prob.cond))

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
    logodds.cond[isDx] <- logodds.cond[isDx] * (1 - beta.diag)

    # Odds, Disclosed
    isDisc <- which(discl == 1)
    logodds.cond[isDisc] <- logodds.cond[isDisc] * (1 - beta.discl)

    # Odds, Tx Full Suppress Type
    isFS <- which(pos.tx == 1 & pos.tt.traj == "YF")
    logodds.cond[isFS] <- logodds.cond[isFS] * (1 + beta.fsupp)

    # Odds, Tx Part Supress Type
    isPS <- which(pos.tx == 1 & pos.tt.traj == "YP")
    logodds.cond[isPS] <- logodds.cond[isPS] * (1 + beta.psupp)

    old.prob.cond <- prob.cond
    prob.cond <- exp(logodds.cond) / (1 + exp(logodds.cond))

    prob.cond[is.na(prob.cond) & old.prob.cond == 0] <- 0
    prob.cond[is.na(prob.cond) & old.prob.cond == 1] <- 1

    uai <- rbinom(length(prob.cond), 1, 1 - prob.cond)


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

