
#' @title Condom Use Module
#'
#' @description Module function stochastically simulates potential condom use
#'              for each act on the discordant edgelist.
#'
#' @inheritParams aging_msm
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
#' @keywords module msm
#' @export
#'
condoms_msm <- function(dat, at) {

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

    # Attributes
    uid <- dat$attr$uid
    diag.status <- dat$attr$diag.status
    race <- dat$attr$race

    # Parameters
    cond.rr.BB <- dat$param$cond.rr.BB
    cond.rr.BW <- dat$param$cond.rr.BW
    cond.rr.WW <- dat$param$cond.rr.WW

    if (type == "main") {
      cond.BB.prob <- dat$param$cond.main.BB.prob
      cond.BW.prob <- dat$param$cond.main.BW.prob
      cond.WW.prob <- dat$param$cond.main.WW.prob
      diag.beta <- dat$param$cond.diag.main.beta
      discl.beta <- dat$param$cond.discl.main.beta
      cond.always <- NULL
      ptype <- 1
    }
    if (type == "pers") {
      cond.BB.prob <- dat$param$cond.pers.BB.prob
      cond.BW.prob <- dat$param$cond.pers.BW.prob
      cond.WW.prob <- dat$param$cond.pers.WW.prob
      diag.beta <- dat$param$cond.diag.pers.beta
      discl.beta <- dat$param$cond.discl.pers.beta
      cond.always <- dat$attr$cond.always.pers
      ptype <- 2
    }
    if (type == "inst") {
      cond.BB.prob <- dat$param$cond.inst.BB.prob
      cond.BW.prob <- dat$param$cond.inst.BW.prob
      cond.WW.prob <- dat$param$cond.inst.WW.prob
      diag.beta <- dat$param$cond.diag.inst.beta
      discl.beta <- dat$param$cond.discl.inst.beta
      cond.always <- dat$attr$cond.always.inst
      ptype <- 3
    }

    el <- dat$temp$el
    elt <- el[el[, "ptype"] == ptype, ]

    ## Process ##

    # Base condom probs
    race.p1 <- race[elt[, 1]]
    race.p2 <- race[elt[, 2]]
    num.B <- (race.p1 == "B") + (race.p2 == "B")
    cond.prob <- (num.B == 2) * (cond.BB.prob * cond.rr.BB) +
                 (num.B == 1) * (cond.BW.prob * cond.rr.BW) +
                 (num.B == 0) * (cond.WW.prob * cond.rr.WW)

    # Transform to UAI logit
    uai.prob <- 1 - cond.prob
    uai.logodds <- log(uai.prob / (1 - uai.prob))

    # Diagnosis modifier
    pos.diag <- diag.status[elt[, 1]]
    isDx <- which(pos.diag == 1)
    uai.logodds[isDx] <- uai.logodds[isDx] + diag.beta

    # Disclosure modifier
    isDiscord <- which((elt[, "st1"] - elt[, "st2"]) == 1)
    delt <- elt[isDiscord, ]
    discl.list <- dat$temp$discl.list
    disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
    delt.cdl <- uid[delt[, 1]] * 1e7 + uid[delt[, 2]]
    discl.disc <- (delt.cdl %in% disclose.cdl)

    discl <- rep(NA, nrow(elt))
    discl[isDiscord] <- discl.disc

    isDisc <- which(discl == 1)
    uai.logodds[isDisc] <- uai.logodds[isDisc] + discl.beta

    # Back transform to prob
    old.uai.prob <- uai.prob
    uai.prob <- exp(uai.logodds) / (1 + exp(uai.logodds))

    uai.prob[is.na(uai.prob) & old.uai.prob == 0] <- 0
    uai.prob[is.na(uai.prob) & old.uai.prob == 1] <- 1

    # UAI group
    if (type %in% c("pers", "inst")) {
      ca1 <- cond.always[elt[, 1]]
      ca2 <- cond.always[elt[, 2]]
      uai.prob <- ifelse(ca1 == 1 | ca2 == 1, 0, uai.prob)
      if (type == "pers") {
        dat$epi$cprob.always.pers[at] <- mean(uai.prob == 0)
      } else {
        dat$epi$cprob.always.inst[at] <- mean(uai.prob == 0)
      }
    }

    ai.vec <- elt[, "ai"]
    pos <- rep(elt[, "p1"], ai.vec)
    neg <- rep(elt[, "p2"], ai.vec)
    ptype <- rep(elt[, "ptype"], ai.vec)

    uai.prob.peract <- rep(uai.prob, ai.vec)
    uai <- rbinom(length(pos), 1, uai.prob.peract)

    if (type == "main") {
      pid <- rep(1:length(ai.vec), ai.vec)
      al <- cbind(pos, neg, ptype, uai, pid)
    } else {
      pid <- rep(max(al[, "pid"]) + (1:length(ai.vec)), ai.vec)
      tmp.al <- cbind(pos, neg, ptype, uai, pid)
      al <- rbind(al, tmp.al)
    }

  } # end ptype loop

  dat$temp$al <- al

  return(dat)
}
