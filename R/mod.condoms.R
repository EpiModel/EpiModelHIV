
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

  # Attributes
  uid <- dat$attr$uid
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  prepStat.la <- dat$attr$prepStat.la
  prepClass.la <- dat$attr$prepClass.la

  # Parameters
  rcomp.prob <- dat$param$rcomp.prob
  rcomp.adh.groups <- dat$param$rcomp.adh.groups
  rcomp.main.only <- dat$param$rcomp.main.only
  rcomp.discl.only <- dat$param$rcomp.discl.only

  el <- dat$temp$el

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

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
    isDiscord <- which((elt[, "st1"] - elt[, "st2"]) == 1) # pull vector of discordant
    pos.diag <- diag.status[elt[, 1]]
    isDx <- which(pos.diag == 1) # pull vector of diagnosis status
    isDiscord.dx <- intersect(isDiscord, isDx)
    uai.logodds[isDiscord.dx] <- uai.logodds[isDiscord.dx] + diag.beta

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
        dat$epi$cprob.always.pers <- NULL
        # dat$epi$cprob.always.pers[at] <- mean(uai.prob == 0)
      } else {
        dat$epi$cprob.always.inst <- NULL
        # dat$epi$cprob.always.inst[at] <- mean(uai.prob == 0)
      }
    }

    # PrEP Status (risk compensation)
    if (rcomp.prob > 0) {

      idsRC <- which((prepStat[elt[, 1]] == 1 & prepClass[elt[, 1]] %in% rcomp.adh.groups) |
                     (prepStat[elt[, 2]] == 1 & prepClass[elt[, 2]] %in% rcomp.adh.groups))
      idsRC.la <- which((prepStat.la[elt[, 1]] == 1 & prepClass.la[elt[, 1]] == 2) |
                        (prepStat.la[elt[, 2]] == 1 & prepClass.la[elt[, 2]] == 2))

      if (rcomp.main.only == TRUE & ptype > 1) {
        idsRC <- NULL
      }
      if (rcomp.discl.only == TRUE) {
        idsRC <- intersect(idsRC, isDisc)
      }
      uai.prob[idsRC] <- 1 - (1 - uai.prob[idsRC]) * (1 - rcomp.prob)
    }

    ai.vec <- elt[, "ai"]
    p1 <- rep(elt[, "p1"], ai.vec)
    p2 <- rep(elt[, "p2"], ai.vec)
    ptype <- rep(elt[, "ptype"], ai.vec)

    uai.prob.peract <- rep(uai.prob, ai.vec)
    uai <- rbinom(length(p1), 1, uai.prob.peract)

    if (type == "main") {
      pid <- rep(1:length(ai.vec), ai.vec)
      al <- cbind(p1, p2, ptype, uai, pid)
    } else {
      pid <- rep(max(al[, "pid"]) + (1:length(ai.vec)), ai.vec)
      tmp.al <- cbind(p1, p2, ptype, uai, pid)
      al <- rbind(al, tmp.al)
    }

  } # end ptype loop

  dat$temp$al <- al

  if (at == 2) {
    dat$epi$ai.events <- rep(NA, 2)
    dat$epi$uai.events <- rep(NA, 2)
  }
  dat$epi$ai.events[at] <- nrow(al)
  dat$epi$uai.events[at] <- sum(al[, "uai"])

  return(dat)
}
