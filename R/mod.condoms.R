
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
    }
    if (type == "pers") {
      cond.BB.prob <- dat$param$cond.pers.BB.prob
      cond.BW.prob <- dat$param$cond.pers.BW.prob
      cond.WW.prob <- dat$param$cond.pers.WW.prob
      diag.beta <- dat$param$cond.diag.pers.beta
      discl.beta <- dat$param$cond.discl.pers.beta
      cond.always <- dat$attr$cond.always.pers
    }
    if (type == "inst") {
      cond.BB.prob <- dat$param$cond.inst.BB.prob
      cond.BW.prob <- dat$param$cond.inst.BW.prob
      cond.WW.prob <- dat$param$cond.inst.WW.prob
      diag.beta <- dat$param$cond.diag.inst.beta
      discl.beta <- dat$param$cond.discl.inst.beta
      cond.always <- dat$attr$cond.always.inst
    }

    el <- dat$temp$el
    elt <- el[el$type == type, ]

    ## Process ##

    # Base condom probs
    race.p1 <- race[elt$p1]
    race.p2 <- race[elt$p2]
    num.B <- (race.p1 == "B") + (race.p2 == "B")
    cond.prob <- (num.B == 2) * (cond.BB.prob * cond.rr.BB) +
                 (num.B == 1) * (cond.BW.prob * cond.rr.BW) +
                 (num.B == 0) * (cond.WW.prob * cond.rr.WW)

    # Transform to UAI logit
    uai.prob <- 1 - cond.prob
    uai.logodds <- log(uai.prob / (1 - uai.prob))

    # Diagnosis modifier
    pos.diag <- diag.status[elt$p1]
    isDx <- which(pos.diag == 1)
    uai.logodds[isDx] <- uai.logodds[isDx] + diag.beta

    # Disclosure modifier
    isDiscord <- which((elt$st1 - elt$st2) == 1)
    delt <- elt[isDiscord, ]
    dlist <- dat$temp$discl.list
    discl.disc <- sapply(1:nrow(delt), function(x) {
      length(intersect(which(uid[delt$p1[x]] == dlist$pos),
                       which(uid[delt$p2[x]] == dlist$neg))) != 0
    })
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
      ca1 <- cond.always[elt$p1]
      ca2 <- cond.always[elt$p2]
      uai.prob <- ifelse(ca1 == 1 | ca2 == 1, 0, uai.prob)
      if (type == "pers") {
        dat$epi$cprob.always.pers[at] <- mean(uai.prob == 0)
      } else {
        dat$epi$cprob.always.inst[at] <- mean(uai.prob == 0)
      }
    }

    elt$uai <- rbinom(nrow(elt), elt$ai, uai.prob)

    ## Output ##
    dat$temp$el[dat$temp$el$type == type, ] <- elt

  } # end ptype loop

  ## Construct discordant act list
  del <- dat$temp$el[which((dat$temp$el$st1 - dat$temp$el$st2) == 1), c(1:2, 5:7)]
  dal <- data.frame(pos = rep(del$p1, del$ai), neg = rep(del$p2, del$ai),
                    type = rep(del$type, del$ai))
  dal$uai <- rep(NA, nrow(dal))
  if (nrow(dal) > 0) {
    uai.vec <- do.call("c", lapply(1:length(del$ai),
                                   function(x) rep(1:0, c(del$uai[x],
                                                          del$ai[x] - del$uai[x]))))
    dal$uai <- uai.vec
  }

  dat$temp$dal <- dal

  return(dat)
}
