
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
#' disease, and full or partial HIV viral suppression
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
  race <- dat$attr$race
  race.flag <- dat$param$netstats$race
  geog.lvl <- dat$param$netstats$race
  age <- dat$attr$age
  diag.status <- dat$attr$diag.status
  prepStat <- dat$attr$prepStat

  # Condom Use Models
  cond.mc.mod <- dat$param$epistats$cond.mc.mod
  cond.oo.mod <- dat$param$epistats$cond.oo.mod

  cond.scale <- dat$param$cond.scale

  # Temp edgelist
  el <- dat$temp$el

  race.combo <- rep(NA, nrow(el))
  race.combo[race[el[, 1]] == 1 & race[el[, 2]] == 1] <- 1
  race.combo[race[el[, 1]] == 1 & race[el[, 2]] %in% 2:3] <- 2
  race.combo[race[el[, 1]] == 2 & race[el[, 2]] %in% c(1, 3)] <- 3
  race.combo[race[el[, 1]] == 2 & race[el[, 2]] == 2] <- 4
  race.combo[race[el[, 1]] == 3 & race[el[, 2]] %in% 1:2] <- 5
  race.combo[race[el[, 1]] == 3 & race[el[, 2]] == 3] <- 6

  comb.age <- age[el[, 1]] + age[el[, 2]]

  hiv.concord.pos <- rep(0, nrow(el))
  cp <- which(diag.status[el[, 1]] == 1 & diag.status[el[, 2]] == 1)
  hiv.concord.pos[cp] <- 1

  any.prep <- as.numeric((prepStat[el[, 1]] + prepStat[el[, 2]]) > 0)

  ## Main/casual partnerships ##
  mc.parts <- which(el[, "ptype"] != 3)
  el.mc <- el[mc.parts, , drop = FALSE]

  if (nrow(el.mc) > 0) {
    if (!is.null(geog.lvl)) {
      if (race.flag == TRUE){
        x <- data.frame(ptype = el.mc[, "ptype"],
                        duration = el.mc[, "durations"],
                        race.combo = race.combo[mc.parts],
                        comb.age = comb.age[mc.parts],
                        hiv.concord.pos = hiv.concord.pos[mc.parts],
                        prep = any.prep[mc.parts],
                        geogYN = 1)
        cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
      } else {
        x <- data.frame(ptype = el.mc[, "ptype"],
                        duration = el.mc[, "durations"],
                        comb.age = comb.age[mc.parts],
                        hiv.concord.pos = hiv.concord.pos[mc.parts],
                        prep = any.prep[mc.parts],
                        geogYN = 1)
        cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
      }
    } else {
      if (race.flag == TRUE){
        x <- data.frame(ptype = el.mc[, "ptype"],
                        duration = el.mc[, "durations"],
                        race.combo = race.combo[mc.parts],
                        comb.age = comb.age[mc.parts],
                        hiv.concord.pos = hiv.concord.pos[mc.parts],
                        prep = any.prep[mc.parts])
        cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
      } else {
        x <- data.frame(ptype = el.mc[, "ptype"],
                        duration = el.mc[, "durations"],
                        comb.age = comb.age[mc.parts],
                        hiv.concord.pos = hiv.concord.pos[mc.parts],
                        prep = any.prep[mc.parts])
        cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
      }
    }
    el.mc <- cbind(el.mc, cond.prob)
  }

  ## One-off partnerships ##
  oo.parts <- which(el[, "ptype"] == 3)
  el.oo <- el[oo.parts, , drop = FALSE]

  if (nrow(el.oo) > 0) {
    if (!is.null(geog.lvl)) {
      if (race.flag == TRUE){
        x <- data.frame(race.combo = race.combo[oo.parts],
                        comb.age = comb.age[oo.parts],
                        hiv.concord.pos = hiv.concord.pos[oo.parts],
                        prep = any.prep[oo.parts],
                        geogYN = 1)
        cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
        el.oo <- cbind(el.oo, cond.prob)
      } else {
        x <- data.frame(comb.age = comb.age[oo.parts],
                        hiv.concord.pos = hiv.concord.pos[oo.parts],
                        prep = any.prep[oo.parts],
                        geogYN = 1)
        cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
        el.oo <- cbind(el.oo, cond.prob)
      }
    } else {
      if (race.flag == TRUE){
        x <- data.frame(race.combo = race.combo[oo.parts],
                        comb.age = comb.age[oo.parts],
                        hiv.concord.pos = hiv.concord.pos[oo.parts],
                        prep = any.prep[oo.parts],
                        geogYN = 1)
        cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
        el.oo <- cbind(el.oo, cond.prob)
      } else {
        x <- data.frame(comb.age = comb.age[oo.parts],
                        hiv.concord.pos = hiv.concord.pos[oo.parts],
                        prep = any.prep[oo.parts],
                        geogYN = 1)
        cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
        el.oo <- cbind(el.oo, cond.prob)
      }
    }

    ## Bind el together
    el <- rbind(el.mc, el.oo)
  } else {
    el <- el.mc
  }

  # Acts
  ai.vec <- el[, "ai"]
  pid <- rep(1:length(ai.vec), ai.vec)
  p1 <- rep(el[, "p1"], ai.vec)
  p2 <- rep(el[, "p2"], ai.vec)
  ptype <- rep(el[, "ptype"], ai.vec)
  cond.prob <- rep(el[, "cond.prob"], ai.vec)

  cond.prob <- cond.prob * cond.scale

  # UAI draw per act
  uai <- rbinom(length(cond.prob), 1, 1 - cond.prob)

  # Act list construction
  al <- cbind(p1, p2, ptype, uai, pid)
  dat$temp$al <- al

  return(dat)
}
