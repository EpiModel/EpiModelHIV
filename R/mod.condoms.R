
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
#' given HIV antiretroviral therapy.
#'
#' @return
#' Updates the discordant edgelist with a `uai` variable indicating whether
#' condoms were used in that act.
#'
#' @export
#'
condoms_msm <- function(dat, at) {

  ## Input
  # Attributes
  race        <- get_attr(dat, "race")
  diag.status <- get_attr(dat, "diag.status")
  prepStat    <- get_attr(dat, "prepStat")
  age         <- get_attr(dat, "age")

  # EDP attributes
  edp.prepStat <- get_attr(dat, "edp.prepStat")

  # Parameters
  netstats   <- get_param(dat, "netstats")
  epistats   <- get_param(dat, "epistats")
  cond.scale <- get_param(dat, "cond.scale")

  race.flag <- netstats[["race"]]
  geog.lvl  <- netstats[["race"]]

  # Condom Use Models
  cond.mc.mod <- epistats[["cond.mc.mod"]]
  cond.oo.mod <- epistats[["cond.oo.mod"]]

  ## Process

  # Temp edgelist
  el <- dat[["temp"]][["el"]]

  race.combo <- get_race_combo(race[el[, 1]], race[el[, 2]])

  comb.age <- age[el[, 1]] + age[el[, 2]]

  hiv.concord.pos <- rep(0, nrow(el))
  cp <- which(diag.status[el[, 1]] == 1 & diag.status[el[, 2]] == 1)
  hiv.concord.pos[cp] <- 1

  any.prep <- as.numeric((prepStat[el[, 1]] + prepStat[el[, 2]] +
                            edp.prepStat[el[, 1]] + edp.prepStat[el[, 2]]) > 0)

  ## Main/casual partnerships ##
  mc.parts <- which(el[, "ptype"] != 3)
  el.mc <- el[mc.parts, ]

  if (!is.null(geog.lvl)) {
    if (race.flag == TRUE) {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = el.mc[, "durations"],
        race.combo = race.combo[mc.parts],
        comb.age = comb.age[mc.parts],
        hiv.concord.pos = hiv.concord.pos[mc.parts],
        prep = any.prep[mc.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
    } else {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = el.mc[, "durations"],
        comb.age = comb.age[mc.parts],
        hiv.concord.pos = hiv.concord.pos[mc.parts],
        prep = any.prep[mc.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
    }
  } else {
    if (race.flag == TRUE) {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = el.mc[, "durations"],
        race.combo = race.combo[mc.parts],
        comb.age = comb.age[mc.parts],
        hiv.concord.pos = hiv.concord.pos[mc.parts],
        prep = any.prep[mc.parts]
      )
      cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
    } else {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = el.mc[, "durations"],
        comb.age = comb.age[mc.parts],
        hiv.concord.pos = hiv.concord.pos[mc.parts],
        prep = any.prep[mc.parts]
      )
      cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
    }
  }
  el.mc <- cbind(el.mc, cond.prob)

  ## One-off partnerships ##
  oo.parts <- which(el[, "ptype"] == 3)
  el.oo <- el[oo.parts, ]

  if (!is.null(geog.lvl)) {
    if (race.flag == TRUE) {
      x <- data.frame(
        race.combo = race.combo[oo.parts],
        comb.age = comb.age[oo.parts],
        hiv.concord.pos = hiv.concord.pos[oo.parts],
        prep = any.prep[oo.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
      el.oo <- cbind(el.oo, cond.prob)
    } else {
      x <- data.frame(
        comb.age = comb.age[oo.parts],
        hiv.concord.pos = hiv.concord.pos[oo.parts],
        prep = any.prep[oo.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
      el.oo <- cbind(el.oo, cond.prob)
    }
  } else {
    if (race.flag == TRUE) {
      x <- data.frame(
        race.combo = race.combo[oo.parts],
        comb.age = comb.age[oo.parts],
        hiv.concord.pos = hiv.concord.pos[oo.parts],
        prep = any.prep[oo.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
      el.oo <- cbind(el.oo, cond.prob)
    } else {
      x <- data.frame(
        comb.age = comb.age[oo.parts],
        hiv.concord.pos = hiv.concord.pos[oo.parts],
        prep = any.prep[oo.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
      el.oo <- cbind(el.oo, cond.prob)
    }
  }

  ## Bind el together
  el <- rbind(el.mc, el.oo)

  # Acts
  ai.vec <- el[, "ai"]
  pid <- rep(seq_along(ai.vec), ai.vec)
  p1 <- rep(el[, "p1"], ai.vec)
  p2 <- rep(el[, "p2"], ai.vec)
  ptype <- rep(el[, "ptype"], ai.vec)
  cond.prob <- rep(el[, "cond.prob"], ai.vec)

  cond.prob <- cond.prob * cond.scale

  # UAI draw per act
  uai <- runif(length(cond.prob)) < 1 - cond.prob

  # Act list construction
  al <- cbind(p1, p2, ptype, uai, pid)
  dat[["temp"]][["al"]] <- al

  ## Output

  return(dat)
}
