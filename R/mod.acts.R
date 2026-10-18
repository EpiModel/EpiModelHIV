
#' @title Sexual Acts Module
#'
#' @description Module function for setting the number of sexual acts on the
#'              discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' The number of acts at each time step is specified as a function of the race of
#' both members in a pair and the expected values within black-black, black-white,
#' and white-white combinations. For one-off partnerships, this is deterministically
#' set at 1, whereas for main and causal partnerships it is a stochastic draw
#' from a Poisson distribution. The number of total acts may further be modified
#' by the level of HIV viral suppression in an infected person.
#'
#' @keywords module msm
#' @export
#'
acts_msm <- function(dat, at) {

  ## Inputs
  # Attributes
  status      <- get_attr(dat, "status")
  diag.status <- get_attr(dat, "diag.status")
  race        <- get_attr(dat, "race")
  age         <- get_attr(dat, "age")
  stage       <- get_attr(dat, "stage")
  vl          <- get_attr(dat, "vl")

  # Parameters
  acts.aids.vl <- get_param(dat, "acts.aids.vl")
  acts.scale   <- get_param(dat, "acts.scale")
  netstats     <- get_param(dat, "netstats")
  epistats     <- get_param(dat, "epistats")
  time.unit    <- get_param(dat, "time.unit")

  geog.lvl     <- netstats[["geog.lvl"]]
  race.flag    <- netstats[["race"]]
  acts.mod     <- epistats[["acts.mod"]]

  el <- get_cumulative_edgelists_df(dat)
  el <- el[is.na(el[["stop"]]), ]
  ptype <- el[["network"]]
  durations <- at - el[["start"]]

  el <- cbind(
    get_posit_ids(dat, el[["head"]]),
    get_posit_ids(dat, el[["tail"]])
  )

  # # Construct edgelist
  st1 <- status[el[, 1]]
  st2 <- status[el[, 2]]

  el <- cbind(el, st1, st2, ptype)
  colnames(el) <- c("p1", "p2", "st1", "st2", "ptype")

  # Subset to main/casual
  el.mc <- el[el[, "ptype"] != 3, ]

  # Base AI rates based on Poisson model for main/casual
  race.combo <- ARTnet::make_race_combo(race[el.mc[, 1]], race[el.mc[, 2]])

  comb.age <- age[el.mc[, 1]] + age[el.mc[, 2]]

  # Current partnership durations (main and casual)
  durations <- durations[el[, "ptype"] != 3]

  # HIV-positive concordant
  hiv.concord.pos <- rep(0, nrow(el.mc))
  cp <- which(diag.status[el.mc[, 1]] == 1 & diag.status[el.mc[, 2]] == 1)
  hiv.concord.pos[cp] <- 1

  # Model predictions
  if (!is.null(geog.lvl)) {
    if (race.flag == TRUE) {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = durations,
        race.combo = race.combo,
        comb.age = comb.age,
        hiv.concord.pos = hiv.concord.pos,
        geogYN = 1
      )
      rates <- unname(predict(acts.mod, newdata = x, type = "response")) / (364 / time.unit)

    } else {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = durations,
        comb.age = comb.age,
        hiv.concord.pos = hiv.concord.pos,
        geogYN = 1
      )
      rates <- unname(predict(acts.mod, newdata = x, type = "response")) / (364 / time.unit)
    }
  } else {
    if (race.flag == TRUE) {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = durations,
        race.combo = race.combo,
        comb.age = comb.age,
        hiv.concord.pos = hiv.concord.pos
      )
      rates <- unname(predict(acts.mod, newdata = x, type = "response")) / (364 / time.unit)
    } else {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = durations,
        comb.age = comb.age,
        hiv.concord.pos = hiv.concord.pos
      )
      rates <- unname(predict(acts.mod, newdata = x, type = "response")) / (364 / time.unit)
    }
  }
  rates <- rates * acts.scale
  ai <- rpois(length(rates), rates)
  el.mc <- cbind(el.mc, durations, ai)

  # Add one-time partnerships
  el.oo <- el[el[, "ptype"] == 3, ]
  ai <- durations <- rep(1, nrow(el.oo))
  el.oo <- cbind(el.oo, durations, ai)

  # Bind el back together
  el <- rbind(el.mc, el.oo)

  # For AIDS cases with VL above acts.aids.vl, reduce their their acts to 0
  p1HIV <- which(el[, "st1"] == 1)
  p1AIDS <- stage[el[p1HIV, "p1"]] == 4 & vl[el[p1HIV, "p1"]] >= acts.aids.vl
  el[p1HIV[p1AIDS == TRUE], "ai"] <- 0

  p2HIV <- which(el[, "st2"] == 1)
  p2AIDS <- stage[el[p2HIV, "p2"]] == 4 & vl[el[p2HIV, "p2"]] >= acts.aids.vl
  el[p2HIV[p2AIDS == TRUE], "ai"] <- 0

  # Flip order of discordant edges
  disc <- abs(el[, "st1"] - el[, "st2"]) == 1
  disc.st2pos <- which(disc == TRUE & el[, "st2"] == 1)
  el[disc.st2pos, 1:4] <- el[disc.st2pos, c(2, 1, 4, 3)]

  # Remove inactive edges from el
  el <- el[-which(el[, "ai"] == 0), ]

  # Save out
  dat[["temp"]][["el"]] <- el

  return(dat)
}
