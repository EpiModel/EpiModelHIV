
#' @title Arrivals Module
#'
#' @description Module function for arrivals into the sexually active
#'              population.
#'
#' @inheritParams aging_msm
#'
#' @details
#' New population members are added based on expected numbers of entries,
#' stochastically determined with draws from Poisson distributions. For each new
#' entry, a set of attributes is added for that node, and the nodes are added onto
#' the network objects. Only attributes that are a part of the network model
#' formulae are updated as vertex attributes on the network objects.
#'
#' @return
#' This function updates the \code{attr} list with new attributes for each new
#' population member, and the \code{nw} objects with new vertices.
#'
#' @keywords module msm
#' @export
#'
arrival_msm <- function(dat, at) {

  ## Input
  # Attributes

  # Parameters
  a.rate   <- get_param(dat, "a.rate")
  netstats <- get_param(dat, "netstats")

  ## Process
  num <- get_epi(dat, "num", at = 1)

  nNew <- rpois(1, a.rate * num)

  ## Update Attr
  if (nNew > 0) {
    dat <- setNewAttr_msm(dat, at, nNew)
  }

  # Update Networks
  if (nNew > 0) {
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::add_vertices(dat$el[[i]], nNew)
    }
  }

  ## Output
  dat <- set_epi(dat, "nNew", at, nNew)

  return(dat)
}


setNewAttr_msm <- function(dat, at, nNew) {

  dat <- append_core_attr(dat, at, nNew)

  ## Inputs
  # Attributes
  active <- get_attr(dat, "active")

  # Parameters
  netstats        <- get_param(dat, "netstats")
  tt.partial.supp <- get_param(dat, "tt.partial.supp")
  tt.full.supp    <- get_param(dat, "tt.full.supp")
  tt.durable.supp <- get_param(dat, "tt.durable.supp")
  circ.prob       <- get_param(dat, "circ.prob")
  arrival.age     <- get_param(dat, "arrival.age")
  stitest.active.int     <- get_param(dat, "stitest.active.int")


  # Set all attributes NA by default
  attrToInit <- names(which(vapply(dat$attr, length, 0) < length(active)))
  for (attr_name in attrToInit) {
    attr_val <- get_attr(dat, attr_name)
    attr_val <- c(attr_val, rep(NA, nNew))
    dat <- set_attr(dat, attr_name, attr_val, override.length.check = TRUE)
  }
  newIds <- which(get_attr(dat, "entrTime") == at)

  race.dist <- prop.table(table(netstats$attr$race))

  race <- sample(sort(unique(get_attr(dat, "race"))), nNew, TRUE, race.dist)
  dat <- set_attr(dat, "race", race, posit_ids = newIds)

  dat <- set_attr(dat, "age", rep(arrival.age, nNew), posit_ids = newIds)
  age.breaks <- netstats$demog$age.breaks
  attr_age.grp <- cut(
    get_attr(dat, "age", posit_ids = newIds),
    age.breaks,
    labels = FALSE,
    right = FALSE
  )
  dat <- set_attr(dat, "age.grp", attr_age.grp, posit_ids = newIds)

  # Disease status and related
  dat <- set_attr(dat, "status", 0, posit_ids = newIds)
  dat <- set_attr(dat, "diag.status", 0, posit_ids = newIds)
  dat <- set_attr(dat, "rGC", 0, posit_ids = newIds)
  dat <- set_attr(dat, "uGC", 0, posit_ids = newIds)
  dat <- set_attr(dat, "rCT", 0, posit_ids = newIds)
  dat <- set_attr(dat, "uCT", 0, posit_ids = newIds)
  dat <- set_attr(dat, "GC.timesInf", 0, posit_ids = newIds)
  dat <- set_attr(dat, "GC.timesInf", 0, posit_ids = newIds)
  dat <- set_attr(dat, "CT.timesInf", 0, posit_ids = newIds)
  dat <- set_attr(dat, "CT.timesInf", 0, posit_ids = newIds)

  # screening-related
  dat <- set_attr(dat, "recentpartners", 0, posit_ids = newIds)
  dat <- set_attr(dat, "stitest.ind.active", 0, posit_ids = newIds)
  dat <- set_attr(dat, "stitest.ind.recentpartners", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.rgc", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.rgc.asympt", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.ugc", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.ugc.asympt", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.gc", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.gc.asympt", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.sti", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.sti.asympt", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.rct", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.rct.asympt", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.uct", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.uct.asympt", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.ct", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.ct.asympt", 0, posit_ids = newIds)

  last.sti.test <- sample(seq_len(stitest.active.int), length(newIds), TRUE)

  dat <- set_attr(dat, "time.since.last.test.rgc", last.sti.test, posit_ids = newIds)
  dat <- set_attr(dat, "time.since.last.test.ugc", last.sti.test, posit_ids = newIds)
  dat <- set_attr(dat, "time.since.last.test.rct", last.sti.test, posit_ids = newIds)
  dat <- set_attr(dat, "time.since.last.test.uct", last.sti.test, posit_ids = newIds)

  dat <- set_attr(dat, "recentpartners", 0, posit_ids = newIds)
  dat <- set_attr(dat, "stitest.ind.active", 0, posit_ids = newIds)
  dat <- set_attr(dat, "stitest.ind.recentpartners", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.sti", 0, posit_ids = newIds)
  dat <- set_attr(dat, "testing.events.sti.asympt", 0, posit_ids = newIds)
  dat <- set_attr(dat, "count.trans", 0, posit_ids = newIds)


  rates <- get_param(dat, "hiv.test.late.prob")[race]
  late.tester <- runif(length(rates)) < rates
  dat <- set_attr(dat, "late.tester", late.tester, posit_ids = newIds)

  race.new <- get_attr(dat, "race", posit_ids = newIds)
  races <- sort(unique(race.new))
  tt.traj <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(get_attr(dat, "race")[newIds] == i)
    tt.traj[ids.race] <- sample(
      1:3, length(ids.race), TRUE,
      c(tt.partial.supp[i], tt.full.supp[i], tt.durable.supp[i])
    )

  }
  dat <- set_attr(dat, "tt.traj", tt.traj, posit_ids = newIds)

  # Circumcision
  circ <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(race.new == i)
    circ[ids.race] <- runif(length(ids.race)) < circ.prob[i]
  }
  dat <- set_attr(dat, "circ", circ, posit_ids = newIds)

  # Role
  ns <- netstats$attr
  role.class <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(race.new == i)
    rc.probs <- prop.table(table(ns$role.class[ns$race == i]))
    role.class[ids.race] <- sample(0:2, length(ids.race), TRUE, rc.probs)
  }
  dat <- set_attr(dat, "role.class", role.class, posit_ids = newIds)

  ins.quot <- rep(NA, nNew)
  role.class.new <- get_attr(dat, "role.class", posit_ids = newIds)
  ins.quot[role.class.new == 0]  <- 1
  ins.quot[role.class.new == 1]  <- 0
  ins.quot[role.class.new == 2]  <- runif(sum(role.class.new == 2))

  dat <- set_attr(dat, "ins.quot", ins.quot, posit_ids = newIds)

  # Degree
  dat <- set_attr(dat, "deg.main", 0, posit_ids = newIds)
  dat <- set_attr(dat, "deg.casl", 0, posit_ids = newIds)
  dat <- set_attr(dat, "deg.tot", 0, posit_ids = newIds)

  # One-off risk group
  dat <- set_attr(dat, "risk.grp", sample(1:5, nNew, TRUE), posit_ids = newIds)

  # PrEP
  dat <- set_attr(dat, "prepStat", 0, posit_ids = newIds)

  return(dat)
}
