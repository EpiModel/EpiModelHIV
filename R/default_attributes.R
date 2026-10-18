#' Get the list the attributes used by the model with their default value
#'
#' @inheritParams aging_msm
#'
#' @details
#' This list must be exhaustive. All attributes should get a default value (even
#' NA) here.
#'
#' @return A named list of all the attributes with default values.
get_default_attrs <- function(dat) {
  list(
    status = 0,
    age = get_param(dat, "arrival.age"),
    sqrt.age = NA,
    active.sex = 1,
    deg.casl = 0,
    deg.main = 0,
    deg.tot = 0,
    diag.status = 0,
    inf.time = NA,
    stage = NA,
    stage.time = NA,
    aids.time = NA,
    diag.stage = NA,
    vl = NA,
    vl.last.usupp = NA,
    vl.last.supp = NA,
    diag.time = NA,
    last.neg.test = NA,
    tx.status = NA,
    cuml.time.on.tx = NA,
    cuml.time.off.tx = NA,
    tx.init.time = NA,
    part.tx.init.time = NA,
    part.tx.reinit.time = NA,
    prep.start.part = NA,
    rGC = 0,
    uGC = 0,
    rGC.sympt = NA,
    uGC.sympt = NA,
    rGC.infTime = NA,
    uGC.infTime = NA,
    rGC.timesInf = 0,
    uGC.timesInf = 0,
    GC.tx = NA,
    diag.status.gc = NA,
    diag.status.gc.prep = NA,
    gc.last.screen.time = NA,
    gc.last.dxpos.time = NA,
    syph = 0,
    syph.timesInf = 0,
    syph.stage = NA,
    syph.stageTime = NA,
    syph.sympt = NA,
    syph.infTime = NA,
    syph.tx = NA,
    diag.status.syph = NA,
    diag.status.syph.prep = NA,
    syph.last.screen.time = 0,
    syph.last.dxpos.time = NA,
    syph.dx.count.window = 0,
    rCT = 0,
    uCT = 0,
    rCT.sympt = NA,
    uCT.sympt = NA,
    rCT.infTime = NA,
    uCT.infTime = NA,
    rCT.timesInf = 0,
    uCT.timesInf = 0,
    CT.tx = NA,
    diag.status.ct = NA,
    diag.status.ct.prep = NA,
    ct.last.screen.time = NA,
    ct.last.dxpos.time = NA,
    gc.dx.count.window = 0,
    ct.dx.count.window = 0,
    sti.dx.count.window = 0,
    prepClass = NA,
    prepElig = 0,
    prepStat = 0,
    prepStartTime = NA,
    prepLastRisk = NA,
    prep.start.counter = NA,
    prep.stop.counter = NA,
    part.scrnd = NA,
    part.ident = NA,
    part.ident.counter = NA,
    prep.ind.uai.mono = NA,
    prep.ind.uai.nmain = NA,
    prep.ind.sti = NA,
    prep.ind.uai.conc = NA,
    # computed
    race = 0,
    role.class = 0,
    risk.grp = 0,
    late.tester = 0,
    circ = 0,
    age.grp = 0,
    ins.quot = 0,
    tt.traj = 0,

    # New attributes
    prepDuration = NA,
    prepStopDuration = NA,
    prepStopTime = NA
  )
}

#' Update the attributes requiring computation for new nodes
#'
#' @inheritParams aging_msm
#' @param n_new The number of new nodes to update
#' @param post_init logical flag, TRUE if not called by the initialization
#' module.
#'
#' @details
#' This function takes care of all the attributes requiring some computations.
#' This includes the random assignment for race or the calculation of the
#' age groups. New attributes that need special assignment should be set here
#' AS WELL AS in the `get_default_attrs` function.
#'
#' @return
#' This function returns the `dat` object with updated attributes for the new
#' nodes.
make_computed_attrs <- function(dat, n_new, post_init) {
  new_nodes_pid <- length(get_attr(dat, "active")) - n_new + seq_len(n_new)
  ns <- get_param(dat, "netstats")
  age_breaks <- ns$demog$age.breaks
  race_dist <- prop.table(table(ns$attr$race))
  races <- as.numeric(names(race_dist))

  n_attr <- list()
  if (post_init) {
    n_attr$race <- sample(races, n_new, TRUE, race_dist)
    n_attr$role.class <- make_role_class(dat, n_attr$race, races)
    n_attr$risk.grp <- sample(seq_len(5), n_new, TRUE)
  } else {
    n_attr$race <- get_attr(dat, "race", posit_ids = new_nodes_pid)
    n_attr$role.class <- get_attr(dat, "role.class", posit_ids = new_nodes_pid)
    n_attr$risk.grp <- get_attr(dat, "risk.grp", posit_ids = new_nodes_pid)
  }

  age <- get_attr(dat, "age", posit_ids = new_nodes_pid)

  n_attr <- c(n_attr, list(
    late.tester   = make_late_tester(dat, n_attr$race),
    circ          = make_circ(dat, n_attr$race, races),
    age.grp       = cut(age, age_breaks, labels = FALSE, right = FALSE),
    ins.quot      = make_ins_quot(n_attr$role.class),
    tt.traj       = make_tt_traj(dat, n_attr$race, races),
    stopper.speed = make_stopper_speed(dat, n_attr$race)
  ))

  for (attr_name in names(n_attr)) {
    dat <- set_attr(dat, attr_name, n_attr[[attr_name]], posit_ids = new_nodes_pid)
  }

  return(dat)
}

# Generate the `late.tester` attributes
make_late_tester <- function(dat, race) {
  rates <- get_param(dat, "hiv.test.late.prob")[race]
  runif(length(rates)) < rates
}

# Generate the `ins.quot` attributes
make_ins_quot <- function(role.class) {
  ins_quot <- numeric(length(role.class))
  ins_quot[role.class == 0]  <- 1
  ins_quot[role.class == 1]  <- 0
  ins_quot[role.class == 2]  <- runif(sum(role.class == 2))
  return(ins_quot)
}

# Generate the `tt.traj` attributes
make_tt_traj <- function(dat, race, races) {
  partial <- get_param(dat, "tt.partial.supp.prob")
  full    <- get_param(dat, "tt.full.supp.prob")
  durable <- get_param(dat, "tt.durable.supp.prob")
  race.new <- race
  tt_traj <- numeric(length(race.new))
  for (r in races) {
    ids.race <- which(race == r)
    tt_traj[ids.race] <- sample(
      seq_len(3), length(ids.race), TRUE,
      c(partial[r], full[r], durable[r])
    )
  }
  return(tt_traj)
}

# Generate the `circ` attributes
make_circ <- function(dat, race, races) {
  hiv.circ.prob <- get_param(dat, "hiv.circ.prob")
  circ <- numeric(length(race))
  for (r in races) {
    ids.race <- which(race == r)
    circ[ids.race] <- runif(length(ids.race)) < hiv.circ.prob[r]
  }
  return(circ)
}

# Generate the `role.class` attributes
make_role_class <- function(dat, race, races) {
  ns <- get_param(dat, "netstats")$attr
  role_class <- numeric(length(race))
  for (r in races) {
    ids.race <- which(race == r)
    rc.probs <- prop.table(table(ns$role.class[ns$race == r]))
    role_class[ids.race] <- sample(0:2, length(ids.race), TRUE, rc.probs)
  }
  return(role_class)
}

# Generate stopper.speed --------------------------------------------------
make_stopper_speed <- function(dat, race) {
  fast.stopper.prob <- get_param(dat, "fast.stopper.prob")
  prob <- fast.stopper.prob[race]
  stopper.speed <- as.numeric(runif(length(race)) < prob)
  return(stopper.speed)
}
