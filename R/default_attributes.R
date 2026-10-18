#' Get the list the attributes used by the model with their default value
#'
#' @inheritParams aging_msm
#'
#' @details
#' This list must be exhaustive. All attributes should get a default value (even
#' NA) here.
#'
#' @return A named list of all the attributes with default values.
#'
#' TODO: document all the attributes
get_default_attrs <- function(dat) {
  list(

    ### OFFSCRIPT additions ----------------------------------------------------
    # Renal function
    gfr = 100,            # set by renal.R
    gfr.prep.baseline = NA,
    gfr.prep.ref = Inf, # reference value for GFR testing on PrEP
    prep.gfr.risk.last = -Inf,
    # HBV
    hbv = 0,
    hbv.flare.last = -Inf,

    tdf.resist = 0,
    ftc.resist = 0,
    ### OFFSCRIPT additions ----------------------------------------------------

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
    tx.status = NA,
    cuml.time.on.tx = NA,
    cuml.time.off.tx = NA,
    tx.init.time = NA,
    part.tx.init.time = NA,
    part.tx.reinit.time = NA,
    ## Gonorrhea
    gono.rect = 0,
    gono.uret = 0,
    gono.rect.sympt = 0,
    gono.uret.sympt = 0,
    gono.rect.inf.last = -Inf,
    gono.uret.inf.last = -Inf,
    gono.rect.inf.count = 0,
    gono.uret.inf.count = 0,
    gono.test.last = -Inf,
    gono.dx = 0,
    gono.dx.last = 0,
    gono.tx = 0,
    ## Chlamidya
    chla.rect = 0,
    chla.uret = 0,
    chla.rect.sympt = 0,
    chla.uret.sympt = 0,
    chla.rect.inf.last = -Inf,
    chla.uret.inf.last = -Inf,
    chla.rect.inf.count = 0,
    chla.uret.inf.count = 0,
    chla.test.last = -Inf,
    chla.dx = 0,
    chla.dx.last = -Inf,
    chla.tx = 0,
    ## Syphilis
    syph.inf = 0,
    syph.sympt = 0,
    syph.inf.last = -Inf,
    syph.inf.count = 0,
    syph.test.last = -Inf,
    syph.dx = 0,
    syph.dx.last = -Inf,
    syph.tx = 0,
    syph.stage = NA,
    syph.stage.last = -Inf,
    ## PrEP
    prep.any = 0,
    prep = 0,
    prep.class = -1, # invalid value, set upon PrEP start, using std or otc
    prep.std.class = 0, # set upon arrival
    prep.start.last = -Inf,
    prep.start.count = 0,
    prep.indic = 0,
    prep.risk.last = -Inf,
    prep.stop.last = -Inf,
    act.last = -Inf,
    uai.last = -Inf,
    usupp.partner.last = -Inf,
    unknown.partner.last = -Inf,
    # PrEP OTC
    prep.otc = 0,
    prep.otc.class = 0, # set upon arrival
    prep.otc.start.last = -Inf,
    prep.otc.start.count = 0,
    prep.otc.indic = 0,
    prep.otc.risk.last = -Inf,
    prep.otc.stop.last = -Inf,
    #
    # computed upon arrival
    race = 0,
    risk.grp = 0,
    late.tester = 0,
    circ = 0,
    age.grp = 0,
    ins.quot = 0,
    tt.traj = 0,
    last.neg.test = -Inf, # entrTime
    role.class = 0
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
  race_lvls <- as.numeric(names(race_dist))

  n_attr <- list()
  if (post_init) {
    n_attr$race <- sample(race_lvls, n_new, TRUE, race_dist)
    n_attr$role.class <- make_role_class(dat, n_attr$race, race_lvls)
    n_attr$risk.grp <- sample(seq_len(5), n_new, TRUE)
  } else {
    n_attr$race <- get_attr(dat, "race", posit_ids = new_nodes_pid)
    n_attr$role.class <- get_attr(dat, "role.class", posit_ids = new_nodes_pid)
    n_attr$risk.grp <- get_attr(dat, "risk.grp", posit_ids = new_nodes_pid)
  }

  age <- get_attr(dat, "age", posit_ids = new_nodes_pid)

  n_attr <- c(n_attr, list(
    late.tester = make_late_tester(dat, n_attr$race),
    circ        = make_circ(dat, n_attr$race, race_lvls),
    age.grp     = cut(age, age_breaks, labels = FALSE, right = FALSE),
    ins.quot    = make_ins_quot(n_attr$role.class),
    tt.traj     = make_tt_traj(dat, n_attr$race, race_lvls),
    prep.std.class   = make_prep_class(dat, length(n_attr$race)),
    prep.otc.class   = make_prep_otc_class(dat, length(n_attr$race)),
    last.neg.test = get_attr(dat, "entrTime", posit_ids = new_nodes_pid)
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
make_tt_traj <- function(dat, race, race_lvls) {
  partial <- get_param(dat, "tt.partial.supp.prob")
  full    <- get_param(dat, "tt.full.supp.prob")
  durable <- get_param(dat, "tt.durable.supp.prob")
  race.new <- race
  tt_traj <- numeric(length(race.new))
  for (r in race_lvls) {
    ids.race <- which(race == r)
    tt_traj[ids.race] <- sample(
      seq_len(3), length(ids.race), TRUE,
      c(partial[r], full[r], durable[r])
    )
  }
  return(tt_traj)
}

# Generate the `circ` attributes
make_circ <- function(dat, race, race_lvls) {
  hiv.circ.prob <- get_param(dat, "hiv.circ.prob")
  circ <- numeric(length(race))
  for (r in race_lvls) {
    ids.race <- which(race == r)
    circ[ids.race] <- runif(length(ids.race)) < hiv.circ.prob[r]
  }
  return(circ)
}

# Generate the `role.class` attributes
make_role_class <- function(dat, race, race_lvls) {
  ns <- get_param(dat, "netstats")$attr
  role_class <- numeric(length(race))
  for (r in race_lvls) {
    ids.race <- which(race == r)
    rc.probs <- prop.table(table(ns$role.class[ns$race == r]))
    role_class[ids.race] <- sample(0:2, length(ids.race), TRUE, rc.probs)
  }
  return(role_class)
}

make_prep_class <- function(dat, n_new) {
  prep.class <- sample(3, n_new, TRUE, get_param(dat, "prep.adhr.dist"))
  return(prep.class)
}

make_prep_otc_class <- function(dat, n_new) {
  prep.class <- sample(3, n_new, TRUE, get_param(dat, "prep.otc.adhr.dist"))
  return(prep.class)
}
