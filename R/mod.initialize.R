
#' @title Initialization Module
#'
#' @description This function initializes the master `dat` object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class [`netest`].
#' @param param An \code{EpiModel} object of class [`param_msm`].
#' @param init An \code{EpiModel} object of class [`init_msm`].
#' @param control An \code{EpiModel} object of class [`control_msm`].
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated `dat` object with the initialized values
#' for demographics and disease-related variables.
#'
#' @export
#'
initialize_msm <- function(x, param, init, control, s) {

  ## Master Data List Setup ##
  dat <- create_dat_object(param, init, control)

  ## Network Setup ##
  # Initial network simulations
  dat[["nw"]] <- list()
  for (i in 1:3) {
    dat[["nw"]][[i]] <- simulate(x[[i]][["formula"]],
      coef = x[[i]][["coef.form.crude"]],
      basis = x[[i]][["newnetwork"]],
      constraints = x[[i]][["constraints"]],
      control = get_control(dat, "set.control.ergm"),
      dynamic = FALSE
    )
  }
  nw <- dat[["nw"]]

  # Pull Network parameters
  dat[["nwparam"]] <- list()
  for (i in 1:3) {
    dat[["nwparam"]][i] <- list(x[[i]][!(names(x[[i]]) %in% c("fit", "newnetwork"))])
  }

  # Convert to tergmLite method
  dat <- init_tergmLite(dat)

  ## Nodal Attributes Setup ##
  dat[["attr"]] <- param[["netstats"]][["attr"]]

  num <- network.size(nw[[1]])
  dat <- append_core_attr(dat, 1, num)

  # Time Unit
  time.unit <- param$epistats$time.unit
  dat <- set_param(dat, "time.unit", time.unit)

  # Circumcision
  hiv.circ.prob <- get_param(dat, "hiv.circ.prob")
  race <- get_attr(dat, "race")
  dat <- set_attr(dat, "circ", runif(num) < hiv.circ.prob[race])

  # Insertivity Quotient
  ins.quot <- rep(NA, num)
  role.class <- get_attr(dat, "role.class")
  ins.quot[role.class == 0]  <- 1
  ins.quot[role.class == 1]  <- 0
  ins.quot[role.class == 2]  <- runif(sum(role.class == 2))
  dat <- set_attr(dat, "ins.quot", ins.quot)

  # HIV-related attributes
  dat <- init_status_msm(dat)

  # STI Status
  dat <- init_sti_msm(dat)

  # PrEP-related attributes
  dat <- set_attr(dat, "prepClass", rep(NA, num))
  dat <- set_attr(dat, "prepElig", rep(NA, num))
  dat <- set_attr(dat, "prepStat", rep(0, num))
  dat <- set_attr(dat, "prepStartTime", rep(NA, num))
  dat <- set_attr(dat, "prepLastRisk", rep(NA, num))
  dat <- set_attr(dat, "prepLastStiScreen", rep(NA, num))
  dat <- set_attr(dat, "prep.start.counter", rep(NA, num))

  # Partner Identification attributes
  dat <- set_attr(dat, "part.scrnd", rep(NA, num))
  dat <- set_attr(dat, "part.ident", rep(NA, num))

  #PSY3 additions
  dat <- set_attr(dat, "foundPrevPos.ident", rep(NA, num))
  dat <- set_attr(dat, "part.art", rep(NA, num))
  dat <- set_attr(dat, "part.ident.counter", rep(NA, num))
  dat <- set_attr(dat, "last.pos.test", rep(NA, num))


  ## Other Setup ##
  dat[["stats"]] <- list()
  dat[["stats"]][["nwstats"]] <- list()

  dat[["temp"]] <- list()
  dat[["epi"]] <- list()

  # Prevalence Tracking
  dat <- set_epi(dat, "num", at = 1,  num)

  # Setup Partner List
  for (n_network in seq_len(3)) {
    dat <- update_cumulative_edgelist(dat, n_network)
  }

  # Network statistics
  if (get_control(dat, "save.nwstats")) {
    for (i in seq_len(3)) {
      new.nwstats <- attributes(dat$nw[[i]])$stats
      keep.cols <- which(!duplicated(colnames(new.nwstats)))
      new.nwstats <- new.nwstats[, keep.cols, drop = FALSE]
      dat$stats$nwstats[[i]] <- new.nwstats
    }
  }

  class(dat) <- "dat"
  return(dat)
}


#' @title Initialize the HIV status of persons in the network
#'
#' @description Sets the initial individual-level disease status of persons
#'              in the network, as well as disease-related attributes for
#'              infected persons.
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#'
init_status_msm <- function(dat) {

  # Attributes
  active <- get_attr(dat, "active")
  # Sub in diag.status from model for status
  status <- get_attr(dat, "diag.status")
  race   <- get_attr(dat, "race")
  age    <- get_attr(dat, "age")

  # Parameters
  netstats             <- get_param(dat, "netstats")
  hiv.test.late.prob   <- get_param(dat, "hiv.test.late.prob")
  hiv.test.rate        <- get_param(dat, "hiv.test.rate")
  tt.partial.supp.prob <- get_param(dat, "tt.partial.supp.prob")
  tt.full.supp.prob    <- get_param(dat, "tt.full.supp.prob")
  tt.durable.supp.prob <- get_param(dat, "tt.durable.supp.prob")
  vl.acute.rise.int    <- get_param(dat, "vl.acute.rise.int")
  vl.acute.fall.int    <- get_param(dat, "vl.acute.fall.int")
  vl.aids.onset.int    <- get_param(dat, "vl.aids.onset.int")
  vl.set.point         <- get_param(dat, "vl.set.point")

  num <- sum(active)
  demog_ages <- netstats[["demog"]][["ages"]]

  # Late (AIDS-stage) tester type
  rates <- hiv.test.late.prob[race]
  dat <- set_attr(dat, "late.tester", runif(num) < rates)

  # Treatment trajectory
  tt.traj <- rep(NA, num)
  races <- sort(unique(race))
  for (i in races) {
    ids.race <- which(race == i)
    tt.traj[ids.race] <- sample(
      x = 1:3,
      size = length(ids.race),
      replace = TRUE,
      prob = c(tt.partial.supp.prob[i], tt.full.supp.prob[i], tt.durable.supp.prob[i])
    )
  }
  dat <- set_attr(dat, "tt.traj", tt.traj)

  ## Infection-related attributes
  dat <- set_attr(dat, "status", status)
  idsInf <- which(status == 1)

  min.ages <- min(demog_ages)
  nweeks <- 365 / 7
  time.sex.active <- pmax(1, round(nweeks * age[idsInf] - nweeks * min.ages, 0))
  min.hiv.time <- round(vl.acute.rise.int + vl.acute.fall.int)
  max.hiv.time <- vl.aids.onset.int - 2

  time.infected <- round(
    pmax(min.hiv.time,
      pmin(time.sex.active,
        sample(min.hiv.time:max.hiv.time, length(idsInf), TRUE)
      )
    )
  )

  dat <- set_attr(dat, "inf.time", rep(NA, num))
  dat <- set_attr(dat, "inf.time", -time.infected, posit_ids = idsInf)

  dat <- set_attr(dat, "stage", rep(NA, num))
  dat <- set_attr(dat, "stage.time", rep(NA, num))
  dat <- set_attr(dat, "aids.time", rep(NA, num))
  dat <- set_attr(dat, "stage", 3, posit_ids = idsInf)

  stage.time <- time.infected - (min.hiv.time - 1)
  dat <- set_attr(dat, "stage.time", stage.time, posit_ids = idsInf)

  dat <- set_attr(dat, "diag.stage", rep(NA, num))

  stage <- get_attr(dat, "stage", posit_ids = idsInf)
  dat <- set_attr(dat, "diag.stage", stage, posit_ids = idsInf)

  dat <- set_attr(dat, "vl", rep(NA, num))
  dat <- set_attr(dat, "vl", vl.set.point, posit_ids = idsInf)
  dat <- set_attr(dat, "vl.last.usupp", rep(NA, num))
  dat <- set_attr(dat, "vl.last.supp", rep(NA, num))

  dat <- set_attr(dat, "diag.time", rep(NA, num))

  inf.time <- get_attr(dat, "inf.time")
  inf.diag.time <- pmin(1, inf.time[idsInf] + round(mean(1 / hiv.test.rate)))
  dat <- set_attr(dat, "diag.time", inf.diag.time, posit_ids = idsInf)

  dat <- set_attr(dat, "last.neg.test", rep(NA, num))

  dat <- set_attr(dat, "tx.status", rep(NA, num))
  dat <- set_attr(dat, "tx.status", 0, posit_ids = idsInf)
  dat <- set_attr(dat, "cuml.time.on.tx", rep(NA, num))
  dat <- set_attr(dat, "cuml.time.on.tx", 0, posit_ids = idsInf)
  dat <- set_attr(dat, "cuml.time.off.tx", rep(NA, num))
  dat <- set_attr(dat, "cuml.time.off.tx", time.infected, posit_ids = idsInf)
  dat <- set_attr(dat, "tx.init.time", rep(NA, num))

  dat <- set_attr(dat, "part.tx.init.time", rep(NA, num))
  dat <- set_attr(dat, "part.tx.reinit.time", rep(NA, num))
  dat <- set_attr(dat, "prep.start.part", rep(NA, num))

  return(dat)
}



#' @title Initialize the STI status of persons in the network
#'
#' @description Sets the initial individual-level disease status of persons
#'              in the network, as well as disease-related attributes for
#'              infected persons.
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#'
init_sti_msm <- function(dat) {

  # Attributes
  active     <- get_attr(dat, "active")
  role.class <- get_attr(dat, "role.class")

  # Parameters
  rgc.sympt.prob <- get_param(dat, "rgc.sympt.prob")
  ugc.sympt.prob <- get_param(dat, "ugc.sympt.prob")
  rct.sympt.prob <- get_param(dat, "rct.sympt.prob")
  uct.sympt.prob <- get_param(dat, "uct.sympt.prob")

  # Init values
  prev.ugc <- get_init(dat, "prev.ugc")
  prev.rgc <- get_init(dat, "prev.rgc")
  prev.uct <- get_init(dat, "prev.uct")
  prev.rct <- get_init(dat, "prev.rct")

  num <- length(active)

  idsUreth <- which(role.class %in% c(0, 2))
  idsRect <- which(role.class %in% c(1, 2))

  uGC <- rep(0, num)
  rGC <- rep(0, num)
  uCT <- rep(0, num)
  rCT <- rep(0, num)

  # Initialize GC infection at both sites
  idsUGC <- sample(idsUreth, size = round(prev.ugc * num), FALSE)
  uGC[idsUGC] <- 1

  idsRGC <- sample(
    setdiff(idsRect, idsUGC),
    size = round(prev.rgc * num),
    FALSE
  )
  rGC[idsRGC] <- 1

  dat <- set_attr(dat, "rGC", rGC)
  dat <- set_attr(dat, "uGC", uGC)

  dat <- set_attr(dat, "rGC.sympt", rep(NA, num))
  dat <- set_attr(dat, "uGC.sympt", rep(NA, num))

  rGC.sympt <- runif(sum(rGC)) < rgc.sympt.prob
  dat <- set_attr(dat, "rGC.sympt", rGC.sympt, posit_ids = as.logical(rGC))

  uGC.sympt <- runif(sum(uGC)) < ugc.sympt.prob
  dat <- set_attr(dat, "uGC.sympt", uGC.sympt, posit_ids = as.logical(uGC))

  dat <- set_attr(dat, "rGC.infTime", rep(NA, num))
  dat <- set_attr(dat, "uGC.infTime", rep(NA, num))
  dat <- set_attr(dat, "rGC.infTime", 1, posit_ids = as.logical(rGC))
  dat <- set_attr(dat, "uGC.infTime", 1, posit_ids = as.logical(uGC))

  dat <- set_attr(dat, "rGC.timesInf", rep(0, num))
  dat <- set_attr(dat, "rGC.timesInf", 1, posit_ids = as.logical(rGC))
  dat <- set_attr(dat, "uGC.timesInf", rep(0, num))
  dat <- set_attr(dat, "uGC.timesInf", 1, posit_ids = as.logical(uGC))

  dat <- set_attr(dat, "rGC.tx", rep(NA, num))
  dat <- set_attr(dat, "uGC.tx", rep(NA, num))
  dat <- set_attr(dat, "rGC.tx.prep", rep(NA, num))
  dat <- set_attr(dat, "uGC.tx.prep", rep(NA, num))

  # Initialize CT infection at both sites
  idsUCT <- sample(idsUreth, size = round(prev.uct * num), FALSE)
  uCT[idsUCT] <- 1

  idsRCT <- sample(
    setdiff(idsRect, idsUCT),
    size = round(prev.rct * num),
    FALSE
  )
  rCT[idsRCT] <- 1

  dat <- set_attr(dat, "rCT", rCT)
  dat <- set_attr(dat, "uCT", uCT)

  dat <- set_attr(dat, "rCT.sympt", rep(NA, num))
  dat <- set_attr(dat, "uCT.sympt", rep(NA, num))

  rCT.sympt <- runif(sum(rCT)) < rct.sympt.prob
  dat <- set_attr(dat, "rCT.sympt", rCT.sympt, posit_ids = as.logical(rCT))

  uCT.sympt <- runif(sum(uCT)) < uct.sympt.prob
  dat <- set_attr(dat, "uCT.sympt", uCT.sympt, posit_ids = as.logical(uCT))

  dat <- set_attr(dat, "rCT.infTime", rep(NA, num))
  dat <- set_attr(dat, "uCT.infTime", rep(NA, num))
  dat <- set_attr(dat, "rCT.infTime", 1, posit_ids = as.logical(rCT))
  dat <- set_attr(dat, "uCT.infTime", 1, posit_ids = as.logical(uCT))

  dat <- set_attr(dat, "rCT.timesInf", rep(0, num))
  dat <- set_attr(dat, "rCT.timesInf", 1, posit_ids = as.logical(rCT))
  dat <- set_attr(dat, "uCT.timesInf", rep(0, num))
  dat <- set_attr(dat, "uCT.timesInf", 1, posit_ids = as.logical(uCT))

  dat <- set_attr(dat, "rCT.tx", rep(NA, num))
  dat <- set_attr(dat, "uCT.tx", rep(NA, num))
  dat <- set_attr(dat, "rCT.tx.prep", rep(NA, num))
  dat <- set_attr(dat, "uCT.tx.prep", rep(NA, num))

  return(dat)

}


#' @title Re-Initialization Module
#'
#' @description This function reinitializes an epidemic model to restart at a
#'              specified time step given an input `netsim` object.
#'
#' @param x An `EpiModel` object of class `netsim`.
#' @inheritParams initialize_msm
#'
#' @details
#' Currently, the necessary components that must be on `x` for a simulation
#' to be restarted must be: param, control, nwparam, epi, attr, temp, el.
#'
#' @return
#' This function resets the data elements on the `dat` master data object
#' in the needed ways for the time loop to function.
#'
#' @export
#'
reinit_msm <- function(x, param, init, control, s) {

  need.for.reinit <- c(
    "param",
    "control",
    "nwparam",
    "epi",
    "attr",
    "temp",
    "el"
  )
  if (!all(need.for.reinit %in% names(x))) {
    stop("x must contain the following elements for restarting: ",
         "param, control, nwparam, epi, attr, temp, el, p",
         call. = FALSE)
  }

  if (length(x[["el"]]) == 1) {
    s <- 1
  }

  dat <- list()

  dat[["param"]] <- param
  dat[["param"]][["modes"]] <- 1

  dat[["control"]] <- control
  # If missing (mcmc.control), `tergmLite::simulate_network` segfault in
  # `tergm_MCMC_slave`
  dat[["control"]][["mcmc.control"]] <- x[["control"]][["mcmc.control"]]
  dat[["control"]][["nwstats.formulas"]] <- x[["control"]][["nwstats.formulas"]]

  dat[["nwparam"]] <- x[["nwparam"]]

  dat[["epi"]] <- sapply(x[["epi"]], function(var) var[s])
  names(dat[["epi"]]) <- names(x[["epi"]])

  dat[["el"]] <- x[["el"]][[s]]
  dat[["el.cuml"]] <- x[["el.cuml"]][[s]]
  dat[["_last_unique_id"]] <- x[["_last_unique_id"]][[s]]

  dat[["attr"]] <- x[["attr"]][[s]]

  if (!is.null(x[["stats"]])) {
    dat[["stats"]] <- list()
    if (!is.null(x[["stats"]][["nwstats"]])) {
      dat[["stats"]][["nwstats"]] <- x[["stats"]][["nwstats"]][[s]]
    }
  }

  dat[["temp"]] <- x[["temp"]][[s]]

  class(dat) <- "dat"

  # Time Unit
  time.unit <- param$epistats$time.unit
  dat <- set_param(dat, "time.unit", time.unit)

  return(dat)
}
