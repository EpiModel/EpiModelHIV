
# MSM -----------------------------------------------------------------

#' @title Initialization Module
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param_msm}}.
#' @param init An \code{EpiModel} object of class \code{\link{init_msm}}.
#' @param control An \code{EpiModel} object of class \code{\link{control_msm}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated \code{dat} object with the initialized values
#' for demographics and disease-related variables.
#'
#' @export
#' @keywords module msm
#'
initialize_msm <- function(x, param, init, control, s) {

  ## Master Data List Setup ##
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control


  ## Network Setup ##
  # Initial network simulations
  dat$nw <- list()
  for (i in 1:3) {
    dat$nw[[i]] <- simulate(x[[i]]$fit, basis = x[[i]]$fit$newnetwork)
  }
  nw <- dat$nw

  # Pull Network parameters
  dat$nwparam <- list()
  for (i in 1:3) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }

  # Convert to tergmLite method
  dat <- init_tergmLite(dat)

  ## Nodal Attributes Setup ##
  dat$attr <- param$netstats$attr

  num <- network.size(nw[[1]])
  dat <- append_core_attr(dat, 1, num)

  # Circumcision
  circ.prob <- get_param(dat, "circ.prob")
  race <- get_attr(dat, "race")
  dat <- set_attr(dat, "circ", runif(num) < circ.prob[race])

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
  dat <- set_attr(dat, "part.art", rep(NA, num))
  dat <- set_attr(dat, "part.ident.counter", rep(NA, num))

  ## Other Setup ##
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()
  dat$epi <- list()

  # Prevalence Tracking
  # dat$temp$max.uid <- num
  dat <- prevalence_msm(dat, at = 1)

  # Setup Partner List
  dat <- init_plist(dat)

  # Network statistics
  if (dat$control$save.nwstats == TRUE) {
    dat <- calc_nwstats(dat, at = 1)
  }

  # dat$param$netstats <- NULL
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
#' @keywords initiation utility msm
#'
init_status_msm <- function(dat) {

  # Attributes
  active <- get_attr(dat, "active")
  # Sub in diag.status from model for status
  status <- get_attr(dat, "diag.status")
  race <- get_attr(dat, "race")
  age <- get_attr(dat, "age")

  # Parameters
  netstats <- get_param(dat, "netstats")
  hiv.test.late.prob <- get_param(dat, "hiv.test.late.prob")
  hiv.test.rate <- get_param(dat, "hiv.test.rate")
  tt.partial.supp <- get_param(dat, "tt.partial.supp")
  tt.full.supp <- get_param(dat, "tt.full.supp")
  tt.durable.supp <- get_param(dat, "tt.durable.supp")
  vl.acute.rise.int <- get_param(dat, "vl.acute.rise.int")
  vl.acute.fall.int <- get_param(dat, "vl.acute.fall.int")
  vl.aids.onset.int <- get_param(dat, "vl.aids.onset.int")
  vl.set.point <- get_param(dat, "vl.set.point")

  num <- sum(active)
  demog_ages <- netstats$demog$ages

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
      prob = c(tt.partial.supp[i], tt.full.supp[i], tt.durable.supp[i])
    )
  }
  dat <- set_attr(dat, "tt.traj", tt.traj)

  ## Infection-related attributes
  dat <- set_attr(dat, "status", status)
  idsInf <- which(status == 1)

  min.ages <- min(demog_ages)
  nweeks <- 365 / 5
  time.sex.active <- pmax(1, round(nweeks * age[idsInf] - nweeks * min.ages, 0))
  min.hiv.time <- round(vl.acute.rise.int + vl.acute.fall.int)
  max.hiv.time <- vl.aids.onset.int

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

  stage.time <- time.infected - min.hiv.time
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

  dat <- set_attr(dat, "inf.role", rep(NA, num))
  dat <- set_attr(dat, "inf.type", rep(NA, num))

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
#' @keywords initiation utility msm
#'
init_sti_msm <- function(dat) {

  # Attributes
  active <- get_attr(dat, "active")
  role.class <- get_attr(dat, "role.class")

  # Parameters
  rgc.sympt.prob     <- get_param(dat, "rgc.sympt.prob")
  ugc.sympt.prob     <- get_param(dat, "ugc.sympt.prob")
  rct.sympt.prob     <- get_param(dat, "rct.sympt.prob")
  uct.sympt.prob     <- get_param(dat, "uct.sympt.prob")
  stitest.active.int <- get_param(dat, "stitest.active.int")

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

  # PAF Specific attributes
  dat <- set_attr(dat, "diag.status.gc", rep(NA, num))
  dat <- set_attr(dat, "diag.status.gc", 0, posit_ids = uGC | rGC)

  dat <- set_attr(dat, "last.diag.time.gc", rep(NA, num))
  dat <- set_attr(dat, "last.neg.test.rgc", rep(NA, num))
  dat <- set_attr(dat, "last.neg.test.ugc", rep(NA, num))
  dat <- set_attr(dat, "last.tx.time.rgc", rep(NA, num))
  dat <- set_attr(dat, "last.tx.time.ugc", rep(NA, num))
  dat <- set_attr(dat, "last.tx.time.rgc.prep", rep(NA, num))
  dat <- set_attr(dat, "last.tx.time.ugc.prep", rep(NA, num))

  dat <- set_attr(dat, "tt.traj.gc.hivpos", rep(NA, num))
  dat <- set_attr(dat, "tt.traj.gc.hivneg", rep(NA, num))
  dat <- set_attr(dat, "testing.events.rgc", rep(0, num))
  dat <- set_attr(dat, "testing.events.rgc.asympt", rep(0, num))
  dat <- set_attr(dat, "testing.events.ugc", rep(0, num))
  dat <- set_attr(dat, "testing.events.ugc.asympt", rep(0, num))
  dat <- set_attr(dat, "testing.events.gc", rep(0, num))
  dat <- set_attr(dat, "testing.events.gc.asympt", rep(0, num))
  dat <- set_attr(dat, "testing.events.sti", rep(0, num))
  dat <- set_attr(dat, "testing.events.sti.asympt", rep(0, num))

  # Initialize CT infection at both sites
  idsUCT <- sample(idsUreth, round(prev.uct * num), FALSE)
  uCT[idsUCT] <- 1

  idsRCT <- sample(setdiff(idsRect, idsUCT), round(prev.rct * num), FALSE)
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

  # Initialize GC infection at both sites
  dat <- set_attr(dat, "diag.status.ct", rep(NA, num))
  dat <- set_attr(dat, "diag.status.ct", 0, posit_ids = uCT | rCT)

  dat <- set_attr(dat, "last.diag.time.ct", rep(NA, num))
  dat <- set_attr(dat, "last.neg.test.rct", rep(NA, num))
  dat <- set_attr(dat, "last.neg.test.uct", rep(NA, num))
  dat <- set_attr(dat, "last.tx.time.rct", rep(NA, num))
  dat <- set_attr(dat, "last.tx.time.uct", rep(NA, num))
  dat <- set_attr(dat, "last.tx.time.rct.prep", rep(NA, num))
  dat <- set_attr(dat, "last.tx.time.uct.prep", rep(NA, num))

  dat <- set_attr(dat, "tt.traj.ct.hivpos", rep(NA, num))
  dat <- set_attr(dat, "tt.traj.ct.hivneg", rep(NA, num))
  dat <- set_attr(dat, "testing.events.rct", rep(0, num))
  dat <- set_attr(dat, "testing.events.rct.asympt", rep(0, num))
  dat <- set_attr(dat, "testing.events.uct", rep(0, num))
  dat <- set_attr(dat, "testing.events.uct.asympt", rep(0, num))
  dat <- set_attr(dat, "testing.events.ct", rep(0, num))
  dat <- set_attr(dat, "testing.events.ct.asympt", rep(0, num))

  # Testing attributes
  last.sti.test <- sample(seq_len(stitest.active.int), num, TRUE)

  dat <- set_attr(dat, "time.since.last.test.rgc", last.sti.test)
  dat <- set_attr(dat, "time.since.last.test.ugc", last.sti.test)
  dat <- set_attr(dat, "time.since.last.test.rct", last.sti.test)
  dat <- set_attr(dat, "time.since.last.test.uct", last.sti.test)

  dat <- set_attr(dat, "time.last.sex", rep(NA, num))
  dat <- set_attr(dat, "recentpartners", rep(0, num))
  dat <- set_attr(dat, "stianntestLastElig", rep(NA, num))
  dat <- set_attr(dat, "stihighrisktestLastElig", rep(NA, num))
  dat <- set_attr(dat, "stitest.ind.active", rep(0, num))
  dat <- set_attr(dat, "stitest.ind.recentpartners", rep(0, num))
  dat <- set_attr(dat, "testing.events.sti", rep(0, num))
  dat <- set_attr(dat, "testing.events.sti.asympt", rep(0, num))
  dat <- set_attr(dat, "count.trans", rep(0, num))

  # EPT variables
  dat <- set_attr(dat, "eptindexElig", rep(NA, num))
  dat <- set_attr(dat, "eptindexStat", rep(NA, num))
  dat <- set_attr(dat, "eptindexEligdate", rep(NA, num))
  dat <- set_attr(dat, "eptpartEligReceive", rep(NA, num))
  dat <- set_attr(dat, "eptpartEligTx_GC", rep(NA, num))
  dat <- set_attr(dat, "eptpartEligTx_CT", rep(NA, num))
  dat <- set_attr(dat, "eptpartEligTxdate", rep(NA, num))
  dat <- set_attr(dat, "eptpartTx", rep(NA, num))
  dat <- set_attr(dat, "rCT.tx.ept", rep(NA, num))
  dat <- set_attr(dat, "uCT.tx.ept", rep(NA, num))
  dat <- set_attr(dat, "rGC.tx.ept", rep(NA, num))
  dat <- set_attr(dat, "uGC.tx.ept", rep(NA, num))

  dat <- set_attr(dat, "tt.traj.syph.hivpos", rep(NA, num))
  dat <- set_attr(dat, "tt.traj.syph.hivneg", rep(NA, num))

  return(dat)

}


#' @title Re-Initialization Module
#'
#' @description This function reinitializes an epidemic model to restart at a
#'              specified time step given an input \code{netsim} object.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}}.
#' @inheritParams initialize_msm
#'
#' @details
#' Currently, the necessary components that must be on \code{x} for a simulation
#' to be restarted must be: param, control, nwparam, epi, attr, temp, el, p.
#' TODO: describe this more.
#'
#' @return
#' This function resets the data elements on the \code{dat} master data object
#' in the needed ways for the time loop to function.
#'
#' @export
#' @keywords module msm
#'
reinit_msm <- function(x, param, init, control, s) {

  need.for.reinit <- c(
    "param",
    "control",
    "nwparam",
    "epi",
    "attr",
    "temp",
    "el",
    "p"
  )
  if (!all(need.for.reinit %in% names(x))) {
    stop("x must contain the following elements for restarting: ",
         "param, control, nwparam, epi, attr, temp, el, p",
         call. = FALSE)
  }

  if (length(x$el) == 1) {
    s <- 1
  }

  dat <- list()

  dat$param <- param
  dat$param$modes <- 1

  dat$control <- control
  dat$nwparam <- x$nwparam

  dat$epi <- sapply(x$epi, function(var) var[s])
  names(dat$epi) <- names(x$epi)

  dat$el <- x$el[[s]]
  dat$p <- x$p[[s]]

  dat$attr <- x$attr[[s]]

  if (!is.null(x$stats)) {
    dat$stats <- list()
    if (!is.null(x$stats$nwstats)) {
      dat$stats$nwstats <- x$stats$nwstats[[s]]
    }
  }

  dat$temp <- x$temp[[s]]

  class(dat) <- "dat"
  return(dat)
}
