
#' @title Initialization Module
#'
#' @description This function initializes the main `netsim_dat` class data
#'              object, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class [`netest`].
#' @param param An \code{EpiModel} object of class [`param_msm`].
#' @param init An \code{EpiModel} object of class [`init_msm`].
#' @param control An \code{EpiModel} object of class [`control_msm`].
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the `netsim_dat` object with the initialized values for
#' demographics and disease-related variables.
#'
#' @export
#'
initialize_msm <- function(x, param, init, control, s) {

  ## Master Data List Setup ##
  dat <- create_dat_object(param, init, control)

  ## network and stats initialization
  dat <- init_nets(dat, x)

  ## Initialize all remaining attributes
  dat <- init_attrs(dat)
  dat <- overwrite_attrs(dat)

  # simulate first time step
  dat <- sim_nets_t1(dat)
  dat <- summary_nets(dat, at = 1L)

  num <- sum(get_attr(dat, "active") == 1)

  # Time Unit
  time.unit <- param$epistats$time.unit
  dat <- set_param(dat, "time.unit", time.unit)

  dat[["temp"]] <- list()

  # Prevalence Tracking
  dat <- set_epi(dat, "num", at = 1,  num)

  return(dat)
}


#' @title Initialize the HIV status of persons in the network
#'
#' @description Sets the initial individual-level disease status of persons
#'              in the network, as well as disease-related attributes for
#'              infected persons.
#'
#' @param dat `netsim_dat` class data object created in initialization module.
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

  dat <- set_attr(dat, "inf.time", -time.infected, posit_ids = idsInf)
  dat <- set_attr(dat, "stage", 3, posit_ids = idsInf)

  stage.time <- time.infected - (min.hiv.time - 1)
  dat <- set_attr(dat, "stage.time", stage.time, posit_ids = idsInf)

  stage <- get_attr(dat, "stage", posit_ids = idsInf)
  dat <- set_attr(dat, "diag.stage", stage, posit_ids = idsInf)

  dat <- set_attr(dat, "vl", vl.set.point, posit_ids = idsInf)

  inf.time <- get_attr(dat, "inf.time")
  inf.diag.time <- pmin(1, inf.time[idsInf] + round(mean(1 / hiv.test.rate)))
  dat <- set_attr(dat, "diag.time", inf.diag.time, posit_ids = idsInf)

  dat <- set_attr(dat, "tx.status", 0, posit_ids = idsInf)
  dat <- set_attr(dat, "cuml.time.on.tx", 0, posit_ids = idsInf)
  dat <- set_attr(dat, "cuml.time.off.tx", time.infected, posit_ids = idsInf)

  return(dat)
}



#' @title Initialize the STI status of persons in the network
#'
#' @description Sets the initial individual-level disease status of persons
#'              in the network, as well as disease-related attributes for
#'              infected persons.
#'
#' @param dat `netsim_dat` class data object created in initialization module.
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
  prev.syph <- get_init(dat, "prev.syph")

  num <- length(active)

  idsUreth <- which(role.class %in% c(0, 2))
  idsRect <- which(role.class %in% c(1, 2))

  uGC <- get_attr(dat, "uGC")
  rGC <- get_attr(dat, "rGC")
  uCT <- get_attr(dat, "uCT")
  rCT <- get_attr(dat, "rCT")
  syph <- get_attr(dat, "syph")

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

  rGC_ids <- which(as.logical(rGC))
  uGC_ids <- which(as.logical(uGC))

  rGC.sympt <- runif(sum(rGC)) < rgc.sympt.prob
  dat <- set_attr(dat, "rGC.sympt", rGC.sympt, posit_ids = rGC_ids)

  uGC.sympt <- runif(sum(uGC)) < ugc.sympt.prob
  dat <- set_attr(dat, "uGC.sympt", uGC.sympt, posit_ids = uGC_ids)

  dat <- set_attr(dat, "rGC.infTime", 1, posit_ids = rGC_ids)
  dat <- set_attr(dat, "uGC.infTime", 1, posit_ids = uGC_ids)

  dat <- set_attr(dat, "rGC.timesInf", 1, posit_ids = rGC_ids)
  dat <- set_attr(dat, "uGC.timesInf", 1, posit_ids = uGC_ids)


  gc.last.screen.time <- 1 - rgeom(num, get_param(dat, "gc.screen.hivneg.rate"))
  dat <- set_attr(dat, "gc.last.screen.time", gc.last.screen.time)

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

  rCT_ids <- which(as.logical(rCT))
  uCT_ids <- which(as.logical(uCT))

  rCT.sympt <- runif(sum(rCT)) < rct.sympt.prob
  dat <- set_attr(dat, "rCT.sympt", rCT.sympt, posit_ids = rCT_ids)

  uCT.sympt <- runif(sum(uCT)) < uct.sympt.prob
  dat <- set_attr(dat, "uCT.sympt", uCT.sympt, posit_ids = uCT_ids)

  dat <- set_attr(dat, "rCT.infTime", 1, posit_ids = rCT_ids)
  dat <- set_attr(dat, "uCT.infTime", 1, posit_ids = uCT_ids)

  dat <- set_attr(dat, "rCT.timesInf", 1, posit_ids = rCT_ids)
  dat <- set_attr(dat, "uCT.timesInf", 1, posit_ids = uCT_ids)


  ct.last.screen.time <- 1 - rgeom(num, get_param(dat, "ct.screen.hivneg.rate"))
  dat <- set_attr(dat, "ct.last.screen.time", ct.last.screen.time)

  # Initialize syphilis infection
  idsSyph <- sample(1:num, size = round(prev.syph * num), FALSE)
  syph[idsSyph] <- 1
  dat <- set_attr(dat, "syph", syph)

  dat <- set_attr(dat, "syph.stage", ifelse(syph == 1, 1, NA))
  dat <- set_attr(dat, "syph.stageTime", ifelse(syph == 1, 1, NA))
  dat <- set_attr(dat, "syph.sympt", ifelse(syph == 1, 0, NA))
  dat <- set_attr(dat, "syph.infTime", ifelse(syph == 1, 1, NA))
  dat <- set_attr(dat, "syph.timesInf", ifelse(syph == 1, 1, 0))

  syph.last.screen.time <- 1 - rgeom(num, get_param(dat, "syph.screen.hivneg.rate"))
  dat <- set_attr(dat, "syph.last.screen.time", syph.last.screen.time)

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
#' @details See \code{\link{initialize.net}}.
#'
#' @return See \code{\link{initialize.net}}.
#'
#' @export
#'
reinit_msm <- initialize.net

init_attrs <- function(dat) {
  n_nodes <- sum(get_attr(dat, "active") == 1)
  def_attrs <- get_default_attrs(dat)
  cur_attrs <- get_attr_list(dat)
  missing <- names(def_attrs)[!names(def_attrs) %in% names(cur_attrs)]
  for (attr_name in missing) {
    dat <- append_attr(dat, attr_name, def_attrs[[attr_name]], n_nodes)
  }

  dat <- make_computed_attrs(dat, n_nodes, post_init = FALSE)

  # HIV-related attributes
  dat <- init_status_msm(dat)

  # STI Status
  dat <- init_sti_msm(dat)

  return(dat)
}


#' @title Helper to use a `data.frame` to initialize some attributes
#'
#' @description Uses `dat$init$init_attr` to overwrite some attributes of the
#' nodes at initialization
#'
#' @details
#' If an `init_attr` `data.frame` is present in `dat$init`, use it to overwrite
#' the attributes it contains.
#' `init_attr` must have a number of rows equal to the number of nodes in the
#' model as the attributes will be overwritten one to one, ensuring the correct
#' ordering.
#' `init_attr` columns MUST have a corresponding attribute already initialized.
#' See "R/default_attributes.R" for adding new attributes to the model.
#' `init_attr` is removed from `dat$init` at the end of the function to free up
#' its memory.
#'
#' @inheritParams aging_msm
#'
overwrite_attrs <- function(dat) {
  init_attr <- get_init(dat, "init_attr", override.null.error = TRUE)
  if (is.null(init_attr)) {
    return(dat)
  }
  message("init_attr used for initialization of attributes")

  status <- get_attr(dat, "status")
  if (nrow(init_attr) != length(status)) {
    stop("init_attr should contains the same number of nodes as the model")
  }

  new_attrs <- setdiff(names(init_attr), names(get_attr_list(dat)))
  if (length(new_attrs) > 0) {
    stop(
      "Some attributes in `init_attr` are not present in `dat`: ",
      paste0(new_attrs, collapse = " ,")
    )
  }

  core_attrs <- c("active", "entrTime", "exitTime", "unique_id")
  for (attr_name in setdiff(names(init_attr), core_attrs)) {
    dat <- set_attr(dat, attr_name, init_attr[[attr_name]])
  }

  dat$init$init_attr <- NULL
  dat
}
