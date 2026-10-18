
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
  dat$run$num <- num
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
  # Sub in diag.status from model for status
  status <- get_attr(dat, "diag.status")
  age    <- get_attr(dat, "age")

  # Parameters
  netstats             <- get_param(dat, "netstats")
  hiv.test.rate        <- get_param(dat, "hiv.test.rate")
  vl.acute.rise.int    <- get_param(dat, "vl.acute.rise.int")
  vl.acute.fall.int    <- get_param(dat, "vl.acute.fall.int")
  vl.aids.onset.int    <- get_param(dat, "vl.aids.onset.int")
  vl.set.point         <- get_param(dat, "vl.set.point")

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
  dat <- gono_chla_init(dat, "gono")
  dat <- gono_chla_init(dat, "chla")
  dat <- syph_init(dat)

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
