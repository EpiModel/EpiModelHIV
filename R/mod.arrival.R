
#' @title Arrivals Module
#'
#' @description Module function for arrivals into the sexually active population.
#'
#' @inheritParams aging_msm
#'
#' @details
#' New population members are added based on expected numbers of entries,
#' stochastically determined with draws from Poisson distributions. For each new
#' entry, a set of attributes is added for that node, and the nodes are added onto
#' the network objects.
#'
#' @return
#' This function updates the `attr` list with new attributes for each new population member.
#'
#' @export
#'
arrival_msm <- function(dat, at) {

  ## Input
  # Attributes

  # Parameters
  a.rate   <- get_param(dat, "a.rate")

  ## Process
  num <- get_epi(dat, "num", at = 1)

  nNew <- rpois(1, a.rate * num)

  ## Update Attr
  if (nNew > 0) {
    dat <- init_new_nodes_attrs(dat, nNew)
  }

  # Update Networks
  dat <- arrive_nodes(dat, nNew)

  ## Output
  dat <- set_epi(dat, "nNew", at, nNew)

  return(dat)
}

#' Set the attributes for new nodes
#'
#' @inheritParams aging_msm
#' @param n_new The number of new nodes to initialize
#'
#' @return
#' This function returns the `dat` object with updated attributes for the new
#' nodes.
init_new_nodes_attrs <- function(dat, n_new) {
  current_timestep <- get_current_timestep(dat)
  # Core attributes (necessary for EpiModel)
  dat <- append_core_attr(dat, current_timestep, n_new)
  def_attrs <- get_default_attrs(dat)
  for (attr_name in names(def_attrs)) {
    dat <- append_attr(dat, attr_name, def_attrs[[attr_name]], n_new)
  }
  dat <- make_computed_attrs(dat, n_new, post_init = TRUE)
  return(dat)
}

