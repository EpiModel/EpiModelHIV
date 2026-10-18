
#' @title Aging Module
#'
#' @description Module for aging over time for active nodes in the population.
#'
#' @param dat Master data list object of class \code{dat} containing networks,
#'        individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @return
#' This function returns \code{dat} after updating the nodal attribute
#' \code{age} and \code{sqrt.age}. The \code{sqrt.age} vertex attribute is also
#' updated on the three networks.
#'
#' @keywords module msm
#' @export
#'
aging_msm <- function(dat, at) {

  age <- dat$attr$age
  alive <- dat$attr$alive
  age.grp <- dat$attr$age.grp

  age[alive == 1] <- age[alive == 1] + 7 / 365

  age.breaks <- dat$param$netstats$demog$age.breaks
  age.grp[alive == 1] <- cut(age[alive == 1], age.breaks, labels = FALSE)

  dat$attr$age.grp <- age.grp
  dat$attr$age <- age

  return(dat)
}


#' @export
#' @rdname aging_msm
aging_het <- function(dat, at) {

  ## Parameters
  time.unit <- dat$param$time.unit

  ## Attributes
  age <- dat$attr$age
  active <- dat$attr$active

  ## Updates
  age[active == 1] <- age[active == 1] + time.unit/365

  ## Save out
  dat$attr$age <- age

  return(dat)
}
