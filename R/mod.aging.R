
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
#' @keywords module
#' @export
#'
aging.mard <- function(dat, at) {

  # system timer
  dat$epi$timer[at] <- proc.time()[3]

  time.unit <- dat$param$time.unit

  age <- dat$attr$age
  active <- dat$attr$active

  age[active == 1] <- age[active == 1] + time.unit / 365

  dat$attr$age <- age
  dat$attr$sqrt.age <- sqrt(age)

  for (i in seq_along(dat$nw)) {
    dat$nw[[i]] <- set.vertex.attribute(dat$nw[[i]], "sqrt.age", sqrt(age))
  }

  return(dat)
}
