#' Module to modify the parameters at specified time steps during the simulation
#'
#' @inheritParams aging_msm
#'
#' @details
#' if a list `dat$param$param_updaters` is present, this module will update the
#' `dat$param` list with new values at given timesteps.
#' `dat$param$param_updaters` is a list containing `updaters`. An updater is a
#' list containing an `at` element telling when the changes will happen and a
#' `param` element which is a named list of the parameters to change.
#'
#' @examples
#'  ## example of a `param_updaters` list
#'  list(
#'    list(
#'      at = 4860,
#'      params = list(
#'        hiv.test.rate = rep(0.0128, 3),
#'        trans.scale = c(1.61, 0.836, 0.622)
#'      )
#'    ),
#'    list(
#'      at = 5160,
#'      params = list(tx.init.prob = c(0.125, 0.158, 0.164))
#'    )
#'  )
#'
param_changer <- function(dat, at) {
  if (is.null(dat$param$param_updaters))
    return(dat)

  param_updaters <- dat$param$param_updaters

  for (i in seq_along(param_updaters)) {
    if (at == param_updaters[[i]][["at"]]) {
      new_params <- param_updaters[[i]][["params"]]

      message(paste0(
        "\n\nAt time step = ", at, " the following parameters were modified: \n",
        paste(names(new_params), collapse = ", ")
      ))

      dat$param[names(new_params)] <- new_params
    }
  }

  return(dat)
}

param_updaters <- list(
  list(at = 4860, params = list(
    hiv.test.rate = rep(0.0128, 3),
    trans.scale = c(1.61, 0.836, 0.622)
  )),
  list(at = 5160, params = list(
    hiv.test.rate = rep(0.012, 3),
    trans.scale = c(1.6, 0.8, 0.6)
  ))
)
