#' Update list `x` using the elements of list `new_x`
#'
#' @param x a list
#' @param new_x a list
#'
#' @return the full `x` list with the modifications added by `new_x`
#'
#' @details
#' This function updates list `x` by name. If `x` and `new_x` elements are not
#' named, the function will not work properly.
#' If a function is provided to replace an element that was originaly not a
#' function, this function will be applied to the original value.
update_list <- function(x, new_x) {
  for (n in names(new_x)) {
    if (is.list(new_x[[n]])) {
      x[[n]] <- update_list(x[[n]], new_x[[n]])
    } else if (is.function(new_x[[n]]) && ! is.function(x[[n]])) {
      x[[n]] <- new_x[[n]](x[[n]])
    } else {
      x[[n]] <- new_x[[n]]
    }
  }

  return(x)
}

#' Module to modify the `param` list at specified time steps during the simulation
#'
#' @inheritParams aging_msm
#'
#' @details
#' if a list `dat$param$param_updaters` is present, this module will update the
#' `param` list with new values at given timesteps.
#' `dat$control$param_updaters` is a list containing `updaters`. An updater is a
#' list containing an `at` element telling when the changes will happend, an
#' optional `verbose` boolean controlling whether to output a message when a
#' change is made (default = TRUE) and a `param` list similar
#' to the simulation's `dat$param` list.
#' if the updater is a function but not the value to replace, the
#' function will be applied to the current element (see example) .
#'
#' @examples
#'  ## example of a `param_updaters` list
#'  list(
#'    list(
#'      at = 4860,
#'      param = list(
#'        hiv.test.rate = rep(0.0128, 3),
#'        trans.scale = c(1.61, 0.836, 0.622)
#'      )
#'    ),
#'    list(
#'      at = 5160,
#'      verbose = FALSE,
#'      param = list(
#'        hiv.test.rate = function(x) x * 3,
#'        trans.scale = function(x) x^2 / 3
#'      )
#'    )
#'  )
#'
param_updater <- function(dat, at) {
  if (is.null(dat$param$param_updaters))
    return(dat)

  param_updaters <- dat$param$param_updaters

  for (i in seq_along(param_updaters)) {
    if (at == param_updaters[[i]][["at"]]) {
      verbose <- param_updaters[[i]][["verbose"]]
      verbose <- if (is.null(verbose)) TRUE else verbose

      new_params <- param_updaters[[i]][["param"]]

      if (verbose) {
        message(paste0(
          "\n\nAt time step = ", at, " the `param` list was modified: \n"))
        message(str(new_params))
      }

      dat$param <- update_list(dat$param, new_params)
    }
  }

  return(dat)
}
