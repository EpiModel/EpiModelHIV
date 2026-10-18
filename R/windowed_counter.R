
#' Update a `counter_name` attribute that counts up events over a time window
#'
#' @inheritParams aging_msm
#'
#' @param posit_ids The positional IDs of the nodes whose counter should be incremented.
#' @param window_size The size of the counting window in time steps (default of `0` results in
#'        no counting, bypassing the function).
#' @param counter_name Character vector containing name of the counter attribute.
#'
#' @details
#' This function is developed to count up events over an arbitrarily defined time window. An example
#' is the number of STI diagnoses over a window of 3 months or 6 months or some other interval. That
#' count could then be used as a indication for interventions. The counting needs to have a vector
#' positional IDS whose counts should be updated (e.g., an ID list of current nodes with a positive
#' STI diagnosis), a window size in time steps, and a name of the count (e.g., `"sti.dx.count"`).
#'
#' @return
#' This function returns `dat` after updating the counter and the history list.
#'
update_windowed_counter <- function(dat, posit_ids, window_size, counter_name) {
  if (window_size == 0)
    return(dat)

  # Initialize the deque and reset the counters to 0 if the window is empty
  if (is.null(dat[["run"]][[counter_name]])) {
    dat[["run"]][[counter_name]] <- collections::deque()
    counter <- rep(0, length(get_attr(dat, "active")))
  } else {
    counter <- get_attr(dat, counter_name)
  }

  counter[posit_ids] <- counter[posit_ids] + 1

  dat[["run"]][[counter_name]]$push(get_unique_ids(dat, posit_ids))

  if (dat[["run"]][[counter_name]]$size() > window_size) {
    # decrement counters when the window "slides"
    uids_old_dx <- dat[["run"]][[counter_name]]$popleft()
    uids_old_dx <- uids_old_dx[is_active_unique_ids(dat, uids_old_dx)]
    pid_old_dx <- get_posit_ids(dat, uids_old_dx)
    counter[pid_old_dx] <- counter[pid_old_dx] - 1
  }

  dat <- set_attr(dat, counter_name, counter)

  return(dat)
}
