# Trackers ---------------------------------------------------------------------
epi_n <- function(races_set) {
  function(dat) {
    needed_attributes <- c("race", "active")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% races_set & active == 1, na.rm = TRUE)
    })
  }
}

# eligible to prep
epi_s_prep_elig <- function(races_set) {
  function(dat) {
    needed_attributes <- c("race", "status", "prepElig")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% races_set & status == 0 & prepElig == 1, na.rm = TRUE)
    })
  }
}

# on prep
epi_s_prep <- function(races_set) {
  function(dat) {
    needed_attributes <- c("race", "status", "prepStat")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% races_set & status == 0 & prepStat == 1, na.rm = TRUE)
    })
  }
}

# HIV infected
epi_i <- function(races_set) {
  function(dat) {
    needed_attributes <- c("race", "status")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% races_set & status == 1, na.rm = TRUE)
    })
  }
}

# HIV infected diagnosed
epi_i_dx <- function(races_set) {
  function(dat) {
    needed_attributes <- c("race", "status", "diag.status")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% races_set & status == 1 & diag.status == 1, na.rm = TRUE)
    })
  }
}

# HIV infected virally suppressed
epi_i_sup <- function(races_set) {
  function(dat) {
    at <- get_current_timestep(dat)
    needed_attributes <- c("race", "status", "vl.last.supp")
    with(get_attr_list(dat, needed_attributes), {
      sum(race %in% races_set & status == 1 & vl.last.supp == at, na.rm = TRUE)
    })
  }
}

# linked in less than `steps` step
epi_linked_time <- function(steps) {
  function(races_set) {
    function(dat) {
      needed_attributes <- c("race", "tx.init.time", "diag.time")
      with(get_attr_list(dat, needed_attributes), {
        sum(
          race %in% races_set &
            tx.init.time - diag.time <= steps,
          na.rm = TRUE
        )
      })
    }
  }
}

epi_med_linked_delay <- function(races_set) {
  function(dat) {
    needed_attributes <-
      c("race", "diag.status", "diag.time", "tx.init.time")
    with(get_attr_list(dat, needed_attributes), {
      pop <- race %in% races_set & diag.status == 1
      tx_delay[pop] <-
        ifelse(is.na(tx.init.time[pop]), Inf, tx.init.time[pop]) - diag.time
      median(tx_delay)
    })
  }
}

# STI trackers
epi_gc_s <- function(hiv_status) {
  function(races_set) {
    function(dat) {
      needed_attributes <- c("race", "rGC", "uGC", "status")
      with(get_attr_list(dat, needed_attributes), {
        sum(
          race %in% races_set &
            status %in% hiv_status &
            (rGC == 0 & uGC == 0),
          na.rm = TRUE
        )
      })
    }
  }
}

epi_ct_s <- function(hiv_status) {
  function(races_set) {
    function(dat) {
      needed_attributes <- c("race", "rCT", "uCT", "status")
      with(get_attr_list(dat, needed_attributes), {
        sum(
          race %in% races_set &
            status %in% hiv_status &
            (rCT == 0 & uCT == 0),
          na.rm = TRUE
        )
      })
    }
  }
}

epi_prep_ret <- function(ret_steps) {
  function(races_set) {
    function(dat) {
      at <- get_current_timestep(dat)
      needed_attributes <- c("race", "prepStartTime")
      with(get_attr_list(dat, needed_attributes), {
        retained <- sum(
          race %in% races_set &
          prepStartTime == at - ret_steps,
          na.rm = TRUE
        )
      })
    }
  }
}

#' Make a List of Trackers for Model Calibration
#'
#' This function produces a list of trackers to be passed to the `.tracker.list`
#' argument of `EpiModel::control.net()`. These trackers are used for
#' calibration purposes
#'
#' @return A List of Trackers
#'
#' @export
make_calibration_trackers <- function() {
  list(
    n__B                = epi_n(1),
    n__H                = epi_n(2),
    n__W                = epi_n(3),
    n__O                = epi_n(4),
    i__B                = epi_i(1),
    i__H                = epi_i(2),
    i__W                = epi_i(3),
    i__O                = epi_i(4),
    i_dx__B             = epi_i_dx(1),
    i_dx__H             = epi_i_dx(2),
    i_dx__W             = epi_i_dx(3),
    i_dx__O             = epi_i_dx(4),
    i_sup__B            = epi_i_sup(1),
    i_sup__H            = epi_i_sup(2),
    i_sup__W            = epi_i_sup(3),
    i_sup__O            = epi_i_sup(4),
    linked1m__B         = epi_linked_time(4)(1),
    linked1m__H         = epi_linked_time(4)(2),
    linked1m__W         = epi_linked_time(4)(3),
    linked1m__O         = epi_linked_time(4)(4),
    gc_s__B             = epi_gc_s(c(0, 1))(1), # gc susceptible HIV+ and -
    gc_s__H             = epi_gc_s(c(0, 1))(2),
    gc_s__W             = epi_gc_s(c(0, 1))(3),
    gc_s__O             = epi_gc_s(c(0, 1))(4),
    ct_s__B             = epi_ct_s(c(0, 1))(1),
    ct_s__H             = epi_ct_s(c(0, 1))(2),
    ct_s__W             = epi_ct_s(c(0, 1))(3),
    ct_s__O             = epi_ct_s(c(0, 1))(4),
    s_prep__B           = epi_s_prep(1),
    s_prep__H           = epi_s_prep(2),
    s_prep__W           = epi_s_prep(3),
    s_prep__O           = epi_s_prep(4),
    s_prep_elig__B      = epi_s_prep_elig(1),
    s_prep_elig__H      = epi_s_prep_elig(2),
    s_prep_elig__W      = epi_s_prep_elig(3),
    s_prep_elig__O      = epi_s_prep_elig(4),
    prep_startat        = epi_prep_ret(0)(1:3),  # n starting at `at`
    prep_ret1y          = epi_prep_ret(52)(1:3), # n starting at `at - 52`
    prep_ret2y          = epi_prep_ret(104)(1:3) # n starting at `at - 104`
  )
}

#' Mutate the Calibration Variables
#'
#' This function is to be used on the `data.frame` version of a `netsim` object.
#' It mutate the variables that are to be matched with calibration targets
#'
#' @param d A `data.frame` representation of a `netsim` object
#' @param year_steps How many time steps constitute a year (default = 52)
#'
#' @return a `tibble` with the calibration variables added
#'
#' @export
mutate_calibration_targets <- function(d, year_steps = 52) {
  dplyr::mutate(d,
    cc.dx.B         = .data$i_dx__B / .data$i__B,
    cc.dx.H         = .data$i_dx__H / .data$i__H,
    cc.dx.W         = .data$i_dx__W / .data$i__W,
    cc.dx.O         = .data$i_dx__O / .data$i__O,
    cc.linked1m.B   = .data$linked1m__B / .data$i_dx__B,
    cc.linked1m.H   = .data$linked1m__H / .data$i_dx__H,
    cc.linked1m.W   = .data$linked1m__W / .data$i_dx__W,
    cc.linked1m.O   = .data$linked1m__O / .data$i_dx__O,
    cc.vsupp.B      = .data$i_sup__B / .data$i_dx__B,
    cc.vsupp.H      = .data$i_sup__H / .data$i_dx__H,
    cc.vsupp.W      = .data$i_sup__W / .data$i_dx__W,
    cc.vsupp.O      = .data$i_sup__O / .data$i_dx__O,
    gc_s            = .data$gc_s__B +
                      .data$gc_s__H +
                      .data$gc_s__W +
                      .data$gc_s__O,
    ir100.gc        = .data$incid.gc / .data$gc_s * year_steps * 100,
    ct_s            = .data$ct_s__B +
                      .data$ct_s__H +
                      .data$ct_s__W +
                      .data$ct_s__O,
    ir100.ct        = .data$incid.ct / .data$ct_s * year_steps * 100,
    i.prev.dx.B     = .data$i_dx__B / .data$n__B,
    i.prev.dx.H     = .data$i_dx__H / .data$n__H,
    i.prev.dx.W     = .data$i_dx__W / .data$n__W,
    i.prev.dx.O     = .data$i_dx__O / .data$n__O,
    cc.prep.B       = .data$s_prep__B / .data$s_prep_elig__B,
    cc.prep.H       = .data$s_prep__H / .data$s_prep_elig__H,
    cc.prep.W       = .data$s_prep__W / .data$s_prep_elig__W,
    cc.prep.O       = .data$s_prep__O / .data$s_prep_elig__O,
    prep_users      = .data$s_prep__B +
                      .data$s_prep__H +
                      .data$s_prep__W +
                      .data$s_prep__O,
    prep_elig       = .data$s_prep_elig__B +
                      .data$s_prep_elig__H +
                      .data$s_prep_elig__W +
                      .data$s_prep_elig__O,
    cc.prep         = .data$prep_users / .data$prep_elig,
    prep_prop_ret1y = .data$prep_ret1y /
                      dplyr::lag(.data$prep_startat, year_steps),
    prep_prop_ret2y = .data$prep_ret2y /
                      dplyr::lag(.data$prep_startat, year_steps * 2)
  )
}

#' Return the calibration targets for the current model
#'
#' @return a named vector of calibration targets
#'
#' @export
get_calibration_targets <- function() {
  c(
    # 1st calibration set (all independant)
    cc.dx.B         = 0.851,
    cc.dx.H         = 0.817,
    cc.dx.W         = 0.863,
    cc.dx.O         = 0.828,
    cc.linked1m.B   = 0.881,
    cc.linked1m.H   = 0.874,
    cc.linked1m.W   = 0.873,
    cc.linked1m.O   = 0.962,
    # CombPrev appendix 8.2.2
    # 2nd calibration set (all independant)
    cc.vsupp.B      = 0.645,
    cc.vsupp.H      = 0.643,
    cc.vsupp.W      = 0.685,
    cc.vsupp.O      = 0.614,
    # STIs
    ir100.gc        = 8.06,
    ir100.ct        = 8.727,
    # 3rd calibration set
    i.prev.dx.B     = 0.215,
    i.prev.dx.H     = 0.260,
    i.prev.dx.W     = 0.108,
    i.prev.dx.O     = 0.083,
    # i.prev.dx       = 0.17, # not used yet but maybe?
    # cc.prep.ind.B   = 0.387,
    # cc.prep.ind.H   = 0.379,
    # cc.prep.ind.W   = 0.407,
    # cc.prep.ind     = 0.402,
    cc.prep.B       = 0.319,
    cc.prep.H       = 0.368,
    cc.prep.W       = 0.516,
    cc.prep.O       = 0.258,
    # cc.prep         = 0.203,
    # prep_prop_ret1y = 0.56, # DOI: 10.1002/jia2.25252 and 10.1093/cid/ciaa037 (54%)
    # prep_prop_ret2y = 0.41,
    disease.mr100   = 0.241
  )
}

#' Mutate the Calibration Distances From the Targets
#'
#' This function is to be used on the `data.frame` version of a `netsim` object.
#' It mutate the variables that are to be matched with calibration targets. It
#' then substract for each of this variable the target value.
#'
#' @param d A `data.frame` representation of a `netsim` object
#' @param year_steps How many time steps constitute a year (default = 52)
#'
#' @return a `tibble` with the calibration distances added
#' @export
mutate_calibration_distances <- function(d, year_steps = 52) {
  d <- mutate_calibration_targets(d, year_steps)
  targets <- get_calibration_targets()
  for (t_name in intersect(names(targets), names(d)))
    d[[t_name]] <- d[[t_name]] - targets[[t_name]]
  return(d)
}
