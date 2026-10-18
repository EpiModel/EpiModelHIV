# The windows counters has been removed as they were not working correctly on
# the latest versions.
# To add them back, see the `update_windowed_counter` function on the
# "R/windowed_counter.R" file.
#
# It requires adding the counters in the `default_attributes`. With a 0 default
# value. The value for the counter for each node can then be retrieve as usual
# with `get_attr`

#' @title STI (Non HIV) Module
#'
#' @description Module function calling the function responsibles for
#'              transmission, progression, screening and treatment of Gonorrhea,
#'              Chlamydia and Syphilis.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This module orchestrate the STI epidemics: Gonorrhea, Chlamydia and
#' Syphilis. Similar to HIV, each disease has a transmission, progression,
#' screening and treatment step. These steps are handled by separate functions
#' called by this module in a specific order. See the documentation of each
#' subfunction for details.
#'
#' @return
#' This function updates the `dat` object for the STIs epidemics (Non HIV).
#'
#' @export
#'
sti_msm <- function(dat, at) {
  # Transmission ---------------------------------------------------------------
  # Gonorrhea
  dat <- gono_chla_transmission(dat, "gono", i_insertive = TRUE)
  dat <- gono_chla_transmission(dat, "gono", i_insertive = FALSE)
  dat <- set_epi(dat, "incid.gono", at, sum(
    get_attr(dat, "gono.rect.inf.last") == at |
      get_attr(dat, "gono.uret.inf.last") == at
  ))
  # Chlamydia
  dat <- gono_chla_transmission(dat, "chla", i_insertive = TRUE)
  dat <- gono_chla_transmission(dat, "chla", i_insertive = FALSE)
  dat <- set_epi(dat, "incid.chla", at, sum(
    get_attr(dat, "chla.rect.inf.last") == at |
      get_attr(dat, "chla.uret.inf.last") == at
  ))
  # Syphilis
  dat <- syphilis_transmission(dat)

  # Progression ----------------------------------------------------------------
  dat <- gono_chla_progression(dat, "gono")
  dat <- gono_chla_progression(dat, "chla")
  dat <- syphilis_progression(dat)
  # Testing --------------------------------------------------------------------
  dat <- gono_chla_testing(dat, "gono")
  dat <- gono_chla_testing(dat, "chla")
  dat <- syphilis_testing(dat)
  # Treatment ------------------------------------------------------------------
  dat <- gono_chla_treatment(dat, "gono")
  dat <- gono_chla_treatment(dat, "chla")
  dat <- syphilis_treatment(dat)

  return(dat)
}
