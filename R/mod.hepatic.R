#' @title Hepatic Function Module
#'
#' @description Module function managing the hepatic function of MSMs in the
#' model. Currently only PrEP stoppage induced flare ups are implemented
#'
#' @inheritParams aging_msm
#'
#' @export
#'
hepatic_msm <- function(dat, at) {
  if (at == get_param(dat, "hbv.start")) {
    dat <- set_hbv_status(dat)
  }

  # Attributes
  hbv <- get_attr(dat, "hbv")
  prep.stop.last <- get_attr(dat, "prep.stop.last")
  prep.otc.stop.last <- get_attr(dat, "prep.otc.stop.last")
  # Parameters
  prep.hbv.flare.prob <- get_param(dat, "prep.hbv.flare.prob")
  prep.otc.hbv.flare.prob <- get_param(dat, "prep.otc.hbv.flare.prob")

  # NOTE: this module runs BEFORE PrEP, so must check for PrEP stop on previous
  # step

  # Pre-subset: all HBV positive
  hbv_ids <- which(hbv == 1)

  # PrEP related HBV flares
  prep_elig_ids <- hbv_ids[prep.stop.last[hbv_ids] == at - 1]
  prep_flare_ids <- which(runif(length(prep_elig_ids)) < prep.hbv.flare.prob)

  # PrEP OTC related HBV flares
  otc_elig_ids <- hbv_ids[prep.otc.stop.last[hbv_ids] == at - 1]
  otc_flare_ids <- which(runif(length(otc_elig_ids)) < prep.otc.hbv.flare.prob)

  hbv_flare_ids <- c(prep_flare_ids, otc_flare_ids)
  dat <- set_attr(dat, "hbv.flare.last", at, posit_ids = hbv_flare_ids)

  if (get_control(dat, "debug")) {
    dat <- set_epi(dat, "dbg_hbv_cur", at, sum(hbv))
    dat <- set_epi(dat, "dbg_stop_std_prep", at, length(prep_elig_ids))
    dat <- set_epi(dat, "dbg_hbv_flares_std", at, length(prep_flare_ids))
    dat <- set_epi(dat, "dbg_stop_otc_prep", at, length(otc_elig_ids))
    dat <- set_epi(dat, "dbg_hbv_flares_otc", at, length(otc_flare_ids))
    # HBV per age grp
    age <- get_attr(dat, "age")
    age_breaks <- c(0, 17, 49, 70)
    age_grps <- cut(age, age_breaks, right = FALSE)
    hbv_grp <- tapply(hbv, age_grps, mean)
    for (i in seq_along(hbv_grp)) {
      dat <- set_epi(
        dat,
        paste0("dbg_hbv_prop_", names(hbv_grp)[i]),
        at, hbv_grp[i]
      )
    }
  }

  return(dat)
}

set_hbv_status <- function(dat) {
  # attributes
  age <- get_attr(dat, "age")
  # parameters
  hbv.init.perc <- get_param(dat, "hbv.init.perc")

  age_breaks <- c(0, 17, 49, 70)
  age_grps <- cut(age, age_breaks, labels = FALSE, right = FALSE)
  probs <- hbv.init.perc[age_grps]

  hbv <- as.numeric(runif(length(probs)) < probs)
  dat <- set_attr(dat, "hbv", hbv)

  return(dat)
}