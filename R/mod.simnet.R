
#' @title Network Resimulation Module
#'
#' @description Module function for resimulating the sexual networks for one
#'              time step.
#'
#' @inheritParams aging_msm
#'
#' @export
#'
simnet_msm <- function(dat, at) {

  ## Parameters
  cumulative.edgelist <- get_control(dat, "cumulative.edgelist")
  truncate.el.cuml <- get_control(dat, "truncate.el.cuml")
  set.control.tergm <- get_control(dat, "set.control.tergm")
  nwstats.formulas <- get_control(dat, "nwstats.formulas")
  tergmLite.track.duration <- get_control(dat, "tergmLite.track.duration")

  ## Edges correction
  dat <- edges_correct_msm(dat, at)

  ## Main network
  nwparam <- EpiModel::get_nwparam(dat, network = 1)

  dat <- set_attr(dat, "deg.casl", EpiModel::get_degree(dat$el[[2]]))

  nwL <- networkLite(dat[["el"]][[1]], dat[["attr"]])
  if (tergmLite.track.duration == TRUE) {
    nwL %n% "time" <- dat$nw[[1]] %n% "time"
    nwL %n% "lasttoggle" <- dat$nw[[1]] %n% "lasttoggle"
  }

  dat[["nw"]][[1]] <- simulate(
    nwL ~ Form(nwparam[["formation"]]) +
          Persist(nwparam[["coef.diss"]][["dissolution"]]),
    coef = c(nwparam[["coef.form"]], nwparam[["coef.diss"]][["coef.adj"]]),
    constraints = nwparam[["constraints"]],
    time.start = at - 1,
    time.slices = 1,
    time.offset = 1,
    monitor = nwstats.formulas[[1]],
    control = set.control.tergm,
    output = "final",
    dynamic = TRUE
  )

  dat[["el"]][[1]] <- as.edgelist(dat[["nw"]][[1]])

  ## Casual network
  nwparam <- EpiModel::get_nwparam(dat, network = 2)

  dat <- set_attr(dat, "deg.main", EpiModel::get_degree(dat$el[[1]]))

  nwL <- networkLite(dat[["el"]][[2]], dat[["attr"]])
  if (tergmLite.track.duration == TRUE) {
    nwL %n% "time" <- dat$nw[[2]] %n% "time"
    nwL %n% "lasttoggle" <- dat$nw[[2]] %n% "lasttoggle"
  }

  dat[["nw"]][[2]] <- simulate(
    nwL ~ Form(nwparam[["formation"]]) +
          Persist(nwparam[["coef.diss"]][["dissolution"]]),
    coef = c(nwparam[["coef.form"]], nwparam[["coef.diss"]][["coef.adj"]]),
    constraints = nwparam[["constraints"]],
    time.start = at - 1,
    time.slices = 1,
    time.offset = 1,
    monitor = nwstats.formulas[[2]],
    control = set.control.tergm,
    output = "final",
    dynamic = TRUE
  )

  dat[["el"]][[2]] <- as.edgelist(dat[["nw"]][[2]])

  ## One-off network
  nwparam <- EpiModel::get_nwparam(dat, network = 3)

  dat <- set_attr(dat, "deg.tot",
    pmin(get_attr(dat, "deg.main") + EpiModel::get_degree(dat[["el"]][[2]]), 3))

  nwL <- networkLite(dat[["el"]][[3]], dat[["attr"]])

  set.control.tergm$MCMC.prop.args <- list(discordance_fraction = 0)
  dat[["nw"]][[3]] <- simulate(
    basis = nwL,
    object = nwparam[["formation"]],
    coef = nwparam[["coef.form"]],
    constraints = nwparam[["constraints"]],
    monitor = nwstats.formulas[[3]],
    control = set.control.tergm,
    time.start = at - 1,
    time.slices = 1,
    time.offset = 1,
    dynamic = TRUE,
    output = "final"
  )

  dat[["el"]][[3]] <- as.edgelist(dat[["nw"]][[3]])

  if (get_control(dat, "save.nwstats")) {
    for (i in seq_len(3)) {
      new.nwstats <- tail(attributes(dat$nw[[i]])$stats, 1)
      keep.cols <- which(!duplicated(colnames(new.nwstats)))
      new.nwstats <- new.nwstats[, keep.cols, drop = FALSE]
      dat$stats$nwstats[[i]] <- rbind(dat$stats$nwstats[[i]], new.nwstats)
    }
  }

  if (cumulative.edgelist) {
    for (n_network in seq_len(3)) {
      dat <- update_cumulative_edgelist(dat, n_network, truncate.el.cuml)
    }
  }

  return(dat)
}

#' @title Adjustment for the Edges Coefficient with Changing Network Size
#'
#' @description Adjusts the edges coefficients in a dynamic network model
#'              to preserve the mean degree.
#'
#' @inheritParams aging_msm
#'
#' @details
#' In HIV/STI modeling, there is typically an assumption that changes in
#' population size do not affect one's number of partners, specified as the
#' mean degree for network models. A person would not have 10 times the number
#' of partners should he move from a city 10 times as large. This module uses
#' the adjustment of Krivitsky et al. to adjust the edges coefficients on the
#' three network models to account for varying population size in order to
#' preserve that mean degree.
#'
#' @return
#' The network model parameters stored in `dat$nwparam` are updated for
#' each of the three network models.
#'
#' @references
#' Krivitsky PN, Handcock MS, and Morris M. "Adjusting for network size and
#' composition effects in exponential-family random graph models." Statistical
#' Methodology. 2011; 8.4: 319-339.
#'
#' @keywords module msm
#'
#' @export
#'
edges_correct_msm <- function(dat, at) {

  old.num <- get_epi(dat, "num", at - 1)
  new.num <- sum(get_attr(dat, "active") == 1, na.rm = TRUE)
  adjust <- log(old.num) - log(new.num)

  coef.form.m <- get_nwparam(dat, network = 1)[["coef.form"]]
  coef.form.m[1] <- coef.form.m[1] + adjust
  dat[["nwparam"]][[1]][["coef.form"]] <- coef.form.m

  coef.form.p <- get_nwparam(dat, network = 2)[["coef.form"]]
  coef.form.p[1] <- coef.form.p[1] + adjust
  dat[["nwparam"]][[2]][["coef.form"]] <- coef.form.p

  coef.form.i <- get_nwparam(dat, network = 3)[["coef.form"]]
  coef.form.i[1] <- coef.form.i[1] + adjust
  dat[["nwparam"]][[3]][["coef.form"]] <- coef.form.i

  return(dat)
}
