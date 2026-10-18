# MSM -----------------------------------------------------------------

#' @title Network Resimulation Module
#'
#' @description Module function for resimulating the sexual networks for one
#'              time step.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
simnet_msm <- function(dat, at) {

  ## Parameters
  truncate.plist <- get_param(dat, "truncate.plist")

  ## Edges correction
  dat <- edges_correct_msm(dat, at)

  ## Main network
  nwparam.m <- EpiModel::get_nwparam(dat, network = 1)

  dat <- set_attr(dat, "deg.casl", EpiModel::get_degree(dat$el[[2]]))
  dat <- tergmLite::updateModelTermInputs(dat, network = 1)

  dat$el[[1]] <- tergmLite::simulate_network(
    p = dat$p[[1]],
    el = dat$el[[1]],
    coef.form = nwparam.m$coef.form,
    coef.diss = nwparam.m$coef.diss$coef.adj,
    save.changes = TRUE
  )

  plist1 <- update_plist(dat, at, ptype = 1)


  ## Casual network
  nwparam.p <- EpiModel::get_nwparam(dat, network = 2)

  dat <- set_attr(dat, "deg.main", EpiModel::get_degree(dat$el[[1]]))
  dat <- tergmLite::updateModelTermInputs(dat, network = 2)

  dat$el[[2]] <- tergmLite::simulate_network(
    p = dat$p[[2]],
    el = dat$el[[2]],
    coef.form = nwparam.p$coef.form,
    coef.diss = nwparam.p$coef.diss$coef.adj,
    save.changes = TRUE
  )

  plist2 <- update_plist(dat, at, ptype = 2)

  ## One-off network
  nwparam.i <- EpiModel::get_nwparam(dat, network = 3)

  dat <- set_attr(dat, "deg.tot",
    pmin(get_attr(dat, "deg.main") + EpiModel::get_degree(dat$el[[2]]), 3))
  dat <- tergmLite::updateModelTermInputs(dat, network = 3)

  dat$el[[3]] <- tergmLite::simulate_ergm(
    p = dat$p[[3]],
    el = dat$el[[3]],
    coef = nwparam.i$coef.form
  )

  plist3 <- update_plist(dat, at, ptype = 3)

  dat$temp$plist <- rbind(plist1, plist2, plist3)
  if (truncate.plist != Inf) {
    to.keep.stop <- at - dat$temp$plist[, "stop"]
    to.keep <- which(is.na(dat$temp$plist[, "stop"]) | to.keep.stop < truncate.plist)
    dat$temp$plist <- dat$temp$plist[to.keep, , drop = FALSE]
  }

  if (get_control(dat, "save.nwstats") == TRUE) {
    dat <- calc_nwstats(dat, at)
  }

  return(dat)
}

calc_nwstats <- function(dat, at) {

  for (nw in 1:3) {
    n <- attr(dat$el[[nw]], "n")
    edges <- nrow(dat$el[[nw]])
    meandeg <- round(edges * (2/n), 3)
    concurrent <- round(mean(EpiModel::get_degree(dat$el[[nw]]) > 1), 3)
    mat <- matrix(c(edges, meandeg, concurrent), ncol = 3, nrow = 1)
    if (at == 1) {
      dat$stats$nwstats[[nw]] <- mat
      colnames(dat$stats$nwstats[[nw]]) <- c("edges", "mdeg", "conc")
    }
    if (at > 1) {
      dat$stats$nwstats[[nw]] <- rbind(dat$stats$nwstats[[nw]], mat)
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
#' The network model parameters stored in \code{dat$nwparam} are updated for
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

  coef.form.m <- get_nwparam(dat, network = 1)$coef.form
  coef.form.m[1] <- coef.form.m[1] + adjust
  dat$nwparam[[1]]$coef.form <- coef.form.m

  coef.form.p <- get_nwparam(dat, network = 2)$coef.form
  coef.form.p[1] <- coef.form.p[1] + adjust
  dat$nwparam[[2]]$coef.form <- coef.form.p

  coef.form.i <- get_nwparam(dat, network = 3)$coef.form
  coef.form.i[1] <- coef.form.i[1] + adjust
  dat$nwparam[[3]]$coef.form <- coef.form.i

  return(dat)
}


#' @title Initialize Cumulative Partner List
#'
#' @description This function initializes the cumulative partners list, which
#'              stores the members, and start/end times, of all partnerships
#'              in a simulation.
#' @inheritParams aging_msm
#'
#' @return An updated \code{dat} object with an initialized partner list.
#' @export
#'
init_plist <- function(dat) {
  plist <- cbind(dat$el[[1]], ptype = 1)
  if (length(dat$el) > 1) {
    for (i in 2:length(dat$el)) {
      plist <- rbind(plist, cbind(dat$el[[i]], ptype = i))
    }
  }

  plist <- cbind(plist, start = 1, stop = NA)
  plist[,"stop"][which(plist[,3] == 3)] <- 1
  colnames(plist)[c(1:2)] <- c("p1", "p2")
  # Save to temporary list location
  dat$temp$plist <- plist

  return(dat)
}

#' @title Update Cumulative Partner List
#'
#' @description This function updates the cumulative partners list, which
#'              stores the members, and start/end times, of all partnerships
#'              in a simulation.
#'
#' @inheritParams aging_msm
#' @param ptype Partnership type to update, as an integer corresponding to
#'              the edgelist position on \code{dat$el}.
#'
#' @return
#' This function returns the updated \code{plist} object with new partnerships
#' and updates to existing partnerships.
#'
#' @export
#'
update_plist <- function(dat, at, ptype) {

  # pull existing partner type specific list
  plist1 <- dat$temp$plist[dat$temp$plist[, "ptype"] == ptype, ]
  active <- get_attr(dat, "active")
  uid <- get_unique_ids(dat)

  # Main/Casl: update stop times, add new edges
  if (ptype %in% 1:2) {
    news <- attr(dat$el[[ptype]], "changes")
    news_uid <- cbind(matrix(uid[news[, 1:2]], ncol = 2), news[, 3])
    news_uid_stop <- news_uid[news_uid[, 3] == 0, , drop = FALSE]
    pid_plist1 <- plist1[, 1]*1e7 + plist1[, 2]
    pid_stop <- news_uid_stop[, 1]*1e7 + news_uid_stop[, 2]
    matches_stop <- match(pid_stop, pid_plist1)
    plist1[matches_stop, "stop"] <- at

    news_uid_start <- news_uid[news_uid[, 3] == 1, , drop = FALSE]
    if (nrow(news_uid_start) > 0) {
      plist1 <- rbind(
        plist1,
        cbind(news_uid_start[, 1:2, drop = FALSE], ptype, at, NA)
      )
    }
  }

  # Inst: add new edges, start=stop times
  if (ptype == 3) {
    news <- dat$el[[ptype]]
    news_uid <- matrix(uid[news[, 1:2]], ncol = 2)
    if (nrow(news_uid) > 0) {
      plist1 <- rbind(
        plist1,
        cbind(news_uid[, 1:2, drop = FALSE], ptype, at, at)
      )
    }
  }

  return(plist1)
}

