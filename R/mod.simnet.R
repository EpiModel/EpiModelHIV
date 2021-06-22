
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

  ## Edges correction
  dat <- edges_correct_msm(dat, at)

  ## Main network
  nwparam.m <- EpiModel::get_nwparam(dat, network = 1)

  nwparam.m$coef.form[1] <-
    nwparam.m$coef.form[1] + log(dat$param$netresim.form.rr[1])
  nwparam.m$coef.diss$coef.adj[1] <-
    nwparam.m$coef.diss$coef.adj[1] - log(dat$param$netresim.disl.rr[1])

  dat$attr$deg.casl <- get_degree(dat$el[[2]])
  dat <- tergmLite::updateModelTermInputs(dat, network = 1)

  dat$el[[1]] <- tergmLite::simulate_network(p = dat$p[[1]],
                                             el = dat$el[[1]],
                                             coef.form = nwparam.m$coef.form,
                                             coef.diss = nwparam.m$coef.diss$coef.adj,
                                             save.changes = TRUE)

  plist1 <- update_plist(dat, at, ptype = 1)

  ## Casual network
  nwparam.p <- EpiModel::get_nwparam(dat, network = 2)

  nwparam.p$coef.form[1] <-
    nwparam.p$coef.form[1] + log(dat$param$netresim.form.rr[2])
  nwparam.p$coef.diss$coef.adj[1] <-
    nwparam.p$coef.diss$coef.adj[1] - log(dat$param$netresim.disl.rr[2])

  dat$attr$deg.main <- get_degree(dat$el[[1]])
  dat <- tergmLite::updateModelTermInputs(dat, network = 2)

  dat$el[[2]] <- tergmLite::simulate_network(p = dat$p[[2]],
                                             el = dat$el[[2]],
                                             coef.form = nwparam.p$coef.form,
                                             coef.diss = nwparam.p$coef.diss$coef.adj,
                                             save.changes = TRUE)

  plist2 <- update_plist(dat, at, ptype = 2)

  dat$temp$plist <- rbind(plist1, plist2)
  if (dat$control$truncate.plist == TRUE) {
    to.keep <- which(is.na(dat$temp$plist[, "stop"]))
    dat$temp$plist <- dat$temp$plist[to.keep, ]
  }

  ## One-off network
  nwparam.i <- EpiModel::get_nwparam(dat, network = 3)
  nwparam.i$coef.form[1] <-
    nwparam.i$coef.form[1] + log(dat$param$netresim.form.rr[3])

  dat$attr$deg.tot <- pmin(dat$attr$deg.main + get_degree(dat$el[[2]]), 3)
  dat <- tergmLite::updateModelTermInputs(dat, network = 3)

  dat$el[[3]] <- tergmLite::simulate_ergm(p = dat$p[[3]],
                                          el = dat$el[[3]],
                                          coef = nwparam.i$coef.form)

  if (dat$control$save.nwstats == TRUE) {
    dat <- calc_nwstats(dat, at)
  }

  return(dat)
}

# updates the partnership list
update_plist <- function(dat, at, ptype) {
  # pull existing partner type specific list
  plist1 <- dat$temp$plist[dat$temp$plist[, "ptype"] == ptype, ]

  # look up dissolutions, update stop time
  uid <- dat$attr$uid
  news <- attr(dat$el[[ptype]], "changes")
  news_uid <- cbind(matrix(uid[news[, 1:2]], ncol = 2), news[, 3])
  news_uid_stop <- news_uid[news_uid[, 3] == 0, , drop = FALSE]
  pid_plist1 <- plist1[, 1]*1e7 + plist1[, 2]
  pid_stop <- news_uid_stop[, 1]*1e7 + news_uid_stop[, 2]
  matches_stop <- match(pid_stop, pid_plist1)
  plist1[matches_stop, "stop"] <- at

  # look up new formations, row bind them
  news_uid_start <- news_uid[news_uid[, 3] == 1, , drop = FALSE]
  if (nrow(news_uid_start) > 0) {
    plist1 <- rbind(plist1, cbind(news_uid_start[, 1:2, drop = FALSE], ptype, at, NA))
  }

  return(plist1)
}


calc_nwstats <- function(dat, at) {

  for (nw in 1:3) {
    n <- attr(dat$el[[nw]], "n")
    edges <- nrow(dat$el[[nw]])
    meandeg <- round(edges * (2/n), 3)
    concurrent <- round(mean(get_degree(dat$el[[nw]]) > 1), 3)
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

  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)
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




# HET -----------------------------------------------------------------


#' @export
#' @rdname simnet_msm
simnet_het <- function(dat, at) {

  # Update edges coefficients
  dat <- edges_correct_het(dat, at)

  # Update internal ergm data
  dat <- tergmLite::updateModelTermInputs(dat, network = 1)

  # Pull network parameters
  nwparam <- get_nwparam(dat, network = 1)

  # Simulate edgelist
  dat$el[[1]] <- tergmLite::simulate_network(p = dat$p[[1]],
                                             el = dat$el[[1]],
                                             coef.form = nwparam$coef.form,
                                             coef.diss = nwparam$coef.diss$coef.adj)

  return(dat)
}


#' @export
#' @rdname edges_correct_msm
edges_correct_het <- function(dat, at) {

  # Popsize
  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)

  # New Coefs
  coef.form <- get_nwparam(dat)$coef.form
  coef.form[1] <- coef.form[1] + log(old.num) - log(new.num)
  dat$nwparam[[1]]$coef.form <- coef.form

  return(dat)
}

