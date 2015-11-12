
#' @title Network Resimulation Module
#'
#' @description Module function for resimulating the main, casual, and one-off
#'              networks for one time step.
#'
#' @inheritParams aging.mard
#'
#' @details
#' The three networks must be simulated for one time step at a time because of
#' the underlying changes to the node set. Nodes have been added and removed due
#' to births and deaths, and vertex attributes of nodes may have changed. For
#' each of the three networks in the epidemic simulation, this module pulls the
#' relevant network model coefficients, resimulates the network with the
#' necessary functions, calculates and stores the network statistics in an
#' external data frame, updates the cross-network degree statistics using
#' \code{\link{update_degree}}, and calculates which edges on the newly
#' simulated network are new (needed for \code{\link{disclose.mard}} module).
#'
#' The main and causal networks are resimulated using the \code{simulate.network}
#' function in the \code{tergm} package, whereas the one-off network is
#' resimulated with the \code{simulate.formula} function in the \code{ergm}
#' package.
#'
#' For the main and casual
#' network, this involves extracting the \code{networkDynamic} object at the
#' current time point after the resimulation has occurre (dead nodes have already
#' been deactivated in \code{\link{deaths.mard}} so that they will not be allowed
#' to form partnerships). For the one-off network, this involves deleting the
#' vertices before the resimulation. This approach for deleting nodes is needed
#' given the book-keeping of calculating the new edges list.
#'
#' @return
#' The three network objects on \code{dat$nw} are returned with newly simulated
#' edges, and potentially reduced node sets. The network statistics on
#' \code{dat$nw} are updated, and the newly formed edge list in
#' \code{dat$temp$new.edges} is also updated.
#'
#' @keywords module
#' @export
#'
simnet.mard <- function(dat, at) {

  ## Main network
  nwparam.m <- EpiModel::get_nwparam(dat, network = 1)
  dat <- updatenwp.mard(dat, network = 1)

  dat$el[[1]] <- tergmLite::simulate_network(p = dat$p[[1]],
                                             el = dat$el[[1]],
                                             coef.form = nwparam.m$coef.form,
                                             coef.diss = nwparam.m$coef.diss$coef.adj,
                                             save.changes = TRUE)


  dat$temp$new.edges <- NULL
  if (at == 2) {
    new.edges.m <- matrix(dat$el[[1]], ncol = 2)
  } else {
    new.edges.m <- attributes(dat$el[[1]])$changes
    new.edges.m <- new.edges.m[new.edges.m[, "to"] == 1, 1:2, drop = FALSE]
  }
  dat$temp$new.edges <- matrix(dat$attr$uid[new.edges.m], ncol = 2)


  ## Casual network
  nwparam.p <- EpiModel::get_nwparam(dat, network = 2)
  dat <- updatenwp.mard(dat, network = 2)

  dat$el[[2]] <- tergmLite::simulate_network(p = dat$p[[2]],
                                             el = dat$el[[2]],
                                             coef.form = nwparam.p$coef.form,
                                             coef.diss = nwparam.p$coef.diss$coef.adj,
                                             save.changes = TRUE)

  if (at == 2) {
    new.edges.p <- matrix(dat$el[[2]], ncol = 2)
  } else {
    new.edges.p <- attributes(dat$el[[2]])$changes
    new.edges.p <- new.edges.p[new.edges.p[, "to"] == 1, 1:2, drop = FALSE]
  }
  dat$temp$new.edges <- rbind(dat$temp$new.edges,
                              matrix(dat$attr$uid[new.edges.p], ncol = 2))


  ## One-off network
  nwparam.i <- EpiModel::get_nwparam(dat, network = 3)
  dat <- updatenwp.mard(dat, network = 3)

  dat$el[[3]] <- tergmLite::simulate_ergm(p = dat$p[[3]],
                                          el = dat$el[[3]],
                                          coef = nwparam.i$coef.form)

  return(dat)
}
