
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
#' The module also deletes inactive nodes on all three networks depending on the
#' value of the \code{delete.nodes} control setting. For the main and casual
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

  resim.int <- dat$control$resim.int

  # Main network ------------------------------------------------------------

  if (at %% resim.int == 0) {
    nwparam.m <- EpiModel::get_nwparam(dat, network = 1)

    suppressWarnings(
      dat$nw$m <- simulate(dat$nw$m,
                           formation = nwparam.m$formation,
                           dissolution = nwparam.m$coef.diss$dissolution,
                           coef.form = nwparam.m$coef.form,
                           coef.diss = nwparam.m$coef.diss$coef.adj,
                           constraints = nwparam.m$constraints,
                           time.start = at,
                           time.slices = 1 * resim.int,
                           time.offset = 0,
                           monitor = "all",
                           output = "networkDynamic"))

    stats <- tail(as.data.frame(attributes(dat$nw$m)$stats), 1 * resim.int)
    stats <- calc_meandegs(dat, at, stats, "main")  ## Fix this for resim.int > 1
    if (at == 2) {
      dat$stats$nwstats$m <- stats
    } else {
      dat$stats$nwstats$m <- rbind(dat$stats$nwstats$m, stats)
    }
  }

  dat$nw$p <- update_degree(dat$nw$p, dat$nw$m, deg.type = "main", at = at)
  dat$nw$i <- update_degree(dat$nw$i, dat$nw$m, deg.type = "main", at = at)

  asdf <- as.data.frame(dat$nw$m)
  dat$temp$new.edges <- NULL
  new.edges.m <- as.matrix(asdf[asdf$onset == at, c("tail", "head")], ncol = 2)[, 2:1]
  dat$temp$new.edges <- matrix(dat$attr$uid[new.edges.m], ncol = 2)


  # Casual network ----------------------------------------------------------
  if (at %% resim.int == 0) {
    nwparam.p <- EpiModel::get_nwparam(dat, network = 2)

    suppressWarnings(
      dat$nw$p <- simulate(dat$nw$p,
                           formation = nwparam.p$formation,
                           dissolution = nwparam.p$coef.diss$dissolution,
                           coef.form = nwparam.p$coef.form,
                           coef.diss = nwparam.p$coef.diss$coef.adj,
                           constraints = nwparam.p$constraints,
                           time.start = at,
                           time.slices = 1 * resim.int,
                           time.offset = 0,
                           monitor = "all",
                           output = "networkDynamic"))

    stats <- tail(as.data.frame(attributes(dat$nw$p)$stats), 1 * resim.int)
    stats <- calc_meandegs(dat, at, stats, "pers") ## Fix this for resim.int > 1
    if (at == 2) {
      dat$stats$nwstats$p <- stats
    } else {
      dat$stats$nwstats$p <- rbind(dat$stats$nwstats$p, stats)
    }
  }

  dat$nw$m <- update_degree(dat$nw$m, dat$nw$p, deg.type = "pers", at = at)
  dat$nw$i <- update_degree(dat$nw$i, dat$nw$p, deg.type = "pers", at = at)

  asdf <- as.data.frame(dat$nw$p)
  new.edges.p <- as.matrix(asdf[asdf$onset == at, c("tail", "head")], ncol = 2)[, 2:1]
  dat$temp$new.edges <- rbind(dat$temp$new.edges,
                              matrix(dat$attr$uid[new.edges.p], ncol = 2))


  # Delete nodes ------------------------------------------------------------
  inactive <- which(dat$attr$active == 0)
  if (dat$control$delete.nodes == TRUE && ((at %% resim.int) == (resim.int - 1))
      && length(inactive) > 0) {
    for (i in 1:2) {
      dat$nw[[i]] <- network.extract(dat$nw[[i]], at = at)
    }
    dat$nw$i <- delete.vertices(dat$nw$i, vid = inactive)
    dat$attr <- EpiModel::deleteAttr(dat$attr, inactive)

    for (i in 1:length(dat$riskh)) {
      dat$riskh[[i]] <- dat$riskh[[i]][-inactive, ]
    }
  }

  # Stop checks for consistency in nw / attr sizes
  stopifnot(length(unique(sapply(dat$nw, function(x) x$gal$n))) == 1)
  stopifnot(length(unique(sapply(dat$attr, length))) == 1)
  stopifnot(length(dat$attr[[1]]) == dat$nw$m$gal$n)


  # Instant network ---------------------------------------------------------
  nwparam.i <- EpiModel::get_nwparam(dat, network = 3)

  inst.formula <- update.formula(nwparam.i$formation, dat$nw$i ~ .)
  environment(inst.formula) <- environment()
  dat$nw$i <- simulate(inst.formula, coef = nwparam.i$coef.form)

  stats <- data.frame(t(summary(inst.formula)))
  stats <- calc_meandegs(dat, at, stats, "inst")
  if (at == 2) {
    dat$stats$nwstats$i <- stats
  } else {
    dat$stats$nwstats$i <- rbind(dat$stats$nwstats$i, stats)
  }

  return(dat)
}


calc_meandegs <- function(dat, at, stats, type) {

  if (type == "main") {
    nf.m <- summary(dat$nw$m ~ nodefactor("race", base = 0), at = at)
    md.MB <- unname(round(nf.m[1] / sum(dat$attr$active == 1 & dat$attr$race == "B"), 3))
    md.MW <- unname(round(nf.m[2] / sum(dat$attr$active == 1 & dat$attr$race == "W"), 3))
    stats$md.MB <- md.MB
    stats$md.MW <- md.MW
  }

  if (type == "pers") {
    nf.p <- summary(dat$nw$p ~ nodefactor("race", base = 0), at = at)
    md.PB <- unname(round(nf.p[1] / sum(dat$attr$active == 1 & dat$attr$race == "B"), 3))
    md.PW <- unname(round(nf.p[2] / sum(dat$attr$active == 1 & dat$attr$race == "W"), 3))
    stats$md.PB <- md.PB
    stats$md.PW <- md.PW
  }

  if (type == "inst") {
    nf.i <- summary(dat$nw$i ~ nodefactor("race", base = 0))
    md.IB <- unname(round(nf.i[1] / sum(dat$attr$active == 1 & dat$attr$race == "B"), 3))
    md.IW <- unname(round(nf.i[2] / sum(dat$attr$active == 1 & dat$attr$race == "W"), 3))
    stats$md.IB <- md.IB
    stats$md.IW <- md.IW
  }

  return(stats)
}
