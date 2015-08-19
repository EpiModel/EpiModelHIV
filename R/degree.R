
#' @title Assign Degree Vertex Attribute on Network Objects
#'
#' @description Assigns the degree vertex attributes on network objects
#'              conditional on their values from the other networks.
#'
#' @param nw Object of class \code{network} that is the target for the vertex
#'        attribute.
#' @param deg.type Type of degree to assign to \code{nw}, with options of
#'        \code{"pers"} to assign casual degree onto main network and
#'        \code{"main"} to assign main degree to casual network.
#' @param nwstats Object of class \code{nwstats}.
#'
#' @details
#' This function assigns the degree of other networks as a vertex attribute on the
#' target network given a bivariate degree mixing matrix of main, casual, and
#' one-partnerships contained in the \code{nwstats} data.
#'
#' @seealso \code{\link{update_degree}}
#' @export
#'
assign_degree <- function(nw, deg.type, nwstats) {

  if (!("network" %in% class(nw))) {
    stop("nw must be of class network")
  }

  if (deg.type == "main") {
    attr.name <- "deg.main"
    dist.B <- rowSums(nwstats$deg.mp.B)
    dist.W <- rowSums(nwstats$deg.mp.W)
  }
  if (deg.type == "pers") {
    attr.name <- "deg.pers"
    dist.B <- colSums(nwstats$deg.mp.B)
    dist.W <- colSums(nwstats$deg.mp.W)
  }

  if (sum(dist.B) != 1 || sum(dist.W) != 1) {
    stop("One of the degree distributions do not sum to 1")
  }

  race <- get.vertex.attribute(nw, "race")
  vB <- which(race == "B")
  vW <- which(race == "W")
  nB <- length(vB)
  nW <- length(vW)

  num.degrees.B <- length(dist.B)
  num.degrees.W <- length(dist.W)

  deg.B <- apportion_lr(nB, 0:(num.degrees.B - 1), dist.B, shuffled = TRUE)
  deg.W <- apportion_lr(nW, 0:(num.degrees.W - 1), dist.W, shuffled = TRUE)

  if (nwstats$method == 2) {
    deg.B <- paste0("B", deg.B)
    deg.W <- paste0("W", deg.W)
  }

  nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.B, v = vB)
  nw <- set.vertex.attribute(nw, attrname = attr.name, value = deg.W, v = vW)

  return(nw)
}


#' @title Update Degree Vertex Attribute on Network Objects
#'
#' @description Updates the degree vertex attributes on network objects
#'              conditional on their values from the other networks.
#'
#' @param nw Object of class \code{network} that is the target for updated
#'        degree vertex attribute.
#' @param nw.source Object of class \code{network} that is the source of the
#'        updated degree vertex attribute.
#' @param deg.type Type of degree to assign to \code{nw}, with options of
#'        \code{"pers"} to assign casual degree onto main network and
#'        \code{"main"} to assign main degree to casual network.
#' @param at Current time step.
#'
#' @details
#' This function queries the degree of each node within a network using the
#' sociality term in a summary function call to the network, and then sets this
#' as a race-based vertex attribute on the output network.
#'
#' @seealso \code{\link{assign_degree}}
#' @export
#'
update_degree <- function(nw, nw.source, deg.type, at) {

  if (!("network" %in% class(nw))) {
    stop("Argument nw must be of class network")
  }
  if (!("network" %in% class(nw.source))) {
    stop("Argument nw.source must be of class network")
  }

  if (deg.type == "main") {
    attr.name <- "deg.main"
  }
  if (deg.type == "pers") {
    attr.name <- "deg.pers"
  }

  deg.dist <- as.numeric(summary(nw.source ~ sociality(base = 0), at = at))
  if (is.character(nw %v% attr.name)) {
    race <- get.vertex.attribute(nw, "race")
    deg.dist <- paste0(race, deg.dist)
  }

  nw <- set.vertex.attribute(nw, attr.name, deg.dist)
  return(nw)
}
