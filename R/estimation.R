
#' @title Calculate Target Statistics for Network Model Estimation
#'
#' @description Calculates the target statistics for the formation and dissolution
#'              components of the network model to be estimated with \code{netest}.
#'
#' @param time.unit Time unit relative to 1 for daily.
#' @param num.B Population size of black MSM.
#' @param num.W Population size of white MSM.
#' @param deg.mp.B Degree distribution matrix for main and casual partners for
#'        black MSM, as a 2 by 3 matrix.
#' @param deg.mp.W Degree distribution matrix for main and causal partners for
#'        white MSM, as a 2 by 3 matrix.
#' @param mdeg.inst.B Mean degree, or rate, of one-off partnerships per time step
#'        for black MSM.
#' @param mdeg.inst.W Mean degree, or rate, of one-off partnerships per time step
#'        for white MSM.
#' @param prop.hom.mpi.B A vector of length 3 for the proportion of main, casual,
#'        and one-off partnerships in same race for black MSM.
#' @param prop.hom.mpi.W A vector of length 3 for the proportion of main, casual,
#'        and one-off partnerships in same race for white MSM.
#' @param balance Method for balancing of edges by race for number of mixed-race
#'        partnerships, with options of \code{"black"} to apply black MSM counts,
#'        \code{"white"} to apply white MSM counts, and \code{"mean"} to take
#'        the average of the two expectations.
#' @param sqrt.adiff.BB Vector of length 3 with the mean absolute differences
#'        in the square root of ages in main, casual, and one-off black-black
#'        partnerships.
#' @param sqrt.adiff.WW Vector of length 3 with the mean absolute differences
#'        in the square root of ages in main, casual, and one-off white-white
#'        partnerships.
#' @param sqrt.adiff.BW Vector of length 3 with the mean absolute differences
#'        in the square root of ages in main, casual, and one-off black-white
#'        partnerships.
#' @param age.method Method for calculating the square root of age differences,
#'        with \code{"heterogeneous"} using the individual \code{sqrt.adiff.XX}
#'        values and \code{"homogeneous"} using a weighted average.
#' @param diss.main Dissolution model formula for main partnerships.
#' @param diss.pers Dissolution model formula for casual partnerships.
#' @param durs.main Vector of length 3 with the duration of BB, BW, and WW main
#'        partnerships in days.
#' @param durs.pers Vector of length 3 with the duration of BB, BW, and WW
#'        casual partnerships in days.
#' @param ages Integer vector of ages in years that defines range of possible
#'        initial ages in the population.
#' @param asm.B.rate Vector of length 40 defining the age-specific mortality rate for
#'        persons within that age slot, for black MSM.
#' @param asm.W.rate Vector of length 40 defining the age-specific mortality rate for
#'        persons within that age slot, for white MSM.
#' @param role.B.prob Vector of length 3 for the probability of sexual role as
#'        insertive, receptive, and versatile, for black MSM.
#' @param role.W.prob Vector of length 3 for the probability of sexual role as
#'        insertive, receptive, and versatile, for white MSM.
#'
#' @details
#' This function performs basic calculations to determine the components of the
#' formationa and dissolution models for the network model estimation to be
#' conducted with \code{\link{netest}}. The inputs inputs for this function are
#' calculated externally to the package in a setup scenario file.
#'
#' @seealso
#' Network statistics calculated here are entered into \code{\link{base_nw.mard}}
#' to construct the base network, and then into the parameters in
#' \code{\link{param.mard}}.
#'
#' @export
#'
calc_nwstats.mard <- function(time.unit = 7,
                              num.B,
                              num.W,
                              deg.mp.B,
                              deg.mp.W,
                              mdeg.inst.B,
                              mdeg.inst.W,
                              prop.hom.mpi.B,
                              prop.hom.mpi.W,
                              balance = "mean",
                              sqrt.adiff.BB,
                              sqrt.adiff.WW,
                              sqrt.adiff.BW,
                              age.method = "heterogeneous",
                              diss.main,
                              diss.pers,
                              durs.main,
                              durs.pers,
                              ages,
                              asm.B.rate,
                              asm.W.rate,
                              role.B.prob,
                              role.W.prob) {

  if (sum(deg.mp.B) != 1) {
    stop("deg.mp.B must sum to 1.")
  }
  if (sum(deg.mp.W) != 1) {
    stop("deg.mp.W must sum to 1.")
  }

  # deg.pers nodal attribute
  deg.pers.B <- apportion.lr(num.B, c("B0", "B1", "B2"), colSums(deg.mp.B))
  deg.pers.W <- apportion.lr(num.W, c("W0", "W1", "W2"), colSums(deg.mp.W))

  # deg main nodal attribute
  deg.main.B <- apportion.lr(num.B, c("B0", "B1"), rowSums(deg.mp.B))
  deg.main.W <- apportion.lr(num.W, c("W0", "W1"), rowSums(deg.mp.W))

  # Main partnerships -------------------------------------------------------

  # Persons in partnerships by casual degree by race
  totdeg.m.by.dp <- c(num.B * deg.mp.B[2, ], num.W * deg.mp.W[2, ])

  # Persons in partnerships by race
  totdeg.m.by.race <- c(sum(totdeg.m.by.dp[1:3]), sum(totdeg.m.by.dp[4:6]))

  # Number of partnerships
  edges.m <- (sum(totdeg.m.by.dp)) / 2

  # Number of mixed-race partnerships, with balancing to decide
  edges.m.B2W <- totdeg.m.by.race[1] * (1 - prop.hom.mpi.B[1])
  edges.m.W2B <- totdeg.m.by.race[2] * (1 - prop.hom.mpi.W[1])
  edges.het.m <- switch(balance,
                        black = edges.m.B2W,
                        white = edges.m.W2B,
                        mean = (edges.m.B2W + edges.m.W2B) / 2)

  # Number of same-race partnerships
  edges.hom.m <- (totdeg.m.by.race - edges.het.m) / 2

  # Nodemix target stat: numer of BB, BW, WW partnerships
  edges.nodemix.m <- c(edges.hom.m[1], edges.het.m, edges.hom.m[2])


  # Sqrt absdiff term for age
  if (!(age.method %in% c("heterogeneous", "homogeneous"))) {
    stop("age.method must be \"heterogeneous\" or \"homogeneous\" ", call. = FALSE)
  }
  if (age.method == "heterogeneous") {
    sqrt.adiff.m <- edges.nodemix.m * c(sqrt.adiff.BB[1],
                                        sqrt.adiff.BW[1],
                                        sqrt.adiff.WW[1])
  }
  if (age.method == "homogeneous") {
    weighted.avg <- sum(edges.nodemix.m * c(sqrt.adiff.BB[1],
                                            sqrt.adiff.BW[1],
                                            sqrt.adiff.WW[1])) /
      sum(edges.nodemix.m)
    sqrt.adiff.m <- edges.nodemix.m * weighted.avg
  }

  # Compile target stats
  stats.m <- c(edges.m,
               edges.nodemix.m[2:3],
               totdeg.m.by.dp[c(2:3, 5:6)],
               sqrt.adiff.m)


  # Dissolution model
  exp.mort <- (mean(asm.B.rate[ages]) + mean(asm.W.rate[ages])) / 2

  coef.diss.m <- dissolution_coefs(dissolution = diss.main,
                                   duration = durs.main,
                                   d.rate = exp.mort)



  # Casual partnerships -----------------------------------------------------

  # Persons in partnerships by main degree by race
  totdeg.p.by.dm <- c(num.B * deg.mp.B[, 2] + num.B * deg.mp.B[, 3] * 2,
                      num.W * deg.mp.W[, 2] + num.W * deg.mp.W[, 3] * 2)

  # Persons in partnerships by race
  totdeg.p.by.race <- c(sum(totdeg.p.by.dm[1:2]), sum(totdeg.p.by.dm[3:4]))

  # Persons concurrent by race
  conc.p.by.race <- c(sum(deg.mp.B[, 3]) * num.B, sum(deg.mp.W[, 3]) * num.W)

  # Number of partnerships
  edges.p <- (sum(totdeg.p.by.dm)) / 2

  # Number of mixed-race partnerships, with balancing to decide
  edges.p.B2W <- totdeg.p.by.race[1] * (1 - prop.hom.mpi.B[2])
  edges.p.W2B <- totdeg.p.by.race[2] * (1 - prop.hom.mpi.W[2])
  edges.het.p <- switch(balance,
                        black = edges.p.B2W,
                        white = edges.p.W2B,
                        mean = (edges.p.B2W + edges.p.W2B) / 2)

  # Number of same-race partnerships
  edges.hom.p <- (totdeg.p.by.race - edges.het.p) / 2

  # Nodemix target stat: number of BB, BW, WW partnerships
  edges.nodemix.p <- c(edges.hom.p[1], edges.het.p, edges.hom.p[2])

  # Sqrt absdiff term for age
  if (age.method == "heterogeneous") {
    sqrt.adiff.p <- edges.nodemix.p * c(sqrt.adiff.BB[2],
                                        sqrt.adiff.BW[2],
                                        sqrt.adiff.WW[2])
  }
  if (age.method == "homogeneous") {
    weighted.avg <- sum(edges.nodemix.p * c(sqrt.adiff.BB[2],
                                            sqrt.adiff.BW[2],
                                            sqrt.adiff.WW[2])) /
      sum(edges.nodemix.p)
    sqrt.adiff.p <- edges.nodemix.p * weighted.avg
  }

  # Compile target statistics
  stats.p <- c(edges.p,
               edges.nodemix.p[2:3],
               totdeg.p.by.dm[c(2, 4)],
               conc.p.by.race,
               sqrt.adiff.p)

  # Dissolution model
  coef.diss.p <- dissolution_coefs(dissolution = diss.pers,
                                   duration = durs.pers,
                                   d.rate = exp.mort)



  # Instant partnerships ----------------------------------------------------

  # Number of instant partnerships per time step, by main and casl degree, for race
  num.inst.B <- num.B * deg.mp.B * mdeg.inst.B
  num.inst.W <- num.W * deg.mp.W * mdeg.inst.W

  # Number of instant partnerships per time step, by race
  totdeg.i.by.race <- c(sum(num.inst.B), sum(num.inst.W))

  # Number of partnerships
  edges.i <- sum(totdeg.i.by.race) / 2

  # Number of mixed-race partnerships, with balancing to decide
  edges.i.B2W <- totdeg.i.by.race[1] * (1 - prop.hom.mpi.B[3])
  edges.i.W2B <- totdeg.i.by.race[2] * (1 - prop.hom.mpi.W[3])
  edges.het.i <- switch(balance,
                        black = edges.i.B2W,
                        white = edges.i.W2B,
                        mean = (edges.i.B2W + edges.i.W2B) / 2)

  # Number of same-race partnerships
  edges.hom.i <- edges.i - edges.het.i

  # Nodemix target stat: number of BB, BW, WW partnerships
  edges.nodemix.i <- c((totdeg.i.by.race[1] - edges.het.i) / 2,
                       edges.het.i,
                       (totdeg.i.by.race[1] - edges.het.i) / 2)

  if (age.method == "heterogeneous") {
    sqrt.adiff.i <- edges.nodemix.i * c(sqrt.adiff.BB[3],
                                        sqrt.adiff.BW[3],
                                        sqrt.adiff.WW[3])
  }
  if (age.method == "homogeneous") {
    weighted.avg <- sum(edges.nodemix.i * c(sqrt.adiff.BB[3],
                                            sqrt.adiff.BW[3],
                                            sqrt.adiff.WW[3])) /
      sum(edges.nodemix.i)
    sqrt.adiff.i <- edges.nodemix.i * weighted.avg
  }

  stats.i <- c(edges.i,
               num.inst.B[-1],
               num.inst.W,
               edges.hom.i,
               sqrt.adiff.i)


  # Compile results ---------------------------------------------------------
  out <- list()
  out$deg.pers <- c(deg.pers.B, deg.pers.W)
  out$deg.main <- c(deg.main.B, deg.main.W)
  out$stats.m <- stats.m
  out$stats.p <- stats.p
  out$stats.i <- stats.i
  out$coef.diss.m <- coef.diss.m
  out$coef.diss.p <- coef.diss.p

  out$ages <- ages
  out$asm.B.rate <- asm.B.rate
  out$asm.W.rate <- asm.W.rate

  out$time.unit <- time.unit
  out$num.B <- num.B
  out$num.W <- num.W

  out$deg.mp.B <- deg.mp.B
  out$deg.mp.W <- deg.mp.W

  out$role.B.prob <- role.B.prob
  out$role.W.prob <- role.W.prob

  class(out) <- "nwstats"
  return(out)
}


#' @title Construct Base Network for Model Estimation and Simulation
#'
#' @description Initializes the base network for model estimation within
#'              \code{netest}.
#'
#' @param nwstats An object of class \code{nwstats}, as output from
#'        \code{\link{calc_nwstats.mard}}.
#'
#' @details
#' This function takes the output of \code{\link{calc_nwstats.mard}} and constructs
#' an empty network with the necessary attributes for race, square root of age,
#' and sexual role class. This base network is used for all three network
#' estimations.
#'
#' @seealso
#' The final vertex attributes on the network for cross-network degree are
#' calculated and set on the network with \code{\link{assign_degree}}.
#'
#' @export
#'
base_nw.mard <- function(nwstats) {

  num.B <- nwstats$num.B
  num.W <- nwstats$num.W

  # Initialize network
  n <- num.B + num.W
  nw <- network::network.initialize(n, directed = FALSE)

  # Calculate attributes
  race <- c(rep("B", num.B), rep("W", num.W))
  race <- sample(race)

  ager <- nwstats$ages
  ages <- seq(min(ager), max(ager) + 1, 1 / (365 / nwstats$time.unit))
  age <- sample(ages, n, TRUE)
  sqrt.age <- sqrt(age)

  role.B <- sample(apportion.lr(num.B, c("I", "R", "V"), nwstats$role.B.prob))
  role.W <- sample(apportion.lr(num.W, c("I", "R", "V"), nwstats$role.W.prob))
  role <- rep(NA, n)
  role[race == "B"] <- role.B
  role[race == "W"] <- role.W

  attr.names <- c("race", "sqrt.age", "role.class")
  attr.values <- list(race, sqrt.age, role)
  nw <- network::set.vertex.attribute(nw, attr.names, attr.values)

  return(nw)
}



