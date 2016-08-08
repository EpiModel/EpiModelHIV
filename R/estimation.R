

# MSM -----------------------------------------------------------------

#' @title Calculate Target Statistics for Network Model Estimation
#'
#' @description Calculates the target statistics for the formation and dissolution
#'              components of the network model to be estimated with \code{netest}.
#'
#' @param time.unit Time unit relative to 1 for daily.
#' @param method Method for calculating target statistics by race, with options of
#'        \code{2} for preserving race-specific statistics and \code{1} for
#'        averaging over the statistics and dropping the race-specific terms.
#' @param num.B Population size of black MSM.
#' @param num.W Population size of white MSM.
#' @param deg.mp.B Degree distribution matrix for main and casual partners for
#'        black MSM, as a 2 by 3 matrix.
#' @param deg.mp.W Degree distribution matrix for main and causal partners for
#'        white MSM, as a 2 by 3 matrix.
#' @param mdeg.inst.B Mean degree, or rate, of one-off partnerships per day
#'        for black MSM.
#' @param mdeg.inst.W Mean degree, or rate, of one-off partnerships per day
#'        for white MSM.
#' @param qnts.B Means of one-off rates split into quintiles for white MSM. Use
#'        \code{NA} to ignore these quantiles in the target statistics.
#' @param qnts.W Means of one-off rates split into quintiles for black MSM. Use
#'        \code{NA} to ignore these quantiles in the target statistics.
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
#' @param diss.main Dissolution model formula for main partnerships.
#' @param diss.pers Dissolution model formula for casual partnerships.
#' @param durs.main Vector of length 3 with the duration of BB, BW, and WW main
#'        partnerships in days.
#' @param durs.pers Vector of length 3 with the duration of BB, BW, and WW
#'        casual partnerships in days.
#' @param ages Integer vector of ages in years that defines range of possible
#'        initial ages in the population.
#' @param asmr.B Vector of length 40 defining the age-specific
#'        mortality rate for persons within that age slot, for black MSM.
#' @param asmr.W Vector of length 40 defining the age-specific
#'        mortality rate for persons within that age slot, for white MSM.
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
#' @keywords msm
#'
#' @seealso
#' Network statistics calculated here are entered into \code{\link{base_nw_msm}}
#' to construct the base network, and then into the parameters in
#' \code{\link{param_msm}}.
#'
#' @export
#'
calc_nwstats_msm <- function(time.unit = 7,
                             method = 2,
                             num.B,
                             num.W,
                             deg.mp.B,
                             deg.mp.W,
                             mdeg.inst.B,
                             mdeg.inst.W,
                             qnts.B,
                             qnts.W,
                             prop.hom.mpi.B,
                             prop.hom.mpi.W,
                             balance = "mean",
                             sqrt.adiff.BB,
                             sqrt.adiff.WW,
                             sqrt.adiff.BW,
                             diss.main,
                             diss.pers,
                             durs.main,
                             durs.pers,
                             ages,
                             asmr.B,
                             asmr.W,
                             role.B.prob,
                             role.W.prob) {

  if (sum(deg.mp.B) != 1) {
    stop("deg.mp.B must sum to 1.")
  }
  if (sum(deg.mp.W) != 1) {
    stop("deg.mp.W must sum to 1.")
  }
  if (!(method %in% 1:2)) {
    stop("method must either be 1 for one-race models or 2 for two-race models", call. = FALSE)
  }

  num <- num.B + num.W

  # deg.pers nodal attribute
  if (method == 2) {
    deg.pers.B <- apportion_lr(num.B, c("B0", "B1", "B2"), colSums(deg.mp.B))
    deg.pers.W <- apportion_lr(num.W, c("W0", "W1", "W2"), colSums(deg.mp.W))
  }
  if (method == 1) {
    deg.pers <- apportion_lr(num, 0:2, colSums(deg.mp.W))
  }

  # deg main nodal attribute
  if (method == 2) {
    deg.main.B <- apportion_lr(num.B, c("B0", "B1"), rowSums(deg.mp.B))
    deg.main.W <- apportion_lr(num.W, c("W0", "W1"), rowSums(deg.mp.W))
  }
  if (method == 1) {
    deg.main <- apportion_lr(num, 0:1, rowSums(deg.mp.W))
  }


  # Main partnerships -------------------------------------------------------

  # Persons in partnerships by casual degree
  if (method == 2) {
    totdeg.m.by.dp <- c(num.B * deg.mp.B[2, ], num.W * deg.mp.W[2, ])
  }
  if (method == 1) {
    totdeg.m.by.dp <- c(num * deg.mp.B[2, ])
  }

  # Persons in partnerships by race
  if (method == 2) {
    totdeg.m.by.race <- c(sum(totdeg.m.by.dp[1:3]), sum(totdeg.m.by.dp[4:6]))
  }

  # Number of partnerships
  edges.m <- (sum(totdeg.m.by.dp)) / 2

  # Mixing
  if (method == 2) {
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
  }

  # Sqrt absdiff term for age
  if (method == 2) {
    sqrt.adiff.m <- edges.nodemix.m * c(sqrt.adiff.BB[1], sqrt.adiff.BW[1], sqrt.adiff.WW[1])
  }
  if (method == 1) {
    sqrt.adiff.m <- edges.m * mean(c(sqrt.adiff.BB[1], sqrt.adiff.BW[1], sqrt.adiff.WW[1]))
  }

  # Compile target stats
  if (method == 2) {
    stats.m <- c(edges.m, edges.nodemix.m[2:3], totdeg.m.by.dp[c(2:3, 5:6)], sqrt.adiff.m)
  }
  if (method == 1) {
    stats.m <- c(edges.m, totdeg.m.by.dp[2:3], sqrt.adiff.m)
  }

  # Dissolution model
  exp.mort <- (mean(asmr.B[ages]) + mean(asmr.W[ages])) / 2

  coef.diss.m <- dissolution_coefs(dissolution = diss.main,
                                   duration = durs.main / time.unit,
                                   d.rate = exp.mort)



  # Casual partnerships -----------------------------------------------------

  # Persons in partnerships by main degree
  if (method == 2) {
    totdeg.p.by.dm <- c(num.B * deg.mp.B[, 2] + num.B * deg.mp.B[, 3] * 2,
                        num.W * deg.mp.W[, 2] + num.W * deg.mp.W[, 3] * 2)
  }
  if (method == 1) {
    totdeg.p.by.dm <- c(num * deg.mp.B[, 2] + num * deg.mp.B[, 3] * 2)
  }

  # Persons in partnerships by race
  if (method == 2) {
    totdeg.p.by.race <- c(sum(totdeg.p.by.dm[1:2]), sum(totdeg.p.by.dm[3:4]))
  }

  # Persons concurrent
  if (method == 2) {
    conc.p.by.race <- c(sum(deg.mp.B[, 3]) * num.B, sum(deg.mp.W[, 3]) * num.W)
  }
  if (method == 1) {
    conc.p <- sum(deg.mp.B[, 3] * num)
  }

  # Number of partnerships
  edges.p <- sum(totdeg.p.by.dm) / 2

  # Mixing
  if (method == 2) {
    # Number of mixed-race partnerships, with balancing to decide
    edges.p.B2W <- totdeg.p.by.race[1] * (1 - prop.hom.mpi.B[2])
    edges.p.W2B <- totdeg.p.by.race[2] * (1 - prop.hom.mpi.W[2])
    edges.het.p <- switch(balance,
                          black = edges.p.B2W, white = edges.p.W2B,
                          mean = (edges.p.B2W + edges.p.W2B) / 2)

    # Number of same-race partnerships
    edges.hom.p <- (totdeg.p.by.race - edges.het.p) / 2

    # Nodemix target stat: number of BB, BW, WW partnerships
    edges.nodemix.p <- c(edges.hom.p[1], edges.het.p, edges.hom.p[2])
  }

  # Sqrt absdiff term for age
  if (method == 2) {
    sqrt.adiff.p <- edges.nodemix.p * c(sqrt.adiff.BB[2], sqrt.adiff.BW[2], sqrt.adiff.WW[2])
  }
  if (method == 1) {
    sqrt.adiff.p <- edges.p * mean(c(sqrt.adiff.BB[2], sqrt.adiff.BW[2], sqrt.adiff.WW[2]))
  }

  # Compile target statistics
  if (method == 2) {
    stats.p <- c(edges.p, edges.nodemix.p[2:3], totdeg.p.by.dm[c(2, 4)],
                 conc.p.by.race, sqrt.adiff.p)
  }
  if (method == 1) {
    stats.p <- c(edges.p, totdeg.p.by.dm[2], conc.p, sqrt.adiff.p)
  }

  # Dissolution model
  coef.diss.p <- dissolution_coefs(dissolution = diss.pers,
                                   duration = durs.pers / time.unit,
                                   d.rate = exp.mort)



  # Instant partnerships ----------------------------------------------------

  # Number of instant partnerships per time step, by main and casl degree
  if (method == 2) {
    num.inst.B <- num.B * deg.mp.B * mdeg.inst.B * time.unit
    num.inst.W <- num.W * deg.mp.W * mdeg.inst.W * time.unit
  }
  if (method == 1) {
    num.inst <- num * deg.mp.W * mdeg.inst.W * time.unit
  }

  # Risk quantiles
  if (!is.na(qnts.B[1]) & !is.na(qnts.W[1])) {
    if (method == 2) {
      num.riskg.B <- (0.2*num.B) * qnts.B * time.unit
      num.riskg.W <- (0.2*num.W) * qnts.W * time.unit
    }
    if (method == 1) {
      num.riskg <- 0.2 * num * qnts.B * time.unit
    }
  }

  # Number of instant partnerships per time step, by race
  if (method == 2) {
    totdeg.i <- c(sum(num.inst.B), sum(num.inst.W))
  }
  if (method == 1) {
    totdeg.i <- sum(num.inst)
  }

  # Number of partnerships
  edges.i <- sum(totdeg.i) / 2

  # Mixing
  if (method == 2) {
    # Number of mixed-race partnerships, with balancing to decide
    edges.i.B2W <- totdeg.i[1] * (1 - prop.hom.mpi.B[3])
    edges.i.W2B <- totdeg.i[2] * (1 - prop.hom.mpi.W[3])
    edges.het.i <- switch(balance,
                          black = edges.i.B2W, white = edges.i.W2B,
                          mean = (edges.i.B2W + edges.i.W2B) / 2)

    # Number of same-race partnerships
    edges.hom.i <- edges.i - edges.het.i

    # Nodemix target stat: number of BB, BW, WW partnerships
    edges.nodemix.i <- c((totdeg.i[1] - edges.het.i) / 2,
                         edges.het.i,
                         (totdeg.i[1] - edges.het.i) / 2)
  }

  if (method == 2) {
    sqrt.adiff.i <- edges.nodemix.i * c(sqrt.adiff.BB[3], sqrt.adiff.BW[3], sqrt.adiff.WW[3])
  }
  if (method == 1) {
    sqrt.adiff.i <- edges.i * mean(c(sqrt.adiff.BB[3], sqrt.adiff.BW[3], sqrt.adiff.WW[3]))
  }

  if (!is.na(qnts.B[1]) & !is.na(qnts.W[1])) {
    if (method == 2) {
      stats.i <- c(edges.i, num.inst.B[-1], num.inst.W,
                   num.riskg.B[-3], num.riskg.W[-3],
                   edges.hom.i, sqrt.adiff.i)
    }
    if (method == 1) {
      stats.i <- c(edges.i, num.inst[-1], num.riskg[-3], sqrt.adiff.i)
    }

  } else {
    if (method == 2) {
      stats.i <- c(edges.i, num.inst.B[-1], num.inst.W, edges.hom.i, sqrt.adiff.i)
    }
    if (method == 1) {
      stats.i <- c(edges.i, num.inst[-1], sqrt.adiff.i)
    }
  }


  # Compile results ---------------------------------------------------------
  out <- list()
  out$method <- method
  if (method == 2) {
    out$deg.pers <- c(deg.pers.B, deg.pers.W)
    out$deg.main <- c(deg.main.B, deg.main.W)
  }
  if (method == 1) {
    out$deg.pers <- deg.pers
    out$deg.main <- deg.main
  }

  out$stats.m <- stats.m
  out$stats.p <- stats.p
  out$stats.i <- stats.i

  out$coef.diss.m <- coef.diss.m
  out$coef.diss.p <- coef.diss.p

  out$ages <- ages
  out$asmr.B <- asmr.B
  out$asmr.W <- asmr.W

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
#'        \code{\link{calc_nwstats_msm}}.
#'
#' @details
#' This function takes the output of \code{\link{calc_nwstats_msm}} and constructs
#' an empty network with the necessary attributes for race, square root of age,
#' and sexual role class. This base network is used for all three network
#' estimations.
#'
#' @seealso
#' The final vertex attributes on the network for cross-network degree are
#' calculated and set on the network with \code{\link{assign_degree}}.
#'
#' @keywords msm
#' @export
#'
base_nw_msm <- function(nwstats) {

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

  role.B <- sample(apportion_lr(num.B, c("I", "R", "V"), nwstats$role.B.prob))
  role.W <- sample(apportion_lr(num.W, c("I", "R", "V"), nwstats$role.W.prob))
  role <- rep(NA, n)
  role[race == "B"] <- role.B
  role[race == "W"] <- role.W

  riskg.B <- sample(apportion_lr(num.B, 1:5, rep(0.2, 5)))
  riskg.W <- sample(apportion_lr(num.W, 1:5, rep(0.2, 5)))
  riskg <- rep(NA, n)
  riskg[race == "B"] <- riskg.B
  riskg[race == "W"] <- riskg.W

  attr.names <- c("race", "riskg", "sqrt.age", "role.class")
  attr.values <- list(race, riskg, sqrt.age, role)
  nw <- network::set.vertex.attribute(nw, attr.names, attr.values)

  return(nw)
}


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
#' @keywords msm
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

  if (!isTRUE(all.equal(sum(colSums(nwstats$deg.mp.B)), 1, tolerance = 5e-6))) {
    stop("B degree distributions do not sum to 1")
  }

  if (!isTRUE(all.equal(sum(colSums(nwstats$deg.mp.W)), 1, tolerance = 5e-6))) {
    stop("W degree distributions do not sum to 1")
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


# Het -----------------------------------------------------------------

#' @title Calculate Network Statistics
#'
#' @description This function calculates the target statistics for the formation
#'              and dissolution models estimated in \code{netest}.
#'
#' @param n Population size.
#' @param meandeg Mean degree.
#' @param prop.male Percent of the population that is male.
#' @param start.prev Starting HIV prevalence in the population.
#' @param part.dur Mean duration of partnerships.
#' @param time.unit Time unit used, relative to days.
#'
#' @keywords het
#' @export
#'
make_nw_het <- function(n = 10000,
                        meandeg = 0.8,
                        prop.male = 0.5,
                        start.prev = 0.05,
                        part.dur = 1000,
                        time.unit = 7) {

  nMale <- round(n * prop.male)

  male <- rep(0, n)
  male[sample(1:n, nMale)] <- 1

  # Set vertex attributes
  nw <- network.initialize(n = n, directed = FALSE)
  nw <- set.vertex.attribute(nw, attrname = "male", value = male)

  # Formation Model
  formation <- ~edges + offset(nodematch("male"))

  # Target stats
  edges.ts <- meandeg * (n/2)

  stats <- edges.ts

  # Dissolution model
  dissolution <- ~offset(edges)
  dur <- part.dur/time.unit
  d.rate <- time.unit * (((1 - start.prev) * 1/(55 - 18)/365) + (start.prev * 1/12/365))
  coef.diss <- dissolution_coefs(dissolution, duration = dur, d.rate = d.rate)

  out <- list()
  out$nw <- nw
  out$time.unit <- time.unit
  out$formation <- formation
  out$stats <- stats
  out$coef.diss <- coef.diss

  return(out)
}
