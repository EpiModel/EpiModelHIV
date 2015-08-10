
#' @title Initialization Module
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param.mard}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.mard}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.mard}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated \code{dat} object with the initialized values
#' for demographics and disease-related variables.
#'
#' @export
#' @keywords module
#'
initialize.mard <- function(x, param, init, control, s) {

  # Master data list --------------------------------------------------------
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control

  dat$attr <- list()
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()
  dat$epi <- list()

  # Network simulation ------------------------------------------------------
  dat$nw <- list()
  for (i in 1:3) {
    dat$nw[[i]] <- simulate(x[[i]]$fit)
    dat$nw[[i]] <- remove_bad_roles(dat$nw[[i]])
    if (i %in% 1:2) {
      dat$nw[[i]] <- activate.vertices(dat$nw[[i]], onset = 1, terminus = Inf)
      if (control$delete.nodes == TRUE) {
        dat$nw[[i]] <- network.extract(dat$nw[[i]], at = 1)
      }
    }
  }
  names(dat$nw) <- c("m", "p", "i")

  # Network parameters
  dat$nwparam <- list()
  for (i in 1:3) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }


  # Nodal attributes --------------------------------------------------------

  # Degree terms
  deg.pers <- x[[1]]$fit$network %v% "deg.pers"
  deg.main <- x[[2]]$fit$network %v% "deg.main"
  dat$nw$m <- set.vertex.attribute(dat$nw$m, "deg.pers", deg.pers)
  dat$nw$p <- set.vertex.attribute(dat$nw$p, "deg.main", deg.main)
  dat$nw$i <- set.vertex.attribute(dat$nw$i, "deg.main", deg.main)
  dat$nw$i <- set.vertex.attribute(dat$nw$i, "deg.pers", deg.pers)


  # Race
  race <- dat$nw[[1]] %v% "race"
  dat$attr$race <- race
  num.B <- dat$init$num.B
  num.W <- dat$init$num.W
  num <- num.B + num.W
  ids.B <- which(dat$attr$race == "B")
  ids.W <- which(dat$attr$race == "W")

  dat$attr$active <- rep(1, num)
  dat$attr$uid <- 1:num
  dat$temp$max.uid <- num

  # Age
  dat$attr$sqrt.age <- dat$nw[[1]] %v% "sqrt.age"
  age <- dat$attr$sqrt.age^2
  dat$attr$age <- age

  # Risk group
  dat$attr$riskg <- dat$nw[[3]] %v% "riskg"

  # Arrival and departure
  dat$attr$arrival.time <- rep(1, num)
  dat$attr$depart.time <- rep(NA, num)

  # Circumcision
  circ <- rep(NA, num)
  circ[ids.B] <- sample(apportion.lr(num.B, 0:1, 1 - param$circ.B.prob))
  circ[ids.W] <- sample(apportion.lr(num.B, 0:1, 1 - param$circ.W.prob))
  dat$attr$circ <- circ

  # PrEP Attributes
  dat$attr$prepClass <- rep(NA, num)
  dat$attr$prepElig <- rep(NA, num)
  dat$attr$prepEligTime <- rep(NA, num)
  dat$attr$prepStat <- rep(NA, num)
  dat$attr$prepStartTime <- rep(NA, num)
  dat$attr$prepEver <- rep(NA, num)

  # One-off AI class
  inst.ai.class <- rep(NA, num)
  ncl <- param$num.inst.ai.classes
  inst.ai.class[ids.B] <- sample(apportion.lr(num.B, 1:ncl, rep(1 / ncl, ncl)))
  inst.ai.class[ids.W] <- sample(apportion.lr(num.W, 1:ncl, rep(1 / ncl, ncl)))
  dat$attr$inst.ai.class <- inst.ai.class

  # Role class
  role.class <- rep(NA, num)
  role.class[ids.B] <- sample(apportion.lr(num.B, c("I", "R", "V"),
                                           param$role.B.prob))
  role.class[ids.W] <- sample(apportion.lr(num.W, c("I", "R", "V"),
                                           param$role.W.prob))
  dat$attr$role.class <- role.class

  # Ins.quot
  ins.quot <- rep(NA, num)
  ins.quot[role.class == "I"]  <- 1
  ins.quot[role.class == "R"]  <- 0
  ins.quot[role.class == "V"]  <- runif(sum(role.class == "V"))
  dat$attr$ins.quot <- ins.quot

  # HIV-related attributes
  dat <- init_status.mard(dat)

  # CCR5
  dat <- init_ccr5(dat)


  # Network statistics ------------------------------------------------------

  dat$stats$nwstats <- list()


  # Prevalence Tracking -----------------------------------------------------

  dat$temp$dal <- list()
  if (dat$control$save.dal == TRUE) {
    dat$temp$dal[[1]] <- list()
  }
  dat$temp$deg.dists <- list()
  dat$temp$discl.list <- as.data.frame(matrix(NA, 0, 4))
  names(dat$temp$discl.list) <- c("pos", "neg", "discl.time", "discl.type")

  dat <- prevalence.mard(dat, at = 1)

  class(dat) <- "dat"
  return(dat)

}


remove_bad_roles <- function(nw) {

  el <- as.edgelist(nw)

  rc <- get.vertex.attribute(nw, "role.class")
  rc.el <- matrix(rc[el], ncol = 2)

  rc.el.bad <- which((rc.el[, 1] == "R" & rc.el[, 2] == "R") |
                     (rc.el[, 1] == "I" & rc.el[, 2] == "I"))

  if (length(rc.el.bad) > 0) {
    el.bad <- el[rc.el.bad, , drop = FALSE]

    eid <- rep(NA, nrow(el.bad))
    for (i in 1:nrow(el.bad)) {
      eid[i] <- get.edgeIDs(nw, v = el.bad[i, 1], alter = el.bad[i, 2])
    }
    nw <- delete.edges(nw, eid)
  }

  return(nw)
}


init_status.mard <- function(dat) {

  num.B <- dat$init$num.B
  num.W <- dat$init$num.W
  num <- num.B + num.W
  ids.B <- which(dat$attr$race == "B")
  ids.W <- which(dat$attr$race == "W")
  age <- dat$attr$age
  race <- dat$attr$race

  # Infection Status
  nInfB <- round(dat$init$prev.B * num.B)
  nInfW <- round(dat$init$prev.W * num.W)

  # Age-based infection probability
  probInfCrB <- age[ids.B] * dat$init$init.prev.age.slope.B
  probInfB <- probInfCrB + (nInfB - sum(probInfCrB)) / num.B

  probInfCrW <- age[ids.W] * dat$init$init.prev.age.slope.W
  probInfW <- probInfCrW + (nInfW - sum(probInfCrW)) / num.W

  if (any(probInfB <= 0) | any(probInfW <= 0)) {
    stop("Slope of initial prevalence by age must be sufficiently low to ",
         "avoid non-positive probabilities.", call. = FALSE)
  }

  # Infection status
  status <- rep(0, num)
  while (sum(status[ids.B]) != nInfB) {
    status[ids.B] <- rbinom(num.B, 1, probInfB)
  }
  while (sum(status[ids.W]) != nInfW) {
    status[ids.W] <- rbinom(num.W, 1, probInfW)
  }
  dat$attr$status <- status


  # Treatment trajectory
  tt.traj <- rep(NA, num)
  tt.traj[ids.B] <- sample(apportion.lr(num.B, c("NN", "YN", "YP", "YF"),
                                        dat$param$tt.traj.B.prob))
  tt.traj[ids.W] <- sample(apportion.lr(num.W, c("NN", "YN", "YP", "YF"),
                                        dat$param$tt.traj.W.prob))
  dat$attr$tt.traj <- tt.traj


  ## Infection-related attributes

  stage <- rep(NA, num)
  stage.time <- rep(NA, num)
  inf.time <- rep(NA, num)
  vl <- rep(NA, num)
  diag.status <- rep(NA, num)
  diag.time <- rep(NA, num)
  last.neg.test <- rep(NA, num)
  tx.status <- rep(NA, num)
  tx.init.time <- rep(NA, num)
  cum.time.on.tx <- rep(NA, num)
  cum.time.off.tx <- rep(NA, num)
  infector <- rep(NA, num)
  inf.role <- rep(NA, num)
  inf.type <- rep(NA, num)
  inf.diag <- rep(NA, num)
  inf.tx <- rep(NA, num)
  inf.stage <- rep(NA, num)

  time.sex.active <- pmax(1,
                          round((365 / dat$param$time.unit) * age - (365 / dat$param$time.unit) *
                                  min(dat$init$ages), 0))

  vlar.int <- dat$param$vl.acute.rise.int
  vlap <- dat$param$vl.acute.peak
  vlaf.int <- dat$param$vl.acute.fall.int
  vlsp <- dat$param$vl.set.point
  vldo.int <- dat$param$vl.aids.onset.int
  vl.aids.int <- dat$param$vl.aids.int
  vlf  <- dat$param$vl.fatal
  vlds <- (vlf - vlsp) / vl.aids.int
  vl.acute.int <- vlar.int + vlaf.int


  ### Non-treater type: tester and non-tester
  selected <- which(status == 1 & tt.traj %in% c("NN", "YN"))
  max.inf.time <- pmin(time.sex.active[selected], vldo.int + vl.aids.int)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  tx.status[selected] <- 0
  cum.time.on.tx[selected] <- 0
  cum.time.off.tx[selected] <- time.since.inf

  stage[selected[time.since.inf <= vlar.int]] <- "AR"
  stage[selected[time.since.inf > vlar.int & time.since.inf <= vl.acute.int]] <- "AF"
  stage[selected[time.since.inf > vl.acute.int & time.since.inf <= vldo.int]] <- "C"
  stage[selected[time.since.inf > vldo.int]] <- "D"

  stage.time[selected][stage[selected] == "AR"] <-
    time.since.inf[stage[selected] == "AR"]
  stage.time[selected][stage[selected] == "AF"] <-
    time.since.inf[stage[selected] == "AF"] - vlar.int
  stage.time[selected][stage[selected] == "C"] <-
    time.since.inf[stage[selected] == "C"] - vl.acute.int
  stage.time[selected][stage[selected] == "D"] <-
    time.since.inf[stage[selected] == "D"] - vldo.int

  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) * (time.since.inf <= vldo.int) * (vlsp) +
                  (time.since.inf > vldo.int) * (vlsp + (time.since.inf - vldo.int) * vlds)

  selected <- which(status == 1 & tt.traj == "NN")
  diag.status[selected] <- 0

  selected <- which(status == 1 & tt.traj == "YN")

  # Time to next test
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                             min = 0,
                             max = dat$param$mean.test.B.int * (race[selected] == "B") +
                                   dat$param$mean.test.W.int * (race[selected] == "W")))
  }
  if (dat$param$testing.pattern == "memoryless") {
    ttntest <- rgeom(length(selected),
                     1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                          dat$param$mean.test.W.int * (race[selected] == "W")))
  }

  twind.int <- dat$param$test.window.int
  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
                           -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]

  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1


  ### Full adherent type

  # Create set of expected values for (cum.time.off.tx, cum.time.on.tx)

  tx.init.time.B <- twind.int + dat$param$last.neg.test.B.int + 1 / dat$param$tx.init.B.prob
  tx.init.time.W <- twind.int + dat$param$last.neg.test.W.int + 1 / dat$param$tx.init.W.prob

  # Stage for Blacks
  prop.time.on.tx.B <- dat$param$tx.reinit.B.prob /
                       (dat$param$tx.halt.B.prob + dat$param$tx.reinit.B.prob)
  offon.B <- matrix(c(1:tx.init.time.B, rep(0, tx.init.time.B)),
                    nrow = tx.init.time.B)
  numsteps.B <- (dat$param$max.time.off.tx.full.int - tx.init.time.B) /
                (1 - prop.time.on.tx.B)
  offon.B <- rbind(offon.B,
                   cbind(tx.init.time.B + (1 - prop.time.on.tx.B) * 1:numsteps.B,
                         prop.time.on.tx.B * 1:numsteps.B))
  offon.B <- round(offon.B)
  exp.dur.chronic.B <- nrow(offon.B) - vl.acute.int
  exp.onset.aids.B <- nrow(offon.B)
  offon.last.B <- offon.B[nrow(offon.B), ]
  offon.B <- rbind(offon.B,
                   matrix(c(offon.last.B[1] + (1:vl.aids.int),
                            rep(offon.last.B[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.B <- nrow(offon.B)
  offon.B[, 2] <- (1:max.possible.inf.time.B) - offon.B[, 1]
  stage.B <- rep(c("AR", "AF", "C", "D"), c(vlar.int, vlaf.int, exp.dur.chronic.B, vl.aids.int))
  stage.time.B <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B, 1:vl.aids.int)

  # Stage for Whites
  prop.time.on.tx.W <- dat$param$tx.reinit.W.prob /
    (dat$param$tx.halt.W.prob + dat$param$tx.reinit.W.prob)
  offon.W <- matrix(c(1:tx.init.time.W, rep(0, tx.init.time.W)),
                    nrow = tx.init.time.W)
  numsteps.W <- (dat$param$max.time.off.tx.full.int - tx.init.time.W) /
    (1 - prop.time.on.tx.W)
  offon.W <- rbind(offon.W,
                   cbind(tx.init.time.W + (1 - prop.time.on.tx.W) * 1:numsteps.W,
                         prop.time.on.tx.W * 1:numsteps.W))
  offon.W <- round(offon.W)
  exp.dur.chronic.W <- nrow(offon.W) - vl.acute.int
  exp.onset.aids.W <- nrow(offon.W)
  offon.last.W <- offon.W[nrow(offon.W), ]
  offon.W <- rbind(offon.W,
                   matrix(c(offon.last.W[1] + (1:vl.aids.int),
                            rep(offon.last.W[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.W <- nrow(offon.W)
  offon.W[, 2] <- (1:max.possible.inf.time.W) - offon.W[, 1]
  stage.W <- rep(c("AR", "AF", "C", "D"), c(vlar.int, vlaf.int, exp.dur.chronic.W, vl.aids.int))
  stage.time.W <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W, 1:vl.aids.int)

  # Vl for Blacks
  selected <- which(status == 1 & tt.traj == "YF" & race == "B")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B[time.since.inf, 1]
  stage[selected] <- stage.B[time.since.inf]
  stage.time[selected] <- stage.time.B[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == "C" & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == "C" & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.B) * (vlsp) +
                  (time.since.inf > exp.onset.aids.B) *
                  (vlsp + (time.since.inf - exp.onset.aids.B) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp

  # VL for Whites
  selected <- which(status == 1 & tt.traj == "YF" & race == "W")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W[time.since.inf, 1]
  stage[selected] <- stage.W[time.since.inf]
  stage.time[selected] <- stage.time.W[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == "C" & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == "C" & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.W) * (vlsp) +
                  (time.since.inf > exp.onset.aids.W) *
                  (vlsp + (time.since.inf - exp.onset.aids.W) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp

  # Diagnosis
  selected <- which(status == 1 & tt.traj == "YF")
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                             min = 0,
                             max = dat$param$mean.test.B.int * (race[selected] == "B") +
                                   dat$param$mean.test.W.int * (race[selected] == "W")))
  }
  if (dat$param$testing.pattern == "memoryless") {
    ttntest <- rgeom(length(selected),
                     1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                          dat$param$mean.test.W.int * (race[selected] == "W")))
  }

  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
                           -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]
  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1
  diag.status[selected][cum.time.on.tx[selected] > 0] <- 1
  last.neg.test[selected][cum.time.on.tx[selected] > 0] <- NA


  ### Part adherent type

  # Create set of expected values for (cum.time.off.tx,cum.time.on.tx)

  prop.time.on.tx.B <- dat$param$tx.reinit.B.prob /
                       (dat$param$tx.halt.B.prob + dat$param$tx.reinit.B.prob)
  offon.B <- matrix(c(1:tx.init.time.B, rep(0, tx.init.time.B)),
                    nrow = tx.init.time.B)
  while (offon.B[nrow(offon.B), 1] / dat$param$max.time.off.tx.part.int +
         offon.B[nrow(offon.B), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.B <- rbind(offon.B,
                     offon.B[nrow(offon.B), ] + c(1 - prop.time.on.tx.B,
                                                      prop.time.on.tx.B))
  }
  offon.B <- round(offon.B)
  exp.dur.chronic.B <- nrow(offon.B) - vl.acute.int
  exp.onset.aids.B <- nrow(offon.B)
  offon.last.B <- offon.B[nrow(offon.B), ]
  offon.B <- rbind(offon.B,
                   matrix(c(offon.last.B[1] + (1:vl.aids.int),
                            rep(offon.last.B[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.B <- nrow(offon.B)
  offon.B[, 2] <- (1:max.possible.inf.time.B) - offon.B[, 1]
  stage.B <- rep(c("AR", "AF", "C", "D"), c(vlar.int, vlaf.int, exp.dur.chronic.B, vl.aids.int))
  stage.time.B <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B, 1:vl.aids.int)

  prop.time.on.tx.W <- dat$param$tx.reinit.W.prob /
                       (dat$param$tx.halt.W.prob + dat$param$tx.reinit.W.prob)
  offon.W <- matrix(c(1:tx.init.time.W, rep(0, tx.init.time.W)),
                    nrow = tx.init.time.W)

  while (offon.W[nrow(offon.W), 1] / dat$param$max.time.off.tx.part.int +
         offon.W[nrow(offon.W), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.W <- rbind(offon.W,
                     offon.W[nrow(offon.W), ] + c(1 - prop.time.on.tx.W,
                                                  prop.time.on.tx.W))
  }
  offon.W <- round(offon.W)
  exp.dur.chronic.W <- nrow(offon.W) - vl.acute.int
  exp.onset.aids.W <- nrow(offon.W)
  offon.last.W <- offon.W[nrow(offon.W), ]
  offon.W <- rbind(offon.W,
                   matrix(c(offon.last.W[1] + (1:vl.aids.int),
                            rep(offon.last.W[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.W <- nrow(offon.W)
  offon.W[, 2] <- (1:max.possible.inf.time.W) - offon.W[, 1]
  stage.W <- rep(c("AR", "AF", "C", "D"), c(vlar.int, vlaf.int, exp.dur.chronic.W, vl.aids.int))
  stage.time.W <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W, 1:vl.aids.int)

  # VL for Blacks
  selected <- which(status == 1 & tt.traj == "YP" & race == "B")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B[time.since.inf, 1]
  stage[selected] <- stage.B[time.since.inf]
  stage.time[selected] <- stage.time.B[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == "C" & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == "C" & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.B) * (vlsp) +
                  (time.since.inf > exp.onset.aids.B) *
                  (vlsp + (time.since.inf - exp.onset.aids.B) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp

  # VL for Whites
  selected <- which(status == 1 & tt.traj == "YP" & race == "W")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W[time.since.inf, 1]
  stage[selected] <- stage.W[time.since.inf]
  stage.time[selected] <- stage.time.W[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == "C" & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == "C" & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.W) * (vlsp) +
                  (time.since.inf > exp.onset.aids.W) *
                  (vlsp + (time.since.inf - exp.onset.aids.W) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp

  # Implement diagnosis for both
  selected <- which(status == 1 & tt.traj == "YP")
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                             min = 0,
                             max = dat$param$mean.test.B.int * (race[selected] == "B") +
                                   dat$param$mean.test.W.int * (race[selected] == "W")))
  }

  if (dat$param$testing.pattern == "memoryless") {
    ttntest <- rgeom(length(selected),
                     1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                          dat$param$mean.test.W.int * (race[selected] == "W")))
  }


  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
    -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]

  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1
  diag.status[selected][cum.time.on.tx[selected] > 0] <- 1
  last.neg.test[selected][cum.time.on.tx[selected] > 0] <- NA


  # Last neg test before present for negatives
  selected <- which(status == 0 & tt.traj %in% c("YN", "YP", "YF"))

  if (dat$param$testing.pattern == "interval") {
    tslt <- ceiling(runif(length(selected),
                          min = 0,
                          max = dat$param$mean.test.B.int * (race[selected] == "B") +
                                dat$param$mean.test.W.int * (race[selected] == "W")))
  }
  if (dat$param$testing.pattern == "memoryless") {
    tslt <- rgeom(length(selected),
                  1 / (dat$param$mean.test.B.int * (race[selected] == "B") +
                       dat$param$mean.test.W.int * (race[selected] == "W")))
  }
  last.neg.test[selected] <- -tslt


  ## Set all onto dat$attr
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time
  dat$attr$inf.time <- inf.time
  dat$attr$vl <- vl
  dat$attr$diag.status <- diag.status
  dat$attr$diag.time <- diag.time
  dat$attr$last.neg.test <- last.neg.test
  dat$attr$tx.status <- tx.status
  dat$attr$tx.init.time <- tx.init.time
  dat$attr$cum.time.on.tx <- cum.time.on.tx
  dat$attr$cum.time.off.tx <- cum.time.off.tx
  dat$attr$infector <- infector
  dat$attr$inf.role <- inf.role
  dat$attr$inf.type <- inf.type
  dat$attr$inf.diag <- inf.diag
  dat$attr$inf.tx <- inf.tx
  dat$attr$inf.stage <- inf.stage

  return(dat)

}


init_ccr5 <- function(dat) {

  num.B <- dat$init$num.B
  num.W <- dat$init$num.W
  num <- num.B + num.W
  ids.B <- which(dat$attr$race == "B")
  ids.W <- which(dat$attr$race == "W")
  race <- dat$attr$race
  status <- dat$attr$status

  nInfB <- sum(race == "B" & status == 1)
  nInfW <- sum(race == "W" & status == 1)

  ##  CCR5 genotype
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr
  ccr5 <- rep("WW", num)

  # homozygotes for deletion
  num.ccr5.DD.B <- dat$param$ccr5.B.prob[1] * num.B
  # heterozygotes
  num.ccr5.DW.B <- dat$param$ccr5.B.prob[2] * num.B
  # homozygotes for deletion
  num.ccr5.WW.B <- num.B - num.ccr5.DD.B - num.ccr5.DW.B
  # DD's can't be infected
  num.uninf.ccr5.DD.B <- round(num.ccr5.DD.B)
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.B <- round(num.ccr5.DW.B * nInfB * ccr5.heteroz.rr /
                             (num.ccr5.WW.B + num.ccr5.DW.B * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.B <- round(num.ccr5.DW.B - num.inf.ccr5.DW.B)
  inf.B <- which(status == 1 & race == "B")
  inf.ccr5.DW.B <- sample(inf.B, num.inf.ccr5.DW.B, replace = FALSE)
  ccr5[inf.ccr5.DW.B] <- "DW"
  uninf.B <- which(status == 0 & race == "B")
  uninf.ccr5.DWDD.B <- sample(uninf.B, num.uninf.ccr5.DW.B + num.uninf.ccr5.DD.B)
  uninf.ccr5.DW.B <- sample(uninf.ccr5.DWDD.B, num.uninf.ccr5.DW.B)
  uninf.ccr5.DD.B <- setdiff(uninf.ccr5.DWDD.B, uninf.ccr5.DW.B)
  ccr5[uninf.ccr5.DW.B] <- "DW"
  ccr5[uninf.ccr5.DD.B] <- "DD"

  num.ccr5.DD.W <- dat$param$ccr5.W.prob[1] * num.W
  num.ccr5.DW.W <- dat$param$ccr5.W.prob[2] * num.W
  num.ccr5.WW.W <- num.W - num.ccr5.DD.W - num.ccr5.DW.W
  num.uninf.ccr5.DD.W <- round(num.ccr5.DD.W)
  num.inf.ccr5.DW.W <- round(num.ccr5.DW.W * nInfW * ccr5.heteroz.rr /
                             (num.ccr5.WW.W + num.ccr5.DW.W * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.W <- round(num.ccr5.DW.W - num.inf.ccr5.DW.W)
  inf.W <- which(status == 1 & race == "W")
  inf.ccr5.DW.W <- sample(inf.W, num.inf.ccr5.DW.W)
  ccr5[inf.ccr5.DW.W] <- "DW"
  uninf.W <- which(status == 0 & race == "W")
  uninf.ccr5.DWDD.W <- sample(uninf.W, num.uninf.ccr5.DW.W + num.uninf.ccr5.DD.W)
  uninf.ccr5.DW.W <- sample(uninf.ccr5.DWDD.W, num.uninf.ccr5.DW.W)
  uninf.ccr5.DD.W <- setdiff(uninf.ccr5.DWDD.W, uninf.ccr5.DW.W)
  ccr5[uninf.ccr5.DW.W] <- "DW"
  ccr5[uninf.ccr5.DD.W] <- "DD"

  dat$attr$ccr5 <- ccr5

  return(dat)
}


#' @title Re-Initialization Module
#'
#' @description This function reinitializes an epidemic model to restart at a
#'              specified time step given an input \code{netsim} object.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}}.
#' @inheritParams initialize.mard
#'
#' @return
#' This function resets the data elements on the \code{dat} master data object
#' in the needed ways for the time loop to function.
#'
#' @export
#' @keywords module
#'
reinit.mard <- function(x, param, init, control, s) {

  if (is.null(x$network)) {
    stop("x must contain network to restart simulation", call. = FALSE)
  }
  if (is.null(x$attr)) {
    stop("x must contain attr to restart simulation", call. = FALSE)
  }
  if (is.null(x$temp)) {
    stop("x must contain temp to restart simulation", call. = FALSE)
  }

  if (!is.null(control$currsim)) {
    s <- control$currsim
  }

  dat <- list()
  dat$nw <- x$network[[s]]
  if (!is.null(x$last.ts)) {
    for (i in 1:2) {
      dat$nw[[i]] <- network.extract(dat$nw[[i]], at = x$last.ts)
    }
  }
  dat$param <- param
  dat$param$modes <- 1
  dat$control <- control
  dat$nwparam <- x$nwparam
  dat$epi <- sapply(x$epi, function(var) var[s])
  names(dat$epi) <- names(x$epi)
  dat$attr <- x$attr[[s]]
  dat$stats <- list()
  dat$stats$nwstats <- x$stats$nwstats[[s]]
  dat$temp <- x$temp[[s]]

  class(dat) <- "dat"

  return(dat)
}
