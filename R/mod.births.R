
#' @title Births Module
#'
#' @description Module function for births or entries into the sexually active
#'              population.
#'
#' @inheritParams aging.mard
#'
#' @details
#' New population members are added based on expected numbers of entries among
#' black and white MSM, stochastically determined with draws from Poisson
#' distributions. For each new entry, a set of attributes is added for that node,
#' and the nodes are added onto the network objects. Only attributes that are
#' a part of the network model formulae are updated as vertex attributes on the
#' network objects.
#'
#' @return
#' This function updates the \code{attr} list with new attributes for each new
#' population member, and the \code{nw} objects with new vertices.
#'
#' @keywords module
#' @export
#'
births.mard <- function(dat, at){

  # Variables ---------------------------------------------------------------
  b.rate.B <- dat$param$b.rate.B
  b.rate.W <- dat$param$b.rate.W

  active <- dat$attr$active
  race <- dat$attr$race
  numB <- sum(race == "B")
  numW <- sum(race == "W")

  currNwSize <- network.size(dat$nw$m)


  # Process -----------------------------------------------------------------
  nBirths.B <- rpois(1, b.rate.B * numB)
  nBirths.W <- rpois(1, b.rate.W * numW)
  nBirths <- nBirths.B + nBirths.W


  # Update Attr -------------------------------------------------------------
  if (nBirths > 0) {
    dat <- setBirthAttr.mard(dat, at, nBirths.B, nBirths.W)
  }


  # Update Networks ---------------------------------------------------------
  newIds <- NULL
  if (nBirths > 0) {
    newIds <- (currNwSize + 1):(currNwSize + nBirths)

    stopifnot(unique(sapply(dat$attr, length)) == (currNwSize + nBirths))

    new.race <- dat$attr$race[newIds]
    new.srage <- dat$attr$sqrt.age[newIds]
    new.rc <- dat$attr$role.class[newIds]
    new.degm <- new.degp <- paste0(new.race, 0)

    for (i in 1:3) {
      dat$nw[[i]] <- add.vertices.active(x = dat$nw[[i]], nv = nBirths,
                                         onset = at, terminus = Inf)
    }

    dat$nw$m <- set.vertex.attribute(x = dat$nw$m,
                                     attrname = c("race", "sqrt.age",
                                                  "role.class", "deg.pers"),
                                     value = list(race = new.race,
                                                  sqrt.age = new.srage,
                                                  role.class = new.rc,
                                                  deg.pers = new.degp),
                                     v = newIds)

    dat$nw$p <- set.vertex.attribute(x = dat$nw$p,
                                     attrname = c("race", "sqrt.age",
                                                  "role.class", "deg.main"),
                                     value = list(race = new.race,
                                                  sqrt.age = new.srage,
                                                  role.class = new.rc,
                                                  deg.main = new.degm),
                                     v = newIds)

    dat$nw$i <- set.vertex.attribute(x = dat$nw$i,
                                     attrname = c("race", "sqrt.age",
                                                  "role.class", "deg.pers",
                                                  "deg.main"),
                                     value = list(race = new.race,
                                                  sqrt.age = new.srage,
                                                  role.class = new.rc,
                                                  deg.pers = new.degp,
                                                  deg.main = new.degm),
                                     v = newIds)
  }


  # Output ------------------------------------------------------------------
  age <- dat$attr$age
  race <- dat$attr$race
  dat$epi$nBirths[at] <- nBirths
  dat$epi$nBirths.B[at] <- nBirths.B
  dat$epi$nBirths.W[at] <- nBirths.W
  dat$epi$nBirths.y[at] <- sum(age[newIds] < 30)
  dat$epi$nBirths.o[at] <- sum(age[newIds] >= 30)
  dat$epi$nBirths.B.y[at] <- sum(race[newIds] == "B" & age[newIds] < 30)
  dat$epi$nBirths.B.o[at] <- sum(race[newIds] == "B" & age[newIds] >= 30)
  dat$epi$nBirths.W.y[at] <- sum(race[newIds] == "W" & age[newIds] < 30)
  dat$epi$nBirths.W.o[at] <- sum(race[newIds] == "W" & age[newIds] >= 30)

  return(dat)
}


#' @title Set Attributes for Incoming Nodes
#'
#' @description Updates the attr list for new entries into the population as
#'              set by the births module.
#'
#' @inheritParams births.mard
#' @param nBirths.B Number of new entries among black MSM.
#' @param nBirths.W Number of new entries among white MSM.
#'
#' @keywords submodule
#' @export
#'
setBirthAttr.mard <- function(dat, at, nBirths.B, nBirths.W) {

  nBirths <- nBirths.B + nBirths.W

  # Set all attributes NA by default
  dat$attr <- lapply(dat$attr, {
    function(x)
      c(x, rep(NA, nBirths))
  })
  newIds <- which(is.na(dat$attr$active))

  # Demographic
  dat$attr$active[newIds] <- rep(1, nBirths)
  dat$attr$uid[newIds] <- dat$temp$max.uid + (1:nBirths)
  dat$temp$max.uid <- dat$temp$max.uid + nBirths

  dat$attr$arrival.time[newIds] <- rep(at, nBirths)

  race <- sample(rep(c("B", "W"), c(nBirths.B, nBirths.W)))
  newB <- which(race == "B")
  newW <- which(race == "W")
  dat$attr$race[newIds] <- race

  dat$attr$age[newIds] <- rep(dat$param$birth.age, nBirths)
  dat$attr$sqrt.age[newIds] <- sqrt(dat$attr$age[newIds])

  # Disease status and related
  dat$attr$status[newIds] <- rep(0, nBirths)

  dat$attr$inst.ai.class[newIds] <- sample(1:dat$param$num.inst.ai.classes,
                                           nBirths, replace = TRUE)

  dat$attr$tt.traj[newIds[newB]] <- sample(c("NN", "YN", "YP", "YF"),
                                           nBirths.B, replace = TRUE,
                                           prob = dat$param$tt.traj.freq.B)
  dat$attr$tt.traj[newIds[newW]] <- sample(c("NN", "YN", "YP", "YF"),
                                           nBirths.W, replace = TRUE,
                                           prob = dat$param$tt.traj.freq.W)

  dat$attr$circ[newIds[newB]] <- rbinom(nBirths.B, 1, dat$param$circ.prev.B)
  dat$attr$circ[newIds[newW]] <- rbinom(nBirths.W, 1, dat$param$circ.prev.W)

  dat$attr$role.class[newIds[newB]] <- sample(c("I", "R", "V"),
                                              nBirths.B, replace = TRUE,
                                              prob = dat$param$role.freq.B)
  dat$attr$role.class[newIds[newW]] <- sample(c("I", "R", "V"),
                                              nBirths.W, replace = TRUE,
                                              prob = dat$param$role.freq.W)

  ins.quot <- rep(NA, nBirths)
  ins.quot[dat$attr$role.class[newIds] == "I"]  <- 1
  ins.quot[dat$attr$role.class[newIds] == "R"]  <- 0
  ins.quot[dat$attr$role.class[newIds] == "V"]  <-
                                  runif(sum(dat$attr$role.class[newIds] == "V"))
  dat$attr$ins.quot[newIds] <- ins.quot

  # CCR5
  ccr5.freq.B <- dat$param$ccr5.freq.B
  ccr5.freq.W <- dat$param$ccr5.freq.W
  dat$attr$ccr5[newIds[newB]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.B, replace = TRUE,
                                        prob = c(1 - sum(ccr5.freq.B),
                                                 ccr5.freq.B[2],
                                                 ccr5.freq.B[1]))
  dat$attr$ccr5[newIds[newW]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.W, replace = TRUE,
                                        prob = c(1 - sum(ccr5.freq.W),
                                                 ccr5.freq.W[2],
                                                 ccr5.freq.W[1]))

  return(dat)

}
