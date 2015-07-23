
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

  # Parameters
  b.B.rate <- dat$param$b.B.rate
  b.W.rate <- dat$param$b.W.rate
  b.method <- dat$param$b.method


  # Process -----------------------------------------------------------------
  if (b.method == "fixed") {
    numB <- dat$epi$num.B[1]
    numW <- dat$epi$num.W[1]
  }
  if (b.method == "varying") {
    numB <- dat$epi$num.B[at - 1]
    numW <- dat$epi$num.W[at - 1]
  }

  nBirths.B <- rpois(1, b.B.rate * numB)
  nBirths.W <- rpois(1, b.W.rate * numW)
  nBirths <- nBirths.B + nBirths.W


  # Update Attr -------------------------------------------------------------
  if (nBirths > 0) {
    dat <- setBirthAttr.mard(dat, at, nBirths.B, nBirths.W)
  }


  # Update Networks ---------------------------------------------------------
  currNwSize <- network.size(dat$nw$m)
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
                                           prob = dat$param$tt.traj.B.prob)
  dat$attr$tt.traj[newIds[newW]] <- sample(c("NN", "YN", "YP", "YF"),
                                           nBirths.W, replace = TRUE,
                                           prob = dat$param$tt.traj.W.prob)

  dat$attr$circ[newIds[newB]] <- rbinom(nBirths.B, 1, dat$param$circ.B.prob)
  dat$attr$circ[newIds[newW]] <- rbinom(nBirths.W, 1, dat$param$circ.W.prob)

  dat$attr$role.class[newIds[newB]] <- sample(c("I", "R", "V"),
                                              nBirths.B, replace = TRUE,
                                              prob = dat$param$role.B.prob)
  dat$attr$role.class[newIds[newW]] <- sample(c("I", "R", "V"),
                                              nBirths.W, replace = TRUE,
                                              prob = dat$param$role.W.prob)

  ins.quot <- rep(NA, nBirths)
  ins.quot[dat$attr$role.class[newIds] == "I"]  <- 1
  ins.quot[dat$attr$role.class[newIds] == "R"]  <- 0
  ins.quot[dat$attr$role.class[newIds] == "V"]  <-
                                  runif(sum(dat$attr$role.class[newIds] == "V"))
  dat$attr$ins.quot[newIds] <- ins.quot

  # CCR5
  ccr5.B.prob <- dat$param$ccr5.B.prob
  ccr5.W.prob <- dat$param$ccr5.W.prob
  dat$attr$ccr5[newIds[newB]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.B, replace = TRUE,
                                        prob = c(1 - sum(ccr5.B.prob),
                                                 ccr5.B.prob[2],
                                                 ccr5.B.prob[1]))
  dat$attr$ccr5[newIds[newW]] <- sample(c("WW", "DW", "DD"),
                                        nBirths.W, replace = TRUE,
                                        prob = c(1 - sum(ccr5.W.prob),
                                                 ccr5.W.prob[2],
                                                 ccr5.W.prob[1]))

  return(dat)

}
