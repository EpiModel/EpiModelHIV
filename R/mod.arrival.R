
#' @title Arrivals Module
#'
#' @description Module function for arrivals into the sexually active
#'              population.
#'
#' @inheritParams aging_msm
#'
#' @details
#' New population members are added based on expected numbers of entries,
#' stochastically determined with draws from Poisson distributions. For each new
#' entry, a set of attributes is added for that node, and the nodes are added onto
#' the network objects. Only attributes that are a part of the network model
#' formulae are updated as vertex attributes on the network objects.
#'
#' @return
#' This function updates the \code{attr} list with new attributes for each new
#' population member, and the \code{nw} objects with new vertices.
#'
#' @keywords module msm
#' @export
#'
arrival_msm <- function(dat, at){

  ## Variables

  # Parameters
  a.rate <- dat$param$a.rate

  ## Process
  num <- dat$epi$num[1]

  nNew <- rpois(1, a.rate * num)

  ## Update Attr
  if (nNew > 0) {
    dat <- setNewAttr_msm(dat, at, nNew)
  }

  # Update Networks
  if (nNew > 0) {
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::add_vertices(dat$el[[i]], nNew)
    }
  }

  ## Output
  dat$epi$nNew[at] <- nNew

  return(dat)
}


setNewAttr_msm <- function(dat, at, nNew) {

  # Set all attributes NA by default
  dat$attr <- lapply(dat$attr, {
    function(x)
      c(x, rep(NA, nNew))
  })
  newIds <- which(is.na(dat$attr$active))

  # Demographic
  dat$attr$active[newIds] <- rep(1, nNew)
  dat$attr$uid[newIds] <- dat$temp$max.uid + (1:nNew)
  dat$temp$max.uid <- dat$temp$max.uid + nNew

  dat$attr$arrival.time[newIds] <- rep(at, nNew)

  race.dist <- prop.table(table(dat$param$netstats$attr$race))

  race <- sample(sort(unique(dat$attr$race)), nNew, TRUE, race.dist)
  dat$attr$race[newIds] <- race

  dat$attr$age[newIds] <- rep(dat$param$arrival.age, nNew)
  age.breaks <- dat$param$netstats$demog$age.breaks
  attr_age.grp <- cut(dat$attr$age[newIds], age.breaks, labels = FALSE)
  dat$attr$age.grp[newIds] <- attr_age.grp

  # Disease status and related
  dat$attr$status[newIds] <- rep(0, nNew)
  dat$attr$diag.status[newIds] <- rep(0, nNew)
  dat$attr$rGC[newIds] <- dat$attr$GC.timesInf[newIds] <- 0
  dat$attr$uGC[newIds] <- dat$attr$GC.timesInf[newIds] <- 0
  dat$attr$rCT[newIds] <- dat$attr$CT.timesInf[newIds] <- 0
  dat$attr$uCT[newIds] <- dat$attr$CT.timesInf[newIds] <- 0

  dat$attr$count.trans[newIds] <- 0

  rates <- dat$param$hiv.test.late.prob[race]
  dat$attr$late.tester[newIds] <- rbinom(length(rates), 1, rates)

  races <- sort(unique(dat$attr$race[newIds]))
  tt.traj <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(dat$attr$race[newIds] == i)
    tt.traj[ids.race] <- sample(1:3, length(ids.race), TRUE,
                                c(dat$param$tt.part.supp[i],
                                  dat$param$tt.full.supp[i],
                                  dat$param$tt.dur.supp[i]))

  }
  dat$attr$tt.traj[newIds] <- tt.traj

  # Circumcision
  circ <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(dat$attr$race[newIds] == i)
    circ[ids.race] <- rbinom(length(ids.race), 1, dat$param$circ.prob[i])
  }
  dat$attr$circ[newIds] <- circ

  # Role
  ns <- dat$param$netstats$attr
  role.class <- rep(NA, nNew)
  for (i in races) {
    ids.race <- which(dat$attr$race[newIds] == i)
    rc.probs <- prop.table(table(ns$role.class[ns$race == i]))
    role.class[ids.race] <- sample(0:2, length(ids.race), TRUE, rc.probs)
  }
  dat$attr$role.class[newIds] <- role.class

  ins.quot <- rep(NA, nNew)
  ins.quot[dat$attr$role.class[newIds] == 0]  <- 1
  ins.quot[dat$attr$role.class[newIds] == 1]  <- 0
  ins.quot[dat$attr$role.class[newIds] == 2]  <-
                                  runif(sum(dat$attr$role.class[newIds] == 2))
  dat$attr$ins.quot[newIds] <- ins.quot

  # Degree
  dat$attr$deg.main[newIds] <- 0
  dat$attr$deg.casl[newIds] <- 0
  dat$attr$deg.tot[newIds] <- 0

  # One-off risk group
  dat$attr$risk.grp[newIds] <- sample(1:5, nNew, TRUE)

  # PrEP
  dat$attr$prepStat[newIds] <- 0
  dat$attr$prepStat.la[newIds] <- 0

  # HIV screening
  dat$attr$num.neg.tests[newIds] <- 0

  # Update clinical history
  if (dat$control$save.clin.hist == TRUE & length(newIds) > 0) {
    m <- dat$temp$clin.hist
    for (i in 1:length(m)) {
      new.m <- array(dim = c(length(newIds), dat$control$nsteps))
      m[[i]] <- rbind(m[[i]], new.m)
    }
    dat$temp$clin.hist <- m
  }

  ## Check attributes written as expected
  # cbind(sapply(dat$attr, function(x) is.na(tail(x, 1))))

  return(dat)
}



#' @export
#' @rdname arrival_msm
births_het <- function(dat, at) {

  # Variables
  b.rate.method <- dat$param$b.rate.method
  b.rate <- dat$param$b.rate
  active <- dat$attr$active


  # Process
  nBirths <- 0
  if (b.rate.method == "stgrowth") {
    exptPopSize <- dat$epi$num[1] * (1 + b.rate*at)
    numNeeded <- exptPopSize - sum(active == 1)
    if (numNeeded > 0) {
      nBirths <- rpois(1, numNeeded)
    }
  }
  if (b.rate.method == "totpop") {
    nElig <- dat$epi$num[at - 1]
    if (nElig > 0) {
      nBirths <- rpois(1, nElig * b.rate)
    }
  }
  if (b.rate.method == "fpop") {
    nElig <- dat$epi$num.feml[at - 1]
    if (nElig > 0) {
      nBirths <- rpois(1, nElig * b.rate)
    }
  }


  # Update Population Structure
  if (nBirths > 0) {
    dat <- setBirthAttr_het(dat, at, nBirths)
    dat$el[[1]] <- tergmLite::add_vertices(dat$el[[1]], nBirths)
  }

  if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
    stop("mismatch between el and attr length in births mod")
  }

  # Output
  dat$epi$b.flow[at] <- nBirths

  return(dat)
}


setBirthAttr_het <- function(dat, at, nBirths) {

  # Set attributes for new births to NA
  dat$attr <- lapply(dat$attr, function(x) c(x, rep(NA, nBirths)))
  newIds <- which(is.na(dat$attr$active))


  # Network Status
  dat$attr$active[newIds] <- rep(1, nBirths)
  dat$attr$entTime[newIds] <- rep(at, nBirths)


  # Demography
  prop.male <- ifelse(is.null(dat$param$b.propmale),
                      dat$epi$propMale[1],
                      dat$param$b.propmale)
  dat$attr$male[newIds] <- rbinom(nBirths, 1, prop.male)

  dat$attr$age[newIds] <- rep(18, nBirths)

  # Circumcision
  entTime <- dat$attr$entTime

  idsNewMale <- which(dat$attr$male == 1 & entTime == at)

  if (length(idsNewMale) > 0) {
    age <- dat$attr$age[idsNewMale]
    newCirc <- rbinom(length(idsNewMale), 1, dat$param$circ.prob.birth)
    isCirc <- which(newCirc == 1)

    newCircTime <- rep(NA, length(idsNewMale))
    newCircTime[isCirc] <- round(-age[isCirc] * (365 / dat$param$time.unit))

    dat$attr$circStat[idsNewMale] <- newCirc
    dat$attr$circTime[idsNewMale] <- newCircTime
  }


  # Epi/Clinical
  dat$attr$status[newIds] <- rep(0, nBirths)

  if (length(unique(sapply(dat$attr, length))) != 1) {
    sapply(dat$attr, length)
    stop("Attribute dimensions not unique")
  }

  return(dat)
}
