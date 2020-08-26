
# MSM -----------------------------------------------------------------

#' @title Initialization Module
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param_msm}}.
#' @param init An \code{EpiModel} object of class \code{\link{init_msm}}.
#' @param control An \code{EpiModel} object of class \code{\link{control_msm}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated \code{dat} object with the initialized values
#' for demographics and disease-related variables.
#'
#' @export
#' @keywords module msm
#'
initialize_msm <- function(x, param, init, control, s) {

  ## Master Data List Setup ##
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control


  ## Network Setup ##
  # Initial network simulations
  dat$nw <- list()
  for (i in 1:3) {
    dat$nw[[i]] <- simulate(x[[i]]$fit, basis = x[[i]]$fit$newnetwork)
  }
  nw <- dat$nw

  # Pull Network parameters
  dat$nwparam <- list()
  for (i in 1:3) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }

  # Convert to tergmLite method
  dat <- init_tergmLite(dat)

  ## Nodal Attributes Setup ##
  dat$attr <- param$netstats$attr

  num <- network.size(nw[[1]])
  dat$attr$active <- rep(1, num)
  dat$attr$arrival.time <- rep(1, num)
  dat$attr$uid <- 1:num

  # Circumcision
  rates <- param$circ.prob[dat$attr$race]
  dat$attr$circ <- rbinom(length(rates), 1, rates)

  # Insertivity Quotient
  ins.quot <- rep(NA, num)
  role.class <- dat$attr$role.class
  ins.quot[role.class == 0]  <- 1
  ins.quot[role.class == 1]  <- 0
  ins.quot[role.class == 2]  <- runif(sum(role.class == 2))
  dat$attr$ins.quot <- ins.quot

  # HIV-related attributes
  dat <- init_status_msm(dat)

  # STI Status
  dat <- init_sti_msm(dat)

  # PrEP-related attributes
  dat$attr$prepClass <- rep(NA, num)
  dat$attr$prepElig <- rep(NA, num)
  dat$attr$prepStat <- rep(0, num)
  dat$attr$prepStartTime <- rep(NA, num)
  dat$attr$prepLastRisk <- rep(NA, num)
  dat$attr$prepLastStiScreen <- rep(NA, num)

  # PrEP LA attributes
  dat$attr$prepStat.la <- rep(0, num)
  dat$attr$prepClass.la <- rep(NA, num)
  dat$attr$prepElig.la <- rep(NA, num)
  dat$attr$prepTimeLastInj <- rep(NA, num)
  dat$attr$prepLA.dlevel <- rep(NA, num)
  dat$attr$prepLA.dlevel.int <- rep(NA, num)





  ## Other Setup ##
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()
  dat$epi <- list()

  # Prevalence Tracking
  dat$temp$max.uid <- num
  dat <- prevalence_msm(dat, at = 1)

  # Setup Partner List
  plist <- cbind(dat$el[[1]], ptype = 1)
  plist <- rbind(plist, cbind(dat$el[[2]], ptype = 2))
  plist <- cbind(plist, start = 1, stop = NA)
  colnames(plist)[1:2] <- c("p1", "p2")
  dat$temp$plist <- plist

  # Clinical history
  if (dat$control$save.clin.hist == TRUE) {
    dat <- save_clin_hist(dat, at = 1)
  }

  # Network statistics
  if (dat$control$save.nwstats == TRUE) {
    dat <- calc_nwstats(dat, at = 1)
  }

  # dat$param$netstats <- NULL
  class(dat) <- "dat"
  return(dat)
}


#' @title Initialize the HIV status of persons in the network
#'
#' @description Sets the initial individual-level disease status of persons
#'              in the network, as well as disease-related attributes for
#'              infected persons.
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#' @keywords initiation utility msm
#'
init_status_msm <- function(dat) {

  num <- sum(dat$attr$active == 1)

  # Sub in diag.status from model for status
  status <- dat$attr$diag.status

  # Late (AIDS-stage) tester type
  rates <- dat$param$hiv.test.late.prob[dat$attr$race]
  dat$attr$late.tester <- rbinom(length(rates), 1, rates)

  # Treatment trajectory
  tt.traj <- rep(NA, num)
  races <- sort(unique(dat$attr$race))
  for (i in races) {
    ids.race <- which(dat$attr$race == i)
    tt.traj[ids.race] <- sample(1:3, length(ids.race), TRUE,
                                 c(dat$param$tt.part.supp[i],
                                   dat$param$tt.full.supp[i],
                                   dat$param$tt.dur.supp[i]))

  }
  dat$attr$tt.traj <- tt.traj

  ## Infection-related attributes
  dat$attr$status <- status
  idsInf <- which(status == 1)

  age <- dat$attr$age
  min.ages <- min(dat$param$netstats$demog$ages)
  time.sex.active <- pmax(1, round((365/7)*age[idsInf] - (365/7)*min.ages, 0))
  min.hiv.time <- round(dat$param$vl.acute.rise.int + dat$param$vl.acute.fall.int)
  max.hiv.time <- dat$param$vl.aids.onset.int

  time.infected <- round(pmax(min.hiv.time,
                            pmin(time.sex.active,
                              sample(min.hiv.time:max.hiv.time, length(idsInf), TRUE))))

  dat$attr$inf.time <- rep(NA, num)
  dat$attr$inf.time[idsInf] <- -time.infected

  dat$attr$stage <- rep(NA, num)
  dat$attr$stage.time <- rep(NA, num)
  dat$attr$aids.time <- rep(NA, num)
  dat$attr$stage[idsInf] <- 3
  dat$attr$stage.time[idsInf] <- time.infected - min.hiv.time

  dat$attr$diag.stage <- rep(NA, num)
  dat$attr$diag.stage[idsInf] <- dat$attr$stage[idsInf]

  dat$attr$vl <- rep(NA, num)
  dat$attr$vl[idsInf] <- dat$param$vl.set.point
  dat$attr$vl.last.usupp <- rep(NA, num)
  dat$attr$vl.last.supp <- rep(NA, num)

  dat$attr$diag.time <- rep(NA, num)
  dat$attr$diag.time[idsInf] <- dat$attr$inf.time[idsInf] + round(mean(1/dat$param$hiv.test.rate))
  dat$attr$last.neg.test <- rep(NA, num)

  dat$attr$tx.status <- rep(NA, num)
  dat$attr$tx.status[idsInf] <- 0
  dat$attr$cuml.time.on.tx <- rep(NA, num)
  dat$attr$cuml.time.on.tx[idsInf] <- 0
  dat$attr$cuml.time.off.tx <- rep(NA, num)
  dat$attr$cuml.time.off.tx[idsInf] <- time.infected
  dat$attr$tx.period.first <- rep(NA, num)
  dat$attr$tx.period.last <- rep(NA, num)
  dat$attr$tx.init.time <- rep(NA, num)

  dat$attr$count.trans <- rep(0, num)
  dat$attr$num.neg.tests <- rep(0, length(status))

  return(dat)
}



#' @title Initialize the STI status of persons in the network
#'
#' @description Sets the initial individual-level disease status of persons
#'              in the network, as well as disease-related attributes for
#'              infected persons.
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#' @keywords initiation utility msm
#'
init_sti_msm <- function(dat) {

  role.class <- dat$attr$role.class
  num <- length(role.class)

  idsUreth <- which(role.class %in% c(0, 2))
  idsRect <- which(role.class %in% c(1, 2))

  uGC <- rGC <- rep(0, num)
  uCT <- rCT <- rep(0, num)

  # Initialize GC infection at both sites
  idsUGC <- sample(idsUreth, size = round(dat$init$prev.ugc * num), FALSE)
  uGC[idsUGC] <- 1

  idsRGC <- sample(setdiff(idsRect, idsUGC), size = round(dat$init$prev.rgc * num), FALSE)
  rGC[idsRGC] <- 1

  dat$attr$rGC <- rGC
  dat$attr$uGC <- uGC

  dat$attr$rGC.sympt <- dat$attr$uGC.sympt <- rep(NA, num)
  dat$attr$rGC.sympt[rGC == 1] <- rbinom(sum(rGC == 1), 1, dat$param$rgc.sympt.prob)
  dat$attr$uGC.sympt[uGC == 1] <- rbinom(sum(uGC == 1), 1, dat$param$ugc.sympt.prob)

  dat$attr$rGC.infTime <- dat$attr$uGC.infTime <- rep(NA, length(dat$attr$active))
  dat$attr$rGC.infTime[rGC == 1] <- 1
  dat$attr$uGC.infTime[uGC == 1] <- 1

  dat$attr$rGC.timesInf <- rep(0, num)
  dat$attr$rGC.timesInf[rGC == 1] <- 1
  dat$attr$uGC.timesInf <- rep(0, num)
  dat$attr$uGC.timesInf[uGC == 1] <- 1

  dat$attr$rGC.tx <- dat$attr$uGC.tx <- rep(NA, num)
  dat$attr$rGC.tx.prep <- dat$attr$uGC.tx.prep <- rep(NA, num)

  # Initialize CT infection at both sites
  idsUCT <- sample(idsUreth, size = round(dat$init$prev.uct * num), FALSE)
  uCT[idsUCT] <- 1

  idsRCT <- sample(setdiff(idsRect, idsUCT), size = round(dat$init$prev.rct * num), FALSE)
  rCT[idsRCT] <- 1

  dat$attr$rCT <- rCT
  dat$attr$uCT <- uCT

  dat$attr$rCT.sympt <- dat$attr$uCT.sympt <- rep(NA, num)
  dat$attr$rCT.sympt[rCT == 1] <- rbinom(sum(rCT == 1), 1, dat$param$rct.sympt.prob)
  dat$attr$uCT.sympt[uCT == 1] <- rbinom(sum(uCT == 1), 1, dat$param$uct.sympt.prob)

  dat$attr$rCT.infTime <- dat$attr$uCT.infTime <- rep(NA, num)
  dat$attr$rCT.infTime[dat$attr$rCT == 1] <- 1
  dat$attr$uCT.infTime[dat$attr$uCT == 1] <- 1

  dat$attr$rCT.timesInf <- rep(0, num)
  dat$attr$rCT.timesInf[rCT == 1] <- 1
  dat$attr$uCT.timesInf <- rep(0, num)
  dat$attr$uCT.timesInf[uCT == 1] <- 1

  dat$attr$rCT.tx <- dat$attr$uCT.tx <- rep(NA, num)
  dat$attr$rCT.tx.prep <- dat$attr$uCT.tx.prep <- rep(NA, num)

  return(dat)

}


#' @title Re-Initialization Module
#'
#' @description This function reinitializes an epidemic model to restart at a
#'              specified time step given an input \code{netsim} object.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}}.
#' @inheritParams initialize_msm
#'
#' @details
#' Currently, the necessary components that must be on \code{x} for a simulation
#' to be restarted must be: param, control, nwparam, epi, attr, temp, el, p.
#' TODO: describe this more.
#'
#' @return
#' This function resets the data elements on the \code{dat} master data object
#' in the needed ways for the time loop to function.
#'
#' @export
#' @keywords module msm
#'
reinit_msm <- function(x, param, init, control, s) {

  need.for.reinit <- c("param", "control", "nwparam", "epi", "attr", "temp", "el", "p")
  if (!all(need.for.reinit %in% names(x))) {
    stop("x must contain the following elements for restarting: ",
         "param, control, nwparam, epi, attr, temp, el, p",
         call. = FALSE)
  }

  if (length(x$el) == 1) {
    s <- 1
  }

  dat <- list()

  dat$param <- param
  dat$param$modes <- 1
  dat$control <- control
  dat$nwparam <- x$nwparam

  dat$epi <- sapply(x$epi, function(var) var[s])
  names(dat$epi) <- names(x$epi)

  dat$el <- x$el[[s]]
  dat$p <- x$p[[s]]

  dat$attr <- x$attr[[s]]

  if (!is.null(x$stats)) {
    dat$stats <- list()
    if (!is.null(x$stats$nwstats)) {
      dat$stats$nwstats <- x$stats$nwstats[[s]]
    }
  }

  dat$temp <- x$temp[[s]]

  class(dat) <- "dat"
  return(dat)
}



# HET -----------------------------------------------------------------


#' @export
#' @rdname initialize_msm
initialize_het <- function(x, param, init, control, s) {

  dat <- list()
  dat$temp <- list()
  nw <- simulate(x$fit, control = control.simulate.ergm(MCMC.burnin = 1e6))

  dat$el <- list()
  dat$el[[1]] <- as.edgelist(nw)
  attributes(dat$el)$vnames <- NULL
  p <- tergmLite::stergm_prep(nw, x$formation, x$coef.diss$dissolution, x$coef.form,
                              x$coef.diss$coef.adj, x$constraints)
  p$model.form$formula <- NULL
  p$model.diss$formula <- NULL
  dat$p <- list()
  dat$p[[1]] <- p

  ## Network Model Parameters
  dat$nwparam <- list(x[-which(names(x) == "fit")])

  ## Simulation Parameters
  dat$param <- param
  dat$param$modes <- 1

  dat$init <- init
  dat$control <- control

  ## Nodal Attributes
  dat$attr <- list()

  dat$attr$male <- get.vertex.attribute(nw, "male")

  n <- network.size(nw)
  dat$attr$active <- rep(1, n)
  dat$attr$entTime <- rep(1, n)

  dat <- initStatus_het(dat)

  age <- rep(NA, n)
  age[dat$attr$male == 0] <- sample(init$ages.feml, sum(dat$attr$male == 0), TRUE)
  age[dat$attr$male == 1] <- sample(init$ages.male, sum(dat$attr$male == 1), TRUE)
  dat$attr$age <- age

  dat <- initInfTime_het(dat)
  dat <- initDx_het(dat)
  dat <- initTx_het(dat)

  # Circumcision
  male <- dat$attr$male
  nMales <- sum(male == 1)
  age <- dat$attr$age

  circStat <- circTime <- rep(NA, n)

  circStat[male == 1] <- rbinom(nMales, 1, dat$param$circ.prob.birth)

  isCirc <- which(circStat == 1)
  circTime[isCirc] <- round(-age[isCirc] * (365 / dat$param$time.unit))

  dat$attr$circStat <- circStat
  dat$attr$circTime <- circTime


  ## Stats List
  dat$stats <- list()

  ## Final steps
  dat$epi <- list()
  dat <- prevalence_het(dat, at = 1)

}


#' @title Reinitialization Module
#'
#' @description This function reinitializes the master \code{dat} object on which
#'              data are stored, simulates the initial state of the network, and
#'              simulates disease status and other attributes.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param_het}}.
#' @param init An \code{EpiModel} object of class \code{\link{init_het}}.
#' @param control An \code{EpiModel} object of class \code{\link{control_het}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated \code{dat} object with the initialized values
#' for demographics and disease-related variables.
#'
#' @keywords module het
#'
#' @export
#'
reinit_het <- function(x, param, init, control, s) {

  need.for.reinit <- c("param", "control", "nwparam", "epi",
                       "attr", "temp", "el", "p")
  if (!all(need.for.reinit %in% names(x))) {
    stop("x must contain the following elements for restarting: ",
         "param, control, nwparam, epi, attr, temp, el, p",
         call. = FALSE)
  }

  if (length(x$el) == 1) {
    s <- 1
  }

  dat <- list()

  dat$param <- param
  dat$param$modes <- 1
  dat$control <- control
  dat$nwparam <- x$nwparam

  dat$epi <- sapply(x$epi, function(var) var[s])
  names(dat$epi) <- names(x$epi)

  dat$el <- x$el[[s]]
  dat$p <- x$p[[s]]

  dat$attr <- x$attr[[s]]

  if (!is.null(x$stats)) {
    dat$stats <- list()
    if (!is.null(x$stats$nwstats)) {
      dat$stats$nwstats <- x$stats$nwstats[[s]]
    }
  }

  dat$temp <- x$temp[[s]]

  class(dat) <- "dat"

  return(dat)
}


initStatus_het <- function(dat) {

  ## Variables
  i.prev.male <- dat$init$i.prev.male
  i.prev.feml <- dat$init$i.prev.feml

  male <- dat$attr$male
  idsMale <- which(male == 1)
  idsFeml <- which(male == 0)
  nMale <- length(idsMale)
  nFeml <- length(idsFeml)
  n <- nMale + nFeml

  ## Process
  status <- rep(0, n)
  status[sample(idsMale, round(i.prev.male * nMale))] <- 1
  status[sample(idsFeml, round(i.prev.feml * nFeml))] <- 1

  dat$attr$status <- status

  return(dat)
}


initInfTime_het <- function(dat) {

  status <- dat$attr$status
  n <- length(status)

  infecteds <- which(status == 1)
  infTime <- rep(NA, n)

  inf.time.dist <- dat$init$inf.time.dist

  if (inf.time.dist == "allacute") {
    max.inf.time <- dat$param$vl.acute.topeak + dat$param$vl.acute.toset
    infTime[infecteds] <- sample(0:(-max.inf.time), length(infecteds), TRUE)
  } else {
    max.inf.time <- dat$init$max.inf.time / dat$param$time.unit
    if (inf.time.dist == "geometric") {
      total.d.rate <- 1/max.inf.time
      infTime[infecteds] <- -rgeom(length(infecteds), total.d.rate)
    }
    if (inf.time.dist == "uniform") {
      infTime[infecteds] <- sample(0:(-max.inf.time), length(infecteds), TRUE)
    }
  }

  ## Enforce that time infected < age
  infTime[infecteds] <- pmax(infTime[infecteds],
                             1 - dat$attr$age[infecteds] * (365 / dat$param$time.unit))

  dat$attr$infTime <- infTime

  timeInf <- 1 - infTime
  dat$attr$ageInf <- pmax(0, dat$attr$age - round(timeInf) * (dat$param$time.unit / 365))

  stopifnot(all(dat$attr$ageInf[infecteds] <= dat$attr$age[infecteds]),
            all(dat$attr$ageInf[infecteds] >= 0))

  return(dat)
}


initDx_het <- function(dat) {

  n <- sum(dat$attr$active == 1)
  status <- dat$attr$status

  dxStat <- rep(NA, n)
  dxStat[status == 1] <- 0

  dxTime <- rep(NA, n)

  dat$attr$dxStat <- dxStat
  dat$attr$dxTime <- dxTime

  return(dat)
}


initTx_het <- function(dat) {

  ## Variables
  status <- dat$attr$status
  n <- sum(dat$attr$active == 1)
  nInf <- sum(status == 1)

  tx.init.cd4.mean <- dat$param$tx.init.cd4.mean
  tx.init.cd4.sd <- dat$param$tx.init.cd4.sd
  tx.elig.cd4 <- dat$param$tx.elig.cd4


  ## Process
  dat$attr$txStat <- rep(NA, n)
  dat$attr$txStartTime <- rep(NA, n)
  dat$attr$txStops <- rep(NA, n)
  dat$attr$txTimeOn <- rep(NA, n)
  dat$attr$txTimeOff <- rep(NA, n)

  txCD4min <- rep(NA, n)
  txCD4min[status == 1] <- pmin(rnbinom(nInf,
                                        size = nbsdtosize(tx.init.cd4.mean,
                                                          tx.init.cd4.sd),
                                        mu = tx.init.cd4.mean), tx.elig.cd4)
  dat$attr$txCD4min <- txCD4min
  dat$attr$txCD4start <- rep(NA, n)
  dat$attr$txType <- rep(NA, n)

  return(dat)
}
