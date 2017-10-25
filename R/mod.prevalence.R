
#' @title Prevalence Calculations within Time Steps
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Summary statistic calculations are of two broad forms: prevalence and
#' incidence. This function establishes the summary statistic vectors for both
#' prevalence and incidence at time 1, and then calculates the prevalence
#' statistics for times 2 onward. Incidence statistics (e.g., number of new
#' infections or deaths) are calculated within the modules as they depend on
#' vectors that are not stored external to the module.
#'
#' @return
#' This function returns the \code{dat} object with an updated summary of current
#' attributes stored in \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
prevalence_msm <- function(dat, at) {

  ## Variables

  # Attributes

  active <- dat$attr$active
  race <- dat$attr$race
  status <- dat$attr$status

  prepAware <- dat$attr$prepAware
  prepAccess <- dat$attr$prepAccess
  prepIndic <- dat$attr$prepIndic
  prepIndic1 <- dat$attr$prepIndic1
  prepIndic2 <- dat$attr$prepIndic2
  prepIndic3 <- dat$attr$prepIndic3
  prepIndic4 <- dat$attr$prepIndic4
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  if (at == 1) {
    dat$epi$num <- dat$epi$num.B <- dat$epi$num.W <- rNA
    dat$epi$s.num <- dat$epi$s.num.B <- dat$epi$s.num.W <- rNA
    dat$epi$i.num <- dat$epi$i.num.B <- dat$epi$i.num.W <- rNA
    dat$epi$i.prev <- dat$epi$i.prev.B <- dat$epi$i.prev.W <- rNA

    dat$epi$incid <- dat$epi$incid.B <- dat$epi$incid.W <- rNA
    dat$epi$ir100 <- dat$epi$ir100.B <- dat$epi$ir100.W <- rNA

    dat$epi$prepAware.B <- dat$epi$prepAware.W <- rNA
    dat$epi$prepAccess.B <- dat$epi$prepAccess.W <- rNA
    dat$epi$prepIndic.B <- dat$epi$prepIndic.W <- rNA

    dat$epi$prepRx.B <- dat$epi$prepRx.W <- rNA
    dat$epi$prepCurr.B <- dat$epi$prepCurr.W <- rNA
    dat$epi$prepHiAdr.B <- dat$epi$prepHiAdr.W <- rNA

    dat$epi$prev.gc <- rNA
    dat$epi$prev.ct <- rNA

    dat$epi$prep.sens <- rNA
    dat$epi$prep.spec <- rNA

    dat$epi$prep.sens.1y <- dat$epi$prep.sens.2y <- dat$epi$prep.sens.3y <-
      dat$epi$prep.sens.4y <- dat$epi$prep.sens.5y <- rNA
    dat$epi$prep.spec.1y <- dat$epi$prep.spec.2y <- dat$epi$prep.spec.3y <-
      dat$epi$prep.spec.4y <- dat$epi$prep.spec.5y <- rNA

    dat$epi$prep.sens.ftime <- rNA
    dat$epi$prep.sens.ltime <- rNA

    dat$epi$incid.gc <- dat$epi$incid.gc.B <- dat$epi$incid.gc.W <- rNA
    dat$epi$incid.ct <- dat$epi$incid.ct.B <- dat$epi$incid.ct.W <- rNA

    dat$epi$ir100.gc <- dat$epi$ir100.gc.B <- dat$epi$ir100.gc.W <- rNA
    dat$epi$ir100.ct <- dat$epi$ir100.ct.B <- dat$epi$ir100.ct.W <- rNA
  }

  if (is.null(dat$epi$prepIndic1.B)) {
    dat$epi$prepIndic1.B <- dat$epi$prepIndic1.W <- rNA
    dat$epi$prepIndic2.B <- dat$epi$prepIndic2.W <- rNA
    dat$epi$prepIndic3.B <- dat$epi$prepIndic3.W <- rNA
    dat$epi$prepIndic4.B <- dat$epi$prepIndic4.W <- rNA
  }

  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$num.B[at] <- sum(race == "B", na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == "W", na.rm = TRUE)

  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)
  dat$epi$s.num.B[at] <- sum(status == 0 & race == "B", na.rm = TRUE)
  dat$epi$s.num.W[at] <- sum(status == 0 & race == "W", na.rm = TRUE)

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$i.num.B[at] <- sum(status == 1 & race == "B", na.rm = TRUE)
  dat$epi$i.num.W[at] <- sum(status == 1 & race == "W", na.rm = TRUE)

  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.B[at] <- dat$epi$i.num.B[at] / dat$epi$num.B[at]
  dat$epi$i.prev.W[at] <- dat$epi$i.num.W[at] / dat$epi$num.W[at]

  dat$epi$ir100[at] <- (dat$epi$incid[at] / sum(status == 0, na.rm = TRUE)) * 5200
  dat$epi$ir100.B[at] <- (dat$epi$incid.B[at] / sum(status == 0 & race == "B", na.rm = TRUE)) * 5200
  dat$epi$ir100.W[at] <- (dat$epi$incid.W[at] / sum(status == 0 & race == "W", na.rm = TRUE)) * 5200

  dat$epi$prepAware.B[at] <- sum(prepAware == 1 & race == "B", na.rm = TRUE) /
                             sum(race == "B", na.rm = TRUE)
  dat$epi$prepAware.W[at] <- sum(prepAware == 1 & race == "W", na.rm = TRUE) /
                             sum(race == "W", na.rm = TRUE)
  dat$epi$prepAccess.B[at] <- sum(prepAccess == 1 & prepAware == 1 & race == "B", na.rm = TRUE) /
                              sum(prepAware == 1 & race == "B", na.rm = TRUE)
  dat$epi$prepAccess.W[at] <- sum(prepAccess == 1 & prepAware == 1 & race == "W", na.rm = TRUE) /
                              sum(prepAware == 1 & race == "W", na.rm = TRUE)
  dat$epi$prepIndic.B[at] <- sum(prepIndic == 1 & race == "B", na.rm = TRUE) /
                             sum(race == "B", na.rm = TRUE)
  dat$epi$prepIndic.W[at] <- sum(prepIndic == 1 & race == "W", na.rm = TRUE) /
                             sum(race == "W", na.rm = TRUE)

  dat$epi$prepIndic1.B[at] <- sum(prepIndic1 == 1 & race == "B", na.rm = TRUE) /
                              sum(race == "B", na.rm = TRUE)
  dat$epi$prepIndic1.W[at] <- sum(prepIndic1 == 1 & race == "W", na.rm = TRUE) /
                              sum(race == "W", na.rm = TRUE)
  dat$epi$prepIndic2.B[at] <- sum(prepIndic2 == 1 & race == "B", na.rm = TRUE) /
                              sum(race == "B", na.rm = TRUE)
  dat$epi$prepIndic2.W[at] <- sum(prepIndic2 == 1 & race == "W", na.rm = TRUE) /
                              sum(race == "W", na.rm = TRUE)
  dat$epi$prepIndic3.B[at] <- sum(prepIndic3 == 1 & race == "B", na.rm = TRUE) /
                              sum(race == "B", na.rm = TRUE)
  dat$epi$prepIndic3.W[at] <- sum(prepIndic3 == 1 & race == "W", na.rm = TRUE) /
                              sum(race == "W", na.rm = TRUE)
  dat$epi$prepIndic4.B[at] <- sum(prepIndic4 == 1 & race == "B", na.rm = TRUE) /
                              sum(race == "B", na.rm = TRUE)
  dat$epi$prepIndic4.W[at] <- sum(prepIndic4 == 1 & race == "W", na.rm = TRUE) /
                              sum(race == "W", na.rm = TRUE)

  dat$epi$prepRx.B[at] <- sum(prepAccess == 1 & prepIndic == 1 & prepStat == 1 & race == "B", na.rm = TRUE) /
                          sum(prepAccess == 1 & prepIndic == 1 & race == "B", na.rm = TRUE)
  dat$epi$prepRx.W[at] <- sum(prepAccess == 1 & prepIndic == 1 & prepStat == 1 & race == "W", na.rm = TRUE) /
                          sum(prepAccess == 1 & prepIndic == 1 & race == "W", na.rm = TRUE)
  dat$epi$prepCurr.B[at] <- sum(prepStat == 1 & race == "B", na.rm = TRUE)
  dat$epi$prepCurr.W[at] <- sum(prepStat == 1 & race == "W", na.rm = TRUE)
  dat$epi$prepHiAdr.B[at] <- sum(prepStat == 1 & prepClass == 3 & race == "B", na.rm = TRUE) /
                             sum(prepStat == 1 & race == "B", na.rm = TRUE)
  dat$epi$prepHiAdr.W[at] <- sum(prepStat == 1 & prepClass == 3 & race == "W", na.rm = TRUE) /
                             sum(prepStat == 1 & race == "W", na.rm = TRUE)

  dat$epi$prev.gc[at] <- sum((rGC == 1 | uGC == 1), na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ct[at] <- sum((rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at]
  dat$epi$ir100.gc[at] <- (dat$epi$incid.gc[at] /
                             (sum(rGC == 0, na.rm = TRUE) +
                                sum(uGC == 0, na.rm = TRUE))) * 5200
  dat$epi$ir100.gc.B[at] <- (dat$epi$incid.gc.B[at] /
                               (sum(rGC == 0 & race == "B", na.rm = TRUE) +
                                  sum(uGC == 0 & race == "B", na.rm = TRUE))) * 5200
  dat$epi$ir100.gc.W[at] <- (dat$epi$incid.gc.W[at] /
                               (sum(rGC == 0 & race == "W", na.rm = TRUE) +
                                  sum(uGC == 0 & race == "W", na.rm = TRUE))) * 5200
  dat$epi$ir100.ct[at] <- (dat$epi$incid.ct[at] /
                             (sum(rCT == 0, na.rm = TRUE) +
                                sum(uCT == 0, na.rm = TRUE))) * 5200

  dat$epi$ir100.ct.B[at] <- (dat$epi$incid.ct.B[at] /
                             (sum(rCT == 0 & race == "B", na.rm = TRUE) +
                                sum(uCT == 0 & race == "B", na.rm = TRUE))) * 5200
  dat$epi$ir100.ct.W[at] <- (dat$epi$incid.ct.W[at] /
                               (sum(rCT == 0 & race == "W", na.rm = TRUE) +
                                  sum(uCT == 0 & race == "W", na.rm = TRUE))) * 5200

  # new stats
  prepElig.ever <- dat$attr$prepElig.ever
  prepElig.first <- dat$attr$prepElig.first
  prepElig.last <- dat$attr$prepElig.last
  inf.time <- dat$attr$inf.time

  dat$epi$prep.sens[at] <- sum(prepElig.ever[status == 1] == 1, na.rm = TRUE) /
                           sum(status == 1, na.rm = TRUE)

  dat$epi$prep.spec[at] <- sum(prepElig.ever[status == 0] == 0, na.rm = TRUE) /
                           sum(status == 0, na.rm = TRUE)

  for (j in 1:5) {
    num <- which(status == 1 & (at - inf.time <= j * 52) & (at - prepElig.last <= j * 52))
    den <- which(status == 1 & (at - inf.time <= j * 52))
    vname <- paste0("prep.sens.", j, "y")
    dat$epi[[vname]][at] <- length(num) / length(den)

    num <- which(dat$attr$status == 0 & (prepElig.ever == 0 | (at - prepElig.last >= j * 52)))
    den <- which(dat$attr$status == 0)
    vname <- paste0("prep.spec.", j, "y")
    dat$epi[[vname]][at] <- length(num) / length(den)
  }

  Infprep <- which(status == 1 & dat$attr$prepElig.ever == 1)
  dat$epi$prep.sens.ftime[at] <- mean(inf.time[Infprep] - prepElig.first[Infprep], na.rm = TRUE)
  dat$epi$prep.sens.ltime[at] <- mean(inf.time[Infprep] - prepElig.last[Infprep], na.rm = TRUE)

  return(dat)
}


#' @title Prevalence Module
#'
#' @description Module function to calculate and store summary statistics for
#'              disease prevalence, demographics, and other epidemiological
#'              outcomes.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
prevalence_het <- function(dat, at) {

  status <- dat$attr$status
  male <- dat$attr$male
  age <- dat$attr$age

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  # Initialize vectors
  if (at == 1) {
    dat$epi$i.num <- rNA
    dat$epi$num <- rNA

    dat$epi$i.num.male <- rNA
    dat$epi$i.num.feml <- rNA
    dat$epi$i.prev.male <- rNA
    dat$epi$i.prev.feml <- rNA

    dat$epi$num.male <- rNA
    dat$epi$num.feml <- rNA
    dat$epi$meanAge <- rNA
    dat$epi$propMale <- rNA

    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA
  }

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$num[at] <- length(status)

  dat$epi$i.num.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.num.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE)
  dat$epi$i.prev.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE) /
    sum(male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE) /
    sum(male == 0, na.rm = TRUE)

  dat$epi$num.male[at] <- sum(male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(male == 0, na.rm = TRUE)
  dat$epi$meanAge[at] <- mean(age, na.rm = TRUE)
  dat$epi$propMale[at] <- mean(male, na.rm = TRUE)

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}
