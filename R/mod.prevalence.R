
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

  active <- dat$attr$active
  status <- dat$attr$status
  race <- dat$attr$race

  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat

  prepElig.la <- dat$attr$prepElig.la
  prepStat.la  <- dat$attr$prepStat.la

  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT

  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$num.B[at] <- sum(race == "B", na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == "W", na.rm = TRUE)

  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)

  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.prep[at] <- sum(status == 1 & (prepStat == 1 | prepStat.la == 1), na.rm = TRUE) /
                             sum(prepStat == 1 | prepStat.la == 1, na.rm = TRUE)
  dat$epi$i.prev.oprep[at] <- sum(status == 1 & (prepStat == 1), na.rm = TRUE) /
                              sum(prepStat == 1, na.rm = TRUE)
  dat$epi$i.prev.iprep[at] <- sum(status == 1 & (prepStat.la == 1), na.rm = TRUE) /
                              sum(prepStat.la == 1, na.rm = TRUE)

  dat$epi$ir100[at] <- (dat$epi$incid[at] / sum(status == 0, dat$epi$incid[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.prep[at] <- (dat$epi$incid.prep[at] / sum(status == 0 &
                                                            (prepStat == 1), na.rm = TRUE)) * 5200

  dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepElig.la[at] <- sum(prepElig.la == 1, na.rm = TRUE)
  dat$epi$prepCurr.la[at] <- sum(prepStat.la == 1, na.rm = TRUE)
  dat$epi$prepCurr.8w.la[at] <- sum(at - dat$attr$prepTimeLastInj <= 8, na.rm = TRUE)


  dat$epi$prev.gc[at] <- sum((rGC == 1 | uGC == 1), na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ct[at] <- sum((rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at]
  ir100.rgc <- (dat$epi$incid.rgc[at]/sum(rGC == 0, dat$epi$incid.rgc[at], na.rm = TRUE))*5200
  ir100.ugc <- (dat$epi$incid.ugc[at]/sum(uGC == 0, dat$epi$incid.ugc[at], na.rm = TRUE))*5200
  dat$epi$ir100.gc[at] <- ir100.rgc + ir100.ugc
  ir100.rct <- (dat$epi$incid.rct[at]/sum(rCT == 0, dat$epi$incid.rct[at], na.rm = TRUE))*5200
  ir100.uct <- (dat$epi$incid.uct[at]/sum(uCT == 0, dat$epi$incid.uct[at], na.rm = TRUE))*5200
  dat$epi$ir100.ct[at] <- ir100.rct + ir100.uct

  if (!is.null(dat$epi$incid.B)) {
    dat$epi$incid.B <- NULL
    dat$epi$incid.W <- NULL
  }

  return(dat)
}



#' @export
#' @rdname prevalence_msm
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

    dat$epi$prep.ind.met <- rNA
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
