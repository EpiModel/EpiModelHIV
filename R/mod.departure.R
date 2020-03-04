
#' @title Depature Module
#'
#' @description Module function for simulting both general and disease-related
#'              departures, including deaths, among population members.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Deaths are divided into two categories: general deaths, for which demographic
#' data on age-specific mortality rates applies; and disease-related diseases,
#' for which the rate of death is a function of progression to end-stage AIDS.
#'
#' @return
#' This function returns the updated \code{dat} object accounting for deaths.
#' The deaths are deactivated from the main and casual networks, as those are in
#' \code{networkDynamic} class objects; dead nodes are not deleted from the
#' instant network until the \code{\link{simnet_msm}} module for bookkeeping
#' purposes.
#'
#' @keywords module msm
#' @export
#'
departure_msm <- function(dat, at) {

  ## General departures
  active <- dat$attr$active
  age <- floor(dat$attr$age)
  race <- dat$attr$race
  status <- dat$attr$status
  stage <- dat$attr$stage
  tx.status <- dat$attr$tx.status

  aids.mr <- dat$param$aids.mr
  asmr <- dat$param$netstats$demog$asmr

  idsElig <- which(active == 1)
  rates <- rep(NA, length(idsElig))

  races <- sort(unique(race))
  for (i in seq_along(races)) {
    ids.race <- which(race == races[i])
    rates[ids.race] <- asmr[age[ids.race], i + 1]
  }
  idsDep <- idsElig[rbinom(length(rates), 1, rates) == 1]

  ## HIV-related deaths
  idsEligAIDS <- which(stage == 4)
  idsDepAIDS <- idsEligAIDS[rbinom(length(idsEligAIDS), 1, aids.mr) == 1]

  idsDepAll <- unique(c(idsDep, idsDepAIDS))
  depHIV <- intersect(idsDepAll, which(status == 1))
  depHIV.old <- intersect(depHIV, which(age >= 65))

  # Cumulative R0 calculations
  # if (at == 2) {
  #   dat$temp$R0 <- NA
  # }
  # if (length(depHIV) > 0) {
  #   newR0 <- dat$attr$count.trans[depHIV]
  #   dat$temp$R0 <- c(dat$temp$R0, newR0)
  # }

  if (length(idsDepAll) > 0) {
    dat$attr$active[idsDepAll] <- 0
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::delete_vertices(dat$el[[i]], idsDepAll)
    }
    dat$attr <- deleteAttr(dat$attr, idsDepAll)
    if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
      stop("mismatch between el and attr length in departures mod")
    }
  }

  # Update clinical history
  if (dat$control$save.clin.hist == TRUE & length(idsDepAll) > 0) {
    m <- dat$temp$clin.hist
    for (i in 1:length(m)) {
      m[[i]] <- m[[i]][-idsDepAll, ]
    }
    dat$temp$clin.hist <- m
  }

  ## Summary Output
  dat$epi$dep.gen[at] <- length(idsDep)
  dat$epi$dep.AIDS[at] <- length(idsDepAIDS)
  dat$epi$dep.HIV[at] <- length(depHIV)
  dat$epi$dep.HIV.old[at] <- length(depHIV.old)

  return(dat)
}


#' @export
#' @rdname departure_msm
deaths_het <- function(dat, at) {

  ### 1. Susceptible Deaths ###

  ## Variables
  male <- dat$attr$male
  age <- dat$attr$age
  cd4Count <- dat$attr$cd4Count

  di.cd4.aids <- dat$param$di.cd4.aids
  ds.exit.age <- dat$param$ds.exit.age

  ## Eligible are: active uninf, pre-death infected, unhealthy old
  idsEligSus <- which((is.na(cd4Count) |
                       cd4Count > di.cd4.aids |
                       (cd4Count <= di.cd4.aids & age > ds.exit.age)))
  nEligSus <- length(idsEligSus)

  # Set age-sex specific rates
  ds.rates <- dat$param$ds.rates
  if (nEligSus > 0) {
    rates <- ds.rates$mrate[100*male[idsEligSus] + age[idsEligSus]]
  }


  ## Process
  nDeathsSus <- 0; idsDeathsSus <- NULL
  if (nEligSus > 0) {
    vecDeathsSus <- which(rbinom(nEligSus, 1, rates) == 1)
    nDeathsSus <- length(vecDeathsSus)
  }


  ## Update Attributes
  if (nDeathsSus > 0) {
    idsDeathsSus <- idsEligSus[vecDeathsSus]
    dat$attr$active[idsDeathsSus] <- 0
  }


  ### 2. Infected Deaths ###

  ## Variables
  active <- dat$attr$active
  di.cd4.rate <- dat$param$di.cd4.rate

  ## Process
  nDeathsInf <- 0; idsDeathsInf <- NULL

  cd4Count <- dat$attr$cd4Count
  stopifnot(length(active) == length(cd4Count))

  idsEligInf <- which(active == 1 & cd4Count <= di.cd4.aids)
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {
    vecDeathsInf <- which(rbinom(nEligInf, 1, di.cd4.rate) == 1)
    if (length(vecDeathsInf) > 0) {
      idsDeathsInf <- idsEligInf[vecDeathsInf]
      nDeathsInf <- length(idsDeathsInf)
    }
  }

  idsDeathsDet <- which(cd4Count <= 0)
  if (length(idsDeathsDet) > 0) {
    idsDeathsInf <- c(idsDeathsInf, idsDeathsDet)
    nDeathsInf <- nDeathsInf + length(idsDeathsDet)
  }


  ### 3. Update Attributes ###
  if (nDeathsInf > 0) {
    dat$attr$active[idsDeathsInf] <- 0
  }

  ## 4. Update Population Structure ##
  inactive <- which(dat$attr$active == 0)
  dat$el[[1]] <- tergmLite::delete_vertices(dat$el[[1]], inactive)
  dat$attr <- deleteAttr(dat$attr, inactive)

  if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
    stop("mismatch between el and attr length in death mod")
  }

  ### 5. Summary Statistics ###
  dat$epi$ds.flow[at] <- nDeathsSus
  dat$epi$di.flow[at] <- nDeathsInf

  return(dat)
}
