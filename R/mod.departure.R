
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

  ## Input
  # Attributes
  active <- get_attr(dat, "active")
  age    <- get_attr(dat, "age")
  race   <- get_attr(dat, "race")
  status <- get_attr(dat, "status")
  stage  <- get_attr(dat, "stage")
  uid    <- get_unique_ids(dat)

  age <- floor(age)

  # Parameters
  aids.mr  <- get_param(dat, "aids.mr")
  netstats <- get_param(dat, "netstats")

  asmr <- netstats$demog$asmr

  idsElig <- which(as.logical(active))
  rates <- rep(NA, length(idsElig))

  races <- sort(unique(race))
  for (i in seq_along(races)) {
    ids.race <- which(race == races[i])
    rates[ids.race] <- asmr[age[ids.race], i + 1]
  }
  idsDep <- idsElig[runif(length(rates)) < rates]

  ## HIV-related deaths
  idsEligAIDS <- which(stage == 4)
  idsDepAIDS <- idsEligAIDS[runif(length(idsEligAIDS)) < aids.mr]

  idsDepAll <- unique(c(idsDep, idsDepAIDS))
  depHIV <- intersect(idsDepAll, which(as.logical(status)))
  depHIV.old <- intersect(depHIV, which(age >= 65))

  # Update partnership list
  if (length(idsDepAll) > 0) {
    plist.dep <- dat$temp$plist
    uid.depart <- uid[idsDepAll]
    depart1 <- which(plist.dep[, "p1"] %in% uid.depart)
    depart2 <- which(plist.dep[, "p2"] %in% uid.depart)
    depart <- unique(depart1, depart2)
    ongoing <- which(is.na(plist.dep[, "stop"]))
    depart.ongoing <- intersect(depart, ongoing)
    plist.dep[depart.ongoing, "stop"] <- at
    dat$temp$plist <- plist.dep
  }

  if (length(idsDepAll) > 0) {
    dat <- set_attr(dat, "active", 0, posit_ids = idsDepAll)
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::delete_vertices(dat$el[[i]], idsDepAll)
    }
    dat$attr <- deleteAttr(dat$attr, idsDepAll)
    if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
      stop("mismatch between el and attr length in departures mod")
    }
  }


  ## Summary Output
  dat <- set_epi(dat, "dep.gen", at, length(idsDep))
  dat <- set_epi(dat, "dep.AIDS", at, length(idsDepAIDS))
  dat <- set_epi(dat, "dep.HIV", at, length(depHIV))
  dat <- set_epi(dat, "dep.HIV.old", at, length(depHIV.old))

  return(dat)
}
