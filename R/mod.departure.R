
#' @title Departure Module
#'
#' @description Module function for simulating departures departures, including deaths,
#'              among population members.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Departures are divided into three categories: general deaths, for which demographic
#' data on age- and race-specific mortality rates apply; AIDS-related deaths, for which the
#' rate of death depends on current ART status; and sexual activity cessation, which may remove
#' persons from the active network for reasons other than mortality.
#'
#' @return
#' This function returns the updated `dat` object accounting for deaths and other departures
#' from the population.
#'
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
  tx.status <- get_attr(dat, "tx.status")

  age <- floor(age)

  # Parameters
  aids.on.tx.mort.rate  <- get_param(dat, "aids.on.tx.mort.rate")
  aids.off.tx.mort.rate <- get_param(dat, "aids.off.tx.mort.rate")
  time.unit <- get_param(dat, "time.unit")
  netstats <- get_param(dat, "netstats")

  asmr <- netstats[["demog"]][["asmr"]]

  ## AIDS-related deaths
  #1. On tx
  idsEligOn <- which(tx.status == 1 & stage == 4)
  idsDepAIDSOn <- idsEligOn[runif(length(idsEligOn)) < aids.on.tx.mort.rate]

  #2. Off tx
  idsEligOff <- which(tx.status == 0 & stage == 4)
  idsDepAIDSOff <- idsEligOff[runif(length(idsEligOff)) < aids.off.tx.mort.rate]

  idsDepAIDS <- union(idsDepAIDSOn, idsDepAIDSOff)

  ## General deaths
  idsElig <- which(as.logical(active))
  idsElig <- setdiff(idsElig, idsDepAIDS)
  rates <- rep(NA, length(idsElig))

  races <- sort(unique(race))
  for (i in seq_along(races)) {
    ids.race <- which(race[idsElig] == races[i])
    rates[ids.race] <- asmr[age[idsElig[ids.race]], i + 1]
  }
  idsDep <- idsElig[runif(length(rates)) < rates]

  idsDepAll <- unique(c(idsDep, idsDepAIDS))

  if (length(idsDepAll) > 0) {
    dat <- set_attr(dat, "active", 0, posit_ids = idsDepAll)
    tergmLite.track.duration <- get_control(dat, "tergmLite.track.duration")
    for (i in 1:3) {
      dat[["el"]][[i]] <- delete_vertices(
        dat[["el"]][[i]],
        idsDepAll
      )

      if (i < 3 && tergmLite.track.duration == TRUE) {
        dat$nw[[i]] %n% "lasttoggle" <-
          delete_vertices(dat$nw[[i]] %n% "lasttoggle", idsDepAll)
      }
    }
    dat[["attr"]] <- deleteAttr(dat[["attr"]], idsDepAll)
    attr.length <- unique(vapply(dat[["attr"]], length, numeric(1)))
    if (attr.length != attributes(dat[["el"]][[1]])[["n"]]) {
      stop("mismatch between el and attr length in departures mod")
    }
  }


  ## Summary Output
  dat <- set_epi(dat, "dep.gen", at, length(idsDep))
  dat <- set_epi(dat, "dep.AIDS.on.tx", at, length(idsDepAIDSOn))
  dat <- set_epi(dat, "dep.AIDS.off.tx", at, length(idsDepAIDSOff))

  disease.mr <- (length(idsDepAIDSOn) + length(idsDepAIDSOff)) /
    sum(status == 1, na.rm = TRUE) * (364 / time.unit)
  dat <- set_epi(dat, "disease.mr100", at, disease.mr * 100)

  return(dat)
}
