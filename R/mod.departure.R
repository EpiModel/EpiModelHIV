
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
  mdd.suitry.status <- get_attr(dat, "mdd.suitry.status")
  mdd.suitry.time   <- get_attr(dat, "mdd.suitry.time")
  mdd.suicompl.status <- get_attr(dat, "mdd.suicompl.status")
  mdd.suicompl.time   <- get_attr(dat, "mdd.suicompl.time")


  age <- floor(age)


  # Parameters
  aids.on.tx.mort.rate  <- get_param(dat, "aids.on.tx.mort.rate")
  aids.off.tx.mort.rate <- get_param(dat, "aids.off.tx.mort.rate")
  time.unit <- get_param(dat, "time.unit")
  netstats <- get_param(dat, "netstats")
  age.sexual.cessation <- get_param(dat, "epistats")[["age.sexual.cessation"]]
  mdd.suicompl.prob <- get_param(dat, "mdd.suicompl.prob")

  asmr <- netstats[["demog"]][["asmr"]]



  ## AIDS-related deaths
  #1. On tx
  idsEligOn <- which(tx.status == 1 & stage == 4)
  idsDepAIDSOn <- idsEligOn[runif(length(idsEligOn)) < aids.on.tx.mort.rate]

  #2. Off tx
  idsEligOff <- which(tx.status == 0 & stage == 4)
  idsDepAIDSOff <- idsEligOff[runif(length(idsEligOff)) < aids.off.tx.mort.rate]

  idsDepAIDS <- union(idsDepAIDSOn, idsDepAIDSOff)



  ##Depression-related deaths
  idsEligSuicompl <- which(
    active == 1 &
      mdd.suitry.status == 1 &
      mdd.suitry.time == at - 1)
  idsEligSuicompl <- setdiff(idsEligSuicompl, idsDepAIDS)
  idsSuicompl <- idsEligSuicompl[runif(length(idsEligSuicompl)) < mdd.suicompl.prob]


  idsAIDSSuiCompl <- c(idsDepAIDS, idsSuicompl)



  ## General deaths
  idsElig <- which(as.logical(active))
  idsElig <- setdiff(idsElig, idsAIDSSuiCompl)
  rates <- rep(NA, length(idsElig))



  races <- sort(unique(race))
  for (i in seq_along(races)) {
    ids.race <- which(race[idsElig] == races[i])
    rates[ids.race] <- asmr[age[idsElig[ids.race]], i + 1]
  }
  idsDep <- idsElig[runif(length(rates)) < rates]



  idsDepAll <- unique(c(idsDep, idsDepAIDS, idsSuicompl))

  if (length(idsDepAll) > 0) {
    dat <- set_attr(dat, "active", 0, posit_ids = idsDepAll)
    dat <- depart_nodes(dat, idsDepAll)
    attr.length <- unique(vapply(dat[["attr"]], length, numeric(1)))
    if (attr.length != attributes(dat[["el"]][[1]])[["n"]]) {
      stop("mismatch between el and attr length in departures mod")
    }
  }


  ## Summary Output
  dat <- set_epi(dat, "dep.gen", at, length(idsDep))
  dat <- set_epi(dat, "dep.AIDS.on.tx", at, length(idsDepAIDSOn))
  dat <- set_epi(dat, "dep.AIDS.off.tx", at, length(idsDepAIDSOff))
  dat <- set_epi(dat, "mdd.eligsuicompl.numall", at, length(idsEligSuicompl))
  dat <- set_epi(dat, "mdd.suicompl.numall", at, length(idsSuicompl))
  dat <- set_epi(dat, "dep.all", at, length(idsDepAll))
  dat <- set_epi(dat, "mdd.suicompl.prpdepall", at,
                 length(idsSuicompl) / length(idsDepAll))


  disease.mr <- (length(idsDepAIDSOn) + length(idsDepAIDSOff)) /
    sum(status == 1, na.rm = TRUE) * (364 / time.unit)
  dat <- set_epi(dat, "disease.mr100", at, disease.mr * 100)

  ## handle sexual cessation age

  # updated age attr with any deleted vertices removed
  age <- get_attr(dat, "age")
  # which nodes are at or above the sexual cessation age?
  nodes_sexual_cessation <- which(age >= age.sexual.cessation)
  # delete any corresponding edges and lasttoggle information
  if (length(nodes_sexual_cessation) > 0) {
    for (i in 1:3) {
      dat[["el"]][[i]] <- delete_edges(dat[["el"]][[i]], nodes_sexual_cessation)

      if (i < 3 && get_network_control(dat, i, "tergmLite.track.duration") == TRUE) {
        dat$net_attr[[i]][["lasttoggle"]] <-
          delete_edges(dat$net_attr[[i]][["lasttoggle"]], nodes_sexual_cessation)
      }
    }
  }

  return(dat)
}
