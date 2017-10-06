
#' @title Death Module for up to 5 race groups heterosexuals and MSM.
#'
#' @description Module function for simulting both general and disease-related
#'              deaths among population members.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Deaths are divided into two categories: general deaths, for which demographic
#' data on age-specific mortality rates applies; and disease-related diseases,
#' for which the rate of death is a function of progression to end-stage AIDS.
#' Which nodes have died is determined stochastically for general deaths using
#' draws from a binomial distribution, and deterministically for disease-related
#' deaths after nodes have reach a maximum viral load value set in the
#' \code{vl.fatal} parameter.
#'
#' @return
#' This function returns the updated \code{dat} object accounting for deaths.
#' The deaths are deactivated from the main and casual networks, as those are in
#' \code{networkDynamic} class objects; dead nodes are not deleted from the
#' instant network until the \code{\link{simnet_shamp}} module for bookkeeping
#' purposes.
#'
#' @keywords module shamp msm het
#' @export
#'
deaths_shamp <- function(dat, at) {

  ## General deaths
  age <- floor(dat$attr$age)
  race <- dat$attr$race
  sex <- dat$attr$sex

  alive.B.f <- which(race == "B" & sex == "F")
  age.B.f <- age[alive.B.f]
  death.B.f.prob <- dat$param$asmr.B.f[age.B.f]
  deaths.B.f <- alive.B.f[rbinom(length(death.B.f.prob), 1, death.B.f.prob) == 1]

  alive.BI.f <- which(race == "BI" & sex == "F")
  age.BI.f <- age[alive.BI.f]
  death.BI.f.prob <- dat$param$asmr.BI.f[age.BI.f]
  deaths.BI.f <- alive.BI.f[rbinom(length(death.BI.f.prob), 1, death.BI.f.prob) == 1]

  alive.H.f <- which(race == "H" & sex == "F")
  age.H.f <- age[alive.H.f]
  death.H.f.prob <- dat$param$asmr.H.f[age.H.f]
  deaths.H.f <- alive.H.f[rbinom(length(death.H.f.prob), 1, death.H.f.prob) == 1]
  
  alive.HI.f <- which(race == "HI" & sex == "F")
  age.HI.f <- age[alive.HI.f]
  death.HI.f.prob <- dat$param$asmr.HI.f[age.HI.f]
  deaths.HI.f <- alive.HI.f[rbinom(length(death.HI.f.prob), 1, death.HI.f.prob) == 1]
  
  alive.W.f <- which(race == "W" & sex == "F")
  age.W.f <- age[alive.W.f]
  death.W.f.prob <- dat$param$asmr.W.f[age.W.f]
  deaths.W.f <- alive.W.f[rbinom(length(death.W.f.prob), 1, death.W.f.prob) == 1]
  
  alive.B.m <- which(race == "B" & sex == "M")
  age.B.m <- age[alive.B.m]
  death.B.m.prob <- dat$param$asmr.B.m[age.B.m]
  deaths.B.m <- alive.B.m[rbinom(length(death.B.m.prob), 1, death.B.m.prob) == 1]
  
  alive.BI.m <- which(race == "BI" & sex == "M")
  age.BI.m <- age[alive.BI.m]
  death.BI.m.prob <- dat$param$asmr.BI.m[age.BI.m]
  deaths.BI.m <- alive.BI.m[rbinom(length(death.BI.m.prob), 1, death.BI.m.prob) == 1]
  
  alive.H.m <- which(race == "H" & sex == "M")
  age.H.m <- age[alive.H.m]
  death.H.m.prob <- dat$param$asmr.H.m[age.H.m]
  deaths.H.m <- alive.H.m[rbinom(length(death.H.m.prob), 1, death.H.m.prob) == 1]
  
  alive.HI.m <- which(race == "HI" & sex == "M")
  age.HI.m <- age[alive.HI.m]
  death.HI.m.prob <- dat$param$asmr.HI.m[age.HI.m]
  deaths.HI.m <- alive.HI.m[rbinom(length(death.HI.m.prob), 1, death.HI.m.prob) == 1]
  
  alive.W.m <- which(race == "W" & sex == "M")
  age.W.m <- age[alive.W.m]
  death.W.m.prob <- dat$param$asmr.W.m[age.W.m]
  deaths.W.m <- alive.W.m[rbinom(length(death.W.m.prob), 1, death.W.m.prob) == 1]
  
  dth.gen <- c(deaths.B.f, deaths.BI.f, deaths.H.f, deaths.HI.f, deaths.W.f,
               deaths.B.m, deaths.BI.m, deaths.H.m, deaths.HI.m, deaths.W.m)


  ## Disease deaths
  dth.dis <- which(dat$attr$stage == 4 &
                   dat$attr$vl >= dat$param$vl.fatal)

  dth.all <- NULL
  dth.all <- unique(c(dth.gen, dth.dis))
  
  ##Those that will age out this time step.
  dth.age<-which(age >= dat$param$exit.age)

  if (length(dth.all) > 0) {
    dat$attr$active[dth.all] <- 0
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::delete_vertices(dat$el[[i]], dth.all)
    }
    dat$attr <- deleteAttr(dat$attr, dth.all)
    if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
      stop("mismatch between el and attr length in death mod")
    }
  }



  ## Summary Output
  dat$epi$dth.gen[at] <- max(0,length(dth.gen)-length(dth.age))
  dat$epi$dth.dis[at] <- max(0,length(dth.dis))
  dat$epi$dth.age[at] <-max(0,length(dth.age))
  
  return(dat)
}

