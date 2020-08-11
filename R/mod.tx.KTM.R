
#' @title Treatment Module for up to 5 race groups heterosexuals and MSM.
#'
#' @description Module function for anti-retroviral treatment initiation and
#'              adherence over time.
#'
#' @inheritParams aging_shamp
#'
#' @details
#' Persons enter into the simulation with one of four ART "patterns": never
#' tested, tested but never treated, treated and achieving partial HIV viral
#' suppression, or treated with full viral suppression (these types are stored
#' as individual-level attributes in \code{tt.traj}). This module initiates ART
#' for treatment naive persons in the latter two types, and then cycles them on
#' and off treatment conditional on empirical race and sexual idenitity adherence rates. ART
#' initiation, non-adherence, and restarting are all stochastically simulated
#' based on binomial statistical models.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{tx.status},
#' \code{tx.init.time}, \code{cum.time.on.tx}, \code{cum.time.off.tx} attributes.
#'
#' @keywords module msm het
#'
#' @export
#'
tx_KTM <- function(dat, at) {

  ## Variables

  # Attributes
  race <- dat$attr$race
  sex <- dat$attr$sex
  sex.ident <- dat$attr$sex.ident
  status <- dat$attr$status
  tx.status <- dat$attr$tx.status
  diag.status <- dat$attr$diag.status
  diag.time <- dat$attr$diag.time
  tt.traj <- dat$attr$tt.traj
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  stage <- dat$attr$stage
  PS.diag.pos.time <- dat$attr$PS.diag.pos.time
  age<-dat$attr$age

  # Parameters
  #The init rate is used for the first 8 weeks to get to 75% starting within 8 weeks
  tx.init.B.f.prob <- dat$param$tx.init.B.f.prob
  tx.init.B.msf.prob <- dat$param$tx.init.B.msf.prob

  #The Halt matches 90% retention at 1 year
  tx.halt.B.f.prob <- dat$param$tx.halt.B.f.prob
  tx.halt.B.msf.prob <- dat$param$tx.halt.B.msf.prob

  #reinit is use for starting after 8 weeks and re-inititation (may be used as a free parameter to hit 75% coberage among PLHIV)
  tx.reinit.B.f.prob <- dat$param$tx.reinit.B.f.prob
  tx.reinit.B.msf.prob <- dat$param$tx.reinit.B.msf.prob
  
  #The start rate for those that are identified via partner services
  tx.init.PS.prob <- dat$param$tx.init.PS.prob

  PS.time <- dat$param$PS.time
  
  ## Initiation
  #Females
  
  #First 8 weeks
  tx.init.elig.B.f <- which(race == "B" & status == 1 & sex == "F" &
                          tx.status == 0 & diag.status == 1 & diag.time >= (at-8) & cum.time.on.tx == 0 & age < 40)
  
  tx.init.B.f <- tx.init.elig.B.f[rbinom(length(tx.init.elig.B.f), 1,
                                     tx.init.B.f.prob) == 1]
  
  #After 8 weeks
  tx.init.elig.B.f <- which(race == "B" & status == 1 & sex == "F" &
                              tx.status == 0 & diag.status == 1 & diag.time < (at-8) & cum.time.on.tx == 0 & age < 40)
  
  tx.init.B.f.late <- tx.init.elig.B.f[rbinom(length(tx.init.elig.B.f), 1,
                                              tx.reinit.B.f.prob) == 1]
  
  #First 8 weeks
  tx.init.elig.B.msf <- which(race == "B" & status == 1 & sex == "M" & sex.ident=="msf" &
                            tx.status == 0 & diag.status == 1 & diag.time >= (at-8) & cum.time.on.tx == 0 & age < 40)
  
  tx.init.B.msf <- tx.init.elig.B.msf[rbinom(length(tx.init.elig.B.msf), 1,
                                     tx.init.B.msf.prob) == 1]
  
  #After 8 weeks
  tx.init.elig.B.msf <- which(race == "B" & status == 1 & sex == "M" & sex.ident=="msf" &
                                tx.status == 0 & diag.status == 1 & diag.time < (at-8) & cum.time.on.tx == 0 & age < 40)
  
  tx.init.B.msf.late <- tx.init.elig.B.msf[rbinom(length(tx.init.elig.B.msf), 1,
                                             tx.reinit.B.msf.prob) == 1]
  
  #Due to partner services
  tx.init.elig.PS <- which(PS.diag.pos.time > at - PS.time & tx.status == 0 & diag.status == 1)
  
  tx.init.PS <- tx.init.elig.PS[rbinom(length(tx.init.elig.PS), 1, tx.init.PS.prob) == 1]
  
  dat$epi$tx.init.ps[at] <- max(0,length(tx.init.PS))
  
  
  tx.init <- c(tx.init.B.f, tx.init.B.f.late, tx.init.B.msf, tx.init.B.msf.late, tx.init.PS)

  dat$attr$tx.status[tx.init] <- 1
  dat$attr$tx.init.time[tx.init] <- at


  ## Halting
  ##females
  tx.halt.elig.B.f <- which(race == "B" & tx.status == 1 & sex=="F" & age < 40)
  tx.halt.B.f <- tx.halt.elig.B.f[rbinom(length(tx.halt.elig.B.f), 1,
                                     tx.halt.B.f.prob) == 1]
  
  ##males HET msf
  tx.halt.elig.B.msf <- which(race == "B" & tx.status == 1 & sex=="M" & sex.ident=="msf" & age < 40)
  tx.halt.B.msf <- tx.halt.elig.B.msf[rbinom(length(tx.halt.elig.B.msf), 1,
                                     tx.halt.B.msf.prob) == 1]
  

  
  tx.halt <- c(tx.halt.B.f, tx.halt.B.msf)
  
  dat$attr$tx.status[tx.halt] <- 0


  ## Restarting
  #females
  tx.reinit.elig.B.f <- which(race == "B" & tx.status == 0 & sex=="F" &
                            cum.time.on.tx > 0 & age < 40)
  tx.reinit.B.f <- tx.reinit.elig.B.f[rbinom(length(tx.reinit.elig.B.f),
                                         1, tx.reinit.B.f.prob) == 1]
  
  #males HET msf
  tx.reinit.elig.B.msf <- which(race == "B" & tx.status == 0 & sex=="M" &
                              cum.time.on.tx > 0 & age < 40)
  tx.reinit.B.msf <- tx.reinit.elig.B.msf[rbinom(length(tx.reinit.elig.B.msf),
                                         1, tx.reinit.B.msf.prob) == 1]
  
 
  tx.reinit <- c(tx.reinit.B.f, tx.reinit.B.msf)
  
  dat$attr$tx.status[tx.reinit] <- 1


  ## Other output
  dat$attr$cum.time.on.tx <- dat$attr$cum.time.on.tx +
                             ((dat$attr$tx.status == 1) %in% TRUE)
  dat$attr$cum.time.off.tx <- dat$attr$cum.time.off.tx +
                              ((dat$attr$tx.status == 0) %in% TRUE)

  ## Summary statistics
  dat$epi$tx.init.inc[at] <- length(tx.init)
  dat$epi$tx.halt.inc[at] <- length(tx.halt)
  dat$epi$tx.resm.inc[at] <- length(tx.reinit)
  dat$epi$tx.cov.diag[at] <-sum(dat$attr$tx.status, na.rm = TRUE) / sum(dat$attr$diag.status, na.rm = TRUE)
  dat$epi$tx.cov.PLHIV[at] <-sum(dat$attr$tx.status, na.rm = TRUE) / sum(dat$attr$status, na.rm = TRUE)
  
  poi<-which(age < 40)
  dat$epi$tx.cov.diag.poi[at] <-sum(dat$attr$tx.status[poi], na.rm = TRUE) / sum(dat$attr$diag.status[poi], na.rm = TRUE)
  dat$epi$tx.cov.PLHIV.poi[at] <-sum(dat$attr$tx.status[poi], na.rm = TRUE) / sum(dat$attr$status[poi], na.rm = TRUE)
  
  
  return(dat)
}


