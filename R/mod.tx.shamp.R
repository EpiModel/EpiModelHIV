
#' @title Treatment Module for up to 5 race groups heterosexuals and MSM.
#'
#' @description Module function for anti-retroviral treatment initiation and
#'              adherence over time.
#'
#' @inheritParams aging_msm
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
tx_shamp <- function(dat, at) {

  ## Variables

  # Attributes
  race <- dat$attr$race
  sex <- dat$attr$sex
  sex.ident <- dat$attr$sex.ident
  immig.loc <- dat$attr$immig.loc
  status <- dat$attr$status
  tx.status <- dat$attr$tx.status
  diag.status <- dat$attr$diag.status
  tt.traj <- dat$attr$tt.traj
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  stage <- dat$attr$stage

  # Parameters
  tx.init.B.f.prob <- dat$param$tx.init.B.f.prob
  tx.init.BI.f.prob <- dat$param$tx.init.BI.f.prob
  tx.init.H.f.prob <- dat$param$tx.init.H.f.prob
  tx.init.HI.f.prob <- dat$param$tx.init.HI.f.prob
  tx.init.W.f.prob <- dat$param$tx.init.W.f.prob
  tx.init.B.msf.prob <- dat$param$tx.init.B.msf.prob
  tx.init.BI.msf.prob <- dat$param$tx.init.BI.msf.prob
  tx.init.H.msf.prob <- dat$param$tx.init.H.msf.prob
  tx.init.HI.msf.prob <- dat$param$tx.init.HI.msf.prob
  tx.init.W.msf.prob <- dat$param$tx.init.W.msf.prob
  tx.init.B.msm.prob <- dat$param$tx.init.B.msm.prob
  tx.init.BI.msm.prob <- dat$param$tx.init.BI.msm.prob
  tx.init.H.msm.prob <- dat$param$tx.init.H.msm.prob
  tx.init.HI.msm.prob <- dat$param$tx.init.HI.msm.prob
  tx.init.W.msm.prob <- dat$param$tx.init.W.msm.prob 
  tx.init.B.msmf.prob <- dat$param$tx.init.B.msmf.prob
  tx.init.BI.msmf.prob <- dat$param$tx.init.BI.msmf.prob
  tx.init.H.msmf.prob <- dat$param$tx.init.H.msmf.prob
  tx.init.HI.msmf.prob <- dat$param$tx.init.HI.msmf.prob
  tx.init.W.msmf.prob <- dat$param$tx.init.W.msmf.prob 
  
  
  
  tx.halt.B.f.prob <- dat$param$tx.halt.B.f.prob
  tx.halt.BI.f.prob <- dat$param$tx.halt.BI.f.prob
  tx.halt.H.f.prob <- dat$param$tx.halt.H.f.prob
  tx.halt.HI.f.prob <- dat$param$tx.halt.HI.f.prob
  tx.halt.W.f.prob <- dat$param$tx.halt.W.f.prob
  tx.halt.B.msf.prob <- dat$param$tx.halt.B.msf.prob
  tx.halt.BI.msf.prob <- dat$param$tx.halt.BI.msf.prob
  tx.halt.H.msf.prob <- dat$param$tx.halt.H.msf.prob
  tx.halt.HI.msf.prob <- dat$param$tx.halt.HI.msf.prob
  tx.halt.W.msf.prob <- dat$param$tx.halt.W.msf.prob
  tx.halt.B.msm.prob <- dat$param$tx.halt.B.msm.prob
  tx.halt.BI.msm.prob <- dat$param$tx.halt.BI.msm.prob
  tx.halt.H.msm.prob <- dat$param$tx.halt.H.msm.prob
  tx.halt.HI.msm.prob <- dat$param$tx.halt.HI.msm.prob
  tx.halt.W.msm.prob <- dat$param$tx.halt.W.msm.prob
  tx.halt.B.msmf.prob <- dat$param$tx.halt.B.msmf.prob
  tx.halt.BI.msmf.prob <- dat$param$tx.halt.BI.msmf.prob
  tx.halt.H.msmf.prob <- dat$param$tx.halt.H.msmf.prob
  tx.halt.HI.msmf.prob <- dat$param$tx.halt.HI.msmf.prob
  tx.halt.W.msmf.prob <- dat$param$tx.halt.W.msmf.prob
  
  tx.reinit.B.f.prob <- dat$param$tx.reinit.B.f.prob
  tx.reinit.BI.f.prob <- dat$param$tx.reinit.BI.f.prob
  tx.reinit.H.f.prob <- dat$param$tx.reinit.H.f.prob
  tx.reinit.HI.f.prob <- dat$param$tx.reinit.HI.f.prob
  tx.reinit.W.f.prob <- dat$param$tx.reinit.W.f.prob
  tx.reinit.B.msf.prob <- dat$param$tx.reinit.B.msf.prob
  tx.reinit.BI.msf.prob <- dat$param$tx.reinit.BI.msf.prob
  tx.reinit.H.msf.prob <- dat$param$tx.reinit.H.msf.prob
  tx.reinit.HI.msf.prob <- dat$param$tx.reinit.HI.msf.prob
  tx.reinit.W.msf.prob <- dat$param$tx.reinit.W.msf.prob
  tx.reinit.B.msm.prob <- dat$param$tx.reinit.B.msm.prob
  tx.reinit.BI.msm.prob <- dat$param$tx.reinit.BI.msm.prob
  tx.reinit.H.msm.prob <- dat$param$tx.reinit.H.msm.prob
  tx.reinit.HI.msm.prob <- dat$param$tx.reinit.HI.msm.prob
  tx.reinit.W.msm.prob <- dat$param$tx.reinit.W.msm.prob
  tx.reinit.B.msmf.prob <- dat$param$tx.reinit.B.msmf.prob
  tx.reinit.BI.msmf.prob <- dat$param$tx.reinit.BI.msmf.prob
  tx.reinit.H.msmf.prob <- dat$param$tx.reinit.H.msmf.prob
  tx.reinit.HI.msmf.prob <- dat$param$tx.reinit.HI.msmf.prob
  tx.reinit.W.msmf.prob <- dat$param$tx.reinit.W.msmf.prob
  
  ## Initiation
  #Females
  tx.init.elig.B.f <- which(race == "B" & status == 1 & sex == "F" & immig.loc== 0 &
                          tx.status == 0 & diag.status == 1 &
                          tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                          stage != 4)
  tx.init.B.f <- tx.init.elig.B.f[rbinom(length(tx.init.elig.B.f), 1,
                                     tx.init.B.f.prob) == 1]
  
  tx.init.elig.BI.f <- which(race == "BI" & status == 1 & sex == "F" & immig.loc== 0 &
                            tx.status == 0 & diag.status == 1 &
                            tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                            stage != 4)
  tx.init.BI.f <- tx.init.elig.BI.f[rbinom(length(tx.init.elig.BI.f), 1,
                                     tx.init.BI.f.prob) == 1]
  
  tx.init.elig.H.f <- which(race == "H" & status == 1 & sex == "F" & immig.loc== 0 &
                            tx.status == 0 & diag.status == 1 &
                            tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                            stage != 4)
  tx.init.H.f <- tx.init.elig.H.f[rbinom(length(tx.init.elig.H.f), 1,
                                     tx.init.H.f.prob) == 1]
  
  tx.init.elig.HI.f <- which(race == "HI" & status == 1 & sex == "F" & immig.loc== 0 &
                             tx.status == 0 & diag.status == 1 &
                             tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                             stage != 4)
  tx.init.HI.f <- tx.init.elig.HI.f[rbinom(length(tx.init.elig.HI.f), 1,
                                       tx.init.HI.f.prob) == 1]

  tx.init.elig.W.f <- which(race == "W" & status == 1 & sex == "F" & immig.loc== 0 &
                          tx.status == 0 & diag.status == 1 &
                          tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                          stage != 4)
  tx.init.W.f <- tx.init.elig.W.f[rbinom(length(tx.init.elig.W.f), 1,
                                     tx.init.W.f.prob) == 1]

  #males HET msf
  
  tx.init.elig.B.msf <- which(race == "B" & status == 1 & sex == "M" & sex.ident=="msf" & immig.loc== 0 &
                            tx.status == 0 & diag.status == 1 &
                            tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                            stage != 4)
  tx.init.B.msf <- tx.init.elig.B.msf[rbinom(length(tx.init.elig.B.msf), 1,
                                     tx.init.B.msf.prob) == 1]
  
  tx.init.elig.BI.msf <- which(race == "BI" & status == 1 & sex == "M" & sex.ident=="msf" & immig.loc== 0 &
                             tx.status == 0 & diag.status == 1 &
                             tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                             stage != 4)
  tx.init.BI.msf <- tx.init.elig.BI.msf[rbinom(length(tx.init.elig.BI.msf), 1,
                                       tx.init.BI.msf.prob) == 1]
  
  tx.init.elig.H.msf <- which(race == "H" & status == 1 & sex == "M" & sex.ident=="msf" & immig.loc== 0 &
                            tx.status == 0 & diag.status == 1 &
                            tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                            stage != 4)
  tx.init.H.msf <- tx.init.elig.H.msf[rbinom(length(tx.init.elig.H.msf), 1,
                                     tx.init.H.msf.prob) == 1]
  
  tx.init.elig.HI.msf <- which(race == "HI" & status == 1 & sex == "M" & sex.ident=="msf" & immig.loc== 0 &
                             tx.status == 0 & diag.status == 1 &
                             tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                             stage != 4)
  tx.init.HI.msf <- tx.init.elig.HI.msf[rbinom(length(tx.init.elig.HI.msf), 1,
                                       tx.init.HI.msf.prob) == 1]
  
  tx.init.elig.W.msf <- which(race == "W" & status == 1 & sex == "M" & sex.ident=="msf" & immig.loc== 0 &
                            tx.status == 0 & diag.status == 1 &
                            tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                            stage != 4)
  tx.init.W.msf <- tx.init.elig.W.msf[rbinom(length(tx.init.elig.W.msf), 1,
                                     tx.init.W.msf.prob) == 1]
 
   #males MSM
  
  tx.init.elig.B.msm <- which(race == "B" & status == 1 & sex == "M" & sex.ident=="msm" & immig.loc== 0 &
                              tx.status == 0 & diag.status == 1 &
                              tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                              stage != 4)
  tx.init.B.msm <- tx.init.elig.B.msm[rbinom(length(tx.init.elig.B.msm), 1,
                                         tx.init.B.msm.prob) == 1]
  
  tx.init.elig.BI.msm <- which(race == "BI" & status == 1 & sex == "M" & sex.ident=="msm" & immig.loc== 0 &
                               tx.status == 0 & diag.status == 1 &
                               tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                               stage != 4)
  tx.init.BI.msm <- tx.init.elig.BI.msm[rbinom(length(tx.init.elig.BI.msm), 1,
                                           tx.init.BI.msm.prob) == 1]
  
  tx.init.elig.H.msm <- which(race == "H" & status == 1 & sex == "M" & sex.ident=="msm" & immig.loc== 0 &
                              tx.status == 0 & diag.status == 1 &
                              tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                              stage != 4)
  tx.init.H.msm <- tx.init.elig.H.msm[rbinom(length(tx.init.elig.H.msm), 1,
                                         tx.init.H.msm.prob) == 1]
  
  tx.init.elig.HI.msm <- which(race == "HI" & status == 1 & sex == "M" & sex.ident=="msm" & immig.loc== 0 &
                               tx.status == 0 & diag.status == 1 &
                               tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                               stage != 4)
  tx.init.HI.msm <- tx.init.elig.HI.msm[rbinom(length(tx.init.elig.HI.msm), 1,
                                           tx.init.HI.msm.prob) == 1]
  
  tx.init.elig.W.msm <- which(race == "W" & status == 1 & sex == "M" & sex.ident=="msm" & immig.loc== 0 &
                              tx.status == 0 & diag.status == 1 &
                              tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                              stage != 4)
  tx.init.W.msm <- tx.init.elig.W.msm[rbinom(length(tx.init.elig.W.msm), 1,
                                         tx.init.W.msm.prob) == 1]
  
  ##Males MSMF
  
  
  tx.init.elig.B.msmf <- which(race == "B" & status == 1 & sex == "M" & sex.ident=="msmf" & immig.loc== 0 &
                                 tx.status == 0 & diag.status == 1 &
                                 tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                                 stage != 4)
  tx.init.B.msmf <- tx.init.elig.B.msmf[rbinom(length(tx.init.elig.B.msmf), 1,
                                               tx.init.B.msmf.prob) == 1]
  
  tx.init.elig.BI.msmf <- which(race == "BI" & status == 1 & sex == "M" & sex.ident=="msmf" & immig.loc== 0 &
                                  tx.status == 0 & diag.status == 1 &
                                  tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                                  stage != 4)
  tx.init.BI.msmf <- tx.init.elig.BI.msmf[rbinom(length(tx.init.elig.BI.msmf), 1,
                                                 tx.init.BI.msmf.prob) == 1]
  
  tx.init.elig.H.msmf <- which(race == "H" & status == 1 & sex == "M" & sex.ident=="msmf" & immig.loc== 0 &
                                 tx.status == 0 & diag.status == 1 &
                                 tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                                 stage != 4)
  tx.init.H.msmf <- tx.init.elig.H.msmf[rbinom(length(tx.init.elig.H.msmf), 1,
                                               tx.init.H.msmf.prob) == 1]
  
  tx.init.elig.HI.msmf <- which(race == "HI" & status == 1 & sex == "M" & sex.ident=="msmf" & immig.loc== 0 &
                                  tx.status == 0 & diag.status == 1 &
                                  tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                                  stage != 4)
  tx.init.HI.msmf <- tx.init.elig.HI.msmf[rbinom(length(tx.init.elig.HI.msmf), 1,
                                                 tx.init.HI.msmf.prob) == 1]
  
  tx.init.elig.W.msmf <- which(race == "W" & status == 1 & sex == "M" & sex.ident=="msmf" & immig.loc== 0 &
                                 tx.status == 0 & diag.status == 1 &
                                 tt.traj %in% c(3, 4) & cum.time.on.tx == 0 &
                                 stage != 4)
  tx.init.W.msmf <- tx.init.elig.W.msmf[rbinom(length(tx.init.elig.W.msmf), 1,
                                               tx.init.W.msmf.prob) == 1]
  
  tx.init <- c(tx.init.B.f, tx.init.BI.f, tx.init.H.f, tx.init.HI.f, tx.init.W.f,
               tx.init.B.msf, tx.init.BI.msf, tx.init.H.msf, tx.init.HI.msf, tx.init.W.msf,
               tx.init.B.msm, tx.init.BI.msm, tx.init.H.msm, tx.init.HI.msm, tx.init.W.msm,
               tx.init.B.msmf, tx.init.BI.msmf, tx.init.H.msmf, tx.init.HI.msmf, tx.init.W.msmf)

  dat$attr$tx.status[tx.init] <- 1
  dat$attr$tx.init.time[tx.init] <- at


  ## Halting
  ##females
  tx.halt.elig.B.f <- which(race == "B" & tx.status == 1 & sex=="F")
  tx.halt.B.f <- tx.halt.elig.B.f[rbinom(length(tx.halt.elig.B.f), 1,
                                     tx.halt.B.f.prob) == 1]
  
  tx.halt.elig.BI.f <- which(race == "BI" & tx.status == 1 & sex=="F")
  tx.halt.BI.f <- tx.halt.elig.BI.f[rbinom(length(tx.halt.elig.BI.f), 1,
                                     tx.halt.BI.f.prob) == 1]

  tx.halt.elig.H.f <- which(race == "H" & tx.status == 1 & sex=="F")
  tx.halt.H.f <- tx.halt.elig.H.f[rbinom(length(tx.halt.elig.H.f), 1,
                                     tx.halt.H.f.prob) == 1]
  
  tx.halt.elig.HI.f <- which(race == "HI" & tx.status == 1 & sex=="F")
  tx.halt.HI.f <- tx.halt.elig.HI.f[rbinom(length(tx.halt.elig.HI.f), 1,
                                       tx.halt.HI.f.prob) == 1]
  
  tx.halt.elig.W.f <- which(race == "W" & tx.status == 1 & sex=="F")
  tx.halt.W.f <- tx.halt.elig.W.f[rbinom(length(tx.halt.elig.W.f),
                                     1, tx.halt.W.f.prob) == 1]
  
  
  ##males HET msf
  tx.halt.elig.B.msf <- which(race == "B" & tx.status == 1 & sex=="M" & sex.ident=="msf")
  tx.halt.B.msf <- tx.halt.elig.B.msf[rbinom(length(tx.halt.elig.B.msf), 1,
                                     tx.halt.B.msf.prob) == 1]
  
  tx.halt.elig.BI.msf <- which(race == "BI" & tx.status == 1 & sex=="M" & sex.ident=="msf")
  tx.halt.BI.msf <- tx.halt.elig.BI.msf[rbinom(length(tx.halt.elig.BI.msf), 1,
                                       tx.halt.BI.msf.prob) == 1]
  
  tx.halt.elig.H.msf <- which(race == "H" & tx.status == 1 & sex=="M" & sex.ident=="msf")
  tx.halt.H.msf <- tx.halt.elig.H.msf[rbinom(length(tx.halt.elig.H.msf), 1,
                                     tx.halt.H.msf.prob) == 1]
  
  tx.halt.elig.HI.msf <- which(race == "HI" & tx.status == 1 & sex=="M" & sex.ident=="msf")
  tx.halt.HI.msf <- tx.halt.elig.HI.msf[rbinom(length(tx.halt.elig.HI.msf), 1,
                                       tx.halt.HI.msf.prob) == 1]
  
  tx.halt.elig.W.msf <- which(race == "W" & tx.status == 1 & sex=="M" & sex.ident=="msf")
  tx.halt.W.msf <- tx.halt.elig.W.msf[rbinom(length(tx.halt.elig.W.msf),
                                     1, tx.halt.W.msf.prob) == 1]
  
  ##males MSM
  tx.halt.elig.B.msm <- which(race == "B" & tx.status == 1 & sex=="M" & sex.ident=="msm")
  tx.halt.B.msm <- tx.halt.elig.B.msm[rbinom(length(tx.halt.elig.B.msm), 1,
                                         tx.halt.B.msm.prob) == 1]
  
  tx.halt.elig.BI.msm <- which(race == "BI" & tx.status == 1 & sex=="M" & sex.ident=="msm")
  tx.halt.BI.msm <- tx.halt.elig.BI.msm[rbinom(length(tx.halt.elig.BI.msm), 1,
                                           tx.halt.BI.msm.prob) == 1]
  
  tx.halt.elig.H.msm <- which(race == "H" & tx.status == 1 & sex=="M" & sex.ident=="msm")
  tx.halt.H.msm <- tx.halt.elig.H.msm[rbinom(length(tx.halt.elig.H.msm), 1,
                                         tx.halt.H.msm.prob) == 1]
  
  tx.halt.elig.HI.msm <- which(race == "HI" & tx.status == 1 & sex=="M" & sex.ident=="msm")
  tx.halt.HI.msm <- tx.halt.elig.HI.msm[rbinom(length(tx.halt.elig.HI.msm), 1,
                                           tx.halt.HI.msm.prob) == 1]
  
  tx.halt.elig.W.msm <- which(race == "W" & tx.status == 1 & sex=="M" & sex.ident=="msm")
  tx.halt.W.msm <- tx.halt.elig.W.msm[rbinom(length(tx.halt.elig.W.msm),
                                         1, tx.halt.W.msm.prob) == 1]
  
  ##males MSMF
  tx.halt.elig.B.msmf <- which(race == "B" & tx.status == 1 & sex=="M" & sex.ident=="msmf")
  tx.halt.B.msmf <- tx.halt.elig.B.msmf[rbinom(length(tx.halt.elig.B.msmf), 1,
                                               tx.halt.B.msmf.prob) == 1]
  
  tx.halt.elig.BI.msmf <- which(race == "BI" & tx.status == 1 & sex=="M" & sex.ident=="msmf")
  tx.halt.BI.msmf <- tx.halt.elig.BI.msmf[rbinom(length(tx.halt.elig.BI.msmf), 1,
                                                 tx.halt.BI.msmf.prob) == 1]
  
  tx.halt.elig.H.msmf <- which(race == "H" & tx.status == 1 & sex=="M" & sex.ident=="msmf")
  tx.halt.H.msmf <- tx.halt.elig.H.msmf[rbinom(length(tx.halt.elig.H.msmf), 1,
                                               tx.halt.H.msmf.prob) == 1]
  
  tx.halt.elig.HI.msmf <- which(race == "HI" & tx.status == 1 & sex=="M" & sex.ident=="msmf")
  tx.halt.HI.msmf <- tx.halt.elig.HI.msmf[rbinom(length(tx.halt.elig.HI.msmf), 1,
                                                 tx.halt.HI.msmf.prob) == 1]
  
  tx.halt.elig.W.msmf <- which(race == "W" & tx.status == 1 & sex=="M" & sex.ident=="msmf")
  tx.halt.W.msmf <- tx.halt.elig.W.msmf[rbinom(length(tx.halt.elig.W.msmf),
                                               1, tx.halt.W.msmf.prob) == 1]
  
  
  tx.halt <- c(tx.halt.B.f, tx.halt.BI.f, tx.halt.H.f, tx.halt.HI.f, tx.halt.W.f,
               tx.halt.B.msf, tx.halt.BI.msf, tx.halt.H.msf, tx.halt.HI.msf, tx.halt.W.msf,
               tx.halt.B.msm, tx.halt.BI.msm, tx.halt.H.msm, tx.halt.HI.msm, tx.halt.W.msm,
               tx.halt.B.msmf, tx.halt.BI.msmf, tx.halt.H.msmf, tx.halt.HI.msmf, tx.halt.W.msmf)
  
  dat$attr$tx.status[tx.halt] <- 0


  ## Restarting
  #females
  tx.reinit.elig.B.f <- which(race == "B" & tx.status == 0 & sex=="F" & immig.loc==0 &
                            cum.time.on.tx > 0 & stage != 4)
  tx.reinit.B.f <- tx.reinit.elig.B.f[rbinom(length(tx.reinit.elig.B.f),
                                         1, tx.reinit.B.f.prob) == 1]
  
  tx.reinit.elig.BI.f <- which(race == "BI" & tx.status == 0 & sex=="F" & immig.loc==0 &
                              cum.time.on.tx > 0 & stage != 4)
  tx.reinit.BI.f <- tx.reinit.elig.BI.f[rbinom(length(tx.reinit.elig.BI.f),
                                         1, tx.reinit.BI.f.prob) == 1]
  
  tx.reinit.elig.H.f <- which(race == "H" & tx.status == 0 & sex=="F" & immig.loc==0 &
                              cum.time.on.tx > 0 & stage != 4)
  tx.reinit.H.f <- tx.reinit.elig.H.f[rbinom(length(tx.reinit.elig.H.f),
                                         1, tx.reinit.H.f.prob) == 1]
  
  tx.reinit.elig.HI.f <- which(race == "HI" & tx.status == 0 & sex=="F" & immig.loc==0 &
                              cum.time.on.tx > 0 & stage != 4)
  tx.reinit.HI.f <- tx.reinit.elig.HI.f[rbinom(length(tx.reinit.elig.HI.f),
                                         1, tx.reinit.HI.f.prob) == 1]

  tx.reinit.elig.W.f <- which(race == "W" & tx.status == 0 & sex=="F" & immig.loc==0 &
                            cum.time.on.tx > 0 & stage != 4)
  tx.reinit.W.f <- tx.reinit.elig.W.f[rbinom(length(tx.reinit.elig.W.f),
                                         1, tx.reinit.W.f.prob) == 1]
  
  #males HET msf
  tx.reinit.elig.B.msf <- which(race == "B" & tx.status == 0 & sex=="M" & sex.ident=="msf" & immig.loc==0 &
                              cum.time.on.tx > 0 & stage != 4)
  tx.reinit.B.msf <- tx.reinit.elig.B.msf[rbinom(length(tx.reinit.elig.B.msf),
                                         1, tx.reinit.B.msf.prob) == 1]
  
  tx.reinit.elig.BI.msf <- which(race == "BI" & tx.status == 0 & sex=="M" & sex.ident=="msf" & immig.loc==0 &
                               cum.time.on.tx > 0 & stage != 4)
  tx.reinit.BI.msf <- tx.reinit.elig.BI.msf[rbinom(length(tx.reinit.elig.BI.msf),
                                           1, tx.reinit.BI.msf.prob) == 1]
  
  tx.reinit.elig.H.msf <- which(race == "H" & tx.status == 0 & sex=="M" & sex.ident=="msf" & immig.loc==0 &
                              cum.time.on.tx > 0 & stage != 4)
  tx.reinit.H.msf <- tx.reinit.elig.H.msf[rbinom(length(tx.reinit.elig.H.msf),
                                         1, tx.reinit.H.msf.prob) == 1]
  
  tx.reinit.elig.HI.msf <- which(race == "HI" & tx.status == 0 & sex=="M" & sex.ident=="msf" & immig.loc==0 &
                               cum.time.on.tx > 0 & stage != 4)
  tx.reinit.HI.msf <- tx.reinit.elig.HI.msf[rbinom(length(tx.reinit.elig.HI.msf),
                                           1, tx.reinit.HI.msf.prob) == 1]
  
  tx.reinit.elig.W.msf <- which(race == "W" & tx.status == 0 & sex=="M" & sex.ident=="msf" & immig.loc==0 &
                              cum.time.on.tx > 0 & stage != 4)
  tx.reinit.W.msf <- tx.reinit.elig.W.msf[rbinom(length(tx.reinit.elig.W.msf),
                                         1, tx.reinit.W.msf.prob) == 1]
  
  #males MSM
  tx.reinit.elig.B.msm <- which(race == "B" & tx.status == 0 & sex=="M" & sex.ident=="msm" &  immig.loc==0 &
                                cum.time.on.tx > 0 & stage != 4)
  tx.reinit.B.msm <- tx.reinit.elig.B.msm[rbinom(length(tx.reinit.elig.B.msm),
                                             1, tx.reinit.B.msm.prob) == 1]
  
  tx.reinit.elig.BI.msm <- which(race == "BI" & tx.status == 0 & sex=="M" & sex.ident=="msm" & immig.loc==0 &
                                 cum.time.on.tx > 0 & stage != 4)
  tx.reinit.BI.msm <- tx.reinit.elig.BI.msm[rbinom(length(tx.reinit.elig.BI.msm),
                                               1, tx.reinit.BI.msm.prob) == 1]
  
  tx.reinit.elig.H.msm <- which(race == "H" & tx.status == 0 & sex=="M" & sex.ident=="msm" & immig.loc==0 &
                                cum.time.on.tx > 0 & stage != 4)
  tx.reinit.H.msm <- tx.reinit.elig.H.msm[rbinom(length(tx.reinit.elig.H.msm),
                                             1, tx.reinit.H.msm.prob) == 1]
  
  tx.reinit.elig.HI.msm <- which(race == "HI" & tx.status == 0 & sex=="M" & sex.ident=="msm" & immig.loc==0 &
                                 cum.time.on.tx > 0 & stage != 4)
  tx.reinit.HI.msm <- tx.reinit.elig.HI.msm[rbinom(length(tx.reinit.elig.HI.msm),
                                               1, tx.reinit.HI.msm.prob) == 1]
  
  tx.reinit.elig.W.msm <- which(race == "W" & tx.status == 0 & sex=="M" & sex.ident=="msm" & immig.loc==0 &
                                cum.time.on.tx > 0 & stage != 4)
  tx.reinit.W.msm <- tx.reinit.elig.W.msm[rbinom(length(tx.reinit.elig.W.msm),
                                             1, tx.reinit.W.msm.prob) == 1]
 
   #males MSMF
  tx.reinit.elig.B.msmf <- which(race == "B" & tx.status == 0 & sex=="M" & sex.ident=="msmf" &  immig.loc==0 &
                                   cum.time.on.tx > 0 & stage != 4)
  tx.reinit.B.msmf <- tx.reinit.elig.B.msmf[rbinom(length(tx.reinit.elig.B.msmf),
                                                   1, tx.reinit.B.msmf.prob) == 1]
  
  tx.reinit.elig.BI.msmf <- which(race == "BI" & tx.status == 0 & sex=="M" & sex.ident=="msmf" & immig.loc==0 &
                                    cum.time.on.tx > 0 & stage != 4)
  tx.reinit.BI.msmf <- tx.reinit.elig.BI.msmf[rbinom(length(tx.reinit.elig.BI.msmf),
                                                     1, tx.reinit.BI.msmf.prob) == 1]
  
  tx.reinit.elig.H.msmf <- which(race == "H" & tx.status == 0 & sex=="M" & sex.ident=="msmf" & immig.loc==0 &
                                   cum.time.on.tx > 0 & stage != 4)
  tx.reinit.H.msmf <- tx.reinit.elig.H.msmf[rbinom(length(tx.reinit.elig.H.msmf),
                                                   1, tx.reinit.H.msmf.prob) == 1]
  
  tx.reinit.elig.HI.msmf <- which(race == "HI" & tx.status == 0 & sex=="M" & sex.ident=="msmf" & immig.loc==0 &
                                    cum.time.on.tx > 0 & stage != 4)
  tx.reinit.HI.msmf <- tx.reinit.elig.HI.msmf[rbinom(length(tx.reinit.elig.HI.msmf),
                                                     1, tx.reinit.HI.msmf.prob) == 1]
  
  tx.reinit.elig.W.msmf <- which(race == "W" & tx.status == 0 & sex=="M" & sex.ident=="msmf" & immig.loc==0 &
                                   cum.time.on.tx > 0 & stage != 4)
  tx.reinit.W.msmf <- tx.reinit.elig.W.msmf[rbinom(length(tx.reinit.elig.W.msmf),
                                                   1, tx.reinit.W.msmf.prob) == 1]
  

  tx.reinit <- c(tx.reinit.B.f, tx.reinit.BI.f, tx.reinit.H.f, tx.reinit.HI.f, tx.reinit.W.f,
                 tx.reinit.B.msf, tx.reinit.BI.msf, tx.reinit.H.msf, tx.reinit.HI.msf, tx.reinit.W.msf,
                 tx.reinit.B.msm, tx.reinit.BI.msm, tx.reinit.H.msm, tx.reinit.HI.msm, tx.reinit.W.msm,
                 tx.reinit.B.msmf, tx.reinit.BI.msmf, tx.reinit.H.msmf, tx.reinit.HI.msmf, tx.reinit.W.msmf)
  
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

  return(dat)
}


