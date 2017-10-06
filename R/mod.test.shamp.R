
#' @title HIV Testing Module for up to 5 race groups heterosexuals and MSM.
#'
#' @description Module function for HIV diagnostic testing of infected persons.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This testing module supports two testing parameterizations, input via the
#' \code{testing.pattern} parameter: memoryless for stochastic and
#' geometrically-distributed waiting times to test (constant hazard); and interval
#' for deterministic tested after defined waiting time intervals.
#' If immigration is used those with dat$attr$immig.loc=1 do not test.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{last.neg.test},
#' \code{diag.status} and \code{diag.time} attributes.
#'
#' @keywords module shamp msm het
#'
#' @export
#'
test_shamp <- function(dat, at) {

  ## Variables

  # Attributes
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  sex <- dat$attr$sex
  race <- dat$attr$race
  sex.ident<- dat$attr$sex.ident
  tt.traj <- dat$attr$tt.traj
  status <- dat$attr$status
  inf.time <- dat$attr$inf.time
  immig.loc <- dat$attr$immig.loc

  prepStat <- dat$attr$prepStat
  prep.tst.int <- dat$param$prep.tst.int

  # Parameters
  testing.pattern <- dat$param$testing.pattern
  mean.test.B.f.int <- dat$param$mean.test.B.f.int
  mean.test.BI.f.int <- dat$param$mean.test.BI.f.int
  mean.test.H.f.int <- dat$param$mean.test.H.f.int
  mean.test.HI.f.int <- dat$param$mean.test.HI.f.int
  mean.test.W.f.int <- dat$param$mean.test.W.f.int
  mean.test.B.msf.int <- dat$param$mean.test.B.msf.int
  mean.test.BI.msf.int <- dat$param$mean.test.BI.msf.int
  mean.test.H.msf.int <- dat$param$mean.test.H.msf.int
  mean.test.HI.msf.int <- dat$param$mean.test.HI.msf.int
  mean.test.W.msf.int <- dat$param$mean.test.W.msf.int
  mean.test.B.msm.int <- dat$param$mean.test.B.msm.int
  mean.test.BI.msm.int <- dat$param$mean.test.BI.msm.int
  mean.test.H.msm.int <- dat$param$mean.test.H.msm.int
  mean.test.HI.msm.int <- dat$param$mean.test.HI.msm.int
  mean.test.W.msm.int <- dat$param$mean.test.W.msmf.int
  mean.test.B.msmf.int <- dat$param$mean.test.B.msmf.int
  mean.test.BI.msmf.int <- dat$param$mean.test.BI.msmf.int
  mean.test.H.msmf.int <- dat$param$mean.test.H.msmf.int
  mean.test.HI.msmf.int <- dat$param$mean.test.HI.msmf.int
  mean.test.W.msmf.int <- dat$param$mean.test.W.msmf.int

  twind.int <- dat$param$test.window.int

  tsincelntst <- at - dat$attr$last.neg.test
  tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]

  ## Process

  if (testing.pattern == "memoryless") {
    
#females
    elig.B.f <- which(race == "B" & sex == "F" & immig.loc == 0 &
                    tt.traj != 1 &
                    (diag.status == 0 | is.na(diag.status)) &
                    prepStat == 0)
    rates.B.f <- rep(1/mean.test.B.f.int, length(elig.B.f))
    tst.B.f <- elig.B.f[rbinom(length(elig.B.f), 1, rates.B.f) == 1]
    
    elig.BI.f <- which(race == "BI" & sex == "F" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        prepStat == 0)
    rates.BI.f <- rep(1/mean.test.BI.f.int, length(elig.BI.f))
    tst.BI.f <- elig.BI.f[rbinom(length(elig.BI.f), 1, rates.BI.f) == 1]
    

    elig.H.f <- which(race == "H" & sex == "F" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        prepStat == 0)
    rates.H.f <- rep(1/mean.test.H.f.int, length(elig.H.f))
    tst.H.f <- elig.H.f[rbinom(length(elig.H.f), 1, rates.H.f) == 1]
    
    elig.HI.f <- which(race == "HI" & sex == "F" & immig.loc == 0 &
                         tt.traj != 1 &
                         (diag.status == 0 | is.na(diag.status)) &
                         prepStat == 0)
    rates.HI.f <- rep(1/mean.test.HI.f.int, length(elig.HI.f))
    tst.HI.f <- elig.HI.f[rbinom(length(elig.HI.f), 1, rates.HI.f) == 1]    

    elig.W.f <- which(race == "W" & sex == "F" & immig.loc == 0 &
                    tt.traj != 1 &
                    (diag.status == 0 | is.na(diag.status)) &
                    prepStat == 0)
    rates.W.f <- rep(1/mean.test.W.f.int, length(elig.W.f))
    tst.W.f <- elig.W.f[rbinom(length(elig.W.f), 1, rates.W.f) == 1]
    
    #males_het MSF
    elig.B.msf <- which(race == "B" & sex.ident == "msf" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        prepStat == 0)
    rates.B.msf <- rep(1/mean.test.B.msf.int, length(elig.B.msf))
    tst.B.msf <- elig.B.msf[rbinom(length(elig.B.msf), 1, rates.B.msf) == 1]
    
    elig.BI.msf <- which(race == "BI" & sex == "M" & sex.ident == "msf" & immig.loc == 0 &
                         tt.traj != 1 &
                         (diag.status == 0 | is.na(diag.status)) &
                         prepStat == 0)
    rates.BI.msf <- rep(1/mean.test.BI.msf.int, length(elig.BI.msf))
    tst.BI.msf <- elig.BI.msf[rbinom(length(elig.BI.msf), 1, rates.BI.msf) == 1]
    
    
    elig.H.msf <- which(race == "H" & sex == "M" & sex.ident == "msf" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        prepStat == 0)
    rates.H.msf <- rep(1/mean.test.H.msf.int, length(elig.H.msf))
    tst.H.msf <- elig.H.msf[rbinom(length(elig.H.msf), 1, rates.H.msf) == 1]
    
    elig.HI.msf <- which(race == "HI" & sex == "M" & sex.ident == "msf" & immig.loc == 0 &
                         tt.traj != 1 &
                         (diag.status == 0 | is.na(diag.status)) &
                         prepStat == 0)
    rates.HI.msf <- rep(1/mean.test.HI.msf.int, length(elig.HI.msf))
    tst.HI.msf <- elig.HI.msf[rbinom(length(elig.HI.msf), 1, rates.HI.msf) == 1]    
    
    elig.W.msf <- which(race == "W" & sex == "M" & sex.ident == "msf" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        prepStat == 0)
    rates.W.msf <- rep(1/mean.test.W.msf.int, length(elig.W.msf))
    tst.W.msf <- elig.W.msf[rbinom(length(elig.W.msf), 1, rates.W.msf) == 1]
    
    #males_msm
    
    
    elig.B.msm <- which(race == "B" & sex == "M" & sex.ident == "msm" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        prepStat == 0)
    rates.B.msm <- rep(1/mean.test.B.msm.int, length(elig.B.msm))
    tst.B.msm <- elig.B.msm[rbinom(length(elig.B.msm), 1, rates.B.msm) == 1]
    
    elig.BI.msm <- which(race == "BI" & sex == "M" & sex.ident == "msm" & immig.loc == 0 &
                         tt.traj != 1 &
                         (diag.status == 0 | is.na(diag.status)) &
                         prepStat == 0)
    rates.BI.msm <- rep(1/mean.test.BI.msm.int, length(elig.BI.msm))
    tst.BI.msm <- elig.BI.msm[rbinom(length(elig.BI.msm), 1, rates.BI.msm) == 1]
    
    
    elig.H.msm <- which(race == "H" & sex == "M" & sex.ident == "msm" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        prepStat == 0)
    rates.H.msm <- rep(1/mean.test.H.msm.int, length(elig.H.msm))
    tst.H.msm <- elig.H.msm[rbinom(length(elig.H.msm), 1, rates.H.msm) == 1]
    
    elig.HI.msm <- which(race == "HI" & sex == "M" & sex.ident == "msm" & immig.loc == 0 &
                         tt.traj != 1 &
                         (diag.status == 0 | is.na(diag.status)) &
                         prepStat == 0)
    rates.HI.msm <- rep(1/mean.test.HI.msm.int, length(elig.HI.msm))
    tst.HI.msm <- elig.HI.msm[rbinom(length(elig.HI.msm), 1, rates.HI.msm) == 1]    
    
    elig.W.msm <- which(race == "W" & sex == "M" & sex.ident == "msm" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        prepStat == 0)
    rates.W.msm <- rep(1/mean.test.W.msm.int, length(elig.W.msm))
    tst.W.msm <- elig.W.msm[rbinom(length(elig.W.msm), 1, rates.W.msm) == 1]
    
    #males_msmf
    
    
    elig.B.msmf <- which(race == "B" & sex == "M" & sex.ident == "msmf" & immig.loc == 0 &
                          tt.traj != 1 &
                          (diag.status == 0 | is.na(diag.status)) &
                          prepStat == 0)
    rates.B.msmf <- rep(1/mean.test.B.msmf.int, length(elig.B.msmf))
    tst.B.msmf <- elig.B.msmf[rbinom(length(elig.B.msmf), 1, rates.B.msmf) == 1]
    
    elig.BI.msmf <- which(race == "BI" & sex == "M" & sex.ident == "msmf" & immig.loc == 0 &
                           tt.traj != 1 &
                           (diag.status == 0 | is.na(diag.status)) &
                           prepStat == 0)
    rates.BI.msmf <- rep(1/mean.test.BI.msmf.int, length(elig.BI.msmf))
    tst.BI.msmf <- elig.BI.msmf[rbinom(length(elig.BI.msmf), 1, rates.BI.msmf) == 1]
    
    
    elig.H.msmf <- which(race == "H" & sex == "M" & sex.ident == "msmf" & immig.loc == 0 &
                          tt.traj != 1 &
                          (diag.status == 0 | is.na(diag.status)) &
                          prepStat == 0)
    rates.H.msmf <- rep(1/mean.test.H.msmf.int, length(elig.H.msmf))
    tst.H.msmf <- elig.H.msmf[rbinom(length(elig.H.msmf), 1, rates.H.msmf) == 1]
    
    elig.HI.msmf <- which(race == "HI" & sex == "M" & sex.ident == "msmf" & immig.loc == 0 &
                           tt.traj != 1 &
                           (diag.status == 0 | is.na(diag.status)) &
                           prepStat == 0)
    rates.HI.msmf <- rep(1/mean.test.HI.msmf.int, length(elig.HI.msmf))
    tst.HI.msmf <- elig.HI.msmf[rbinom(length(elig.HI.msmf), 1, rates.HI.msmf) == 1]    
    
    elig.W.msmf <- which(race == "W" & sex == "M" & sex.ident == "msmf" & immig.loc == 0 &
                          tt.traj != 1 &
                          (diag.status == 0 | is.na(diag.status)) &
                          prepStat == 0)
    rates.W.msmf <- rep(1/mean.test.W.msmf.int, length(elig.W.msmf))
    tst.W.msmf <- elig.W.msmf[rbinom(length(elig.W.msmf), 1, rates.W.msmf) == 1]
    
    tst.nprep <- c(tst.B.f, tst.BI.f, tst.H.f, tst.HI.f, tst.W.f, 
                   tst.B.msf, tst.BI.msf, tst.H.msf, tst.HI.msf, tst.W.msf,
                   tst.B.msm, tst.BI.msm, tst.H.msm, tst.HI.msm, tst.W.msm,
                   tst.B.msmf, tst.BI.msmf, tst.H.msmf, tst.HI.msmf, tst.W.msmf)
  }


  if (testing.pattern == "interval") {
    
    #females
    tst.B.f <- which(race == "B" & sex=="F" & immig.loc == 0 &
                   tt.traj != 1 &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= 2*(mean.test.B.f.int) &
                   prepStat == 0)
    
    tst.BI.f <- which(race == "BI" & sex=="F" & immig.loc == 0 &
                       tt.traj != 1 &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= 2*(mean.test.BI.f.int) &
                       prepStat == 0)
    
    tst.H.f <- which(race == "H" &  sex=="F" & immig.loc == 0 &
                       tt.traj != 1 &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= 2*(mean.test.H.f.int) &
                       prepStat == 0)
    
    tst.HI.f <- which(race == "HI" & sex=="F" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        tsincelntst >= 2*(mean.test.HI.f.int) &
                        prepStat == 0)

    tst.W.f <- which(race == "W" & sex=="F" & immig.loc == 0 &
                   tt.traj != 1 &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= 2*(mean.test.W.f.int) &
                   prepStat == 0)
    ##males Het msf
    
    tst.B.msf <- which(race == "B" &  sex == "M" & sex.ident == "msf" & immig.loc == 0 &
                       tt.traj != 1 &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= 2*(mean.test.B.msf.int) &
                       prepStat == 0)
    
    tst.BI.msf <- which(race == "BI" &  sex == "M" & sex.ident == "msf" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        tsincelntst >= 2*(mean.test.BI.msf.int) &
                        prepStat == 0)
    
    tst.H.msf <- which(race == "H" &  sex == "M" & sex.ident == "msf" & immig.loc == 0 &
                       tt.traj != 1 &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= 2*(mean.test.H.msf.int) &
                       prepStat == 0)
    
    tst.HI.msf <- which(race == "HI" &  sex == "M" & sex.ident == "msf" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        tsincelntst >= 2*(mean.test.HI.msf.int) &
                        prepStat == 0)
    
    tst.W.msf <- which(race == "W" &  sex == "M" & sex.ident == "msf" & immig.loc == 0 &
                       tt.traj != 1 &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= 2*(mean.test.W.msf.int) &
                       prepStat == 0)
    
    ##males MSM
    
    tst.B.msm <- which(race == "B" &  sex == "M" & sex.ident == "msm" & immig.loc == 0 &
                       tt.traj != 1 &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= 2*(mean.test.B.msm.int) &
                       prepStat == 0)
    
    tst.BI.msm <- which(race == "BI" &  sex == "M" & sex.ident == "msm" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        tsincelntst >= 2*(mean.test.BI.msm.int) &
                        prepStat == 0)
    
    tst.H.msm <- which(race == "H" &  sex == "M" & sex.ident == "msm" & immig.loc == 0 &
                       tt.traj != 1 &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= 2*(mean.test.H.msm.int) &
                       prepStat == 0)
    
    tst.HI.msm <- which(race == "HI" &  sex == "M" & sex.ident == "msm" & immig.loc == 0 &
                        tt.traj != 1 &
                        (diag.status == 0 | is.na(diag.status)) &
                        tsincelntst >= 2*(mean.test.HI.msm.int) &
                        prepStat == 0)
    
    tst.W.msm <- which(race == "W" &  sex == "M" & sex.ident == "msm" & immig.loc == 0 &
                       tt.traj != 1 &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= 2*(mean.test.W.msm.int) &
                       prepStat == 0)
    
    ##males MSM
    
    tst.B.msmf <- which(race == "B" &  sex == "M" & sex.ident == "msmf" & immig.loc == 0 &
                         tt.traj != 1 &
                         (diag.status == 0 | is.na(diag.status)) &
                         tsincelntst >= 2*(mean.test.B.msmf.int) &
                         prepStat == 0)
    
    tst.BI.msmf <- which(race == "BI" &  sex == "M" & sex.ident == "msmf" & immig.loc == 0 &
                          tt.traj != 1 &
                          (diag.status == 0 | is.na(diag.status)) &
                          tsincelntst >= 2*(mean.test.BI.msmf.int) &
                          prepStat == 0)
    
    tst.H.msmf <- which(race == "H" &  sex == "M" & sex.ident == "msmf" & immig.loc == 0 &
                         tt.traj != 1 &
                         (diag.status == 0 | is.na(diag.status)) &
                         tsincelntst >= 2*(mean.test.H.msmf.int) &
                         prepStat == 0)
    
    tst.HI.msmf <- which(race == "HI" &  sex == "M" & sex.ident == "msmf" & immig.loc == 0 &
                          tt.traj != 1 &
                          (diag.status == 0 | is.na(diag.status)) &
                          tsincelntst >= 2*(mean.test.HI.msmf.int) &
                          prepStat == 0)
    
    tst.W.msmf <- which(race == "W" &  sex == "M" & sex.ident == "msmf" & immig.loc == 0 &
                         tt.traj != 1 &
                         (diag.status == 0 | is.na(diag.status)) &
                         tsincelntst >= 2*(mean.test.W.msmf.int) &
                         prepStat == 0)
    
    tst.nprep <- c(tst.B.f, tst.BI.f, tst.H.f, tst.HI.f, tst.W.f,
                   tst.B.msf, tst.BI.msf, tst.H.msf, tst.HI.msf, tst.W.msf,
                   tst.B.msm, tst.BI.msm, tst.H.msm, tst.HI.msm, tst.W.msm,
                   tst.B.msmf, tst.BI.msmf, tst.H.msmf, tst.HI.msmf, tst.W.msmf)
  }

  # PrEP testing
  tst.prep <- which((diag.status == 0 | is.na(diag.status)) &
                    prepStat == 1 & immig.loc == 0 &
                    tsincelntst >= prep.tst.int)

  tst.all <- c(tst.nprep, tst.prep)

  tst.pos <- tst.all[status[tst.all] == 1 & inf.time[tst.all] <= at - twind.int]
  tst.neg <- setdiff(tst.all, tst.pos)

  # Attributes
  dat$attr$last.neg.test[tst.neg] <- at
  dat$attr$diag.status[tst.pos] <- 1
  dat$attr$diag.time[tst.pos] <- at
  dat$attr$evertest[tst.all] <- 1

  return(dat)
}

