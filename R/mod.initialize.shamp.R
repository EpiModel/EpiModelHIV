
# SHAMP ego -----------------------------------------------------------------

#' @title Initialization Module for up to 5 race groups heterosexuals and MSM.
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, initial state of the network is simulated using ergm.ego, and
#'              simulates disease status and other attributes are initialized.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param_shamp}}.
#' @param init An \code{EpiModel} object of class \code{\link{init_shamp}}.
#' @param control An \code{EpiModel} object of class \code{\link{control_shamp}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated \code{dat} object with the initialized values
#' for demographics and disease-related variables.
#'
#' @export
#' @keywords module HET MSM ego 
#'
initialize_shamp <- function(x, param, init, control, s) {

  # Master data list
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control

  dat$attr <- list()
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()
  dat$epi <- list()

  ## Network simulation ##
  
  nw <- list()
  for (i in 1:3) {

   # nw[[i]] <- simulate(x[[i]]$fit, popsize=sim.size, 
    #                control = control.simulate.ergm.ego(simulate.control = control.simulate(MCMC.burnin = 1e6)))
    
    nw[[i]] <- simulate(x[[i]]$fit)
   }


  ## ergm_prep here
  dat$el <- list()
  dat$p <- list()
  for (i in 1:2) {
    dat$el[[i]] <- as.edgelist(nw[[i]])
    attributes(dat$el[[i]])$vnames <- NULL
    p <- tergmLite::stergm_prep(nw[[i]], x[[i]]$formation, x[[i]]$coef.diss$dissolution,
                                x[[i]]$coef.form, x[[i]]$coef.diss$coef.adj, x[[i]]$constraints)
    p$model.form$formula <- NULL
    p$model.diss$formula <- NULL
    dat$p[[i]] <- p
  }
  dat$el[[3]] <- as.edgelist(nw[[3]])
  attributes(dat$el[[3]])$vnames <- NULL
  p <- tergmLite::ergm_prep(nw[[3]], x[[3]]$formation, x[[3]]$coef.form, x[[3]]$constraints)
  p$model.form$formula <- NULL
  dat$p[[3]] <- p
  

  # Network parameters
  dat$nwparam <- list()
  for (i in 1:3) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
    dat$nwparam[[i]]$formation[2]<-dat$nwparam[[i]]$formation[3]
    dat$nwparam[[i]]$formation[3]<-NULL
  }


  ## Nodal attributes ##

  # Degree terms
  dat$attr$deg.pers <- get.vertex.attribute(x[[1]]$fit$network, "deg.pers")
  dat$attr$deg.cohab <- get.vertex.attribute(x[[2]]$fit$network, "deg.cohab")
  
  dat$attr$deg.pers.c <- get.vertex.attribute(x[[1]]$fit$network, "deg.pers.c")
  dat$attr$deg.cohab.c <- get.vertex.attribute(x[[2]]$fit$network, "deg.cohab.c")
  

  # Race
  dat$attr$race <- get.vertex.attribute(nw[[2]], "race")
  num.B <-sum(dat$attr$race == "B")
  num.BI <-sum(dat$attr$race == "BI")
  num.H <-sum(dat$attr$race == "H")
  num.HI <-sum(dat$attr$race == "HI")
  num.W <-sum(dat$attr$race == "W")
  
  num <-num.B + num.BI + num.H + num.HI + num.W
  ids.B <- which(dat$attr$race == "B")
  ids.BI <- which(dat$attr$race == "BI")
  ids.H <- which(dat$attr$race == "H")
  ids.HI <- which(dat$attr$race == "HI")
  ids.W <- which(dat$attr$race == "W")
  
  # Sex
  dat$attr$sex <- get.vertex.attribute(nw[[1]], "sex")
  num.male <-sum(dat$attr$sex == "M")
  num.female <-sum(dat$attr$sex == "F")
  ids.M<-which(dat$attr$sex=="M")
  ids.F<-which(dat$attr$sex=="F")
  
  #sex by race
  ids.B.f<-na.omit(ids.F[ids.B])
  ids.BI.f<-na.omit(ids.F[ids.BI])
  ids.H.f<-na.omit(ids.F[ids.H])
  ids.HI.f<-na.omit(ids.F[ids.HI])
  ids.W.f<-na.omit(ids.F[ids.W])
  ids.B.m<-na.omit(ids.M[ids.B])
  ids.BI.m<-na.omit(ids.M[ids.BI])
  ids.H.m<-na.omit(ids.M[ids.H])
  ids.HI.m<-na.omit(ids.M[ids.HI])
  ids.W.m<-na.omit(ids.M[ids.W])
  
  #race.sex.cohab.
  dat$attr$race.sex.cohab<-get.vertex.attribute(nw[[1]], "race.sex.cohab")
  
   #race.sex.pers.
  dat$attr$race.sex.pers<-get.vertex.attribute(nw[[1]], "race.sex.pers")
  
  
  # Sex Identity
  dat$attr$sex.ident <- get.vertex.attribute(nw[[1]], "sex.ident")
  num.msf <-sum(dat$attr$sex.ident == "msf")
  num.msm <-sum(dat$attr$sex.ident == "msm")
  num.msmf <-sum(dat$attr$sex.ident == "msmf")
  num.f <-sum(dat$attr$sex.ident == "f")
  
  dat$attr$active <- rep(1, num)
  dat$attr$uid <- 1:num
  dat$temp$max.uid <- num

  # Age
  dat$attr$age <- get.vertex.attribute(nw[[1]], "age")
  partial<-(0:51)* (time.unit / 365)
  partial<-sample(partial,length(dat$attr$age),replace=TRUE)
  
  dat$attr$age<-dat$attr$age+partial
  
  dat$attr$sqrt.age <- get.vertex.attribute(nw[[1]], "sqrt.age")
  
  dat$attr$agecat <- get.vertex.attribute(nw[[1]], "agecat")
  dat$attr$sqrt.age.adj <- get.vertex.attribute(nw[[1]], "sqrt.age.adj")
  
  #demographic catagory.
  dat$attr$demog.cat<-rep(NA,length(dat$attr$sex))
  sex.groups<-sort(unique(dat$attr$sex))
  for (i in 1:(length(sex.groups))){
    dat$attr$demog.cat<-ifelse(dat$attr$sex==sex.groups[i],i*1000,dat$attr$demog.cat)      
  }
  
  race.groups<-sort(unique(dat$attr$race))
  for (i in 1:(length(race.groups))){
    dat$attr$demog.cat<-ifelse(dat$attr$race==race.groups[i],dat$attr$demog.cat+(i*100),dat$attr$demog.cat)      
  }
  
  age.groups<-sort(unique(floor(dat$attr$age)))
  age.temp<-floor(dat$attr$age)
  for (i in 1:(length(age.groups))){   
    dat$attr$demog.cat<-ifelse (age.temp==age.groups[i],dat$attr$demog.cat+(age.groups[i]),dat$attr$demog.cat)      
  }
  
  
  # Risk group
  dat$attr$riskg <- get.vertex.attribute(nw[[1]], "riskg")
  
  # Immigrant status
  dat$attr$immig <- get.vertex.attribute(nw[[1]], "immig") 
  dat$attr$immig.loc <- rep(0,length(dat$attr$age))

  # UAI group
  p1.msm <- dat$param$cond.pers.always.prob.msm
  p2.msm <- dat$param$cond.inst.always.prob.msm
  rho.msm <- dat$param$cond.always.prob.corr.msm
  uai.always.msm <- bindata::rmvbin(num, c(p1.msm, p2.msm), bincorr = (1 - rho.msm) * diag(2) + rho.msm)
  dat$attr$cond.always.pers.msm <- uai.always.msm[, 1]
  dat$attr$cond.always.inst.msm <- uai.always.msm[, 2]
  
  # UVI group
  p1.het <- dat$param$cond.pers.always.prob.het
  p2.het <- dat$param$cond.inst.always.prob.het
  rho.het <- dat$param$cond.always.prob.corr.het
  uvi.always.het <- bindata::rmvbin(num, c(p1.het, p2.het), bincorr = (1 - rho.het) * diag(2) + rho.het)
  dat$attr$cond.always.pers.het <- uvi.always.het[, 1]
  dat$attr$cond.always.inst.het <- uvi.always.het[, 2]
  

  # Arrival and departure
  dat$attr$arrival.time <- rep(1, num)

  # Circumcision
  circ <- rep(NA, num)

  if (length(ids.B.m) > 0) {
  circ[ids.B.m] <- sample(apportion_lr(length(ids.B.m), 0:1, 1 - param$circ.B.prob))}
  if (length(ids.BI.m) > 0) { 
  circ[ids.BI.m] <- sample(apportion_lr(length(ids.BI.m), 0:1, 1 - param$circ.BI.prob))}
  if (length(ids.H.m) > 0) { 
  circ[ids.H.m] <- sample(apportion_lr(length(ids.H.m), 0:1, 1 - param$circ.H.prob))}
  if (length(ids.HI.m) > 0) { 
  circ[ids.HI.m] <- sample(apportion_lr(length(ids.HI.m), 0:1, 1 - param$circ.HI.prob))}
  if (length(ids.W.m) > 0) { 
  circ[ids.W.m] <- sample(apportion_lr(length(ids.W.m), 0:1, 1 - param$circ.W.prob))}

  dat$attr$circ <- circ

  # PrEP Attributes
  dat$attr$prepClass <- rep(NA, num)
  dat$attr$prepElig <- rep(NA, num)
  dat$attr$prepStat <- rep(0, num)

  # One-off AI class
#  inst.ai.class <- rep(NA, num)
#  ncl <- param$num.inst.ai.classes
#  inst.ai.class[ids.B] <- sample(apportion_lr(num.B, 1:ncl, rep(1 / ncl, ncl)))
#  inst.ai.class[ids.W] <- sample(apportion_lr(num.W, 1:ncl, rep(1 / ncl, ncl)))
#  dat$attr$inst.ai.class <- inst.ai.class

  # Role class

  role.class <- get.vertex.attribute(nw[[1]], "role.class")
  dat$attr$role.class <- role.class
  
  # Ins.quot
  ins.quot <- rep(NA, num)
  ins.quot[role.class == "I"]  <- 1
  ins.quot[role.class == "R"]  <- 0
  ins.quot[role.class == "V"]  <- runif(sum(role.class == "V"))
  dat$attr$ins.quot <- ins.quot

  
  # HIV-related attributes
  dat <- init_status_shamp(dat)

  # CCR5
  dat <- init_ccr5_shamp(dat)


  # Network statistics
  #dat$stats$nwstats <- list()


  # Prevalence Tracking
  dat$temp$deg.dists <- list()
  dat$temp$discl.list <- matrix(NA, nrow = 0, ncol = 3)
  colnames(dat$temp$discl.list) <- c("pos", "neg", "discl.time")

  if (control$save.nwstats == TRUE) {
    dat$stats <- list()
    dat$stats$nwstats <- list()

  }

  dat <- prevalence_shamp(dat, at = 1)

  class(dat) <- "dat"
  return(dat)
}



#' @title Initialize the HIV status of persons in the network
#'
#' @description Sets the initial individual-level disease status of persons
#'              in the network, as well as disease-related attributes for
#'              infected persons.
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#' @keywords initiation utility shamp
#'
init_status_shamp <- function(dat) {
  
  age <- dat$attr$age
  race <- dat$attr$race
  sex<- dat$attr$sex
  sex.ident<-dat$attr$sex.ident
  
  num.msm<-max(0,sum(dat$attr$sex.ident=="msm"))
  num.het<-max(0,sum(dat$attr$sex.ident=="f")+sum(dat$attr$sex.ident=="msf"))
  num.bi<-max(0,sum(dat$attr$sex.ident=="msmf"))
                            
  num.B.f <- max(0,length(which(dat$attr$race == "B" & dat$attr$sex.ident=="f")))
  num.BI.f <- max(0,length(which(dat$attr$race == "BI" & dat$attr$sex.ident=="f")))
  num.H.f <- max(0,length(which(dat$attr$race == "H" & dat$attr$sex.ident=="f")))
  num.HI.f <- max(0,length(which(dat$attr$race == "HI" & dat$attr$sex.ident=="f")))
  num.W.f <- max(0,length(which(dat$attr$race == "W" & dat$attr$sex.ident=="f")))
  
  num.B.msf <- max(0,length(which(dat$attr$race == "B" & dat$attr$sex.ident=="msf")))
  num.BI.msf <- max(0,length(which(dat$attr$race == "BI" & dat$attr$sex.ident=="msf")))
  num.H.msf <- max(0,length(which(dat$attr$race == "H" & dat$attr$sex.ident=="msf")))
  num.HI.msf <- max(0,length(which(dat$attr$race == "HI" & dat$attr$sex.ident=="msf")))
  num.W.msf <- max(0,length(which(dat$attr$race == "W" & dat$attr$sex.ident=="msf")))
  
  num.B.msm <- max(0,length(which(dat$attr$race == "B" & dat$attr$sex.ident=="msm")))
  num.BI.msm <- max(0,length(which(dat$attr$race == "BI" & dat$attr$sex.ident=="msm")))
  num.H.msm <- max(0,length(which(dat$attr$race == "H" & dat$attr$sex.ident=="msm")))
  num.HI.msm <- max(0,length(which(dat$attr$race == "HI" & dat$attr$sex.ident=="msm")))
  num.W.msm <- max(0,length(which(dat$attr$race == "W" & dat$attr$sex.ident=="msm")))
  
  num.B.msmf <- max(0,length(which(dat$attr$race == "B" & dat$attr$sex.ident=="msmf")))
  num.BI.msmf <- max(0,length(which(dat$attr$race == "BI" & dat$attr$sex.ident=="msmf")))
  num.H.msmf <- max(0,length(which(dat$attr$race == "H" & dat$attr$sex.ident=="msmf")))
  num.HI.msmf <- max(0,length(which(dat$attr$race == "HI" & dat$attr$sex.ident=="msmf")))
  num.W.msmf <- max(0,length(which(dat$attr$race == "W" & dat$attr$sex.ident=="msmf")))
  
  num <- num.B.f + num.BI.f + num.H.f + num.HI.f + num.W.f + 
         num.B.msf + num.BI.msf + num.H.msf + num.HI.msf + num.W.msf +
         num.B.msm + num.BI.msm + num.H.msm + num.HI.msm + num.W.msm +
         num.B.msmf + num.BI.msmf + num.H.msmf + num.HI.msmf + num.W.msmf  
  
  ##Females
  ids.B.f <- which(dat$attr$race == "B" & dat$attr$sex=="F")
  ids.BI.f <- which(dat$attr$race == "BI" & dat$attr$sex=="F")
  ids.H.f <- which(dat$attr$race == "H" & dat$attr$sex=="F")
  ids.HI.f <- which(dat$attr$race == "HI" & dat$attr$sex=="F")
  ids.W.f <- which(dat$attr$race == "W" & dat$attr$sex=="F")
  

  ##Males MSF.
  ids.B.msf <- which(dat$attr$race == "B" & dat$attr$sex.ident=="msf")
  ids.BI.msf <- which(dat$attr$race == "BI" & dat$attr$sex.ident=="msf")
  ids.H.msf <- which(dat$attr$race == "H" & dat$attr$sex.ident=="msf")
  ids.HI.msf <- which(dat$attr$race == "HI" & dat$attr$sex.ident=="msf")
  ids.W.msf <- which(dat$attr$race == "W" & dat$attr$sex.ident=="msf")
  
  ##Males MSM.
  ids.B.msm <- which(dat$attr$race == "B" & dat$attr$sex.ident=="msm")
  ids.BI.msm <- which(dat$attr$race == "BI" & dat$attr$sex.ident=="msm")
  ids.H.msm <- which(dat$attr$race == "H" & dat$attr$sex.ident=="msm")
  ids.HI.msm <- which(dat$attr$race == "HI" & dat$attr$sex.ident=="msm")
  ids.W.msm <- which(dat$attr$race == "W" & dat$attr$sex.ident=="msm")

  
  ##Males MSMF.
  ids.B.msmf <- which(dat$attr$race == "B" & dat$attr$sex.ident=="msmf")
  ids.BI.msmf <- which(dat$attr$race == "BI" & dat$attr$sex.ident=="msmf")
  ids.H.msmf <- which(dat$attr$race == "H" & dat$attr$sex.ident=="msmf")
  ids.HI.msmf <- which(dat$attr$race == "HI" & dat$attr$sex.ident=="msmf")
  ids.W.msmf <- which(dat$attr$race == "W" & dat$attr$sex.ident=="msmf")
  

  # Infection Status
  nInfB.f <- ifelse(num.B.f > 0 & dat$init$prev.B.f > 0 , max(1,round(dat$init$prev.B.f * num.B.f)),round(dat$init$prev.B.f * num.B.f))
  nInfBI.f <- ifelse(num.BI.f > 0 & dat$init$prev.BI.f > 0 , max(1,round(dat$init$prev.BI.f * num.BI.f)),round(dat$init$prev.BI.f * num.BI.f))
  nInfH.f <- ifelse(num.H.f > 0 & dat$init$prev.H.f > 0 , max(1,round(dat$init$prev.H.f * num.H.f)),round(dat$init$prev.H.f * num.H.f))
  nInfHI.f <- ifelse(num.HI.f > 0 & dat$init$prev.HI.f > 0 , max(1,round(dat$init$prev.HI.f * num.HI.f)),round(dat$init$prev.HI.f * num.HI.f))
  nInfW.f <- ifelse(num.W.f > 0 & dat$init$prev.W.f > 0 , max(1,round(dat$init$prev.W.f * num.W.f)),round(dat$init$prev.W.f * num.W.f))
  
  nInfB.msf <- ifelse(num.B.msf > 0 & dat$init$prev.B.msf > 0 , max(1,round(dat$init$prev.B.msf * num.B.msf)),round(dat$init$prev.B.msf * num.B.msf))
  nInfBI.msf <- ifelse(num.BI.msf > 0 & dat$init$prev.BI.msf > 0 , max(1,round(dat$init$prev.BI.msf * num.BI.msf)),round(dat$init$prev.BI.msf * num.BI.msf))
  nInfH.msf <- ifelse(num.H.msf > 0 & dat$init$prev.H.msf > 0 , max(1,round(dat$init$prev.H.msf * num.H.msf)),round(dat$init$prev.H.msf * num.H.msf))
  nInfHI.msf <- ifelse(num.HI.msf > 0 & dat$init$prev.HI.msf > 0 , max(1,round(dat$init$prev.HI.msf * num.HI.msf)),round(dat$init$prev.HI.msf * num.HI.msf))
  nInfW.msf <- ifelse(num.W.msf > 0 & dat$init$prev.W.msf > 0 , max(1,round(dat$init$prev.W.msf * num.W.msf)),round(dat$init$prev.W.msf * num.W.msf))
  
  nInfB.msm <- ifelse(num.B.msm > 0 & dat$init$prev.B.msm > 0 , max(1,round(dat$init$prev.B.msm * num.B.msm)),round(dat$init$prev.B.msm * num.B.msm))
  nInfBI.msm <- ifelse(num.BI.msm > 0 & dat$init$prev.BI.msm > 0 , max(1,round(dat$init$prev.BI.msm * num.BI.msm)),round(dat$init$prev.BI.msm * num.BI.msm))
  nInfH.msm <- ifelse(num.H.msm > 0 & dat$init$prev.H.msm > 0 , max(1,round(dat$init$prev.H.msm * num.H.msm)),round(dat$init$prev.H.msm * num.H.msm))
  nInfHI.msm <- ifelse(num.HI.msm > 0 & dat$init$prev.HI.msm > 0 , max(1,round(dat$init$prev.HI.msm * num.HI.msm)),round(dat$init$prev.HI.msm * num.HI.msm))
  nInfW.msm <- ifelse(num.W.msm > 0 & dat$init$prev.W.msm > 0 , max(1,round(dat$init$prev.W.msm * num.W.msm)),round(dat$init$prev.W.msm * num.W.msm))
  
  nInfB.msmf <- ifelse(num.B.msmf > 0 & dat$init$prev.B.msmf > 0 , max(1,round(dat$init$prev.B.msmf * num.B.msmf)),round(dat$init$prev.B.msmf * num.B.msmf))
  nInfBI.msmf <- ifelse(num.BI.msmf > 0 & dat$init$prev.BI.msmf > 0 , max(1,round(dat$init$prev.BI.msmf * num.BI.msmf)),round(dat$init$prev.BI.msmf * num.BI.msmf))
  nInfH.msmf <- ifelse(num.H.msmf > 0 & dat$init$prev.H.msmf > 0 , max(1,round(dat$init$prev.H.msmf * num.H.msmf)),round(dat$init$prev.H.msmf * num.H.msmf))
  nInfHI.msmf <- ifelse(num.HI.msmf > 0 & dat$init$prev.HI.msmf > 0 , max(1,round(dat$init$prev.HI.msmf * num.HI.msmf)),round(dat$init$prev.HI.msmf * num.HI.msmf))
  nInfW.msmf <- ifelse(num.W.msmf > 0 & dat$init$prev.W.msmf > 0 , max(1,round(dat$init$prev.W.msmf * num.W.msmf)),round(dat$init$prev.W.msmf * num.W.msmf))

  # Age-based infection probability
  probInfCrB.f <- age[ids.B.f] * dat$init$init.prev.age.slope.B.f
  probInfB.f <- probInfCrB.f + (nInfB.f - sum(probInfCrB.f)) / num.B.f
  
  probInfCrBI.f <- age[ids.BI.f] * dat$init$init.prev.age.slope.BI.f
  probInfBI.f <- probInfCrBI.f + (nInfBI.f - sum(probInfCrBI.f)) / num.BI.f
  
  probInfCrH.f <- age[ids.H.f] * dat$init$init.prev.age.slope.H.f
  probInfH.f <- probInfCrH.f + (nInfH.f - sum(probInfCrH.f)) / num.H.f
  
  probInfCrHI.f <- age[ids.HI.f] * dat$init$init.prev.age.slope.HI.f
  probInfHI.f <- probInfCrHI.f + (nInfHI.f - sum(probInfCrHI.f)) / num.HI.f
  
  probInfCrW.f <- age[ids.W.f] * dat$init$init.prev.age.slope.W.f
  probInfW.f <- probInfCrW.f + (nInfW.f - sum(probInfCrW.f)) / num.W.f
  
  ##Male msf.
  probInfCrB.msf <- age[ids.B.msf] * dat$init$init.prev.age.slope.B.msf
  probInfB.msf <- probInfCrB.msf + (nInfB.msf - sum(probInfCrB.msf)) / num.B.msf
  
  probInfCrBI.msf <- age[ids.BI.msf] * dat$init$init.prev.age.slope.BI.msf
  probInfBI.msf <- probInfCrBI.msf + (nInfBI.msf - sum(probInfCrBI.msf)) / num.BI.msf
  
  probInfCrH.msf <- age[ids.H.msf] * dat$init$init.prev.age.slope.H.msf
  probInfH.msf <- probInfCrH.msf + (nInfH.msf - sum(probInfCrH.msf)) / num.H.msf
  
  probInfCrHI.msf <- age[ids.HI.msf] * dat$init$init.prev.age.slope.HI.msf
  probInfHI.msf <- probInfCrHI.msf + (nInfHI.msf - sum(probInfCrHI.msf)) / num.HI.msf
  
  probInfCrW.msf <- age[ids.W.msf] * dat$init$init.prev.age.slope.W.msf
  probInfW.msf <- probInfCrW.msf + (nInfW.msf - sum(probInfCrW.msf)) / num.W.msf
  
  ##Male msm
  probInfCrB.msm <- age[ids.B.msm] * dat$init$init.prev.age.slope.B.msm
  probInfB.msm <- probInfCrB.msm + (nInfB.msm - sum(probInfCrB.msm)) / num.B.msm
  
  probInfCrBI.msm <- age[ids.BI.msm] * dat$init$init.prev.age.slope.BI.msm
  probInfBI.msm <- probInfCrBI.msm + (nInfBI.msm - sum(probInfCrBI.msm)) / num.BI.msm
  
  probInfCrH.msm <- age[ids.H.msm] * dat$init$init.prev.age.slope.H.msm
  probInfH.msm <- probInfCrH.msm + (nInfH.msm - sum(probInfCrH.msm)) / num.H.msm
  
  probInfCrHI.msm <- age[ids.HI.msm] * dat$init$init.prev.age.slope.HI.msm
  probInfHI.msm <- probInfCrHI.msm + (nInfHI.msm - sum(probInfCrHI.msm)) / num.HI.msm
  
  probInfCrW.msm <- age[ids.W.msm] * dat$init$init.prev.age.slope.W.msm
  probInfW.msm <- probInfCrW.msm + (nInfW.msm - sum(probInfCrW.msm)) / num.W.msm
  
  ##Male msmf
  probInfCrB.msmf <- age[ids.B.msmf] * dat$init$init.prev.age.slope.B.msmf
  probInfB.msmf <- probInfCrB.msmf + (nInfB.msmf - sum(probInfCrB.msmf)) / num.B.msmf
  
  probInfCrBI.msmf <- age[ids.BI.msmf] * dat$init$init.prev.age.slope.BI.msmf
  probInfBI.msmf <- probInfCrBI.msmf + (nInfBI.msmf - sum(probInfCrBI.msmf)) / num.BI.msmf
  
  probInfCrH.msmf <- age[ids.H.msmf] * dat$init$init.prev.age.slope.H.msmf
  probInfH.msmf <- probInfCrH.msmf + (nInfH.msmf - sum(probInfCrH.msmf)) / num.H.msmf
  
  probInfCrHI.msmf <- age[ids.HI.msmf] * dat$init$init.prev.age.slope.HI.msmf
  probInfHI.msmf <- probInfCrHI.msmf + (nInfHI.msmf - sum(probInfCrHI.msmf)) / num.HI.msmf
  
  probInfCrW.msmf <- age[ids.W.msmf] * dat$init$init.prev.age.slope.W.msmf
  probInfW.msmf <- probInfCrW.msmf + (nInfW.msmf - sum(probInfCrW.msmf)) / num.W.msmf
  

  if(num.het>0){
  if (any(probInfB.f < 0 , probInfBI.f < 0 , probInfH.f < 0 , probInfHI.f < 0 , probInfW.f < 0 , 
          probInfB.msf < 0 , probInfBI.msf < 0 , probInfH.msf < 0 , probInfHI.msf < 0 , probInfW.msf < 0)) {
    stop("Slope of initial prevalence by age must be sufficiently low to ",
         "avoid non-positive probabilities.", call. = FALSE)
  }
}

  if(num.msm>0){
  if (any(probInfB.msm < 0 , probInfBI.msm < 0 , probInfH.msm < 0 , probInfHI.msm < 0 , probInfW.msm < 0)) {
    stop("Slope of initial prevalence by age must be sufficiently low to ",
         "avoid non-positive probabilities.", call. = FALSE)
  }
  }  

  if(num.bi>0){
    if (any(probInfB.msmf < 0 , probInfBI.msmf < 0 , probInfH.msmf < 0 , probInfHI.msmf < 0 , probInfW.msmf < 0)) {
      stop("Slope of initial prevalence by age must be sufficiently low to ",
           "avoid non-positive probabilities.", call. = FALSE)
    }
  }  
  
 
   # Infection status
   status <- rep(0, num)
   
   #Females.
    while (sum(status[ids.B.f]) != nInfB.f) {
      status[ids.B.f] <- rbinom(num.B.f, 1, probInfB.f)
    }
   while (sum(status[ids.BI.f]) != nInfBI.f) {
     status[ids.BI.f] <- rbinom(num.BI.f, 1, probInfBI.f)
   }
   while (sum(status[ids.H.f]) != nInfH.f) {
     status[ids.H.f] <- rbinom(num.H.f, 1, probInfH.f)
   }
   while (sum(status[ids.HI.f]) != nInfHI.f) {
     status[ids.HI.f] <- rbinom(num.HI.f, 1, probInfHI.f)
   }
   while (sum(status[ids.W.f]) != nInfW.f) {
     status[ids.W.f] <- rbinom(num.W.f, 1, probInfW.f)
   }
   
   #Males msf.
   while (sum(status[ids.B.msf]) != nInfB.msf) {
     status[ids.B.msf] <- rbinom(num.B.msf, 1, probInfB.msf)
   }
   while (sum(status[ids.BI.msf]) != nInfBI.msf) {
     status[ids.BI.msf] <- rbinom(num.BI.msf, 1, probInfBI.msf)
   }
   while (sum(status[ids.H.msf]) != nInfH.msf) {
     status[ids.H.msf] <- rbinom(num.H.msf, 1, probInfH.msf)
   }
   while (sum(status[ids.HI.msf]) != nInfHI.msf) {
     status[ids.HI.msf] <- rbinom(num.HI.msf, 1, probInfHI.msf)
   }
   while (sum(status[ids.W.msf]) != nInfW.msf) {
     status[ids.W.msf] <- rbinom(num.W.msf, 1, probInfW.msf)
   }

   #Males MSM
   while (sum(status[ids.B.msm]) != nInfB.msm) {
     status[ids.B.msm] <- rbinom(num.B.msm, 1, probInfB.msm)
   }
   while (sum(status[ids.BI.msm]) != nInfBI.msm) {
     status[ids.BI.msm] <- rbinom(num.BI.msm, 1, probInfBI.msm)
   }
   while (sum(status[ids.H.msm]) != nInfH.msm) {
     status[ids.H.msm] <- rbinom(num.H.msm, 1, probInfH.msm)
   }
   while (sum(status[ids.HI.msm]) != nInfHI.msm) {
     status[ids.HI.msm] <- rbinom(num.HI.msm, 1, probInfHI.msm)
   }
   while (sum(status[ids.W.msm]) != nInfW.msm) {
     status[ids.W.msm] <- rbinom(num.W.msm, 1, probInfW.msm)
   }

   #Males MSMF
   while (sum(status[ids.B.msmf]) != nInfB.msmf) {
     status[ids.B.msmf] <- rbinom(num.B.msmf, 1, probInfB.msmf)
   }
   while (sum(status[ids.BI.msmf]) != nInfBI.msmf) {
     status[ids.BI.msmf] <- rbinom(num.BI.msmf, 1, probInfBI.msmf)
   }
   while (sum(status[ids.H.msmf]) != nInfH.msmf) {
     status[ids.H.msmf] <- rbinom(num.H.msmf, 1, probInfH.msmf)
   }
   while (sum(status[ids.HI.msmf]) != nInfHI.msmf) {
     status[ids.HI.msmf] <- rbinom(num.HI.msmf, 1, probInfHI.msmf)
   }
   while (sum(status[ids.W.msmf]) != nInfW.msmf) {
     status[ids.W.msmf] <- rbinom(num.W.msmf, 1, probInfW.msmf)
   }
  
    dat$attr$status <- status
    dat$attr$inf.class[dat$attr$status==1] <- "L"

  # Treatment trajectory
  tt.traj <- rep(NA, num)

  #female
  if (length(ids.B.f) > 0){ 
  tt.traj[ids.B.f] <- sample(apportion_lr(num.B.f, c(1, 2, 3, 4),
                                        dat$param$tt.traj.B.f.prob))}
  if (length(ids.BI.f) > 0){ 
  tt.traj[ids.BI.f] <- sample(apportion_lr(num.BI.f, c(1, 2, 3, 4),
                                          dat$param$tt.traj.BI.f.prob))}
  if (length(ids.H.f) > 0){ 
  tt.traj[ids.H.f] <- sample(apportion_lr(num.H.f, c(1, 2, 3, 4),
                                          dat$param$tt.traj.H.f.prob))}
  if (length(ids.HI.f) > 0){ 
  tt.traj[ids.HI.f] <- sample(apportion_lr(num.HI.f, c(1, 2, 3, 4),
                                           dat$param$tt.traj.HI.f.prob))}
  if (length(ids.W.f) > 0){ 
  tt.traj[ids.W.f] <- sample(apportion_lr(num.W.f, c(1, 2, 3, 4),
                                          dat$param$tt.traj.W.f.prob))}
  
  #Males MSF. 
  if (length(ids.B.msf) > 0){ 
  tt.traj[ids.B.msf] <- sample(apportion_lr(num.B.msf, c(1, 2, 3, 4),
                                          dat$param$tt.traj.B.msf.prob))}
  if (length(ids.BI.msf) > 0){ 
  tt.traj[ids.BI.msf] <- sample(apportion_lr(num.BI.msf, c(1, 2, 3, 4),
                                           dat$param$tt.traj.BI.msf.prob))}
  if (length(ids.H.msf) > 0){ 
  tt.traj[ids.H.msf] <- sample(apportion_lr(num.H.msf, c(1, 2, 3, 4),
                                          dat$param$tt.traj.H.msf.prob))}
  if (length(ids.HI.msf) > 0){
  tt.traj[ids.HI.msf] <- sample(apportion_lr(num.HI.msf, c(1, 2, 3, 4),
                                           dat$param$tt.traj.HI.msf.prob))}
  if (length(ids.W.msf) > 0){
  tt.traj[ids.W.msf] <- sample(apportion_lr(num.W.msf, c(1, 2, 3, 4),
                                          dat$param$tt.traj.W.msf.prob))}
  #Males MSM.
  if (length(ids.B.msm) > 0){ 
    tt.traj[ids.B.msm] <- sample(apportion_lr(num.B.msm, c(1, 2, 3, 4),
                                            dat$param$tt.traj.B.msm.prob))}
  if (length(ids.BI.msm) > 0){ 
    tt.traj[ids.BI.msm] <- sample(apportion_lr(num.BI.msm, c(1, 2, 3, 4),
                                             dat$param$tt.traj.BI.msm.prob))}
  if (length(ids.H.msm) > 0){ 
    tt.traj[ids.H.msm] <- sample(apportion_lr(num.H.msm, c(1, 2, 3, 4),
                                            dat$param$tt.traj.H.msm.prob))}
  if (length(ids.HI.msm) > 0){
    tt.traj[ids.HI.msm] <- sample(apportion_lr(num.HI.msm, c(1, 2, 3, 4),
                                             dat$param$tt.traj.HI.msm.prob))}
  if (length(ids.W.msm) > 0){
    tt.traj[ids.W.msm] <- sample(apportion_lr(num.W.msm, c(1, 2, 3, 4),
                                            dat$param$tt.traj.W.msm.prob))}
  
  #Males MSM.
  if (length(ids.B.msmf) > 0){ 
    tt.traj[ids.B.msmf] <- sample(apportion_lr(num.B.msmf, c(1, 2, 3, 4),
                                              dat$param$tt.traj.B.msmf.prob))}
  if (length(ids.BI.msmf) > 0){ 
    tt.traj[ids.BI.msmf] <- sample(apportion_lr(num.BI.msmf, c(1, 2, 3, 4),
                                               dat$param$tt.traj.BI.msmf.prob))}
  if (length(ids.H.msmf) > 0){ 
    tt.traj[ids.H.msmf] <- sample(apportion_lr(num.H.msmf, c(1, 2, 3, 4),
                                              dat$param$tt.traj.H.msmf.prob))}
  if (length(ids.HI.msmf) > 0){
    tt.traj[ids.HI.msmf] <- sample(apportion_lr(num.HI.msmf, c(1, 2, 3, 4),
                                               dat$param$tt.traj.HI.msmf.prob))}
  if (length(ids.W.msmf) > 0){
    tt.traj[ids.W.msmf] <- sample(apportion_lr(num.W.msmf, c(1, 2, 3, 4),
                                              dat$param$tt.traj.W.msmf.prob))}
  
  dat$attr$tt.traj <- tt.traj



  ## Infection-related attributes

  stage <- rep(NA, num)
  stage.time <- rep(NA, num)
  inf.time <- rep(NA, num)
  vl <- rep(NA, num)
  diag.status <- rep(NA, num)
  diag.time <- rep(NA, num)
  last.neg.test <- rep(NA, num)
  tx.status <- rep(NA, num)
  tx.init.time <- rep(NA, num)
  cum.time.on.tx <- rep(NA, num)
  cum.time.off.tx <- rep(NA, num)
  infector <- rep(NA, num)
  inf.role <- rep(NA, num)
  inf.type <- rep(NA, num)
  inf.diag <- rep(NA, num)
  inf.tx <- rep(NA, num)
  inf.stage <- rep(NA, num)
  inf.class <- rep(NA, num)

  time.sex.active <- pmax(1,
                          round((365 / dat$param$time.unit) * age - (365 / dat$param$time.unit) *
                                  min(age), 0))

  vlar.int <- dat$param$vl.acute.rise.int
  vlap <- dat$param$vl.acute.peak
  vlaf.int <- dat$param$vl.acute.fall.int
  vlsp <- dat$param$vl.set.point
  vldo.int <- dat$param$vl.aids.onset.int
  vl.aids.int <- dat$param$vl.aids.int
  vlf  <- dat$param$vl.fatal
  vlds <- (vlf - vlsp) / vl.aids.int
  vl.acute.int <- vlar.int + vlaf.int


  ### Non-treater type: tester and non-tester
  selected <- which(status == 1 & tt.traj %in% c(1, 2))
  max.inf.time <- pmin(time.sex.active[selected], vldo.int + vl.aids.int)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  tx.status[selected] <- 0
  cum.time.on.tx[selected] <- 0
  cum.time.off.tx[selected] <- time.since.inf

  stage[selected[time.since.inf <= vlar.int]] <- 1
  stage[selected[time.since.inf > vlar.int & time.since.inf <= vl.acute.int]] <- 2
  stage[selected[time.since.inf > vl.acute.int & time.since.inf <= vldo.int]] <- 3
  stage[selected[time.since.inf > vldo.int]] <- 4

  stage.time[selected][stage[selected] == 1] <- time.since.inf[stage[selected] == 1]
  stage.time[selected][stage[selected] == 2] <- time.since.inf[stage[selected] == 2] -
                                                   vlar.int
  stage.time[selected][stage[selected] == 3] <- time.since.inf[stage[selected] == 3] -
                                                  vl.acute.int
  stage.time[selected][stage[selected] == 4] <- time.since.inf[stage[selected] == 4] -
                                                  vldo.int

  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) * (time.since.inf <= vldo.int) * (vlsp) +
                  (time.since.inf > vldo.int) * (vlsp + (time.since.inf - vldo.int) * vlds)

  selected <- which(status == 1 & tt.traj == 1)
  diag.status[selected] <- 0

  selected <- which(status == 1 & tt.traj == 2)

  # Time to next test
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                             min = 0,
                             max = dat$param$mean.test.B.f.int * (race[selected] == "B" & sex[selected] == "F") +
                                   dat$param$mean.test.BI.f.int * (race[selected] == "BI" & sex[selected] == "F") +
                                   dat$param$mean.test.H.f.int * (race[selected] == "H" & sex[selected] == "F") +
                                   dat$param$mean.test.HI.f.int * (race[selected] == "HI" & sex[selected] == "F") +
                                   dat$param$mean.test.W.f.int * (race[selected] == "W" & sex[selected] == "F") +
                                   dat$param$mean.test.B.msf.int * (race[selected] == "B" & sex.ident[selected] == "msf") +
                                   dat$param$mean.test.BI.msf.int * (race[selected] == "BI" & sex.ident[selected] == "msf") +
                                   dat$param$mean.test.H.msf.int * (race[selected] == "H" & sex.ident[selected] == "msf") +
                                   dat$param$mean.test.HI.msf.int * (race[selected] == "HI" & sex.ident[selected] == "msf") +
                                   dat$param$mean.test.W.msf.int * (race[selected] == "W" & sex.ident[selected] == "msf") +
                                   dat$param$mean.test.B.msm.int * (race[selected] == "B" & sex.ident[selected] == "msm") +
                                   dat$param$mean.test.BI.msm.int * (race[selected] == "BI" & sex.ident[selected] == "msm") +
                                   dat$param$mean.test.H.msm.int * (race[selected] == "H" & sex.ident[selected] == "msm") +
                                   dat$param$mean.test.HI.msm.int * (race[selected] == "HI" & sex.ident[selected] == "msm") +
                                   dat$param$mean.test.W.msm.int * (race[selected] == "W" & sex.ident[selected] == "msm") +
                                   dat$param$mean.test.B.msmf.int * (race[selected] == "B" & sex.ident[selected] == "msmf") +
                                   dat$param$mean.test.BI.msmf.int * (race[selected] == "BI" & sex.ident[selected] == "msmf") +
                                   dat$param$mean.test.H.msmf.int * (race[selected] == "H" & sex.ident[selected] == "msmf") +
                                   dat$param$mean.test.HI.msmf.int * (race[selected] == "HI" & sex.ident[selected] == "msmf") +
                                   dat$param$mean.test.W.msmf.int * (race[selected] == "W" & sex.ident[selected] == "msmf")))
  }
  
  if (dat$param$testing.pattern == "memoryless") {
    ttntest <- rgeom(length(selected),
                     1 / (dat$param$mean.test.B.f.int * (race[selected] == "B" & sex[selected] == "F") +
                          dat$param$mean.test.BI.f.int * (race[selected] == "BI" & sex[selected] == "F") +
                          dat$param$mean.test.H.f.int * (race[selected] == "H" & sex[selected] == "F") +
                          dat$param$mean.test.HI.f.int * (race[selected] == "HI" & sex[selected] == "F") +
                          dat$param$mean.test.W.f.int * (race[selected] == "W" & sex[selected] == "F") +
                          dat$param$mean.test.B.msf.int * (race[selected] == "B" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.BI.msf.int * (race[selected] == "BI" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.H.msf.int * (race[selected] == "H" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.HI.msf.int * (race[selected] == "HI" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.W.msf.int * (race[selected] == "W" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.B.msm.int * (race[selected] == "B" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.BI.msm.int * (race[selected] == "BI" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.H.msm.int * (race[selected] == "H" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.HI.msm.int * (race[selected] == "HI" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.W.msm.int * (race[selected] == "W" & sex.ident[selected] == "msm") +
                            dat$param$mean.test.B.msmf.int * (race[selected] == "B" & sex.ident[selected] == "msmf") +
                            dat$param$mean.test.BI.msmf.int * (race[selected] == "BI" & sex.ident[selected] == "msmf") +
                            dat$param$mean.test.H.msmf.int * (race[selected] == "H" & sex.ident[selected] == "msmf") +
                            dat$param$mean.test.HI.msmf.int * (race[selected] == "HI" & sex.ident[selected] == "msmf") +
                            dat$param$mean.test.W.msmf.int * (race[selected] == "W" & sex.ident[selected] == "msmf")
                          ))
  }
  



  twind.int <- dat$param$test.window.int
  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
                           -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]

  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1


  ### Full adherent type

  # Create set of expected values for (cum.time.off.tx, cum.time.on.tx)

  tx.init.time.B.f <- twind.int + dat$param$last.neg.test.B.f.int + 1 / dat$param$tx.init.B.f.prob
  tx.init.time.BI.f <- twind.int + dat$param$last.neg.test.BI.f.int + 1 / dat$param$tx.init.BI.f.prob
  tx.init.time.H.f <- twind.int + dat$param$last.neg.test.H.f.int + 1 / dat$param$tx.init.H.f.prob
  tx.init.time.HI.f <- twind.int + dat$param$last.neg.test.HI.f.int + 1 / dat$param$tx.init.HI.f.prob
  tx.init.time.W.f <- twind.int + dat$param$last.neg.test.W.f.int + 1 / dat$param$tx.init.W.f.prob
  
  tx.init.time.B.msf <- twind.int + dat$param$last.neg.test.B.msf.int + 1 / dat$param$tx.init.B.msf.prob
  tx.init.time.BI.msf <- twind.int + dat$param$last.neg.test.BI.msf.int + 1 / dat$param$tx.init.BI.msf.prob
  tx.init.time.H.msf <- twind.int + dat$param$last.neg.test.H.msf.int + 1 / dat$param$tx.init.H.msf.prob
  tx.init.time.HI.msf <- twind.int + dat$param$last.neg.test.HI.msf.int + 1 / dat$param$tx.init.HI.msf.prob
  tx.init.time.W.msf <- twind.int + dat$param$last.neg.test.W.msf.int + 1 / dat$param$tx.init.W.msf.prob
  
  tx.init.time.B.msm <- twind.int + dat$param$last.neg.test.B.msm.int + 1 / dat$param$tx.init.B.msm.prob
  tx.init.time.BI.msm <- twind.int + dat$param$last.neg.test.BI.msm.int + 1 / dat$param$tx.init.BI.msm.prob
  tx.init.time.H.msm <- twind.int + dat$param$last.neg.test.H.msm.int + 1 / dat$param$tx.init.H.msm.prob
  tx.init.time.HI.msm <- twind.int + dat$param$last.neg.test.HI.msm.int + 1 / dat$param$tx.init.HI.msm.prob
  tx.init.time.W.msm <- twind.int + dat$param$last.neg.test.W.msm.int + 1 / dat$param$tx.init.W.msm.prob
  
  tx.init.time.B.msmf <- twind.int + dat$param$last.neg.test.B.msmf.int + 1 / dat$param$tx.init.B.msmf.prob
  tx.init.time.BI.msmf <- twind.int + dat$param$last.neg.test.BI.msmf.int + 1 / dat$param$tx.init.BI.msmf.prob
  tx.init.time.H.msmf <- twind.int + dat$param$last.neg.test.H.msmf.int + 1 / dat$param$tx.init.H.msmf.prob
  tx.init.time.HI.msmf <- twind.int + dat$param$last.neg.test.HI.msmf.int + 1 / dat$param$tx.init.HI.msmf.prob
  tx.init.time.W.msmf <- twind.int + dat$param$last.neg.test.W.msmf.int + 1 / dat$param$tx.init.W.msmf.prob

  #Stage Females 
  # Stage for Black females
  prop.time.on.tx.B.f <- dat$param$tx.reinit.B.f.prob /
                       (dat$param$tx.halt.B.f.prob + dat$param$tx.reinit.B.f.prob)
  offon.B.f <- matrix(c(1:tx.init.time.B.f, rep(0, tx.init.time.B.f)),
                    nrow = tx.init.time.B.f)
  numsteps.B.f <- (dat$param$max.time.off.tx.full.int - tx.init.time.B.f) /
                (1 - prop.time.on.tx.B.f)
  offon.B.f <- rbind(offon.B.f,
                   cbind(tx.init.time.B.f + (1 - prop.time.on.tx.B.f) * 1:numsteps.B.f,
                         prop.time.on.tx.B.f * 1:numsteps.B.f))
  offon.B.f <- round(offon.B.f)
  exp.dur.chronic.B.f <- nrow(offon.B.f) - vl.acute.int
  exp.onset.aids.B.f <- nrow(offon.B.f)
  offon.last.B.f <- offon.B.f[nrow(offon.B.f), ]
  offon.B.f <- rbind(offon.B.f,
                   matrix(c(offon.last.B.f[1] + (1:vl.aids.int),
                            rep(offon.last.B.f[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.B.f <- nrow(offon.B.f)
  offon.B.f[, 2] <- (1:max.possible.inf.time.B.f) - offon.B.f[, 1]
  stage.B.f <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B.f, vl.aids.int))
  stage.time.B.f <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B.f, 1:vl.aids.int)
  
  # Stage for Black Imigrant females
  prop.time.on.tx.BI.f <- dat$param$tx.reinit.BI.f.prob /
    (dat$param$tx.halt.BI.f.prob + dat$param$tx.reinit.BI.f.prob)
  offon.BI.f <- matrix(c(1:tx.init.time.BI.f, rep(0, tx.init.time.BI.f)),
                      nrow = tx.init.time.BI.f)
  numsteps.BI.f <- (dat$param$max.time.off.tx.full.int - tx.init.time.BI.f) /
    (1 - prop.time.on.tx.BI.f)
  offon.BI.f <- rbind(offon.BI.f,
                     cbind(tx.init.time.BI.f + (1 - prop.time.on.tx.BI.f) * 1:numsteps.BI.f,
                           prop.time.on.tx.BI.f * 1:numsteps.BI.f))
  offon.BI.f <- round(offon.BI.f)
  exp.dur.chronic.BI.f <- nrow(offon.BI.f) - vl.acute.int
  exp.onset.aids.BI.f <- nrow(offon.BI.f)
  offon.last.BI.f <- offon.BI.f[nrow(offon.BI.f), ]
  offon.BI.f <- rbind(offon.BI.f,
                     matrix(c(offon.last.BI.f[1] + (1:vl.aids.int),
                              rep(offon.last.BI.f[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.BI.f <- nrow(offon.BI.f)
  offon.BI.f[, 2] <- (1:max.possible.inf.time.BI.f) - offon.BI.f[, 1]
  stage.BI.f <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.BI.f, vl.aids.int))
  stage.time.BI.f <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.BI.f, 1:vl.aids.int)
  
  
  # Stage for Hispanic females
  prop.time.on.tx.H.f <- dat$param$tx.reinit.H.f.prob /
    (dat$param$tx.halt.H.f.prob + dat$param$tx.reinit.H.f.prob)
  offon.H.f <- matrix(c(1:tx.init.time.H.f, rep(0, tx.init.time.H.f)),
                      nrow = tx.init.time.H.f)
  numsteps.H.f <- (dat$param$max.time.off.tx.full.int - tx.init.time.H.f) /
    (1 - prop.time.on.tx.H.f)
  offon.H.f <- rbind(offon.H.f,
                     cbind(tx.init.time.H.f + (1 - prop.time.on.tx.H.f) * 1:numsteps.H.f,
                           prop.time.on.tx.H.f * 1:numsteps.H.f))
  offon.H.f <- round(offon.H.f)
  exp.dur.chronic.H.f <- nrow(offon.H.f) - vl.acute.int
  exp.onset.aids.H.f <- nrow(offon.H.f)
  offon.last.H.f <- offon.H.f[nrow(offon.H.f), ]
  offon.H.f <- rbind(offon.H.f,
                     matrix(c(offon.last.H.f[1] + (1:vl.aids.int),
                              rep(offon.last.H.f[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.H.f <- nrow(offon.H.f)
  offon.H.f[, 2] <- (1:max.possible.inf.time.H.f) - offon.H.f[, 1]
  stage.H.f <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.H.f, vl.aids.int))
  stage.time.H.f <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.H.f, 1:vl.aids.int)
  
  # Stage for Hisapnic Imigrant females
  prop.time.on.tx.HI.f <- dat$param$tx.reinit.HI.f.prob /
    (dat$param$tx.halt.HI.f.prob + dat$param$tx.reinit.HI.f.prob)
  offon.HI.f <- matrix(c(1:tx.init.time.HI.f, rep(0, tx.init.time.HI.f)),
                       nrow = tx.init.time.HI.f)
  numsteps.HI.f <- (dat$param$max.time.off.tx.full.int - tx.init.time.HI.f) /
    (1 - prop.time.on.tx.HI.f)
  offon.HI.f <- rbind(offon.HI.f,
                      cbind(tx.init.time.HI.f + (1 - prop.time.on.tx.HI.f) * 1:numsteps.HI.f,
                            prop.time.on.tx.HI.f * 1:numsteps.HI.f))
  offon.HI.f <- round(offon.HI.f)
  exp.dur.chronic.HI.f <- nrow(offon.HI.f) - vl.acute.int
  exp.onset.aids.HI.f <- nrow(offon.HI.f)
  offon.last.HI.f <- offon.HI.f[nrow(offon.HI.f), ]
  offon.HI.f <- rbind(offon.HI.f,
                      matrix(c(offon.last.HI.f[1] + (1:vl.aids.int),
                               rep(offon.last.HI.f[2], vl.aids.int)),
                             ncol = 2))
  max.possible.inf.time.HI.f <- nrow(offon.HI.f)
  offon.HI.f[, 2] <- (1:max.possible.inf.time.HI.f) - offon.HI.f[, 1]
  stage.HI.f <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.HI.f, vl.aids.int))
  stage.time.HI.f <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.HI.f, 1:vl.aids.int) 

  # Stage for White Females
  prop.time.on.tx.W.f <- dat$param$tx.reinit.W.f.prob /
    (dat$param$tx.halt.W.f.prob + dat$param$tx.reinit.W.f.prob)
  offon.W.f <- matrix(c(1:tx.init.time.W.f, rep(0, tx.init.time.W.f)),
                    nrow = tx.init.time.W.f)
  numsteps.W.f <- (dat$param$max.time.off.tx.full.int - tx.init.time.W.f) /
    (1 - prop.time.on.tx.W.f)
  offon.W.f <- rbind(offon.W.f,
                   cbind(tx.init.time.W.f + (1 - prop.time.on.tx.W.f) * 1:numsteps.W.f,
                         prop.time.on.tx.W.f * 1:numsteps.W.f))
  offon.W.f <- round(offon.W.f)
  exp.dur.chronic.W.f <- nrow(offon.W.f) - vl.acute.int
  exp.onset.aids.W.f <- nrow(offon.W.f)
  offon.last.W.f <- offon.W.f[nrow(offon.W.f), ]
  offon.W.f <- rbind(offon.W.f,
                   matrix(c(offon.last.W.f[1] + (1:vl.aids.int),
                            rep(offon.last.W.f[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.W.f <- nrow(offon.W.f)
  offon.W.f[, 2] <- (1:max.possible.inf.time.W.f) - offon.W.f[, 1]
  stage.W.f <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W.f, vl.aids.int))
  stage.time.W.f <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W.f, 1:vl.aids.int)

  
  #Stage for Males HET  MSF
  # Stage for Black males
  prop.time.on.tx.B.msf <- dat$param$tx.reinit.B.msf.prob /
    (dat$param$tx.halt.B.msf.prob + dat$param$tx.reinit.B.msf.prob)
  offon.B.msf <- matrix(c(1:tx.init.time.B.msf, rep(0, tx.init.time.B.msf)),
                      nrow = tx.init.time.B.msf)
  numsteps.B.msf <- (dat$param$max.time.off.tx.full.int - tx.init.time.B.msf) /
    (1 - prop.time.on.tx.B.msf)
  offon.B.msf <- rbind(offon.B.msf,
                     cbind(tx.init.time.B.msf + (1 - prop.time.on.tx.B.msf) * 1:numsteps.B.msf,
                           prop.time.on.tx.B.msf * 1:numsteps.B.msf))
  offon.B.msf <- round(offon.B.msf)
  exp.dur.chronic.B.msf <- nrow(offon.B.msf) - vl.acute.int
  exp.onset.aids.B.msf <- nrow(offon.B.msf)
  offon.last.B.msf <- offon.B.msf[nrow(offon.B.msf), ]
  offon.B.msf <- rbind(offon.B.msf,
                     matrix(c(offon.last.B.msf[1] + (1:vl.aids.int),
                              rep(offon.last.B.msf[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.B.msf <- nrow(offon.B.msf)
  offon.B.msf[, 2] <- (1:max.possible.inf.time.B.msf) - offon.B.msf[, 1]
  stage.B.msf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B.msf, vl.aids.int))
  stage.time.B.msf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B.msf, 1:vl.aids.int)
  
  # Stage for Black Immigrant males MSF
  prop.time.on.tx.BI.msf <- dat$param$tx.reinit.BI.msf.prob /
    (dat$param$tx.halt.BI.msf.prob + dat$param$tx.reinit.BI.msf.prob)
  offon.BI.msf <- matrix(c(1:tx.init.time.BI.msf, rep(0, tx.init.time.BI.msf)),
                       nrow = tx.init.time.BI.msf)
  numsteps.BI.msf <- (dat$param$max.time.off.tx.full.int - tx.init.time.BI.msf) /
    (1 - prop.time.on.tx.BI.msf)
  offon.BI.msf <- rbind(offon.BI.msf,
                      cbind(tx.init.time.BI.msf + (1 - prop.time.on.tx.BI.msf) * 1:numsteps.BI.msf,
                            prop.time.on.tx.BI.msf * 1:numsteps.BI.msf))
  offon.BI.msf <- round(offon.BI.msf)
  exp.dur.chronic.BI.msf <- nrow(offon.BI.msf) - vl.acute.int
  exp.onset.aids.BI.msf <- nrow(offon.BI.msf)
  offon.last.BI.msf <- offon.BI.msf[nrow(offon.BI.msf), ]
  offon.BI.msf <- rbind(offon.BI.msf,
                      matrix(c(offon.last.BI.msf[1] + (1:vl.aids.int),
                               rep(offon.last.BI.msf[2], vl.aids.int)),
                             ncol = 2))
  max.possible.inf.time.BI.msf <- nrow(offon.BI.msf)
  offon.BI.msf[, 2] <- (1:max.possible.inf.time.BI.msf) - offon.BI.msf[, 1]
  stage.BI.msf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.BI.msf, vl.aids.int))
  stage.time.BI.msf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.BI.msf, 1:vl.aids.int)
  
  
  # Stage for Hispanic males MSF
  prop.time.on.tx.H.msf <- dat$param$tx.reinit.H.msf.prob /
    (dat$param$tx.halt.H.msf.prob + dat$param$tx.reinit.H.msf.prob)
  offon.H.msf <- matrix(c(1:tx.init.time.H.msf, rep(0, tx.init.time.H.msf)),
                      nrow = tx.init.time.H.msf)
  numsteps.H.msf <- (dat$param$max.time.off.tx.full.int - tx.init.time.H.msf) /
    (1 - prop.time.on.tx.H.msf)
  offon.H.msf <- rbind(offon.H.msf,
                     cbind(tx.init.time.H.msf + (1 - prop.time.on.tx.H.msf) * 1:numsteps.H.msf,
                           prop.time.on.tx.H.msf * 1:numsteps.H.msf))
  offon.H.msf <- round(offon.H.msf)
  exp.dur.chronic.H.msf <- nrow(offon.H.msf) - vl.acute.int
  exp.onset.aids.H.msf <- nrow(offon.H.msf)
  offon.last.H.msf <- offon.H.msf[nrow(offon.H.msf), ]
  offon.H.msf <- rbind(offon.H.msf,
                     matrix(c(offon.last.H.msf[1] + (1:vl.aids.int),
                              rep(offon.last.H.msf[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.H.msf <- nrow(offon.H.msf)
  offon.H.msf[, 2] <- (1:max.possible.inf.time.H.msf) - offon.H.msf[, 1]
  stage.H.msf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.H.msf, vl.aids.int))
  stage.time.H.msf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.H.msf, 1:vl.aids.int)
  
  # Stage for Hispanic Immigrant males MSF
  prop.time.on.tx.HI.msf <- dat$param$tx.reinit.HI.msf.prob /
    (dat$param$tx.halt.HI.msf.prob + dat$param$tx.reinit.HI.msf.prob)
  offon.HI.msf <- matrix(c(1:tx.init.time.HI.msf, rep(0, tx.init.time.HI.msf)),
                       nrow = tx.init.time.HI.msf)
  numsteps.HI.msf <- (dat$param$max.time.off.tx.full.int - tx.init.time.HI.msf) /
    (1 - prop.time.on.tx.HI.msf)
  offon.HI.msf <- rbind(offon.HI.msf,
                      cbind(tx.init.time.HI.msf + (1 - prop.time.on.tx.HI.msf) * 1:numsteps.HI.msf,
                            prop.time.on.tx.HI.msf * 1:numsteps.HI.msf))
  offon.HI.msf <- round(offon.HI.msf)
  exp.dur.chronic.HI.msf <- nrow(offon.HI.msf) - vl.acute.int
  exp.onset.aids.HI.msf <- nrow(offon.HI.msf)
  offon.last.HI.msf <- offon.HI.msf[nrow(offon.HI.msf), ]
  offon.HI.msf <- rbind(offon.HI.msf,
                      matrix(c(offon.last.HI.msf[1] + (1:vl.aids.int),
                               rep(offon.last.HI.msf[2], vl.aids.int)),
                             ncol = 2))
  max.possible.inf.time.HI.msf <- nrow(offon.HI.msf)
  offon.HI.msf[, 2] <- (1:max.possible.inf.time.HI.msf) - offon.HI.msf[, 1]
  stage.HI.msf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.HI.msf, vl.aids.int))
  stage.time.HI.msf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.HI.msf, 1:vl.aids.int) 
  
  # Stage for White males MSF
  prop.time.on.tx.W.msf <- dat$param$tx.reinit.W.msf.prob /
    (dat$param$tx.halt.W.msf.prob + dat$param$tx.reinit.W.msf.prob)
  offon.W.msf <- matrix(c(1:tx.init.time.W.msf, rep(0, tx.init.time.W.msf)),
                      nrow = tx.init.time.W.msf)
  numsteps.W.msf <- (dat$param$max.time.off.tx.full.int - tx.init.time.W.msf) /
    (1 - prop.time.on.tx.W.msf)
  offon.W.msf <- rbind(offon.W.msf,
                     cbind(tx.init.time.W.msf + (1 - prop.time.on.tx.W.msf) * 1:numsteps.W.msf,
                           prop.time.on.tx.W.msf * 1:numsteps.W.msf))
  offon.W.msf <- round(offon.W.msf)
  exp.dur.chronic.W.msf <- nrow(offon.W.msf) - vl.acute.int
  exp.onset.aids.W.msf <- nrow(offon.W.msf)
  offon.last.W.msf <- offon.W.msf[nrow(offon.W.msf), ]
  offon.W.msf <- rbind(offon.W.msf,
                     matrix(c(offon.last.W.msf[1] + (1:vl.aids.int),
                              rep(offon.last.W.msf[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.W.msf <- nrow(offon.W.msf)
  offon.W.msf[, 2] <- (1:max.possible.inf.time.W.msf) - offon.W.msf[, 1]
  stage.W.msf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W.msf, vl.aids.int))
  stage.time.W.msf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W.msf, 1:vl.aids.int)
  
  #Stage for Males MSM
  # Stage for Black MSM
  prop.time.on.tx.B.msm <- dat$param$tx.reinit.B.msm.prob /
    (dat$param$tx.halt.B.msm.prob + dat$param$tx.reinit.B.msm.prob)
  offon.B.msm <- matrix(c(1:tx.init.time.B.msm, rep(0, tx.init.time.B.msm)),
                      nrow = tx.init.time.B.msm)
  numsteps.B.msm <- (dat$param$max.time.off.tx.full.int - tx.init.time.B.msm) /
    (1 - prop.time.on.tx.B.msm)
  offon.B.msm <- rbind(offon.B.msm,
                     cbind(tx.init.time.B.msm + (1 - prop.time.on.tx.B.msm) * 1:numsteps.B.msm,
                           prop.time.on.tx.B.msm * 1:numsteps.B.msm))
  offon.B.msm <- round(offon.B.msm)
  exp.dur.chronic.B.msm <- nrow(offon.B.msm) - vl.acute.int
  exp.onset.aids.B.msm <- nrow(offon.B.msm)
  offon.last.B.msm <- offon.B.msm[nrow(offon.B.msm), ]
  offon.B.msm <- rbind(offon.B.msm,
                     matrix(c(offon.last.B.msm[1] + (1:vl.aids.int),
                              rep(offon.last.B.msm[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.B.msm <- nrow(offon.B.msm)
  offon.B.msm[, 2] <- (1:max.possible.inf.time.B.msm) - offon.B.msm[, 1]
  stage.B.msm <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B.msm, vl.aids.int))
  stage.time.B.msm <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B.msm, 1:vl.aids.int)
  
  # Stage for Black Immigrant MSM
  prop.time.on.tx.BI.msm <- dat$param$tx.reinit.BI.msm.prob /
    (dat$param$tx.halt.BI.msm.prob + dat$param$tx.reinit.BI.msm.prob)
  offon.BI.msm <- matrix(c(1:tx.init.time.BI.msm, rep(0, tx.init.time.BI.msm)),
                       nrow = tx.init.time.BI.msm)
  numsteps.BI.msm <- (dat$param$max.time.off.tx.full.int - tx.init.time.BI.msm) /
    (1 - prop.time.on.tx.BI.msm)
  offon.BI.msm <- rbind(offon.BI.msm,
                      cbind(tx.init.time.BI.msm + (1 - prop.time.on.tx.BI.msm) * 1:numsteps.BI.msm,
                            prop.time.on.tx.BI.msm * 1:numsteps.BI.msm))
  offon.BI.msm <- round(offon.BI.msm)
  exp.dur.chronic.BI.msm <- nrow(offon.BI.msm) - vl.acute.int
  exp.onset.aids.BI.msm <- nrow(offon.BI.msm)
  offon.last.BI.msm <- offon.BI.msm[nrow(offon.BI.msm), ]
  offon.BI.msm <- rbind(offon.BI.msm,
                      matrix(c(offon.last.BI.msm[1] + (1:vl.aids.int),
                               rep(offon.last.BI.msm[2], vl.aids.int)),
                             ncol = 2))
  max.possible.inf.time.BI.msm <- nrow(offon.BI.msm)
  offon.BI.msm[, 2] <- (1:max.possible.inf.time.BI.msm) - offon.BI.msm[, 1]
  stage.BI.msm <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.BI.msm, vl.aids.int))
  stage.time.BI.msm <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.BI.msm, 1:vl.aids.int)
  
  
  # Stage for Hispanic MSM
  prop.time.on.tx.H.msm <- dat$param$tx.reinit.H.msm.prob /
    (dat$param$tx.halt.H.msm.prob + dat$param$tx.reinit.H.msm.prob)
  offon.H.msm <- matrix(c(1:tx.init.time.H.msm, rep(0, tx.init.time.H.msm)),
                      nrow = tx.init.time.H.msm)
  numsteps.H.msm <- (dat$param$max.time.off.tx.full.int - tx.init.time.H.msm) /
    (1 - prop.time.on.tx.H.msm)
  offon.H.msm <- rbind(offon.H.msm,
                     cbind(tx.init.time.H.msm + (1 - prop.time.on.tx.H.msm) * 1:numsteps.H.msm,
                           prop.time.on.tx.H.msm * 1:numsteps.H.msm))
  offon.H.msm <- round(offon.H.msm)
  exp.dur.chronic.H.msm <- nrow(offon.H.msm) - vl.acute.int
  exp.onset.aids.H.msm <- nrow(offon.H.msm)
  offon.last.H.msm <- offon.H.msm[nrow(offon.H.msm), ]
  offon.H.msm <- rbind(offon.H.msm,
                     matrix(c(offon.last.H.msm[1] + (1:vl.aids.int),
                              rep(offon.last.H.msm[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.H.msm <- nrow(offon.H.msm)
  offon.H.msm[, 2] <- (1:max.possible.inf.time.H.msm) - offon.H.msm[, 1]
  stage.H.msm <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.H.msm, vl.aids.int))
  stage.time.H.msm <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.H.msm, 1:vl.aids.int)
  
  # Stage for Hispanic Immigrant MSM
  prop.time.on.tx.HI.msm <- dat$param$tx.reinit.HI.msm.prob /
    (dat$param$tx.halt.HI.msm.prob + dat$param$tx.reinit.HI.msm.prob)
  offon.HI.msm <- matrix(c(1:tx.init.time.HI.msm, rep(0, tx.init.time.HI.msm)),
                       nrow = tx.init.time.HI.msm)
  numsteps.HI.msm <- (dat$param$max.time.off.tx.full.int - tx.init.time.HI.msm) /
    (1 - prop.time.on.tx.HI.msm)
  offon.HI.msm <- rbind(offon.HI.msm,
                      cbind(tx.init.time.HI.msm + (1 - prop.time.on.tx.HI.msm) * 1:numsteps.HI.msm,
                            prop.time.on.tx.HI.msm * 1:numsteps.HI.msm))
  offon.HI.msm <- round(offon.HI.msm)
  exp.dur.chronic.HI.msm <- nrow(offon.HI.msm) - vl.acute.int
  exp.onset.aids.HI.msm <- nrow(offon.HI.msm)
  offon.last.HI.msm <- offon.HI.msm[nrow(offon.HI.msm), ]
  offon.HI.msm <- rbind(offon.HI.msm,
                      matrix(c(offon.last.HI.msm[1] + (1:vl.aids.int),
                               rep(offon.last.HI.msm[2], vl.aids.int)),
                             ncol = 2))
  max.possible.inf.time.HI.msm <- nrow(offon.HI.msm)
  offon.HI.msm[, 2] <- (1:max.possible.inf.time.HI.msm) - offon.HI.msm[, 1]
  stage.HI.msm <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.HI.msm, vl.aids.int))
  stage.time.HI.msm <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.HI.msm, 1:vl.aids.int) 
  
  # Stage for White MSM
  prop.time.on.tx.W.msm <- dat$param$tx.reinit.W.msm.prob /
    (dat$param$tx.halt.W.msm.prob + dat$param$tx.reinit.W.msm.prob)
  offon.W.msm <- matrix(c(1:tx.init.time.W.msm, rep(0, tx.init.time.W.msm)),
                      nrow = tx.init.time.W.msm)
  numsteps.W.msm <- (dat$param$max.time.off.tx.full.int - tx.init.time.W.msm) /
    (1 - prop.time.on.tx.W.msm)
  offon.W.msm <- rbind(offon.W.msm,
                     cbind(tx.init.time.W.msm + (1 - prop.time.on.tx.W.msm) * 1:numsteps.W.msm,
                           prop.time.on.tx.W.msm * 1:numsteps.W.msm))
  offon.W.msm <- round(offon.W.msm)
  exp.dur.chronic.W.msm <- nrow(offon.W.msm) - vl.acute.int
  exp.onset.aids.W.msm <- nrow(offon.W.msm)
  offon.last.W.msm <- offon.W.msm[nrow(offon.W.msm), ]
  offon.W.msm <- rbind(offon.W.msm,
                     matrix(c(offon.last.W.msm[1] + (1:vl.aids.int),
                              rep(offon.last.W.msm[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.W.msm <- nrow(offon.W.msm)
  offon.W.msm[, 2] <- (1:max.possible.inf.time.W.msm) - offon.W.msm[, 1]
  stage.W.msm <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W.msm, vl.aids.int))
  stage.time.W.msm <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W.msm, 1:vl.aids.int)
  
  #Stage for Males MSMF
  # Stage for Black MSMF
  prop.time.on.tx.B.msmf <- dat$param$tx.reinit.B.msmf.prob /
    (dat$param$tx.halt.B.msmf.prob + dat$param$tx.reinit.B.msmf.prob)
  offon.B.msmf <- matrix(c(1:tx.init.time.B.msmf, rep(0, tx.init.time.B.msmf)),
                         nrow = tx.init.time.B.msmf)
  numsteps.B.msmf <- (dat$param$max.time.off.tx.full.int - tx.init.time.B.msmf) /
    (1 - prop.time.on.tx.B.msmf)
  offon.B.msmf <- rbind(offon.B.msmf,
                        cbind(tx.init.time.B.msmf + (1 - prop.time.on.tx.B.msmf) * 1:numsteps.B.msmf,
                              prop.time.on.tx.B.msmf * 1:numsteps.B.msmf))
  offon.B.msmf <- round(offon.B.msmf)
  exp.dur.chronic.B.msmf <- nrow(offon.B.msmf) - vl.acute.int
  exp.onset.aids.B.msmf <- nrow(offon.B.msmf)
  offon.last.B.msmf <- offon.B.msmf[nrow(offon.B.msmf), ]
  offon.B.msmf <- rbind(offon.B.msmf,
                        matrix(c(offon.last.B.msmf[1] + (1:vl.aids.int),
                                 rep(offon.last.B.msmf[2], vl.aids.int)),
                               ncol = 2))
  max.possible.inf.time.B.msmf <- nrow(offon.B.msmf)
  offon.B.msmf[, 2] <- (1:max.possible.inf.time.B.msmf) - offon.B.msmf[, 1]
  stage.B.msmf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B.msmf, vl.aids.int))
  stage.time.B.msmf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B.msmf, 1:vl.aids.int)
  
  # Stage for Black Immigrant MSMF
  prop.time.on.tx.BI.msmf <- dat$param$tx.reinit.BI.msmf.prob /
    (dat$param$tx.halt.BI.msmf.prob + dat$param$tx.reinit.BI.msmf.prob)
  offon.BI.msmf <- matrix(c(1:tx.init.time.BI.msmf, rep(0, tx.init.time.BI.msmf)),
                          nrow = tx.init.time.BI.msmf)
  numsteps.BI.msmf <- (dat$param$max.time.off.tx.full.int - tx.init.time.BI.msmf) /
    (1 - prop.time.on.tx.BI.msmf)
  offon.BI.msmf <- rbind(offon.BI.msmf,
                         cbind(tx.init.time.BI.msmf + (1 - prop.time.on.tx.BI.msmf) * 1:numsteps.BI.msmf,
                               prop.time.on.tx.BI.msmf * 1:numsteps.BI.msmf))
  offon.BI.msmf <- round(offon.BI.msmf)
  exp.dur.chronic.BI.msmf <- nrow(offon.BI.msmf) - vl.acute.int
  exp.onset.aids.BI.msmf <- nrow(offon.BI.msmf)
  offon.last.BI.msmf <- offon.BI.msmf[nrow(offon.BI.msmf), ]
  offon.BI.msmf <- rbind(offon.BI.msmf,
                         matrix(c(offon.last.BI.msmf[1] + (1:vl.aids.int),
                                  rep(offon.last.BI.msmf[2], vl.aids.int)),
                                ncol = 2))
  max.possible.inf.time.BI.msmf <- nrow(offon.BI.msmf)
  offon.BI.msmf[, 2] <- (1:max.possible.inf.time.BI.msmf) - offon.BI.msmf[, 1]
  stage.BI.msmf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.BI.msmf, vl.aids.int))
  stage.time.BI.msmf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.BI.msmf, 1:vl.aids.int)
  
  
  # Stage for Hispanic MSMF
  prop.time.on.tx.H.msmf <- dat$param$tx.reinit.H.msmf.prob /
    (dat$param$tx.halt.H.msmf.prob + dat$param$tx.reinit.H.msmf.prob)
  offon.H.msmf <- matrix(c(1:tx.init.time.H.msmf, rep(0, tx.init.time.H.msmf)),
                         nrow = tx.init.time.H.msmf)
  numsteps.H.msmf <- (dat$param$max.time.off.tx.full.int - tx.init.time.H.msmf) /
    (1 - prop.time.on.tx.H.msmf)
  offon.H.msmf <- rbind(offon.H.msmf,
                        cbind(tx.init.time.H.msmf + (1 - prop.time.on.tx.H.msmf) * 1:numsteps.H.msmf,
                              prop.time.on.tx.H.msmf * 1:numsteps.H.msmf))
  offon.H.msmf <- round(offon.H.msmf)
  exp.dur.chronic.H.msmf <- nrow(offon.H.msmf) - vl.acute.int
  exp.onset.aids.H.msmf <- nrow(offon.H.msmf)
  offon.last.H.msmf <- offon.H.msmf[nrow(offon.H.msmf), ]
  offon.H.msmf <- rbind(offon.H.msmf,
                        matrix(c(offon.last.H.msmf[1] + (1:vl.aids.int),
                                 rep(offon.last.H.msmf[2], vl.aids.int)),
                               ncol = 2))
  max.possible.inf.time.H.msmf <- nrow(offon.H.msmf)
  offon.H.msmf[, 2] <- (1:max.possible.inf.time.H.msmf) - offon.H.msmf[, 1]
  stage.H.msmf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.H.msmf, vl.aids.int))
  stage.time.H.msmf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.H.msmf, 1:vl.aids.int)
  
  # Stage for Hispanic Immigrant MSMF
  prop.time.on.tx.HI.msmf <- dat$param$tx.reinit.HI.msmf.prob /
    (dat$param$tx.halt.HI.msmf.prob + dat$param$tx.reinit.HI.msmf.prob)
  offon.HI.msmf <- matrix(c(1:tx.init.time.HI.msmf, rep(0, tx.init.time.HI.msmf)),
                          nrow = tx.init.time.HI.msmf)
  numsteps.HI.msmf <- (dat$param$max.time.off.tx.full.int - tx.init.time.HI.msmf) /
    (1 - prop.time.on.tx.HI.msmf)
  offon.HI.msmf <- rbind(offon.HI.msmf,
                         cbind(tx.init.time.HI.msmf + (1 - prop.time.on.tx.HI.msmf) * 1:numsteps.HI.msmf,
                               prop.time.on.tx.HI.msmf * 1:numsteps.HI.msmf))
  offon.HI.msmf <- round(offon.HI.msmf)
  exp.dur.chronic.HI.msmf <- nrow(offon.HI.msmf) - vl.acute.int
  exp.onset.aids.HI.msmf <- nrow(offon.HI.msmf)
  offon.last.HI.msmf <- offon.HI.msmf[nrow(offon.HI.msmf), ]
  offon.HI.msmf <- rbind(offon.HI.msmf,
                         matrix(c(offon.last.HI.msmf[1] + (1:vl.aids.int),
                                  rep(offon.last.HI.msmf[2], vl.aids.int)),
                                ncol = 2))
  max.possible.inf.time.HI.msmf <- nrow(offon.HI.msmf)
  offon.HI.msmf[, 2] <- (1:max.possible.inf.time.HI.msmf) - offon.HI.msmf[, 1]
  stage.HI.msmf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.HI.msmf, vl.aids.int))
  stage.time.HI.msmf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.HI.msmf, 1:vl.aids.int) 
  
  # Stage for White MSMF
  prop.time.on.tx.W.msmf <- dat$param$tx.reinit.W.msmf.prob /
    (dat$param$tx.halt.W.msmf.prob + dat$param$tx.reinit.W.msmf.prob)
  offon.W.msmf <- matrix(c(1:tx.init.time.W.msmf, rep(0, tx.init.time.W.msmf)),
                         nrow = tx.init.time.W.msmf)
  numsteps.W.msmf <- (dat$param$max.time.off.tx.full.int - tx.init.time.W.msmf) /
    (1 - prop.time.on.tx.W.msmf)
  offon.W.msmf <- rbind(offon.W.msmf,
                        cbind(tx.init.time.W.msmf + (1 - prop.time.on.tx.W.msmf) * 1:numsteps.W.msmf,
                              prop.time.on.tx.W.msmf * 1:numsteps.W.msmf))
  offon.W.msmf <- round(offon.W.msmf)
  exp.dur.chronic.W.msmf <- nrow(offon.W.msmf) - vl.acute.int
  exp.onset.aids.W.msmf <- nrow(offon.W.msmf)
  offon.last.W.msmf <- offon.W.msmf[nrow(offon.W.msmf), ]
  offon.W.msmf <- rbind(offon.W.msmf,
                        matrix(c(offon.last.W.msmf[1] + (1:vl.aids.int),
                                 rep(offon.last.W.msmf[2], vl.aids.int)),
                               ncol = 2))
  max.possible.inf.time.W.msmf <- nrow(offon.W.msmf)
  offon.W.msmf[, 2] <- (1:max.possible.inf.time.W.msmf) - offon.W.msmf[, 1]
  stage.W.msmf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W.msmf, vl.aids.int))
  stage.time.W.msmf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W.msmf, 1:vl.aids.int)
  
  
    ##Viral load females
  # Vl for Black females
  selected <- which(status == 1 & tt.traj == 4 & race == "B" & sex== "F")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B.f)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B.f[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B.f[time.since.inf, 1]
  stage[selected] <- stage.B.f[time.since.inf]
  stage.time[selected] <- stage.time.B.f[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B.f)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.B.f) * (vlsp) +
                  (time.since.inf > exp.onset.aids.B.f) *
                  (vlsp + (time.since.inf - exp.onset.aids.B.f) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  # Vl for Black Immigrant females
  selected <- which(status == 1 & tt.traj == 4 & race == "BI" & sex== "F")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.BI.f)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.BI.f[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.BI.f[time.since.inf, 1]
  stage[selected] <- stage.BI.f[time.since.inf]
  stage.time[selected] <- stage.time.BI.f[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.BI.f)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.BI.f) * (vlsp) +
    (time.since.inf > exp.onset.aids.BI.f) *
    (vlsp + (time.since.inf - exp.onset.aids.BI.f) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  # Vl for Hispanic females
  selected <- which(status == 1 & tt.traj == 4 & race == "H" & sex== "F")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.H.f)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.H.f[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.H.f[time.since.inf, 1]
  stage[selected] <- stage.H.f[time.since.inf]
  stage.time[selected] <- stage.time.H.f[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.H.f)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.H.f) * (vlsp) +
    (time.since.inf > exp.onset.aids.H.f) *
    (vlsp + (time.since.inf - exp.onset.aids.H.f) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  # Vl for Hispanic Immigrant females
  selected <- which(status == 1 & tt.traj == 4 & race == "HI" & sex== "F")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.HI.f)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.HI.f[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.HI.f[time.since.inf, 1]
  stage[selected] <- stage.HI.f[time.since.inf]
  stage.time[selected] <- stage.time.HI.f[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.HI.f)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.HI.f) * (vlsp) +
    (time.since.inf > exp.onset.aids.HI.f) *
    (vlsp + (time.since.inf - exp.onset.aids.HI.f) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  

  # VL for White females
  selected <- which(status == 1 & tt.traj == 4 & race == "W" & sex== "F")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W.f)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W.f[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W.f[time.since.inf, 1]
  stage[selected] <- stage.W.f[time.since.inf]
  stage.time[selected] <- stage.time.W.f[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W.f)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.W.f) * (vlsp) +
                  (time.since.inf > exp.onset.aids.W.f) *
                  (vlsp + (time.since.inf - exp.onset.aids.W.f) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp

  #Viral Load Males Het MSF.
  # Vl for Black males MSF
  selected <- which(status == 1 & tt.traj == 4 & race == "B" & sex.ident=="msf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B.msf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B.msf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B.msf[time.since.inf, 1]
  stage[selected] <- stage.B.msf[time.since.inf]
  stage.time[selected] <- stage.time.B.msf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B.msf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.B.msf) * (vlsp) +
    (time.since.inf > exp.onset.aids.B.msf) *
    (vlsp + (time.since.inf - exp.onset.aids.B.msf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
   # Vl for Black Immigrant males MSF
  selected <- which(status == 1 & tt.traj == 4 & race == "BI" & sex.ident=="msf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.BI.msf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.BI.msf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.BI.msf[time.since.inf, 1]
  stage[selected] <- stage.BI.msf[time.since.inf]
  stage.time[selected] <- stage.time.BI.msf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.BI.msf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.BI.msf) * (vlsp) +
    (time.since.inf > exp.onset.aids.BI.msf) *
    (vlsp + (time.since.inf - exp.onset.aids.BI.msf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  # Vl for Hispanic males MSF
  selected <- which(status == 1 & tt.traj == 4 & race == "H" & sex.ident=="msf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.H.msf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.H.msf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.H.msf[time.since.inf, 1]
  stage[selected] <- stage.H.msf[time.since.inf]
  stage.time[selected] <- stage.time.H.msf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.H.msf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.H.msf) * (vlsp) +
    (time.since.inf > exp.onset.aids.H.msf) *
    (vlsp + (time.since.inf - exp.onset.aids.H.msf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  # Vl for Hispanic Immigrant males
  selected <- which(status == 1 & tt.traj == 4 & race == "HI" & sex.ident=="msf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.HI.msf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.HI.msf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.HI.msf[time.since.inf, 1]
  stage[selected] <- stage.HI.msf[time.since.inf]
  stage.time[selected] <- stage.time.HI.msf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.HI.msf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.HI.msf) * (vlsp) +
    (time.since.inf > exp.onset.aids.HI.msf) *
    (vlsp + (time.since.inf - exp.onset.aids.HI.msf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  
  # VL for White males
  selected <- which(status == 1 & tt.traj == 4 & race == "W" & sex.ident=="msf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W.msf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W.msf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W.msf[time.since.inf, 1]
  stage[selected] <- stage.W.msf[time.since.inf]
  stage.time[selected] <- stage.time.W.msf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W.msf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.W.msf) * (vlsp) +
    (time.since.inf > exp.onset.aids.W.msf) *
    (vlsp + (time.since.inf - exp.onset.aids.W.msf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  
  
  #Viral Load Males MSM
  # Vl for Black MSM
  selected <- which(status == 1 & tt.traj == 4 & race == "B" & sex.ident=="msm")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B.msm)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B.msm[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B.msm[time.since.inf, 1]
  stage[selected] <- stage.B.msm[time.since.inf]
  stage.time[selected] <- stage.time.B.msm[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B.msm)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.B.msm) * (vlsp) +
    (time.since.inf > exp.onset.aids.B.msm) *
    (vlsp + (time.since.inf - exp.onset.aids.B.msm) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  # Vl for Black Immigrant MSM
  selected <- which(status == 1 & tt.traj == 4 & race == "BI" & sex.ident=="msm")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.BI.msm)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.BI.msm[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.BI.msm[time.since.inf, 1]
  stage[selected] <- stage.BI.msm[time.since.inf]
  stage.time[selected] <- stage.time.BI.msm[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.BI.msm)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.BI.msm) * (vlsp) +
    (time.since.inf > exp.onset.aids.BI.msm) *
    (vlsp + (time.since.inf - exp.onset.aids.BI.msm) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  # Vl for Hispanic MSM
  selected <- which(status == 1 & tt.traj == 4 & race == "H" & sex.ident=="msm")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.H.msm)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.H.msm[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.H.msm[time.since.inf, 1]
  stage[selected] <- stage.H.msm[time.since.inf]
  stage.time[selected] <- stage.time.H.msm[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.H.msm)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.H.msm) * (vlsp) +
    (time.since.inf > exp.onset.aids.H.msm) *
    (vlsp + (time.since.inf - exp.onset.aids.H.msm) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  # Vl for Hispanic Immigrant MSM
  selected <- which(status == 1 & tt.traj == 4 & race == "HI" & sex.ident=="msm")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.HI.msm)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.HI.msm[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.HI.msm[time.since.inf, 1]
  stage[selected] <- stage.HI.msm[time.since.inf]
  stage.time[selected] <- stage.time.HI.msm[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.HI.msm)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.HI.msm) * (vlsp) +
    (time.since.inf > exp.onset.aids.HI.msm) *
    (vlsp + (time.since.inf - exp.onset.aids.HI.msm) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  
  # VL for White MSM
  selected <- which(status == 1 & tt.traj == 4 & race == "W" & sex.ident=="msm")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W.msm)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W.msm[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W.msm[time.since.inf, 1]
  stage[selected] <- stage.W.msm[time.since.inf]
  stage.time[selected] <- stage.time.W.msm[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W.msm)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.W.msm) * (vlsp) +
    (time.since.inf > exp.onset.aids.W.msm) *
    (vlsp + (time.since.inf - exp.onset.aids.W.msm) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  
  #Viral Load Males MSMF
  # Vl for Black MSMF
  selected <- which(status == 1 & tt.traj == 4 & race == "B" & sex.ident=="msmf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B.msmf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B.msmf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B.msmf[time.since.inf, 1]
  stage[selected] <- stage.B.msmf[time.since.inf]
  stage.time[selected] <- stage.time.B.msmf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B.msmf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.B.msmf) * (vlsp) +
    (time.since.inf > exp.onset.aids.B.msmf) *
    (vlsp + (time.since.inf - exp.onset.aids.B.msmf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  # Vl for Black Immigrant MSMF
  selected <- which(status == 1 & tt.traj == 4 & race == "BI" & sex.ident=="msmf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.BI.msmf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.BI.msmf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.BI.msmf[time.since.inf, 1]
  stage[selected] <- stage.BI.msmf[time.since.inf]
  stage.time[selected] <- stage.time.BI.msmf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.BI.msmf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.BI.msmf) * (vlsp) +
    (time.since.inf > exp.onset.aids.BI.msmf) *
    (vlsp + (time.since.inf - exp.onset.aids.BI.msmf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  # Vl for Hispanic MSMF
  selected <- which(status == 1 & tt.traj == 4 & race == "H" & sex.ident=="msmf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.H.msmf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.H.msmf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.H.msmf[time.since.inf, 1]
  stage[selected] <- stage.H.msmf[time.since.inf]
  stage.time[selected] <- stage.time.H.msmf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.H.msmf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.H.msmf) * (vlsp) +
    (time.since.inf > exp.onset.aids.H.msmf) *
    (vlsp + (time.since.inf - exp.onset.aids.H.msmf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  # Vl for Hispanic Immigrant MSMF
  selected <- which(status == 1 & tt.traj == 4 & race == "HI" & sex.ident=="msmf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.HI.msmf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.HI.msmf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.HI.msmf[time.since.inf, 1]
  stage[selected] <- stage.HI.msmf[time.since.inf]
  stage.time[selected] <- stage.time.HI.msmf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.HI.msmf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.HI.msmf) * (vlsp) +
    (time.since.inf > exp.onset.aids.HI.msmf) *
    (vlsp + (time.since.inf - exp.onset.aids.HI.msmf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp
  
  
  # VL for White MSMF
  selected <- which(status == 1 & tt.traj == 4 & race == "W" & sex.ident=="msmf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W.msmf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W.msmf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W.msmf[time.since.inf, 1]
  stage[selected] <- stage.W.msmf[time.since.inf]
  stage.time[selected] <- stage.time.W.msmf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W.msmf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.W.msmf) * (vlsp) +
    (time.since.inf > exp.onset.aids.W.msmf) *
    (vlsp + (time.since.inf - exp.onset.aids.W.msmf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.full.supp 
  
  # Diagnosis
  selected <- which(status == 1 & tt.traj == 4)
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                      min = 0,
                      max = dat$param$mean.test.B.f.int * (race[selected] == "B" & sex[selected] == "F") +
                           dat$param$mean.test.BI.f.int * (race[selected] == "BI" & sex[selected] == "F") +
                           dat$param$mean.test.H.f.int * (race[selected] == "H" & sex[selected] == "F") +
                           dat$param$mean.test.HI.f.int * (race[selected] == "HI" & sex[selected] == "F") +
                           dat$param$mean.test.W.f.int * (race[selected] == "W" & sex[selected] == "F") +
                           dat$param$mean.test.B.msf.int * (race[selected] == "B" & sex.ident[selected] == "msf") +
                           dat$param$mean.test.BI.msf.int * (race[selected] == "BI" & sex.ident[selected] == "msf") +
                           dat$param$mean.test.H.msf.int * (race[selected] == "H" & sex.ident[selected] == "msf") +
                           dat$param$mean.test.HI.msf.int * (race[selected] == "HI" & sex.ident[selected] == "msf") +
                           dat$param$mean.test.W.msf.int * (race[selected] == "W" & sex.ident[selected] == "msf")+
                           dat$param$mean.test.B.msm.int * (race[selected] == "B" & sex.ident[selected] == "msm") +
                           dat$param$mean.test.BI.msm.int * (race[selected] == "BI" & sex.ident[selected] == "msm") +
                           dat$param$mean.test.H.msm.int * (race[selected] == "H" & sex.ident[selected] == "msm") +
                           dat$param$mean.test.HI.msm.int * (race[selected] == "HI" & sex.ident[selected] == "msm") +
                           dat$param$mean.test.W.msm.int * (race[selected] == "W" & sex.ident[selected] == "msm") +
                           dat$param$mean.test.B.msmf.int * (race[selected] == "B" & sex.ident[selected] == "msmf") +
                           dat$param$mean.test.BI.msmf.int * (race[selected] == "BI" & sex.ident[selected] == "msmf") +
                           dat$param$mean.test.H.msmf.int * (race[selected] == "H" & sex.ident[selected] == "msmf") +
                           dat$param$mean.test.HI.msmf.int * (race[selected] == "HI" & sex.ident[selected] == "msmf") +
                           dat$param$mean.test.W.msmf.int * (race[selected] == "W" & sex.ident[selected] == "msmf")))
  }
  if (dat$param$testing.pattern == "memoryless") {
    ttntest <- rgeom(length(selected),
                     1 / (dat$param$mean.test.B.f.int * (race[selected] == "B" & sex[selected] == "F") +
                          dat$param$mean.test.BI.f.int * (race[selected] == "BI" & sex[selected] == "F") +
                          dat$param$mean.test.H.f.int * (race[selected] == "H" & sex[selected] == "F") +
                          dat$param$mean.test.HI.f.int * (race[selected] == "HI" & sex[selected] == "F") +
                          dat$param$mean.test.W.f.int * (race[selected] == "W" & sex[selected] == "F") +
                          dat$param$mean.test.B.msf.int * (race[selected] == "B" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.BI.msf.int * (race[selected] == "BI" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.H.msf.int * (race[selected] == "H" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.HI.msf.int * (race[selected] == "HI" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.W.msf.int * (race[selected] == "W" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.B.msm.int * (race[selected] == "B" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.BI.msm.int * (race[selected] == "BI" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.H.msm.int * (race[selected] == "H" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.HI.msm.int * (race[selected] == "HI" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.W.msm.int * (race[selected] == "W" & sex.ident[selected] == "msm") + 
                          dat$param$mean.test.B.msmf.int * (race[selected] == "B" & sex.ident[selected] == "msmf") +
                          dat$param$mean.test.BI.msmf.int * (race[selected] == "BI" & sex.ident[selected] == "msmf") +
                          dat$param$mean.test.H.msmf.int * (race[selected] == "H" & sex.ident[selected] == "msmf") +
                          dat$param$mean.test.HI.msmf.int * (race[selected] == "HI" & sex.ident[selected] == "msmf") +
                          dat$param$mean.test.W.msmf.int * (race[selected] == "W" & sex.ident[selected] == "msmf")))
  }

  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
                           -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]
  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1
  diag.status[selected][cum.time.on.tx[selected] > 0] <- 1
  last.neg.test[selected][cum.time.on.tx[selected] > 0] <- NA


  ### Part adherent type

  # Create set of expected values for (cum.time.off.tx,cum.time.on.tx)

  #Black females.
  prop.time.on.tx.B.f <- dat$param$tx.reinit.B.f.prob /
                       (dat$param$tx.halt.B.f.prob + dat$param$tx.reinit.B.f.prob)
  offon.B.f <- matrix(c(1:tx.init.time.B.f, rep(0, tx.init.time.B.f)),
                    nrow = tx.init.time.B.f)
  while (offon.B.f[nrow(offon.B.f), 1] / dat$param$max.time.off.tx.part.int +
         offon.B.f[nrow(offon.B.f), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.B.f <- rbind(offon.B.f,
                     offon.B.f[nrow(offon.B.f), ] + c(1 - prop.time.on.tx.B.f,
                                                      prop.time.on.tx.B.f))
  }
  offon.B.f <- round(offon.B.f)
  exp.dur.chronic.B.f <- nrow(offon.B.f) - vl.acute.int
  exp.onset.aids.B.f <- nrow(offon.B.f)
  offon.last.B.f <- offon.B.f[nrow(offon.B.f), ]
  offon.B.f <- rbind(offon.B.f,
                   matrix(c(offon.last.B.f[1] + (1:vl.aids.int),
                            rep(offon.last.B.f[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.B.f <- nrow(offon.B.f)
  offon.B.f[, 2] <- (1:max.possible.inf.time.B.f) - offon.B.f[, 1]
  stage.B.f <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B.f, vl.aids.int))
  stage.time.B.f <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B.f, 1:vl.aids.int)

  #Black Immigrant females.
  prop.time.on.tx.BI.f <- dat$param$tx.reinit.BI.f.prob /
    (dat$param$tx.halt.BI.f.prob + dat$param$tx.reinit.BI.f.prob)
  offon.BI.f <- matrix(c(1:tx.init.time.BI.f, rep(0, tx.init.time.BI.f)),
                      nrow = tx.init.time.BI.f)
  while (offon.BI.f[nrow(offon.BI.f), 1] / dat$param$max.time.off.tx.part.int +
         offon.BI.f[nrow(offon.BI.f), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.BI.f <- rbind(offon.BI.f,
                       offon.BI.f[nrow(offon.BI.f), ] + c(1 - prop.time.on.tx.BI.f,
                                                        prop.time.on.tx.BI.f))
  }
  offon.BI.f <- round(offon.BI.f)
  exp.dur.chronic.BI.f <- nrow(offon.BI.f) - vl.acute.int
  exp.onset.aids.BI.f <- nrow(offon.BI.f)
  offon.last.BI.f <- offon.BI.f[nrow(offon.BI.f), ]
  offon.BI.f <- rbind(offon.BI.f,
                     matrix(c(offon.last.BI.f[1] + (1:vl.aids.int),
                              rep(offon.last.BI.f[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.BI.f <- nrow(offon.BI.f)
  offon.BI.f[, 2] <- (1:max.possible.inf.time.BI.f) - offon.BI.f[, 1]
  stage.BI.f <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.BI.f, vl.aids.int))
  stage.time.BI.f <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.BI.f, 1:vl.aids.int)
  
  #Hispanic females.
  prop.time.on.tx.H.f <- dat$param$tx.reinit.H.f.prob /
    (dat$param$tx.halt.H.f.prob + dat$param$tx.reinit.H.f.prob)
  offon.H.f <- matrix(c(1:tx.init.time.H.f, rep(0, tx.init.time.H.f)),
                      nrow = tx.init.time.H.f)
  while (offon.H.f[nrow(offon.H.f), 1] / dat$param$max.time.off.tx.part.int +
         offon.H.f[nrow(offon.H.f), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.H.f <- rbind(offon.H.f,
                       offon.H.f[nrow(offon.H.f), ] + c(1 - prop.time.on.tx.H.f,
                                                        prop.time.on.tx.H.f))
  }
  offon.H.f <- round(offon.H.f)
  exp.dur.chronic.H.f <- nrow(offon.H.f) - vl.acute.int
  exp.onset.aids.H.f <- nrow(offon.H.f)
  offon.last.H.f <- offon.H.f[nrow(offon.H.f), ]
  offon.H.f <- rbind(offon.H.f,
                     matrix(c(offon.last.H.f[1] + (1:vl.aids.int),
                              rep(offon.last.H.f[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.H.f <- nrow(offon.H.f)
  offon.H.f[, 2] <- (1:max.possible.inf.time.H.f) - offon.H.f[, 1]
  stage.H.f <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.H.f, vl.aids.int))
  stage.time.H.f <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.H.f, 1:vl.aids.int)
  
  #Hispanic Immigrant females.
  prop.time.on.tx.HI.f <- dat$param$tx.reinit.HI.f.prob /
    (dat$param$tx.halt.HI.f.prob + dat$param$tx.reinit.HI.f.prob)
  offon.HI.f <- matrix(c(1:tx.init.time.HI.f, rep(0, tx.init.time.HI.f)),
                       nrow = tx.init.time.HI.f)
  while (offon.HI.f[nrow(offon.HI.f), 1] / dat$param$max.time.off.tx.part.int +
         offon.HI.f[nrow(offon.HI.f), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.HI.f <- rbind(offon.HI.f,
                        offon.HI.f[nrow(offon.HI.f), ] + c(1 - prop.time.on.tx.HI.f,
                                                           prop.time.on.tx.HI.f))
  }
  offon.HI.f <- round(offon.HI.f)
  exp.dur.chronic.HI.f <- nrow(offon.HI.f) - vl.acute.int
  exp.onset.aids.HI.f <- nrow(offon.HI.f)
  offon.last.HI.f <- offon.HI.f[nrow(offon.HI.f), ]
  offon.HI.f <- rbind(offon.HI.f,
                      matrix(c(offon.last.HI.f[1] + (1:vl.aids.int),
                               rep(offon.last.HI.f[2], vl.aids.int)),
                             ncol = 2))
  max.possible.inf.time.HI.f <- nrow(offon.HI.f)
  offon.HI.f[, 2] <- (1:max.possible.inf.time.HI.f) - offon.HI.f[, 1]
  stage.HI.f <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.HI.f, vl.aids.int))
  stage.time.HI.f <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.HI.f, 1:vl.aids.int)
  
  
  #White females
  prop.time.on.tx.W.f <- dat$param$tx.reinit.W.f.prob /
                       (dat$param$tx.halt.W.f.prob + dat$param$tx.reinit.W.f.prob)
  offon.W.f <- matrix(c(1:tx.init.time.W.f, rep(0, tx.init.time.W.f)),
                    nrow = tx.init.time.W.f)

  while (offon.W.f[nrow(offon.W.f), 1] / dat$param$max.time.off.tx.part.int +
         offon.W.f[nrow(offon.W.f), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.W.f <- rbind(offon.W.f,
                     offon.W.f[nrow(offon.W.f), ] + c(1 - prop.time.on.tx.W.f,
                                                  prop.time.on.tx.W.f))
  }
  offon.W.f <- round(offon.W.f)
  exp.dur.chronic.W.f <- nrow(offon.W.f) - vl.acute.int
  exp.onset.aids.W.f <- nrow(offon.W.f)
  offon.last.W.f <- offon.W.f[nrow(offon.W.f), ]
  offon.W.f <- rbind(offon.W.f,
                   matrix(c(offon.last.W.f[1] + (1:vl.aids.int),
                            rep(offon.last.W.f[2], vl.aids.int)),
                          ncol = 2))
  max.possible.inf.time.W.f <- nrow(offon.W.f)
  offon.W.f[, 2] <- (1:max.possible.inf.time.W.f) - offon.W.f[, 1]
  stage.W.f <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W.f, vl.aids.int))
  stage.time.W.f <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W.f, 1:vl.aids.int)

  
  
  #Black males MSF.
  prop.time.on.tx.B.msf <- dat$param$tx.reinit.B.msf.prob /
    (dat$param$tx.halt.B.msf.prob + dat$param$tx.reinit.B.msf.prob)
  offon.B.msf <- matrix(c(1:tx.init.time.B.msf, rep(0, tx.init.time.B.msf)),
                      nrow = tx.init.time.B.msf)
  while (offon.B.msf[nrow(offon.B.msf), 1] / dat$param$max.time.off.tx.part.int +
         offon.B.msf[nrow(offon.B.msf), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.B.msf <- rbind(offon.B.msf,
                       offon.B.msf[nrow(offon.B.msf), ] + c(1 - prop.time.on.tx.B.msf,
                                                        prop.time.on.tx.B.msf))
  }
  offon.B.msf <- round(offon.B.msf)
  exp.dur.chronic.B.msf <- nrow(offon.B.msf) - vl.acute.int
  exp.onset.aids.B.msf <- nrow(offon.B.msf)
  offon.last.B.msf <- offon.B.msf[nrow(offon.B.msf), ]
  offon.B.msf <- rbind(offon.B.msf,
                     matrix(c(offon.last.B.msf[1] + (1:vl.aids.int),
                              rep(offon.last.B.msf[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.B.msf <- nrow(offon.B.msf)
  offon.B.msf[, 2] <- (1:max.possible.inf.time.B.msf) - offon.B.msf[, 1]
  stage.B.msf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B.msf, vl.aids.int))
  stage.time.B.msf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B.msf, 1:vl.aids.int)
  
  #Black Immigrant males.
  prop.time.on.tx.BI.msf <- dat$param$tx.reinit.BI.msf.prob /
    (dat$param$tx.halt.BI.msf.prob + dat$param$tx.reinit.BI.msf.prob)
  offon.BI.msf <- matrix(c(1:tx.init.time.BI.msf, rep(0, tx.init.time.BI.msf)),
                       nrow = tx.init.time.BI.msf)
  while (offon.BI.msf[nrow(offon.BI.msf), 1] / dat$param$max.time.off.tx.part.int +
         offon.BI.msf[nrow(offon.BI.msf), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.BI.msf <- rbind(offon.BI.msf,
                        offon.BI.msf[nrow(offon.BI.msf), ] + c(1 - prop.time.on.tx.BI.msf,
                                                           prop.time.on.tx.BI.msf))
  }
  offon.BI.msf <- round(offon.BI.msf)
  exp.dur.chronic.BI.msf <- nrow(offon.BI.msf) - vl.acute.int
  exp.onset.aids.BI.msf <- nrow(offon.BI.msf)
  offon.last.BI.msf <- offon.BI.msf[nrow(offon.BI.msf), ]
  offon.BI.msf <- rbind(offon.BI.msf,
                      matrix(c(offon.last.BI.msf[1] + (1:vl.aids.int),
                               rep(offon.last.BI.msf[2], vl.aids.int)),
                             ncol = 2))
  max.possible.inf.time.BI.msf <- nrow(offon.BI.msf)
  offon.BI.msf[, 2] <- (1:max.possible.inf.time.BI.msf) - offon.BI.msf[, 1]
  stage.BI.msf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.BI.msf, vl.aids.int))
  stage.time.BI.msf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.BI.msf, 1:vl.aids.int)
  
  #Hispanic males.
  prop.time.on.tx.H.msf <- dat$param$tx.reinit.H.msf.prob /
    (dat$param$tx.halt.H.msf.prob + dat$param$tx.reinit.H.msf.prob)
  offon.H.msf <- matrix(c(1:tx.init.time.H.msf, rep(0, tx.init.time.H.msf)),
                      nrow = tx.init.time.H.msf)
  while (offon.H.msf[nrow(offon.H.msf), 1] / dat$param$max.time.off.tx.part.int +
         offon.H.msf[nrow(offon.H.msf), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.H.msf <- rbind(offon.H.msf,
                       offon.H.msf[nrow(offon.H.msf), ] + c(1 - prop.time.on.tx.H.msf,
                                                        prop.time.on.tx.H.msf))
  }
  offon.H.msf <- round(offon.H.msf)
  exp.dur.chronic.H.msf <- nrow(offon.H.msf) - vl.acute.int
  exp.onset.aids.H.msf <- nrow(offon.H.msf)
  offon.last.H.msf <- offon.H.msf[nrow(offon.H.msf), ]
  offon.H.msf <- rbind(offon.H.msf,
                     matrix(c(offon.last.H.msf[1] + (1:vl.aids.int),
                              rep(offon.last.H.msf[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.H.msf <- nrow(offon.H.msf)
  offon.H.msf[, 2] <- (1:max.possible.inf.time.H.msf) - offon.H.msf[, 1]
  stage.H.msf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.H.msf, vl.aids.int))
  stage.time.H.msf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.H.msf, 1:vl.aids.int)
  
  #Hispanic Immigrant males.
  prop.time.on.tx.HI.msf <- dat$param$tx.reinit.HI.msf.prob /
    (dat$param$tx.halt.HI.msf.prob + dat$param$tx.reinit.HI.msf.prob)
  offon.HI.msf <- matrix(c(1:tx.init.time.HI.msf, rep(0, tx.init.time.HI.msf)),
                       nrow = tx.init.time.HI.msf)
  while (offon.HI.msf[nrow(offon.HI.msf), 1] / dat$param$max.time.off.tx.part.int +
         offon.HI.msf[nrow(offon.HI.msf), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.HI.msf <- rbind(offon.HI.msf,
                        offon.HI.msf[nrow(offon.HI.msf), ] + c(1 - prop.time.on.tx.HI.msf,
                                                           prop.time.on.tx.HI.msf))
  }
  offon.HI.msf <- round(offon.HI.msf)
  exp.dur.chronic.HI.msf <- nrow(offon.HI.msf) - vl.acute.int
  exp.onset.aids.HI.msf <- nrow(offon.HI.msf)
  offon.last.HI.msf <- offon.HI.msf[nrow(offon.HI.msf), ]
  offon.HI.msf <- rbind(offon.HI.msf,
                      matrix(c(offon.last.HI.msf[1] + (1:vl.aids.int),
                               rep(offon.last.HI.msf[2], vl.aids.int)),
                             ncol = 2))
  max.possible.inf.time.HI.msf <- nrow(offon.HI.msf)
  offon.HI.msf[, 2] <- (1:max.possible.inf.time.HI.msf) - offon.HI.msf[, 1]
  stage.HI.msf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.HI.msf, vl.aids.int))
  stage.time.HI.msf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.HI.msf, 1:vl.aids.int)
  
  
  #White males
  prop.time.on.tx.W.msf <- dat$param$tx.reinit.W.msf.prob /
    (dat$param$tx.halt.W.msf.prob + dat$param$tx.reinit.W.msf.prob)
  offon.W.msf <- matrix(c(1:tx.init.time.W.msf, rep(0, tx.init.time.W.msf)),
                      nrow = tx.init.time.W.msf)
  
  while (offon.W.msf[nrow(offon.W.msf), 1] / dat$param$max.time.off.tx.part.int +
         offon.W.msf[nrow(offon.W.msf), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.W.msf <- rbind(offon.W.msf,
                       offon.W.msf[nrow(offon.W.msf), ] + c(1 - prop.time.on.tx.W.msf,
                                                        prop.time.on.tx.W.msf))
  }
  offon.W.msf <- round(offon.W.msf)
  exp.dur.chronic.W.msf <- nrow(offon.W.msf) - vl.acute.int
  exp.onset.aids.W.msf <- nrow(offon.W.msf)
  offon.last.W.msf <- offon.W.msf[nrow(offon.W.msf), ]
  offon.W.msf <- rbind(offon.W.msf,
                     matrix(c(offon.last.W.msf[1] + (1:vl.aids.int),
                              rep(offon.last.W.msf[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.W.msf <- nrow(offon.W.msf)
  offon.W.msf[, 2] <- (1:max.possible.inf.time.W.msf) - offon.W.msf[, 1]
  stage.W.msf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W.msf, vl.aids.int))
  stage.time.W.msf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W.msf, 1:vl.aids.int)
  
 
  #Black males MSM.
  prop.time.on.tx.B.msm <- dat$param$tx.reinit.B.msm.prob /
    (dat$param$tx.halt.B.msm.prob + dat$param$tx.reinit.B.msm.prob)
  offon.B.msm <- matrix(c(1:tx.init.time.B.msm, rep(0, tx.init.time.B.msm)),
                      nrow = tx.init.time.B.msm)
  while (offon.B.msm[nrow(offon.B.msm), 1] / dat$param$max.time.off.tx.part.int +
         offon.B.msm[nrow(offon.B.msm), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.B.msm <- rbind(offon.B.msm,
                       offon.B.msm[nrow(offon.B.msm), ] + c(1 - prop.time.on.tx.B.msm,
                                                        prop.time.on.tx.B.msm))
  }
  offon.B.msm <- round(offon.B.msm)
  exp.dur.chronic.B.msm <- nrow(offon.B.msm) - vl.acute.int
  exp.onset.aids.B.msm <- nrow(offon.B.msm)
  offon.last.B.msm <- offon.B.msm[nrow(offon.B.msm), ]
  offon.B.msm <- rbind(offon.B.msm,
                     matrix(c(offon.last.B.msm[1] + (1:vl.aids.int),
                              rep(offon.last.B.msm[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.B.msm <- nrow(offon.B.msm)
  offon.B.msm[, 2] <- (1:max.possible.inf.time.B.msm) - offon.B.msm[, 1]
  stage.B.msm <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B.msm, vl.aids.int))
  stage.time.B.msm <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B.msm, 1:vl.aids.int)
  
  #Black Immigrant males MSM=.
  prop.time.on.tx.BI.msm <- dat$param$tx.reinit.BI.msm.prob /
    (dat$param$tx.halt.BI.msm.prob + dat$param$tx.reinit.BI.msm.prob)
  offon.BI.msm <- matrix(c(1:tx.init.time.BI.msm, rep(0, tx.init.time.BI.msm)),
                       nrow = tx.init.time.BI.msm)
  while (offon.BI.msm[nrow(offon.BI.msm), 1] / dat$param$max.time.off.tx.part.int +
         offon.BI.msm[nrow(offon.BI.msm), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.BI.msm <- rbind(offon.BI.msm,
                        offon.BI.msm[nrow(offon.BI.msm), ] + c(1 - prop.time.on.tx.BI.msm,
                                                           prop.time.on.tx.BI.msm))
  }
  offon.BI.msm <- round(offon.BI.msm)
  exp.dur.chronic.BI.msm <- nrow(offon.BI.msm) - vl.acute.int
  exp.onset.aids.BI.msm <- nrow(offon.BI.msm)
  offon.last.BI.msm <- offon.BI.msm[nrow(offon.BI.msm), ]
  offon.BI.msm <- rbind(offon.BI.msm,
                      matrix(c(offon.last.BI.msm[1] + (1:vl.aids.int),
                               rep(offon.last.BI.msm[2], vl.aids.int)),
                             ncol = 2))
  max.possible.inf.time.BI.msm <- nrow(offon.BI.msm)
  offon.BI.msm[, 2] <- (1:max.possible.inf.time.BI.msm) - offon.BI.msm[, 1]
  stage.BI.msm <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.BI.msm, vl.aids.int))
  stage.time.BI.msm <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.BI.msm, 1:vl.aids.int)
  
  #Hispanic males MSM.
  prop.time.on.tx.H.msm <- dat$param$tx.reinit.H.msm.prob /
    (dat$param$tx.halt.H.msm.prob + dat$param$tx.reinit.H.msm.prob)
  offon.H.msm <- matrix(c(1:tx.init.time.H.msm, rep(0, tx.init.time.H.msm)),
                      nrow = tx.init.time.H.msm)
  while (offon.H.msm[nrow(offon.H.msm), 1] / dat$param$max.time.off.tx.part.int +
         offon.H.msm[nrow(offon.H.msm), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.H.msm <- rbind(offon.H.msm,
                       offon.H.msm[nrow(offon.H.msm), ] + c(1 - prop.time.on.tx.H.msm,
                                                        prop.time.on.tx.H.msm))
  }
  offon.H.msm <- round(offon.H.msm)
  exp.dur.chronic.H.msm <- nrow(offon.H.msm) - vl.acute.int
  exp.onset.aids.H.msm <- nrow(offon.H.msm)
  offon.last.H.msm <- offon.H.msm[nrow(offon.H.msm), ]
  offon.H.msm <- rbind(offon.H.msm,
                     matrix(c(offon.last.H.msm[1] + (1:vl.aids.int),
                              rep(offon.last.H.msm[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.H.msm <- nrow(offon.H.msm)
  offon.H.msm[, 2] <- (1:max.possible.inf.time.H.msm) - offon.H.msm[, 1]
  stage.H.msm <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.H.msm, vl.aids.int))
  stage.time.H.msm <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.H.msm, 1:vl.aids.int)
  
  #Hispanic Immigrant males MSM.
  prop.time.on.tx.HI.msm <- dat$param$tx.reinit.HI.msm.prob /
    (dat$param$tx.halt.HI.msm.prob + dat$param$tx.reinit.HI.msm.prob)
  offon.HI.msm <- matrix(c(1:tx.init.time.HI.msm, rep(0, tx.init.time.HI.msm)),
                       nrow = tx.init.time.HI.msm)
  while (offon.HI.msm[nrow(offon.HI.msm), 1] / dat$param$max.time.off.tx.part.int +
         offon.HI.msm[nrow(offon.HI.msm), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.HI.msm <- rbind(offon.HI.msm,
                        offon.HI.msm[nrow(offon.HI.msm), ] + c(1 - prop.time.on.tx.HI.msm,
                                                           prop.time.on.tx.HI.msm))
  }
  offon.HI.msm <- round(offon.HI.msm)
  exp.dur.chronic.HI.msm <- nrow(offon.HI.msm) - vl.acute.int
  exp.onset.aids.HI.msm <- nrow(offon.HI.msm)
  offon.last.HI.msm <- offon.HI.msm[nrow(offon.HI.msm), ]
  offon.HI.msm <- rbind(offon.HI.msm,
                      matrix(c(offon.last.HI.msm[1] + (1:vl.aids.int),
                               rep(offon.last.HI.msm[2], vl.aids.int)),
                             ncol = 2))
  max.possible.inf.time.HI.msm <- nrow(offon.HI.msm)
  offon.HI.msm[, 2] <- (1:max.possible.inf.time.HI.msm) - offon.HI.msm[, 1]
  stage.HI.msm <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.HI.msm, vl.aids.int))
  stage.time.HI.msm <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.HI.msm, 1:vl.aids.int)
  
  
  #White males MSM
  prop.time.on.tx.W.msm <- dat$param$tx.reinit.W.msm.prob /
    (dat$param$tx.halt.W.msm.prob + dat$param$tx.reinit.W.msm.prob)
  offon.W.msm <- matrix(c(1:tx.init.time.W.msm, rep(0, tx.init.time.W.msm)),
                      nrow = tx.init.time.W.msm)
  
  while (offon.W.msm[nrow(offon.W.msm), 1] / dat$param$max.time.off.tx.part.int +
         offon.W.msm[nrow(offon.W.msm), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.W.msm <- rbind(offon.W.msm,
                       offon.W.msm[nrow(offon.W.msm), ] + c(1 - prop.time.on.tx.W.msm,
                                                        prop.time.on.tx.W.msm))
  }
  offon.W.msm <- round(offon.W.msm)
  exp.dur.chronic.W.msm <- nrow(offon.W.msm) - vl.acute.int
  exp.onset.aids.W.msm <- nrow(offon.W.msm)
  offon.last.W.msm <- offon.W.msm[nrow(offon.W.msm), ]
  offon.W.msm <- rbind(offon.W.msm,
                     matrix(c(offon.last.W.msm[1] + (1:vl.aids.int),
                              rep(offon.last.W.msm[2], vl.aids.int)),
                            ncol = 2))
  max.possible.inf.time.W.msm <- nrow(offon.W.msm)
  offon.W.msm[, 2] <- (1:max.possible.inf.time.W.msm) - offon.W.msm[, 1]
  stage.W.msm <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W.msm, vl.aids.int))
  stage.time.W.msm <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W.msm, 1:vl.aids.int)
 
  
  
  #Black males MSMF.
  prop.time.on.tx.B.msmf <- dat$param$tx.reinit.B.msmf.prob /
    (dat$param$tx.halt.B.msmf.prob + dat$param$tx.reinit.B.msmf.prob)
  offon.B.msmf <- matrix(c(1:tx.init.time.B.msmf, rep(0, tx.init.time.B.msmf)),
                         nrow = tx.init.time.B.msmf)
  while (offon.B.msmf[nrow(offon.B.msmf), 1] / dat$param$max.time.off.tx.part.int +
         offon.B.msmf[nrow(offon.B.msmf), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.B.msmf <- rbind(offon.B.msmf,
                          offon.B.msmf[nrow(offon.B.msmf), ] + c(1 - prop.time.on.tx.B.msmf,
                                                                 prop.time.on.tx.B.msmf))
  }
  offon.B.msmf <- round(offon.B.msmf)
  exp.dur.chronic.B.msmf <- nrow(offon.B.msmf) - vl.acute.int
  exp.onset.aids.B.msmf <- nrow(offon.B.msmf)
  offon.last.B.msmf <- offon.B.msmf[nrow(offon.B.msmf), ]
  offon.B.msmf <- rbind(offon.B.msmf,
                        matrix(c(offon.last.B.msmf[1] + (1:vl.aids.int),
                                 rep(offon.last.B.msmf[2], vl.aids.int)),
                               ncol = 2))
  max.possible.inf.time.B.msmf <- nrow(offon.B.msmf)
  offon.B.msmf[, 2] <- (1:max.possible.inf.time.B.msmf) - offon.B.msmf[, 1]
  stage.B.msmf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.B.msmf, vl.aids.int))
  stage.time.B.msmf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.B.msmf, 1:vl.aids.int)
  
  #Black Immigrant males MSMF=.
  prop.time.on.tx.BI.msmf <- dat$param$tx.reinit.BI.msmf.prob /
    (dat$param$tx.halt.BI.msmf.prob + dat$param$tx.reinit.BI.msmf.prob)
  offon.BI.msmf <- matrix(c(1:tx.init.time.BI.msmf, rep(0, tx.init.time.BI.msmf)),
                          nrow = tx.init.time.BI.msmf)
  while (offon.BI.msmf[nrow(offon.BI.msmf), 1] / dat$param$max.time.off.tx.part.int +
         offon.BI.msmf[nrow(offon.BI.msmf), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.BI.msmf <- rbind(offon.BI.msmf,
                           offon.BI.msmf[nrow(offon.BI.msmf), ] + c(1 - prop.time.on.tx.BI.msmf,
                                                                    prop.time.on.tx.BI.msmf))
  }
  offon.BI.msmf <- round(offon.BI.msmf)
  exp.dur.chronic.BI.msmf <- nrow(offon.BI.msmf) - vl.acute.int
  exp.onset.aids.BI.msmf <- nrow(offon.BI.msmf)
  offon.last.BI.msmf <- offon.BI.msmf[nrow(offon.BI.msmf), ]
  offon.BI.msmf <- rbind(offon.BI.msmf,
                         matrix(c(offon.last.BI.msmf[1] + (1:vl.aids.int),
                                  rep(offon.last.BI.msmf[2], vl.aids.int)),
                                ncol = 2))
  max.possible.inf.time.BI.msmf <- nrow(offon.BI.msmf)
  offon.BI.msmf[, 2] <- (1:max.possible.inf.time.BI.msmf) - offon.BI.msmf[, 1]
  stage.BI.msmf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.BI.msmf, vl.aids.int))
  stage.time.BI.msmf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.BI.msmf, 1:vl.aids.int)
  
  #Hispanic males MSMF.
  prop.time.on.tx.H.msmf <- dat$param$tx.reinit.H.msmf.prob /
    (dat$param$tx.halt.H.msmf.prob + dat$param$tx.reinit.H.msmf.prob)
  offon.H.msmf <- matrix(c(1:tx.init.time.H.msmf, rep(0, tx.init.time.H.msmf)),
                         nrow = tx.init.time.H.msmf)
  while (offon.H.msmf[nrow(offon.H.msmf), 1] / dat$param$max.time.off.tx.part.int +
         offon.H.msmf[nrow(offon.H.msmf), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.H.msmf <- rbind(offon.H.msmf,
                          offon.H.msmf[nrow(offon.H.msmf), ] + c(1 - prop.time.on.tx.H.msmf,
                                                                 prop.time.on.tx.H.msmf))
  }
  offon.H.msmf <- round(offon.H.msmf)
  exp.dur.chronic.H.msmf <- nrow(offon.H.msmf) - vl.acute.int
  exp.onset.aids.H.msmf <- nrow(offon.H.msmf)
  offon.last.H.msmf <- offon.H.msmf[nrow(offon.H.msmf), ]
  offon.H.msmf <- rbind(offon.H.msmf,
                        matrix(c(offon.last.H.msmf[1] + (1:vl.aids.int),
                                 rep(offon.last.H.msmf[2], vl.aids.int)),
                               ncol = 2))
  max.possible.inf.time.H.msmf <- nrow(offon.H.msmf)
  offon.H.msmf[, 2] <- (1:max.possible.inf.time.H.msmf) - offon.H.msmf[, 1]
  stage.H.msmf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.H.msmf, vl.aids.int))
  stage.time.H.msmf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.H.msmf, 1:vl.aids.int)
  
  #Hispanic Immigrant males MSMF.
  prop.time.on.tx.HI.msmf <- dat$param$tx.reinit.HI.msmf.prob /
    (dat$param$tx.halt.HI.msmf.prob + dat$param$tx.reinit.HI.msmf.prob)
  offon.HI.msmf <- matrix(c(1:tx.init.time.HI.msmf, rep(0, tx.init.time.HI.msmf)),
                          nrow = tx.init.time.HI.msmf)
  while (offon.HI.msmf[nrow(offon.HI.msmf), 1] / dat$param$max.time.off.tx.part.int +
         offon.HI.msmf[nrow(offon.HI.msmf), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.HI.msmf <- rbind(offon.HI.msmf,
                           offon.HI.msmf[nrow(offon.HI.msmf), ] + c(1 - prop.time.on.tx.HI.msmf,
                                                                    prop.time.on.tx.HI.msmf))
  }
  offon.HI.msmf <- round(offon.HI.msmf)
  exp.dur.chronic.HI.msmf <- nrow(offon.HI.msmf) - vl.acute.int
  exp.onset.aids.HI.msmf <- nrow(offon.HI.msmf)
  offon.last.HI.msmf <- offon.HI.msmf[nrow(offon.HI.msmf), ]
  offon.HI.msmf <- rbind(offon.HI.msmf,
                         matrix(c(offon.last.HI.msmf[1] + (1:vl.aids.int),
                                  rep(offon.last.HI.msmf[2], vl.aids.int)),
                                ncol = 2))
  max.possible.inf.time.HI.msmf <- nrow(offon.HI.msmf)
  offon.HI.msmf[, 2] <- (1:max.possible.inf.time.HI.msmf) - offon.HI.msmf[, 1]
  stage.HI.msmf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.HI.msmf, vl.aids.int))
  stage.time.HI.msmf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.HI.msmf, 1:vl.aids.int)
  
  
  #White males MSMF
  prop.time.on.tx.W.msmf <- dat$param$tx.reinit.W.msmf.prob /
    (dat$param$tx.halt.W.msmf.prob + dat$param$tx.reinit.W.msmf.prob)
  offon.W.msmf <- matrix(c(1:tx.init.time.W.msmf, rep(0, tx.init.time.W.msmf)),
                         nrow = tx.init.time.W.msmf)
  
  while (offon.W.msmf[nrow(offon.W.msmf), 1] / dat$param$max.time.off.tx.part.int +
         offon.W.msmf[nrow(offon.W.msmf), 2] / dat$param$max.time.on.tx.part.int < 1) {
    offon.W.msmf <- rbind(offon.W.msmf,
                          offon.W.msmf[nrow(offon.W.msmf), ] + c(1 - prop.time.on.tx.W.msmf,
                                                                 prop.time.on.tx.W.msmf))
  }
  offon.W.msmf <- round(offon.W.msmf)
  exp.dur.chronic.W.msmf <- nrow(offon.W.msmf) - vl.acute.int
  exp.onset.aids.W.msmf <- nrow(offon.W.msmf)
  offon.last.W.msmf <- offon.W.msmf[nrow(offon.W.msmf), ]
  offon.W.msmf <- rbind(offon.W.msmf,
                        matrix(c(offon.last.W.msmf[1] + (1:vl.aids.int),
                                 rep(offon.last.W.msmf[2], vl.aids.int)),
                               ncol = 2))
  max.possible.inf.time.W.msmf <- nrow(offon.W.msmf)
  offon.W.msmf[, 2] <- (1:max.possible.inf.time.W.msmf) - offon.W.msmf[, 1]
  stage.W.msmf <- rep(c(1, 2, 3, 4), c(vlar.int, vlaf.int, exp.dur.chronic.W.msmf, vl.aids.int))
  stage.time.W.msmf <- c(1:vlar.int, 1:vlaf.int, 1:exp.dur.chronic.W.msmf, 1:vl.aids.int)
  
  
  
  # VL for Black females
  selected <- which(status == 1 & tt.traj == 3 & race == "B" & sex == "F")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B.f)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B.f[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B.f[time.since.inf, 1]
  stage[selected] <- stage.B.f[time.since.inf]
  stage.time[selected] <- stage.time.B.f[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B.f)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.B.f) * (vlsp) +
                  (time.since.inf > exp.onset.aids.B.f) *
                  (vlsp + (time.since.inf - exp.onset.aids.B.f) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp

  # VL for Black Immigrant females
  selected <- which(status == 1 & tt.traj == 3 & race == "BI" & sex == "F")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.BI.f)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.BI.f[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.BI.f[time.since.inf, 1]
  stage[selected] <- stage.BI.f[time.since.inf]
  stage.time[selected] <- stage.time.BI.f[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.BI.f)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.BI.f) * (vlsp) +
    (time.since.inf > exp.onset.aids.BI.f) *
    (vlsp + (time.since.inf - exp.onset.aids.BI.f) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for Hispanic females
  selected <- which(status == 1 & tt.traj == 3 & race == "H" & sex == "F")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.H.f)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.H.f[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.H.f[time.since.inf, 1]
  stage[selected] <- stage.H.f[time.since.inf]
  stage.time[selected] <- stage.time.H.f[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.H.f)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.H.f) * (vlsp) +
    (time.since.inf > exp.onset.aids.H.f) *
    (vlsp + (time.since.inf - exp.onset.aids.H.f) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for Hispanic Immigrant females
  selected <- which(status == 1 & tt.traj == 3 & race == "HI" & sex == "F")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.HI.f)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.HI.f[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.HI.f[time.since.inf, 1]
  stage[selected] <- stage.HI.f[time.since.inf]
  stage.time[selected] <- stage.time.HI.f[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.HI.f)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.HI.f) * (vlsp) +
    (time.since.inf > exp.onset.aids.HI.f) *
    (vlsp + (time.since.inf - exp.onset.aids.HI.f) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for White females
  selected <- which(status == 1 & tt.traj == 3 & race == "W" & sex == "F")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W.f)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W.f[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W.f[time.since.inf, 1]
  stage[selected] <- stage.W.f[time.since.inf]
  stage.time[selected] <- stage.time.W.f[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W.f)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
                  (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
                     ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
                  (time.since.inf > vlar.int + vlaf.int) *
                  (time.since.inf <= exp.onset.aids.W.f) * (vlsp) +
                  (time.since.inf > exp.onset.aids.W.f) *
                  (vlsp + (time.since.inf - exp.onset.aids.W.f) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp

  #Males HET MSF.
  # VL for Black males Het
  selected <- which(status == 1 & tt.traj == 3 & race == "B" & sex.ident == "msf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B.msf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B.msf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B.msf[time.since.inf, 1]
  stage[selected] <- stage.B.msf[time.since.inf]
  stage.time[selected] <- stage.time.B.msf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B.msf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.B.msf) * (vlsp) +
    (time.since.inf > exp.onset.aids.B.msf) *
    (vlsp + (time.since.inf - exp.onset.aids.B.msf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for Black Immigrant males MSF
  selected <- which(status == 1 & tt.traj == 3 & race == "BI" & sex.ident == "msf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.BI.msf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.BI.msf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.BI.msf[time.since.inf, 1]
  stage[selected] <- stage.BI.msf[time.since.inf]
  stage.time[selected] <- stage.time.BI.msf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.BI.msf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.BI.msf) * (vlsp) +
    (time.since.inf > exp.onset.aids.BI.msf) *
    (vlsp + (time.since.inf - exp.onset.aids.BI.msf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for Hispanic males MSF
  selected <- which(status == 1 & tt.traj == 3 & race == "H" & sex.ident == "msf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.H.msf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.H.msf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.H.msf[time.since.inf, 1]
  stage[selected] <- stage.H.msf[time.since.inf]
  stage.time[selected] <- stage.time.H.msf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.H.msf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.H.msf) * (vlsp) +
    (time.since.inf > exp.onset.aids.H.msf) *
    (vlsp + (time.since.inf - exp.onset.aids.H.msf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for Hispanic Immigrant males MSF
  selected <- which(status == 1 & tt.traj == 3 & race == "HI" & sex.ident == "msf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.HI.msf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.HI.msf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.HI.msf[time.since.inf, 1]
  stage[selected] <- stage.HI.msf[time.since.inf]
  stage.time[selected] <- stage.time.HI.msf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.HI.msf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.HI.msf) * (vlsp) +
    (time.since.inf > exp.onset.aids.HI.msf) *
    (vlsp + (time.since.inf - exp.onset.aids.HI.msf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for White males MSF
  selected <- which(status == 1 & tt.traj == 3 & race == "W" & sex.ident == "msf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W.msf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W.msf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W.msf[time.since.inf, 1]
  stage[selected] <- stage.W.msf[time.since.inf]
  stage.time[selected] <- stage.time.W.msf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W.msf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.W.msf) * (vlsp) +
    (time.since.inf > exp.onset.aids.W.msf) *
    (vlsp + (time.since.inf - exp.onset.aids.W.msf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
 
  
  #Males MSM
  # VL for Black males msm
  selected <- which(status == 1 & tt.traj == 3 & race == "B" & sex.ident == "msm")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B.msm)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B.msm[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B.msm[time.since.inf, 1]
  stage[selected] <- stage.B.msm[time.since.inf]
  stage.time[selected] <- stage.time.B.msm[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B.msm)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.B.msm) * (vlsp) +
    (time.since.inf > exp.onset.aids.B.msm) *
    (vlsp + (time.since.inf - exp.onset.aids.B.msm) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for Black Immigrant males MSM
  selected <- which(status == 1 & tt.traj == 3 & race == "BI" & sex.ident == "msm")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.BI.msm)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.BI.msm[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.BI.msm[time.since.inf, 1]
  stage[selected] <- stage.BI.msm[time.since.inf]
  stage.time[selected] <- stage.time.BI.msm[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.BI.msm)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.BI.msm) * (vlsp) +
    (time.since.inf > exp.onset.aids.BI.msm) *
    (vlsp + (time.since.inf - exp.onset.aids.BI.msm) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for Hispanic males MSM
  selected <- which(status == 1 & tt.traj == 3 & race == "H" & sex.ident == "msm")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.H.msm)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.H.msm[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.H.msm[time.since.inf, 1]
  stage[selected] <- stage.H.msm[time.since.inf]
  stage.time[selected] <- stage.time.H.msm[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.H.msm)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.H.msm) * (vlsp) +
    (time.since.inf > exp.onset.aids.H.msm) *
    (vlsp + (time.since.inf - exp.onset.aids.H.msm) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for Hispanic Immigrant males MSM
  selected <- which(status == 1 & tt.traj == 3 & race == "HI" & sex.ident == "msm")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.HI.msm)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.HI.msm[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.HI.msm[time.since.inf, 1]
  stage[selected] <- stage.HI.msm[time.since.inf]
  stage.time[selected] <- stage.time.HI.msm[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.HI.msm)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.HI.msm) * (vlsp) +
    (time.since.inf > exp.onset.aids.HI.msm) *
    (vlsp + (time.since.inf - exp.onset.aids.HI.msm) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for White males MSM
  selected <- which(status == 1 & tt.traj == 3 & race == "W" & sex.ident == "msm")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W.msm)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W.msm[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W.msm[time.since.inf, 1]
  stage[selected] <- stage.W.msm[time.since.inf]
  stage.time[selected] <- stage.time.W.msm[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W.msm)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.W.msm) * (vlsp) +
    (time.since.inf > exp.onset.aids.W.msm) *
    (vlsp + (time.since.inf - exp.onset.aids.W.msm) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  #Males MSMF
  # VL for Black males msmf
  selected <- which(status == 1 & tt.traj == 3 & race == "B" & sex.ident == "msmf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.B.msmf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.B.msmf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.B.msmf[time.since.inf, 1]
  stage[selected] <- stage.B.msmf[time.since.inf]
  stage.time[selected] <- stage.time.B.msmf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.B.msmf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.B.msmf) * (vlsp) +
    (time.since.inf > exp.onset.aids.B.msmf) *
    (vlsp + (time.since.inf - exp.onset.aids.B.msmf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for Black Immigrant males MSMF
  selected <- which(status == 1 & tt.traj == 3 & race == "BI" & sex.ident == "msmf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.BI.msmf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.BI.msmf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.BI.msmf[time.since.inf, 1]
  stage[selected] <- stage.BI.msmf[time.since.inf]
  stage.time[selected] <- stage.time.BI.msmf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.BI.msmf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.BI.msmf) * (vlsp) +
    (time.since.inf > exp.onset.aids.BI.msmf) *
    (vlsp + (time.since.inf - exp.onset.aids.BI.msmf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for Hispanic males MSMF
  selected <- which(status == 1 & tt.traj == 3 & race == "H" & sex.ident == "msmf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.H.msmf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.H.msmf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.H.msmf[time.since.inf, 1]
  stage[selected] <- stage.H.msmf[time.since.inf]
  stage.time[selected] <- stage.time.H.msmf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.H.msmf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.H.msmf) * (vlsp) +
    (time.since.inf > exp.onset.aids.H.msmf) *
    (vlsp + (time.since.inf - exp.onset.aids.H.msmf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for Hispanic Immigrant males MSMF
  selected <- which(status == 1 & tt.traj == 3 & race == "HI" & sex.ident == "msmf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.HI.msmf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.HI.msmf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.HI.msmf[time.since.inf, 1]
  stage[selected] <- stage.HI.msmf[time.since.inf]
  stage.time[selected] <- stage.time.HI.msmf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.HI.msmf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.HI.msmf) * (vlsp) +
    (time.since.inf > exp.onset.aids.HI.msmf) *
    (vlsp + (time.since.inf - exp.onset.aids.HI.msmf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  # VL for White males MSMF
  selected <- which(status == 1 & tt.traj == 3 & race == "W" & sex.ident == "msmf")
  max.inf.time <- pmin(time.sex.active[selected], max.possible.inf.time.W.msmf)
  time.since.inf <- ceiling(runif(length(selected), max = max.inf.time))
  inf.time[selected] <- 1 - time.since.inf
  cum.time.on.tx[selected] <- offon.W.msmf[time.since.inf, 2]
  cum.time.off.tx[selected] <- offon.W.msmf[time.since.inf, 1]
  stage[selected] <- stage.W.msmf[time.since.inf]
  stage.time[selected] <- stage.time.W.msmf[time.since.inf]
  tx.status[selected] <- 0
  tx.status[selected][stage[selected] == 3 & cum.time.on.tx[selected] > 0] <-
    rbinom(sum(stage[selected] == 3 & cum.time.on.tx[selected] > 0),
           1, prop.time.on.tx.W.msmf)
  vl[selected] <- (time.since.inf <= vlar.int) * (vlap * time.since.inf / vlar.int) +
    (time.since.inf > vlar.int) * (time.since.inf <= vlar.int + vlaf.int) *
    ((vlsp - vlap) * (time.since.inf - vlar.int) / vlaf.int + vlap) +
    (time.since.inf > vlar.int + vlaf.int) *
    (time.since.inf <= exp.onset.aids.W.msmf) * (vlsp) +
    (time.since.inf > exp.onset.aids.W.msmf) *
    (vlsp + (time.since.inf - exp.onset.aids.W.msmf) * vlds)
  vl[selected][tx.status[selected] == 1] <- dat$param$vl.part.supp
  
  
  # Implement diagnosis for both
  selected <- which(status == 1 & tt.traj == 3)
  if (dat$param$testing.pattern == "interval") {
    ttntest <- ceiling(runif(length(selected),
                       min = 0,
                       max = dat$param$mean.test.B.f.int * (race[selected] == "B" & sex[selected] == "F") +
                              dat$param$mean.test.BI.f.int * (race[selected] == "BI" & sex[selected] == "F") +
                              dat$param$mean.test.H.f.int * (race[selected] == "H" & sex[selected] == "F") +
                              dat$param$mean.test.HI.f.int * (race[selected] == "HI" & sex[selected] == "F") +
                              dat$param$mean.test.W.f.int * (race[selected] == "W" & sex[selected] == "F") +
                              dat$param$mean.test.B.msf.int * (race[selected] == "B" & sex.ident[selected] == "msf") +
                              dat$param$mean.test.BI.msf.int * (race[selected] == "BI" & sex.ident[selected] == "msf") +
                              dat$param$mean.test.H.msf.int * (race[selected] == "H" & sex.ident[selected] == "msf") +
                              dat$param$mean.test.HI.msf.int * (race[selected] == "HI" & sex.ident[selected] == "msf") +
                              dat$param$mean.test.W.msf.int * (race[selected] == "W" & sex.ident[selected] == "msf") +
                              dat$param$mean.test.B.msm.int * (race[selected] == "B" & sex.ident[selected] == "msm") +
                              dat$param$mean.test.BI.msm.int * (race[selected] == "BI" & sex.ident[selected] == "msm") +
                              dat$param$mean.test.H.msm.int * (race[selected] == "H" & sex.ident[selected] == "msm") +
                              dat$param$mean.test.HI.msm.int * (race[selected] == "HI" & sex.ident[selected] == "msm") +
                              dat$param$mean.test.W.msm.int * (race[selected] == "W" & sex.ident[selected] == "msm") +
                              dat$param$mean.test.B.msmf.int * (race[selected] == "B" & sex.ident[selected] == "msmf") +
                              dat$param$mean.test.BI.msmf.int * (race[selected] == "BI" & sex.ident[selected] == "msmf") +
                              dat$param$mean.test.H.msmf.int * (race[selected] == "H" & sex.ident[selected] == "msmf") +
                              dat$param$mean.test.HI.msmf.int * (race[selected] == "HI" & sex.ident[selected] == "msmf") +
                              dat$param$mean.test.W.msmf.int * (race[selected] == "W" & sex.ident[selected] == "msmf")))
  }

  if (dat$param$testing.pattern == "memoryless") {
    ttntest <- rgeom(length(selected),
                     1 / (dat$param$mean.test.B.f.int * (race[selected] == "B" & sex[selected] == "F") +
                          dat$param$mean.test.BI.f.int * (race[selected] == "BI" & sex[selected] == "F") +
                          dat$param$mean.test.H.f.int * (race[selected] == "H" & sex[selected] == "F") +
                          dat$param$mean.test.HI.f.int * (race[selected] == "HI" & sex[selected] == "F") +
                          dat$param$mean.test.W.f.int * (race[selected] == "W" & sex[selected] == "F") +
                          dat$param$mean.test.B.msf.int * (race[selected] == "B" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.BI.msf.int * (race[selected] == "BI" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.H.msf.int * (race[selected] == "H" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.HI.msf.int * (race[selected] == "HI" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.W.msf.int * (race[selected] == "W" & sex.ident[selected] == "msf") +
                          dat$param$mean.test.B.msm.int * (race[selected] == "B" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.BI.msm.int * (race[selected] == "BI" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.H.msm.int * (race[selected] == "H" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.HI.msm.int * (race[selected] == "HI" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.W.msm.int * (race[selected] == "W" & sex.ident[selected] == "msm") +
                          dat$param$mean.test.B.msmf.int * (race[selected] == "B" & sex.ident[selected] == "msmf") +
                          dat$param$mean.test.BI.msmf.int * (race[selected] == "BI" & sex.ident[selected] == "msmf") +
                          dat$param$mean.test.H.msmf.int * (race[selected] == "H" & sex.ident[selected] == "msmf") +
                          dat$param$mean.test.HI.msmf.int * (race[selected] == "HI" & sex.ident[selected] == "msmf") +
                          dat$param$mean.test.W.msmf.int * (race[selected] == "W" & sex.ident[selected] == "msmf")))
  }


  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
    -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]

  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1
  diag.status[selected][cum.time.on.tx[selected] > 0] <- 1
  last.neg.test[selected][cum.time.on.tx[selected] > 0] <- NA


  # Last neg test before present for negatives
  selected <- which(status == 0 & tt.traj %in% c(2, 3, 4))

  if (dat$param$testing.pattern == "interval") {
    tslt <- ceiling(runif(length(selected),
                      min = 0,
                      max = dat$param$mean.test.B.f.int * (race[selected] == "B" & sex[selected] == "F") +
                            dat$param$mean.test.BI.f.int * (race[selected] == "BI" & sex[selected] == "F") +
                            dat$param$mean.test.H.f.int * (race[selected] == "H" & sex[selected] == "F") +
                            dat$param$mean.test.HI.f.int * (race[selected] == "HI" & sex[selected] == "F") +
                            dat$param$mean.test.W.f.int * (race[selected] == "W" & sex[selected] == "F") +
                            dat$param$mean.test.B.msf.int * (race[selected] == "B" & sex.ident[selected] == "msf") +
                            dat$param$mean.test.BI.msf.int * (race[selected] == "BI" & sex.ident[selected] == "msf") +
                            dat$param$mean.test.H.msf.int * (race[selected] == "H" & sex.ident[selected] == "msf") +
                            dat$param$mean.test.HI.msf.int * (race[selected] == "HI" & sex.ident[selected] == "msf") +
                            dat$param$mean.test.W.msf.int * (race[selected] == "W" & sex.ident[selected] == "msf") +
                            dat$param$mean.test.B.msm.int * (race[selected] == "B" & sex.ident[selected] == "msm") +
                            dat$param$mean.test.BI.msm.int * (race[selected] == "BI" & sex.ident[selected] == "msm") +
                            dat$param$mean.test.H.msm.int * (race[selected] == "H" & sex.ident[selected] == "msm") +
                            dat$param$mean.test.HI.msm.int * (race[selected] == "HI" & sex.ident[selected] == "msm") +
                            dat$param$mean.test.W.msm.int * (race[selected] == "W" & sex.ident[selected] == "msm") +
                            dat$param$mean.test.B.msmf.int * (race[selected] == "B" & sex.ident[selected] == "msmf") +
                            dat$param$mean.test.BI.msmf.int * (race[selected] == "BI" & sex.ident[selected] == "msmf") +
                            dat$param$mean.test.H.msmf.int * (race[selected] == "H" & sex.ident[selected] == "msmf") +
                            dat$param$mean.test.HI.msmf.int * (race[selected] == "HI" & sex.ident[selected] == "msmf") +
                            dat$param$mean.test.W.msmf.int * (race[selected] == "W" & sex.ident[selected] == "msmf")
                      ))
  }
  if (dat$param$testing.pattern == "memoryless") {
    tslt <- rgeom(length(selected),
                  1 / (dat$param$mean.test.B.f.int * (race[selected] == "B" & sex[selected] == "F") +
                       dat$param$mean.test.BI.f.int * (race[selected] == "BI" & sex[selected] == "F") +
                       dat$param$mean.test.H.f.int * (race[selected] == "H" & sex[selected] == "F") +
                       dat$param$mean.test.HI.f.int * (race[selected] == "HI" & sex[selected] == "F") +
                       dat$param$mean.test.W.f.int * (race[selected] == "W" & sex[selected] == "F") +
                       dat$param$mean.test.B.msf.int * (race[selected] == "B" & sex.ident[selected] == "msf") +
                       dat$param$mean.test.BI.msf.int * (race[selected] == "BI" & sex.ident[selected] == "msf") +
                       dat$param$mean.test.H.msf.int * (race[selected] == "H" & sex.ident[selected] == "msf") +
                       dat$param$mean.test.HI.msf.int * (race[selected] == "HI" & sex.ident[selected] == "msf") +
                       dat$param$mean.test.W.msf.int * (race[selected] == "W" & sex.ident[selected] == "msf") +
                       dat$param$mean.test.B.msm.int * (race[selected] == "B" & sex.ident[selected] == "msm") +
                       dat$param$mean.test.BI.msm.int * (race[selected] == "BI" & sex.ident[selected] == "msm") +
                       dat$param$mean.test.H.msm.int * (race[selected] == "H" & sex.ident[selected] == "msm") +
                       dat$param$mean.test.HI.msm.int * (race[selected] == "HI" & sex.ident[selected] == "msm") +
                       dat$param$mean.test.W.msm.int * (race[selected] == "W" & sex.ident[selected] == "msm") +
                       dat$param$mean.test.B.msmf.int * (race[selected] == "B" & sex.ident[selected] == "msmf") +
                       dat$param$mean.test.BI.msmf.int * (race[selected] == "BI" & sex.ident[selected] == "msmf") +
                       dat$param$mean.test.H.msmf.int * (race[selected] == "H" & sex.ident[selected] == "msmf") +
                       dat$param$mean.test.HI.msmf.int * (race[selected] == "HI" & sex.ident[selected] == "msmf") +
                       dat$param$mean.test.W.msmf.int * (race[selected] == "W" & sex.ident[selected] == "msmf")))
  }
  last.neg.test[selected] <- -tslt

  
  ##Infection Class
  selected <- which(status == 1)
  inf.class[selected] <- "Lhet"

  ## Set all onto dat$attr
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time
  dat$attr$inf.time <- inf.time
  dat$attr$vl <- vl
  dat$attr$diag.status <- diag.status
  dat$attr$diag.time <- diag.time
  dat$attr$last.neg.test <- last.neg.test
  dat$attr$evertest<-rep(NA, num)
  dat$attr$evertest<-ifelse(dat$attr$last.neg.test<=0 | dat$attr$diag.status==1,1,0)
  dat$attr$tx.status <- tx.status
  dat$attr$tx.init.time <- tx.init.time
  dat$attr$cum.time.on.tx <- cum.time.on.tx
  dat$attr$cum.time.off.tx <- cum.time.off.tx
  dat$attr$infector <- infector
  dat$attr$inf.role <- inf.role
  dat$attr$inf.type <- inf.type
  dat$attr$inf.diag <- inf.diag
  dat$attr$inf.tx <- inf.tx
  dat$attr$inf.stage <- inf.stage
  dat$attr$inf.class <- inf.class

  return(dat)

}


#' @title Sets the CCR5 genetic status of persons
#'
#' @description Initializes the CCR5-delta-32 genetic allele of the men in the
#'              population, based on parameters defining the probability
#'              distribution.
#'
#' @param dat Data object created in initialization module.
#'
#' @export
#' @keywords initiation utility msm
#'
init_ccr5_shamp <- function(dat) {
  
  sex<-dat$attr$sex
  sex.ident<-dat$attr$sex.ident
  race<-dat$attr$race

  num.B.f <- length(which(dat$attr$race == "B" & dat$attr$sex=="F"))
  num.BI.f <- length(which(dat$attr$race == "BI" & dat$attr$sex=="F"))
  num.H.f <- length(which(dat$attr$race == "H" & dat$attr$sex=="F"))
  num.HI.f <- length(which(dat$attr$race == "HI" & dat$attr$sex=="F"))
  num.W.f <- length(which(dat$attr$race == "W" & dat$attr$sex=="F"))
  
  num.B.m <- length(which(dat$attr$race == "B" & dat$attr$sex=="M"))
  num.BI.m <- length(which(dat$attr$race == "BI" & dat$attr$sex=="M"))
  num.H.m <- length(which(dat$attr$race == "H" & dat$attr$sex=="M"))
  num.HI.m <- length(which(dat$attr$race == "HI" & dat$attr$sex=="M"))
  num.W.m <- length(which(dat$attr$race == "W" & dat$attr$sex=="M"))
  num <- num.B.f + num.BI.f + num.H.f + num.HI.f + num.W.f + num.B.m + num.BI.m + num.H.m + num.HI.m + num.W.m
  race <- dat$attr$race
  status <- dat$attr$status

  nInfB.f <- sum(race == "B" & status == 1 & sex=="F")
  nInfBI.f <- sum(race == "BI" & status == 1 & sex=="F")
  nInfH.f <- sum(race == "H" & status == 1 & sex=="F")
  nInfHI.f <- sum(race == "HI" & status == 1 & sex=="F")
  nInfW.f <- sum(race == "W" & status == 1 & sex=="F")
  
  nInfB.m <- sum(race == "B" & status == 1 & sex=="M")
  nInfBI.m <- sum(race == "BI" & status == 1 & sex=="M")
  nInfH.m <- sum(race == "H" & status == 1 & sex=="M")
  nInfHI.m <- sum(race == "HI" & status == 1 & sex=="M")
  nInfW.m <- sum(race == "W" & status == 1 & sex=="M")
  
  ##  CCR5 genotype
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr
  ccr5 <- rep("WW", num)
  
  ##Black females
if (num.B.f > 0){
  # homozygotes for deletion
  num.ccr5.DD.B.f <- dat$param$ccr5.B.f.prob[1] * num.B.f
  # heterozygotes
  num.ccr5.DW.B.f <- dat$param$ccr5.B.f.prob[2] * num.B.f
  # homozygotes for deletion
  num.ccr5.WW.B.f <- num.B.f - num.ccr5.DD.B.f - num.ccr5.DW.B.f
  # DD's can't be infected
  num.uninf.ccr5.DD.B.f <- round(num.ccr5.DD.B.f)
  
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.B.f <- round(num.ccr5.DW.B.f * nInfB.f * ccr5.heteroz.rr /
                             (num.ccr5.WW.B.f + num.ccr5.DW.B.f * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.B.f <- round(num.ccr5.DW.B.f - num.inf.ccr5.DW.B.f)
  inf.B.f <- which(status == 1 & race == "B" & sex == "F")
  inf.ccr5.DW.B.f <- sample(inf.B.f, num.inf.ccr5.DW.B.f, replace = FALSE)
  ccr5[inf.ccr5.DW.B.f] <- "DW"
  uninf.B.f <- which(status == 0 & race == "B" & sex == "F")
  uninf.ccr5.DWDD.B.f <- sample(uninf.B.f, num.uninf.ccr5.DW.B.f + num.uninf.ccr5.DD.B.f)
  uninf.ccr5.DW.B.f <- sample(uninf.ccr5.DWDD.B.f, num.uninf.ccr5.DW.B.f)
  uninf.ccr5.DD.B.f <- setdiff(uninf.ccr5.DWDD.B.f, uninf.ccr5.DW.B.f)
  ccr5[uninf.ccr5.DW.B.f] <- "DW"
  ccr5[uninf.ccr5.DD.B.f] <- "DD"
  }

  ##Black immigrant females
  if (num.BI.f > 0){
  # homozygotes for deletion
  num.ccr5.DD.BI.f <- dat$param$ccr5.BI.f.prob[1] * num.BI.f
  # heterozygotes
  num.ccr5.DW.BI.f <- dat$param$ccr5.BI.f.prob[2] * num.BI.f
  # homozygotes for deletion
  num.ccr5.WW.BI.f <- num.BI.f - num.ccr5.DD.BI.f - num.ccr5.DW.BI.f
  # DD's can't be infected
  num.uninf.ccr5.DD.BI.f <- round(num.ccr5.DD.BI.f)
  
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.BI.f <- round(num.ccr5.DW.BI.f * nInfBI.f * ccr5.heteroz.rr /
                                 (num.ccr5.WW.BI.f + num.ccr5.DW.BI.f * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.BI.f <- round(num.ccr5.DW.BI.f - num.inf.ccr5.DW.BI.f)
  inf.BI.f <- which(status == 1 & race == "BI" & sex == "F")
  inf.ccr5.DW.BI.f <- sample(inf.BI.f, num.inf.ccr5.DW.BI.f, replace = FALSE)
  ccr5[inf.ccr5.DW.BI.f] <- "DW"
  uninf.BI.f <- which(status == 0 & race == "BI" & sex == "F")
  uninf.ccr5.DWDD.BI.f <- sample(uninf.BI.f, num.uninf.ccr5.DW.BI.f + num.uninf.ccr5.DD.BI.f)
  uninf.ccr5.DW.BI.f <- sample(uninf.ccr5.DWDD.BI.f, num.uninf.ccr5.DW.BI.f)
  uninf.ccr5.DD.BI.f <- setdiff(uninf.ccr5.DWDD.BI.f, uninf.ccr5.DW.BI.f)
  ccr5[uninf.ccr5.DW.BI.f] <- "DW"
  ccr5[uninf.ccr5.DD.BI.f] <- "DD"
  }
  
  ##Hisapnic females
  if (num.H.f > 0){
  # homozygotes for deletion
  num.ccr5.DD.H.f <- dat$param$ccr5.H.f.prob[1] * num.H.f
  # heterozygotes
  num.ccr5.DW.H.f <- dat$param$ccr5.H.f.prob[2] * num.H.f
  # homozygotes for deletion
  num.ccr5.WW.H.f <- num.H.f - num.ccr5.DD.H.f - num.ccr5.DW.H.f
  # DD's can't be infected
  num.uninf.ccr5.DD.H.f <- round(num.ccr5.DD.H.f)
  
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.H.f <- round(num.ccr5.DW.H.f * nInfH.f * ccr5.heteroz.rr /
                                 (num.ccr5.WW.H.f + num.ccr5.DW.H.f * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.H.f <- round(num.ccr5.DW.H.f - num.inf.ccr5.DW.H.f)
  inf.H.f <- which(status == 1 & race == "H" & sex == "F")
  inf.ccr5.DW.H.f <- sample(inf.H.f, num.inf.ccr5.DW.H.f, replace = FALSE)
  ccr5[inf.ccr5.DW.H.f] <- "DW"
  uninf.H.f <- which(status == 0 & race == "H" & sex == "F")
  uninf.ccr5.DWDD.H.f <- sample(uninf.H.f, num.uninf.ccr5.DW.H.f + num.uninf.ccr5.DD.H.f)
  uninf.ccr5.DW.H.f <- sample(uninf.ccr5.DWDD.H.f, num.uninf.ccr5.DW.H.f)
  uninf.ccr5.DD.H.f <- setdiff(uninf.ccr5.DWDD.H.f, uninf.ccr5.DW.H.f)
  ccr5[uninf.ccr5.DW.H.f] <- "DW"
  ccr5[uninf.ccr5.DD.H.f] <- "DD"
  }
  
  ##Hispanic immigrant females
 
   if (num.HI.f > 0){
    
  # homozygotes for deletion
  num.ccr5.DD.HI.f <- dat$param$ccr5.HI.f.prob[1] * num.HI.f
  # heterozygotes
  num.ccr5.DW.HI.f <- dat$param$ccr5.HI.f.prob[2] * num.HI.f
  # homozygotes for deletion
  num.ccr5.WW.HI.f <- num.HI.f - num.ccr5.DD.HI.f - num.ccr5.DW.HI.f
  # DD's can't be infected
  num.uninf.ccr5.DD.HI.f <- round(num.ccr5.DD.HI.f)
  
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.HI.f <- round(num.ccr5.DW.HI.f * nInfHI.f * ccr5.heteroz.rr /
                                  (num.ccr5.WW.HI.f + num.ccr5.DW.HI.f * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.HI.f <- round(num.ccr5.DW.HI.f - num.inf.ccr5.DW.HI.f)
  inf.HI.f <- which(status == 1 & race == "HI" & sex == "F")
  inf.ccr5.DW.HI.f <- sample(inf.HI.f, num.inf.ccr5.DW.HI.f, replace = FALSE)
  ccr5[inf.ccr5.DW.HI.f] <- "DW"
  uninf.HI.f <- which(status == 0 & race == "HI" & sex == "F")
  uninf.ccr5.DWDD.HI.f <- sample(uninf.HI.f, num.uninf.ccr5.DW.HI.f + num.uninf.ccr5.DD.HI.f)
  uninf.ccr5.DW.HI.f <- sample(uninf.ccr5.DWDD.HI.f, num.uninf.ccr5.DW.HI.f)
  uninf.ccr5.DD.HI.f <- setdiff(uninf.ccr5.DWDD.HI.f, uninf.ccr5.DW.HI.f)
  ccr5[uninf.ccr5.DW.HI.f] <- "DW"
  ccr5[uninf.ccr5.DD.HI.f] <- "DD"
  }
  
  ## White females
  if (num.W.f > 0){
  # homozygotes for deletion
  num.ccr5.DD.W.f <- dat$param$ccr5.W.f.prob[1] * num.W.f
  # heterozygotes
  num.ccr5.DW.W.f <- dat$param$ccr5.W.f.prob[2] * num.W.f
  # homozygotes for deletion
  num.ccr5.WW.W.f <- num.W.f - num.ccr5.DD.W.f - num.ccr5.DW.W.f
  # DD's can't be infected
  num.uninf.ccr5.DD.W.f <- round(num.ccr5.DD.W.f)
  
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.W.f <- round(num.ccr5.DW.W.f * nInfW.f * ccr5.heteroz.rr /
                             (num.ccr5.WW.W.f + num.ccr5.DW.W.f * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.W.f <- round(num.ccr5.DW.W.f - num.inf.ccr5.DW.W.f)
  inf.W.f <- which(status == 1 & race == "W" & sex == "F")
  inf.ccr5.DW.W.f <- sample(inf.W.f, num.inf.ccr5.DW.W.f)
  ccr5[inf.ccr5.DW.W.f] <- "DW"
  uninf.W.f <- which(status == 0 & race == "W" & sex == "F")
  uninf.ccr5.DWDD.W.f <- sample(uninf.W.f, num.uninf.ccr5.DW.W.f + num.uninf.ccr5.DD.W.f)
  uninf.ccr5.DW.W.f <- sample(uninf.ccr5.DWDD.W.f, num.uninf.ccr5.DW.W.f)
  uninf.ccr5.DD.W.f <- setdiff(uninf.ccr5.DWDD.W.f, uninf.ccr5.DW.W.f)
  ccr5[uninf.ccr5.DW.W.f] <- "DW"
  ccr5[uninf.ccr5.DD.W.f] <- "DD"
  }
  
  ##Black Males
  if (num.B.m > 0){
  # homozygotes for deletion
  num.ccr5.DD.B.m <- dat$param$ccr5.B.m.prob[1] * num.B.m
  # heterozygotes
  num.ccr5.DW.B.m <- dat$param$ccr5.B.m.prob[2] * num.B.m
  # homozygotes for deletion
  num.ccr5.WW.B.m <- num.B.m - num.ccr5.DD.B.m - num.ccr5.DW.B.m
  # DD's can't be infected
  num.uninf.ccr5.DD.B.m <- round(num.ccr5.DD.B.m)
  
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.B.m <- round(num.ccr5.DW.B.m * nInfB.m * ccr5.heteroz.rr /
                                 (num.ccr5.WW.B.m + num.ccr5.DW.B.m * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.B.m <- round(num.ccr5.DW.B.m - num.inf.ccr5.DW.B.m)
  inf.B.m <- which(status == 1 & race == "B" & sex == "M")
  inf.ccr5.DW.B.m <- sample(inf.B.m, num.inf.ccr5.DW.B.m, replace = FALSE)
  ccr5[inf.ccr5.DW.B.m] <- "DW"
  uninf.B.m <- which(status == 0 & race == "B" & sex == "M")
  uninf.ccr5.DWDD.B.m <- sample(uninf.B.m, num.uninf.ccr5.DW.B.m + num.uninf.ccr5.DD.B.m)
  uninf.ccr5.DW.B.m <- sample(uninf.ccr5.DWDD.B.m, num.uninf.ccr5.DW.B.m)
  uninf.ccr5.DD.B.m <- setdiff(uninf.ccr5.DWDD.B.m, uninf.ccr5.DW.B.m)
  ccr5[uninf.ccr5.DW.B.m] <- "DW"
  ccr5[uninf.ccr5.DD.B.m] <- "DD"
  }
  
  ##Black immigrant Males
  if (num.BI.m > 0){
  # homozygotes for deletion
  num.ccr5.DD.BI.m <- dat$param$ccr5.BI.m.prob[1] * num.BI.m
  # heterozygotes
  num.ccr5.DW.BI.m <- dat$param$ccr5.BI.m.prob[2] * num.BI.m
  # homozygotes for deletion
  num.ccr5.WW.BI.m <- num.BI.m - num.ccr5.DD.BI.m - num.ccr5.DW.BI.m
  # DD's can't be infected
  num.uninf.ccr5.DD.BI.m <- round(num.ccr5.DD.BI.m)
  
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.BI.m <- round(num.ccr5.DW.BI.m * nInfBI.m * ccr5.heteroz.rr /
                                  (num.ccr5.WW.BI.m + num.ccr5.DW.BI.m * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.BI.m <- round(num.ccr5.DW.BI.m - num.inf.ccr5.DW.BI.m)
  inf.BI.m <- which(status == 1 & race == "BI" & sex == "M")
  inf.ccr5.DW.BI.m <- sample(inf.BI.m, num.inf.ccr5.DW.BI.m, replace = FALSE)
  ccr5[inf.ccr5.DW.BI.m] <- "DW"
  uninf.BI.m <- which(status == 0 & race == "BI" & sex == "M")
  uninf.ccr5.DWDD.BI.m <- sample(uninf.BI.m, num.uninf.ccr5.DW.BI.m + num.uninf.ccr5.DD.BI.m)
  uninf.ccr5.DW.BI.m <- sample(uninf.ccr5.DWDD.BI.m, num.uninf.ccr5.DW.BI.m)
  uninf.ccr5.DD.BI.m <- setdiff(uninf.ccr5.DWDD.BI.m, uninf.ccr5.DW.BI.m)
  ccr5[uninf.ccr5.DW.BI.m] <- "DW"
  ccr5[uninf.ccr5.DD.BI.m] <- "DD"
  }
  
  ##Hisapnic Males
  if (num.H.m > 0){
  # homozygotes for deletion
  num.ccr5.DD.H.m <- dat$param$ccr5.H.m.prob[1] * num.H.m
  # heterozygotes
  num.ccr5.DW.H.m <- dat$param$ccr5.H.m.prob[2] * num.H.m
  # homozygotes for deletion
  num.ccr5.WW.H.m <- num.H.m - num.ccr5.DD.H.m - num.ccr5.DW.H.m
  # DD's can't be infected
  num.uninf.ccr5.DD.H.m <- round(num.ccr5.DD.H.m)
  
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.H.m <- round(num.ccr5.DW.H.m * nInfH.m * ccr5.heteroz.rr /
                                 (num.ccr5.WW.H.m + num.ccr5.DW.H.m * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.H.m <- round(num.ccr5.DW.H.m - num.inf.ccr5.DW.H.m)
  inf.H.m <- which(status == 1 & race == "H" & sex == "M")
  inf.ccr5.DW.H.m <- sample(inf.H.m, num.inf.ccr5.DW.H.m, replace = FALSE)
  ccr5[inf.ccr5.DW.H.m] <- "DW"
  uninf.H.m <- which(status == 0 & race == "H" & sex == "M")
  uninf.ccr5.DWDD.H.m <- sample(uninf.H.m, num.uninf.ccr5.DW.H.m + num.uninf.ccr5.DD.H.m)
  uninf.ccr5.DW.H.m <- sample(uninf.ccr5.DWDD.H.m, num.uninf.ccr5.DW.H.m)
  uninf.ccr5.DD.H.m <- setdiff(uninf.ccr5.DWDD.H.m, uninf.ccr5.DW.H.m)
  ccr5[uninf.ccr5.DW.H.m] <- "DW"
  ccr5[uninf.ccr5.DD.H.m] <- "DD"
  }
  
  ##Hispanic immigrant Males
  if (num.HI.m > 0){
  # homozygotes for deletion
  num.ccr5.DD.HI.m <- dat$param$ccr5.HI.m.prob[1] * num.HI.m
  # heterozygotes
  num.ccr5.DW.HI.m <- dat$param$ccr5.HI.m.prob[2] * num.HI.m
  # homozygotes for deletion
  num.ccr5.WW.HI.m <- num.HI.m - num.ccr5.DD.HI.m - num.ccr5.DW.HI.m
  # DD's can't be infected
  num.uninf.ccr5.DD.HI.m <- round(num.ccr5.DD.HI.m)
  
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.HI.m <- round(num.ccr5.DW.HI.m * nInfHI.m * ccr5.heteroz.rr /
                                  (num.ccr5.WW.HI.m + num.ccr5.DW.HI.m * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.HI.m <- round(num.ccr5.DW.HI.m - num.inf.ccr5.DW.HI.m)
  inf.HI.m <- which(status == 1 & race == "HI" & sex == "M")
  inf.ccr5.DW.HI.m <- sample(inf.HI.m, num.inf.ccr5.DW.HI.m, replace = FALSE)
  ccr5[inf.ccr5.DW.HI.m] <- "DW"
  uninf.HI.m <- which(status == 0 & race == "HI" & sex == "M")
  uninf.ccr5.DWDD.HI.m <- sample(uninf.HI.m, num.uninf.ccr5.DW.HI.m + num.uninf.ccr5.DD.HI.m)
  uninf.ccr5.DW.HI.m <- sample(uninf.ccr5.DWDD.HI.m, num.uninf.ccr5.DW.HI.m)
  uninf.ccr5.DD.HI.m <- setdiff(uninf.ccr5.DWDD.HI.m, uninf.ccr5.DW.HI.m)
  ccr5[uninf.ccr5.DW.HI.m] <- "DW"
  ccr5[uninf.ccr5.DD.HI.m] <- "DD"
  }
  
  ## White Males
  if (num.W.m > 0){
  # homozygotes for deletion
  num.ccr5.DD.W.m <- dat$param$ccr5.W.m.prob[1] * num.W.m
  # heterozygotes
  num.ccr5.DW.W.m <- dat$param$ccr5.W.m.prob[2] * num.W.m
  # homozygotes for deletion
  num.ccr5.WW.W.m <- num.W.m - num.ccr5.DD.W.m - num.ccr5.DW.W.m
  # DD's can't be infected
  num.uninf.ccr5.DD.W.m <- round(num.ccr5.DD.W.m)
  
  # Unique solution to get relative risk right in init pop
  num.inf.ccr5.DW.W.m <- round(num.ccr5.DW.W.m * nInfW.m * ccr5.heteroz.rr /
                                 (num.ccr5.WW.W.m + num.ccr5.DW.W.m * ccr5.heteroz.rr))
  num.uninf.ccr5.DW.W.m <- round(num.ccr5.DW.W.m - num.inf.ccr5.DW.W.m)
  inf.W.m <- which(status == 1 & race == "W" & sex == "M")
  inf.ccr5.DW.W.m <- sample(inf.W.m, num.inf.ccr5.DW.W.m)
  ccr5[inf.ccr5.DW.W.m] <- "DW"
  uninf.W.m <- which(status == 0 & race == "W" & sex == "M")
  uninf.ccr5.DWDD.W.m <- sample(uninf.W.m, num.uninf.ccr5.DW.W.m + num.uninf.ccr5.DD.W.m)
  uninf.ccr5.DW.W.m <- sample(uninf.ccr5.DWDD.W.m, num.uninf.ccr5.DW.W.m)
  uninf.ccr5.DD.W.m <- setdiff(uninf.ccr5.DWDD.W.m, uninf.ccr5.DW.W.m)
  ccr5[uninf.ccr5.DW.W.m] <- "DW"
  ccr5[uninf.ccr5.DD.W.m] <- "DD"
  }
  
  dat$attr$ccr5 <- ccr5

  return(dat)
}


#' @title Re-Initialization Module
#'
#' @description This function reinitializes an epidemic model to restart at a
#'              specified time step given an input \code{netsim} object.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netsim}}.
#' @inheritParams initialize_msm
#'
#' @details
#' Currently, the necessary components that must be on \code{x} for a simulation
#' to be restarted must be: param, control, nwparam, epi, attr, temp, el, p.
#' TODO: describe this more.
#'
#' @return
#' This function resets the data elements on the \code{dat} master data object
#' in the needed ways for the time loop to function.
#'
#' @export
#' @keywords module shamp
#'
reinit_shamp <- function(x, param, init, control, s) {

  need.for.reinit <- c("param", "control", "nwparam", "epi", "attr", "temp", "el", "p")
  if (!all(need.for.reinit %in% names(x))) {
    stop("x must contain the following elements for restarting: ",
         "param, control, nwparam, epi, attr, temp, el, p",
         call. = FALSE)
  }

  if (length(x$el) == 1) {
    s <- 1
  }

  dat <- list()

  dat$param <- param
  dat$param$modes <- 1
  dat$control <- control
  dat$nwparam <- x$nwparam

  dat$epi <- sapply(x$epi, function(var) var[s])
  names(dat$epi) <- names(x$epi)

  dat$el <- x$el[[s]]
  dat$p <- x$p[[s]]

  dat$attr <- x$attr[[s]]

  if (!is.null(x$stats)) {
    dat$stats <- list()
    if (!is.null(x$stats$nwstats)) {
      dat$stats$nwstats <- x$stats$nwstats[[s]]
    }
  }

  dat$temp <- x$temp[[s]]

  class(dat) <- "dat"
  return(dat)
}

