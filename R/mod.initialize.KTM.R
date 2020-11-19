
# SHAMP ego -----------------------------------------------------------------

#' @title Initialization Module for up to 5 race groups heterosexuals and MSM.
#'
#' @description This function initializes the master \code{dat} object on which
#'              data are stored, initial state of the network is simulated using ergm.ego, and
#'              simulates disease status and other attributes are initialized.
#'
#' @param x An \code{EpiModel} object of class \code{\link{netest}}.
#' @param param An \code{EpiModel} object of class \code{\link{param_KTM}}.
#' @param init An \code{EpiModel} object of class \code{\link{init_KTM}}.
#' @param control An \code{EpiModel} object of class \code{\link{control_KTM}}.
#' @param s Simulation number, used for restarting dependent simulations.
#'
#' @return
#' This function returns the updated \code{dat} object with the initialized values
#' for demographics and disease-related variables.
#'
#' @export
#' @keywords module HET MSM ego 
#'


initialize_KTM <- function(x, param, init, control, s) {

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
  dat$trans.el <- list()
  dat$death.stats <- list()
  
  ## Network simulation ##
 # Do we need to create the nw in each environment?
  #nw<-x[[1]]$fit$network
  #count <- network.edgecount(nw)
  #delete.edges(nw,1:count)
  #environment(x[[2]]$fit$formula) <- environment()
  #environment(x[[3]]$fit$formula) <- environment()
  #x[[1]]$fit$network<-NULL
  
  
  nwl <- list()
  for (i in 1:3) {

        nwl[[i]] <- simulate(x[[i]]$fit)
  }
  
  
  nw<-nwl
  nwl<-NULL

  ## ergm_prep here
  dat$el <- list()
  dat$p <- list()
  for (i in 1:2) {
    dat$el[[i]] <- as.edgelist(nw[[i]])
    
       attributes(dat$el[[i]])$vnames <- NULL
    p <- stergm_prep(nw[[i]], x[[i]]$formation, x[[i]]$coef.diss$dissolution,
                                x[[i]]$coef.form, x[[i]]$coef.diss$coef.adj, x[[i]]$constraints)
    p$model.form$formula <- NULL
    p$model.diss$formula <- NULL
    dat$p[[i]] <- p
  }
  dat$el[[3]] <- as.edgelist(nw[[3]])
  attributes(dat$el[[3]])$vnames <- NULL
  p <- ergm_prep(nw[[3]], x[[3]]$formation, x[[3]]$coef.form, x[[3]]$constraints)
  p$model.form$formula <- NULL
  dat$p[[3]] <- p
  

  # Network parameters
  dat$nwparam <- list()
  for (i in 1:3) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }

  ## Nodal attributes ##

  # Degree terms
  dat$attr$deg.pers <- get_degree(dat$el[[2]])
  dat$attr$deg.cohab <- get_degree(dat$el[[1]])
  
  dat$attr$deg.cohab.c <- ifelse(dat$attr$deg.cohab > 0,1,dat$attr$deg.cohab)
  dat$attr$deg.pers.c <- ifelse(dat$attr$deg.pers > 0,1,dat$attr$deg.pers)
  
  dat$attr$deg.cohab.c_ee <- dat$attr$deg.cohab.c 
  dat$attr$deg.pers.c_ee <- dat$attr$deg.pers.c
  
  dat$attr$deg.inst <- get_degree(dat$el[[3]])
  dat$attr$deg.tot<-dat$attr$deg.cohab + dat$attr$deg.pers + dat$attr$deg.inst
  
  # Race
  dat$attr$race <- rep("B", length(dat$attr$deg.pers.c))
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
 


  dat$attr$cohab.lt <- rep(0, num)
  dat$attr$pers.lt <- rep(0, num)
  dat$attr$onetime.lt <-  rep(0, num)
  
  # Sex Identity
  dat$attr$sex.ident <- get.vertex.attribute(nw[[1]], "sex")
  dat$attr$sex.ident[dat$attr$sex.ident=="M"] <- "msf"
  dat$attr$sex.ident[dat$attr$sex.ident=="F"] <- "f"
  dat$attr$msmf <- ifelse(dat$attr$sex.ident == "msmf" ,1 ,0)
  num.msf <-sum(dat$attr$sex.ident == "msf")
  num.msm <-sum(dat$attr$sex.ident == "msm")
  num.msmf <-sum(dat$attr$sex.ident == "msmf")
  num.f <-sum(dat$attr$sex.ident == "f")
  
  dat$attr$active <- rep(1, num)
  dat$attr$uid <- 1:num
  dat$temp$max.uid <- num

  # Population of interest
  dat$attr$poi <- get.vertex.attribute(nw[[1]], "poi")
  
  # Age
  dat$attr$age.adj <- get.vertex.attribute(nw[[1]], "age.adj")
  dat$attr$age <- get.vertex.attribute(nw[[1]], "est_age")
  
  partial<-(0:51)* (dat$param$time.unit / 365)
  partial<-sample(partial,length(dat$attr$age),replace=TRUE)
  
  dat$attr$age<-dat$attr$age+partial
  dat$attr$age.adj<-dat$attr$age.adj+partial
  dat$attr$sqrt.age <- sqrt(dat$attr$age)
  dat$attr$agesq <- dat$attr$age^2
  
  dat$attr$age.group <- get.vertex.attribute(nw[[1]], "age_grp")
  #dat$attr$age.adj.t<-ifelse(dat$attr$sex=="M",dat$attr$age,
  #                              ifelse(dat$attr$sex=="F",dat$attr$age + dat$param$age.adj, dat$attr$age))
  
  dat$attr$age.inf <-rep(NA, num)
  dat$attr$age.diag <-rep(NA, num)
  
  #sex.age.group.
  #dat$attr$sex.age.group <- get.vertex.attribute(nw[[1]], "sex.age.group")
  
  # Risk group
  dat$attr$riskg <- get.vertex.attribute(nw[[1]], "riskg")
  
  # Immigrant status
  # dat$attr$immig.loc <- rep(0,length(dat$attr$age))
  #dat$attr$arv.BI.pos <- rep(0,length(dat$attr$age))
  #dat$attr$arv.HI.pos <- rep(0,length(dat$attr$age))
  
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
  dat$attr$cond.always.main.het <- rbinom(num,1,dat$param$cond.main.always.prob.het)

  # Arrival and departure
  dat$attr$arrival.time <- rep(1, num)

  # Circumcision
  circ <- rep(0, num)

  if (length(ids.B.m) > 0) {
  circ[ids.B.m] <- sample(apportion_lr(length(ids.B.m), 0:1, 1 - param$circ.B.prob))}


  dat$attr$circ <- circ

  # PrEP Attributes
  dat$attr$prepClass <- rep(NA, num)
  dat$attr$prepElig <- rep(NA, num)
  dat$attr$prepStat <- rep(0, num)
  dat$attr$prep.elig.time <- rep(NA, num) 

  # One-off AI class
#  inst.ai.class <- rep(NA, num)
#  ncl <- param$num.inst.ai.classes
#  inst.ai.class[ids.B] <- sample(apportion_lr(num.B, 1:ncl, rep(1 / ncl, ncl)))
#  inst.ai.class[ids.W] <- sample(apportion_lr(num.W, 1:ncl, rep(1 / ncl, ncl)))
#  dat$attr$inst.ai.class <- inst.ai.class

  # Role class

  dat$attr$role.class <- get.vertex.attribute(nw[[1]], "sex")
  dat$attr$role.class[dat$attr$role.class=="M"] <- "I"
  dat$attr$role.class[dat$attr$role.class=="F"] <- "R"
  role.class <- dat$attr$role.class
  
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
  
  #Kenya TM intervention
  dat$attr$partner.serv.part <- rep(0,num)
  dat$attr$partner.serv.part.time <- rep(0,num)
  dat$attr$PS.index.acute <- rep(0,num)
  dat$attr$PS.index.prev <- rep(0,num)
  dat$attr$PS.diag.neg <- rep(0,num)
  dat$attr$PS.diag.pos.time <- rep(0,num)
  
  dat$attr$cat.5.status  <- rep(0,num)
  dat$attr$testclin  <- rep(0,num)
  
  # Prevalence Tracking
  dat$temp$deg.dists <- list()
  dat$temp$discl.list <- matrix(NA, nrow = 0, ncol = 3)
  colnames(dat$temp$discl.list) <- c("pos", "neg", "discl.time")

  if (control$save.nwstats == TRUE) {
    dat$stats <- list()
    dat$stats$nwstats <- list()

  }

  dat <- prevalence_KTM(dat, at = 1)

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
  nInfB.msf <- ifelse(num.B.msf > 0 & dat$init$prev.B.msf > 0 , max(1,round(dat$init$prev.B.msf * num.B.msf)),round(dat$init$prev.B.msf * num.B.msf))


  # Age-based infection probability
  probInfCrB.f <- age[ids.B.f] * dat$init$init.prev.age.slope.B.f
  probInfB.f <- probInfCrB.f + (nInfB.f - sum(probInfCrB.f)) / num.B.f
  

  
  ##Male msf.
  probInfCrB.msf <- age[ids.B.msf] * dat$init$init.prev.age.slope.B.msf
  probInfB.msf <- probInfCrB.msf + (nInfB.msf - sum(probInfCrB.msf)) / num.B.msf
  


  if(num.het>0){
    if (any(probInfB.f < 0, probInfB.msf < 0)) {
      stop("Slope of initial prevalence by age must be sufficiently low to ",
           "avoid non-positive probabilities.", call. = FALSE)
    }
  }
  
  
#  if(num.het>0){
#  if (any(probInfB.f < 0 , probInfBI.f < 0 , probInfH.f < 0 , probInfHI.f < 0 , probInfW.f < 0 , 
#          probInfB.msf < 0 , probInfBI.msf < 0 , probInfH.msf < 0 , probInfHI.msf < 0 , probInfW.msf < 0)) {
#    stop("Slope of initial prevalence by age must be sufficiently low to ",
#         "avoid non-positive probabilities.", call. = FALSE)
#  }
#  }
  
#  if(num.msm>0){
#  if (any(probInfB.msm < 0 , probInfBI.msm < 0 , probInfH.msm < 0 , probInfHI.msm < 0 , probInfW.msm < 0)) {
#    stop("Slope of initial prevalence by age must be sufficiently low to ",
#         "avoid non-positive probabilities.", call. = FALSE)
#  }
#  }  


#  if(num.bi>0){
#    if (any(probInfB.msmf < 0 , probInfBI.msmf < 0 , probInfH.msmf < 0 , probInfHI.msmf < 0 , probInfW.msmf < 0)) {
#      stop("Slope of initial prevalence by age must be sufficiently low to ",
#           "avoid non-positive probabilities.", call. = FALSE)
#    }
#  }  
  
 
   # Infection status
   status <- rep(0, num)
   
   #Females.
    while (sum(status[ids.B.f]) != nInfB.f) {
      status[ids.B.f] <- rbinom(num.B.f, 1, probInfB.f)
    }

   #Males msf.
   while (sum(status[ids.B.msf]) != nInfB.msf) {
     status[ids.B.msf] <- rbinom(num.B.msf, 1, probInfB.msf)
   }



  dat$attr$status <- status
  dat$attr$infected.gen <- rep(NA,length(status))
  infected<-which(dat$attr$status==1)
  dat$attr$infected.gen[infected] <-0
  # Treatment trajectory
  tt.traj <- rep(NA, num)

  #female
  if (length(ids.B.f) > 0){ 
  tt.traj[ids.B.f] <- sample(apportion_lr(num.B.f, c(1, 2, 3, 4),
                                        dat$param$tt.traj.B.f.prob))}

  #Males MSF. 
  if (length(ids.B.msf) > 0){ 
  tt.traj[ids.B.msf] <- sample(apportion_lr(num.B.msf, c(1, 2, 3, 4),
                                          dat$param$tt.traj.B.msf.prob))}


  
  dat$attr$tt.traj <- tt.traj



  ## Infection-related attributes

  ##Edge list of each transmission event.
  ##Infector, infected, time)
  dat$trans.el$infector <- dat$trans.el$infector.sex <- dat$trans.el$infector.race <- dat$trans.el$infector.age <- NULL
  dat$trans.el$infector.sex.ident <- dat$trans.el$infector.deg.tot <- dat$trans.el$infector.gen <-NULL
  dat$trans.el$infected <- dat$trans.el$infected.sex <- dat$trans.el$infected.race <- dat$trans.el$infected.age <- NULL 
  dat$trans.el$infected.sex.ident <- dat$trans.el$infected.deg.tot <- dat$trans.el$infected.gen <-NULL 
  dat$trans.el$time <-NULL
  
  stage <- rep(NA, num)
  stage.time <- rep(NA, num)
  stage4.onset.time <- rep(NA, num)
  inf.time <- rep(NA, num)
  vl <- rep(NA, num)
  diag.status <- rep(NA, num)
  diag.time <- rep(NA, num)
  PS.diag.pos.time <- rep(NA, num)
  PS.diag.neg <- rep(NA, num)
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


  



  twind.int <- dat$param$test.window.int.ab*7
 
  ### Full adherent type

  # Create set of expected values for (cum.time.off.tx, cum.time.on.tx)

  tx.init.time.B.f <- twind.int + dat$param$last.neg.test.B.f.int + 1 / dat$param$tx.init.B.f.prob
  tx.init.time.B.msf <- twind.int + dat$param$last.neg.test.B.msf.int + 1 / dat$param$tx.init.B.msf.prob



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
  
  
  # Diagnosis
  selected <- which(status == 1 & tt.traj == 4)
    ttntest <- rgeom(length(selected),
                     1 / (dat$param$mean.test.B.f.int * (race[selected] == "B" & sex[selected] == "F") +
                          dat$param$mean.test.B.msf.int * (race[selected] == "B" & sex.ident[selected] == "msf")))

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

  # Implement diagnosis for both
  selected <- which(status == 1 & tt.traj == 3)

    ttntest <- rgeom(length(selected),
                     1 / (dat$param$mean.test.B.f.int * (race[selected] == "B" & sex[selected] == "F") +
                          dat$param$mean.test.B.msf.int * (race[selected] == "B" & sex.ident[selected] == "msf")))


  diag.status[selected][ttntest > cum.time.off.tx[selected] - twind.int] <- 0
  last.neg.test[selected][ttntest > cum.time.off.tx[selected] - twind.int] <-
    -ttntest[ttntest > cum.time.off.tx[selected] - twind.int]

  diag.status[selected][ttntest <= cum.time.off.tx[selected] - twind.int] <- 1
  diag.status[selected][cum.time.on.tx[selected] > 0] <- 1
  last.neg.test[selected][cum.time.on.tx[selected] > 0] <- NA


  # Last neg test before present for negatives
  selected <- which(status == 0 & tt.traj %in% c(2, 3, 4))


    tslt <- rgeom(length(selected),
                  1 / (dat$param$mean.test.B.f.int * (race[selected] == "B" & sex[selected] == "F") +
                       dat$param$mean.test.B.msf.int * (race[selected] == "B" & sex.ident[selected] == "msf")))
  last.neg.test[selected] <- -tslt

  
  ##Infection Class
  selected <- which(status == 1)
  inf.class[selected] <- "Lhet"

  ## Set all onto dat$attr
  dat$attr$stage <- stage
  dat$attr$stage.time <- stage.time
  dat$attr$stage4.onset.time <- stage4.onset.time 
  dat$attr$inf.time <- inf.time
  dat$attr$vl <- vl
  dat$attr$diag.status <- diag.status
  dat$attr$diag.time <- diag.time
  dat$attr$last.neg.test <- last.neg.test
  dat$attr$evertest <-rep(0, num)
  dat$attr$evertest[dat$attr$diag.status==1] <- 1
  #dat$attr$evertest<-ifelse(dat$attr$last.neg.test<=0 | dat$attr$diag.status==1,1,0)
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
  
 
    dat$trans.el$infected <- c(dat$trans.el$infected, dat$attr$uid[dat$attr$status == 1]) 
    dat$trans.el$infected.gen <- c(dat$trans.el$infected.gen, rep(1, sum(dat$attr$status == 1)))
    dat$trans.el$infected.sex <- c(dat$trans.el$infected.sex, dat$attr$sex[dat$attr$status == 1])
    dat$trans.el$infected.race <- c(dat$trans.el$infected.race, dat$attr$race[dat$attr$status == 1])
    dat$trans.el$infected.age <- c(dat$trans.el$infected.age, dat$attr$age[dat$attr$status == 1])
    dat$trans.el$infected.sex.ident <- c(dat$trans.el$infected.sex.ident, dat$attr$sex.ident[dat$attr$status == 1])
    
    dat$trans.el$infector <- c(dat$trans.el$infector, rep("SEED", sum(dat$attr$status, na.rm = TRUE)))
    dat$trans.el$infector.gen <- c(dat$trans.el$infector.gen, rep(0, sum(dat$attr$status, na.rm = TRUE)))
    dat$trans.el$infector.sex <- c(dat$trans.el$infector.sex, rep("SEED", sum(dat$attr$status, na.rm = TRUE)))
    dat$trans.el$infector.race <- c(dat$trans.el$infector.race, rep("SEED", sum(dat$attr$status, na.rm = TRUE)))
    dat$trans.el$infector.age <- c(dat$trans.el$infector.age, rep("SEED", sum(dat$attr$status, na.rm = TRUE)))
    dat$trans.el$infector.sex.ident <- c(dat$trans.el$infector.sex.ident, rep("SEED", sum(dat$attr$status, na.rm = TRUE)))  
    
    dat$trans.el$time <- c(dat$trans.el$time, rep(0, sum(dat$attr$status, na.rm = TRUE)))
    


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
  status<-dat$attr$status

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
#' @inheritParams initialize_KTM
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
#' @keywords module KTM
#'
reinit_shamp <- function(x, param, init, control, s) {

  need.for.reinit <- c("param", "control", "nwparam", "epi", "attr", "temp", "el", "p", "trans.el")
  if (!all(need.for.reinit %in% names(x))) {
    stop("x must contain the following elements for restarting: ",
         "param, control, nwparam, epi, attr, temp, el, p, trans.el",
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
  dat$trans.el <- x$trans.el[[s]]

  if (!is.null(x$stats)) {
    dat$stats <- list()
    if (!is.null(x$stats$nwstats)) {
      dat$stats$nwstats <- x$stats$nwstats[[s]]
    }
  }
  
  if (!is.null(x$cel.temp)) {
    dat$cel.temp <- x$cel.temp[[s]]
  }
  
  if (!is.null(x$cel.complete)) {
    dat$cel.complete <- x$cel.complete[[s]]
  }
  
  if (!is.null(x$death.stats)) {
    dat$death.stats <- x$death.stats[[s]]
  }

  dat$temp <- x$temp[[s]]

  class(dat) <- "dat"
  return(dat)
}

