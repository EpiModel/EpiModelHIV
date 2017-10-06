
#' @title Births Module for up to 5 race groups heterosexuals and MSM for preserving demographic distribution.
#'
#' @description Module function for births or entries into the sexually active
#'              population.
#'
#' @inheritParams aging_msm
#'
#' @details
#' New population members are added based on expected numbers of entries among
#' all five race/immigrant groups and two sexes, stochastically determined with draws from Poisson
#' distributions. The proportion of men who are MSM are determined by \code{msm.frac}. 
#' The proportion of men who are MSMF are determined by \code{msmf.frac}. 
#' For each new entry, a set of attributes is added for that node,
#' and the nodes are added onto the network objects. Only attributes that are
#' a part of the network model formula are updated as vertex attributes on the
#' network objects.
#'
#' @return
#' This function updates the \code{attr} list with new attributes for each new
#' population member, and the \code{nw} objects with new vertices.
#'
#' @keywords module msm het
#' @export
#'
births_shamp2 <- function(dat, at){

  if(at > 1){
    
  ## Variables

  # Parameters
  demog.list <- dat$param$demog.list
  demog.dist <- dat$param$demog.dist



  ## Process

  nBirths.gen <- dat$epi$dth.gen[at] 
  nBirths.dis <- dat$epi$dth.dis[at]
  nBirths.age <- dat$epi$dth.age[at]
  

  ##For now we will not replace disease deaths.
  nBirths <- nBirths.age + nBirths.gen

  ## Update Attr

  if (nBirths.gen > 0 | nBirths.dis > 0 | nBirths.age > 0) {
    dat <- setBirthAttr_shamp(dat, at, nBirths.gen, nBirths.age, nBirths.dis)
  }

 
  
  # Update Networks
  if (nBirths > 0) {
    for (i in 1:3) {
      dat$el[[i]] <- add_vertices(dat$el[[i]], nBirths)
    }
  }


  ## Output
  dat$epi$nBirths[at] <- nBirths

  return(dat)
  }
}


setBirthAttr_shamp <- function(dat, at, nBirths.gen, nBirths.age, nBirths.dis) {

  # Parameters
  demog.list <- dat$param$demog.list
  demog.dist <- dat$param$demog.dist
  sex.groups <- dat$param$sex.groups
  race.groups <- dat$param$race.groups
  age.groups <- dat$param$age.groups
  
  ##Not replcing disease deathe at this time
  nBirths <- nBirths.gen + nBirths.age

  # Set all attributes NA by default
  dat$attr <- lapply(dat$attr, {
    function(x)
      c(x, rep(NA, nBirths))
  })
  
  newIds <- which(is.na(dat$attr$active))

  # Demographic
  dat$attr$active[newIds] <- rep(1, nBirths)
  dat$attr$uid[newIds] <- dat$temp$max.uid + (1:nBirths)
  dat$temp$max.uid <- dat$temp$max.uid + nBirths
  dat$attr$sex.ident[newIds] <- rep("NA", nBirths)
  dat$attr$immig.loc[newIds] <- rep(0, nBirths)
  dat$attr$sex[newIds] <- rep(0, nBirths)
  dat$attr$race[newIds] <- rep(0, nBirths)
  dat$attr$age[newIds] <- rep(0, nBirths)
  dat$attr$demog.cat[newIds] <- rep(0, nBirths)
  dat$attr$evertest[newIds] <- rep(0, nBirths)

  dat$attr$arrival.time[newIds] <- rep(at, nBirths)

 
  demog.cat <- sample(demog.list,size=nBirths,prob=demog.dist,replace=TRUE)
  dat$attr$demog.cat[newIds]<-demog.cat

  sex<-rep(NA,nBirths)
  race<-rep(NA,nBirths)
  age<-rep(NA,nBirths)
  
  sex.temp<-as.numeric(substr(demog.cat,1,1))
  race.temp<-as.numeric(substr(demog.cat,2,2))
  age.temp<-as.numeric(substr(demog.cat,3,4))-(dat$param$birth.age-1)
  
  for (i in 1:length(demog.cat)){
    sex[i]<-sex.groups[sex.temp[i]]
    race[i]<-race.groups[race.temp[i]]
    age[i]<-age.groups[age.temp[i]]
  }
 
  
  #set age for nBirths.age to 18.
  ids.birth<-sample(1:length(age),nBirths.age,replace=FALSE)
  age[ids.birth]<-dat$param$birth.age
  

  dat$attr$sex[newIds] <- sex
  dat$attr$race[newIds] <- race
  dat$attr$age[newIds] <- age
  dat$attr$sqrt.age[newIds] <- sqrt(dat$attr$age[newIds]) 
  
  dat$attr$sqrt.age.adj<-ifelse(dat$attr$sex=="M",dat$attr$sqrt.age,
                                ifelse(dat$attr$sex=="F",dat$attr$sqrt.age + dat$param$age.adj,dat$attr$sqrt.age))
  
  #Set agecat
  
  dat$attr$agecat[newIds]<-ifelse(dat$att$age[newIds] < 26 ,"18-25",
                          ifelse(dat$attr$age[newIds] > 25 & dat$attr$age[newIds] < 36,"26-35",
                                 ifelse(dat$attr$age[newIds] > 35 & dat$attr$age[newIds] < 46,"36-45",
                                        ifelse(dat$attr$age[newIds] > 45,"46-60",dat$attr$agecat[newIds]))))
                      
  
  newF<-which(sex=="F")
  newM<-which(sex=="M")

  
  newB<-which(race=="B")
  newBI<-which(race=="BI")
  newH<-which(race=="H")
  newHI<-which(race=="HI")
  newW<-which(race=="W")
  
  sex.ident<-rep(NA,nBirths)
  sex.ident[newF]<-"f"
  sex.ident[newM]<-sample(c("msm","msmf","msf"),length(newM),
                                    prob=c(dat$param$msm.frac,dat$param$msmf.frac,(1-(dat$param$msm.frac+dat$param$msmf.frac))),
                                           replace=TRUE)
  dat$attr$sex.ident[newIds] <- sex.ident
  new.msm<-which(sex.ident=="msm")
  new.msmf<-which(sex.ident=="msmf")
  new.msf<-which(sex.ident=="msf")
 

 

  dat$attr$status[newIds] <- rep(0, nBirths)
  dat$attr$aq.class[newIds] <- rep(NA, nBirths)

 
  # dat$attr$inst.ai.class[newIds] <- sample(1:dat$param$num.inst.ai.classes,
  #                                         nBirths, replace = TRUE)
 
   #### females f. 
  
  if (length(newF) > 0) {
  dat$attr$tt.traj[newIds[intersect(newB,newF)]] <- sample(c(1, 2, 3, 4),
                                           length(intersect(newB,newF)), replace = TRUE,
                                           prob = dat$param$tt.traj.B.f.prob)
  dat$attr$tt.traj[newIds[intersect(newBI,newF)]] <- sample(c(1, 2, 3, 4),
                                           length(intersect(newBI,newF)), replace = TRUE,
                                           prob = dat$param$tt.traj.BI.f.prob)
  dat$attr$tt.traj[newIds[intersect(newH,newF)]] <- sample(c(1, 2, 3, 4),
                                           length(intersect(newH,newF)), replace = TRUE,
                                           prob = dat$param$tt.traj.H.f.prob)
  dat$attr$tt.traj[newIds[intersect(newHI,newF)]] <- sample(c(1, 2, 3, 4),
                                           length(intersect(newHI,newF)), replace = TRUE,
                                           prob = dat$param$tt.traj.HI.f.prob)
  dat$attr$tt.traj[newIds[intersect(newW,newF)]] <- sample(c(1, 2, 3, 4),
                                           length(intersect(newW,newF)), replace = TRUE,
                                           prob = dat$param$tt.traj.W.f.prob)
  }
  
  #### Males msf. 
  if (length(new.msf) > 0) {
  dat$attr$tt.traj[newIds[intersect(newB,new.msf)]] <- sample(c(1, 2, 3, 4),
                                                  length(intersect(newB,new.msf)), replace = TRUE,
                                                 prob = dat$param$tt.traj.B.msf.prob)
  dat$attr$tt.traj[newIds[intersect(newBI,new.msf)]] <- sample(c(1, 2, 3, 4),
                                                  length(intersect(newBI,new.msf)), replace = TRUE,
                                                  prob = dat$param$tt.traj.BI.msf.prob)
  dat$attr$tt.traj[newIds[intersect(newH,new.msf)]] <- sample(c(1, 2, 3, 4),
                                                  length(intersect(newH,new.msf)), replace = TRUE,
                                                 prob = dat$param$tt.traj.H.msf.prob)
  dat$attr$tt.traj[newIds[intersect(newHI,new.msf)]] <- sample(c(1, 2, 3, 4),
                                                  length(intersect(newHI,new.msf)), replace = TRUE,
                                                  prob = dat$param$tt.traj.HI.msf.prob)
  dat$attr$tt.traj[newIds[intersect(newW,new.msf)]] <- sample(c(1, 2, 3, 4),
                                                  length(intersect(newW,new.msf)), replace = TRUE,
                                                 prob = dat$param$tt.traj.W.msf.prob)
  }
  
  #### Males msm. 
  if(length(new.msm) > 0){
  dat$attr$tt.traj[newIds[intersect(newB,new.msm)]] <- sample(c(1, 2, 3, 4),
                                                      length(intersect(newB,new.msm)), replace = TRUE,
                                                      prob = dat$param$tt.traj.B.msm.prob)
  dat$attr$tt.traj[newIds[intersect(newBI,new.msm)]] <- sample(c(1, 2, 3, 4),
                                                      length(intersect(newBI,new.msm)), replace = TRUE,
                                                       prob = dat$param$tt.traj.BI.msm.prob)
  dat$attr$tt.traj[newIds[intersect(newH,new.msm)]] <- sample(c(1, 2, 3, 4),
                                                      length(intersect(newH,new.msm)), replace = TRUE,
                                                      prob = dat$param$tt.traj.H.msm.prob)
  dat$attr$tt.traj[newIds[intersect(newHI,new.msm)]] <- sample(c(1, 2, 3, 4),
                                                      length(intersect(newHI,new.msm)), replace = TRUE,
                                                       prob = dat$param$tt.traj.HI.msm.prob)
  dat$attr$tt.traj[newIds[intersect(newW,new.msm)]] <- sample(c(1, 2, 3, 4),
                                                      length(intersect(newW,new.msm)), replace = TRUE,
                                                      prob = dat$param$tt.traj.W.msm.prob)
  }
  
  #### Males msmf. 
  if(length(new.msmf) > 0){
    dat$attr$tt.traj[newIds[intersect(newB,new.msmf)]] <- sample(c(1, 2, 3, 4),
                                                            length(intersect(newB,new.msmf)), replace = TRUE,
                                                            prob = dat$param$tt.traj.B.msmf.prob)
    dat$attr$tt.traj[newIds[intersect(newBI,new.msmf)]] <- sample(c(1, 2, 3, 4),
                                                             length(intersect(newBI,new.msmf)), replace = TRUE,
                                                             prob = dat$param$tt.traj.BI.msmf.prob)
    dat$attr$tt.traj[newIds[intersect(newH,new.msmf)]] <- sample(c(1, 2, 3, 4),
                                                            length(intersect(newH,new.msmf)), replace = TRUE,
                                                            prob = dat$param$tt.traj.H.msmf.prob)
    dat$attr$tt.traj[newIds[intersect(newHI,new.msmf)]] <- sample(c(1, 2, 3, 4),
                                                             length(intersect(newHI,new.msmf)), replace = TRUE,
                                                             prob = dat$param$tt.traj.HI.msmf.prob)
    dat$attr$tt.traj[newIds[intersect(newW,new.msmf)]] <- sample(c(1, 2, 3, 4),
                                                            length(intersect(newW,new.msmf)), replace = TRUE,
                                                            prob = dat$param$tt.traj.W.msmf.prob)
  }
  
  
  # Circumcision
  dat$attr$circ[newIds[intersect(newB,newM)]] <- rbinom(length(intersect(newB,newM)), 1, dat$param$circ.B.prob)
  dat$attr$circ[newIds[intersect(newBI,newM)]] <- rbinom(length(intersect(newBI,newM)), 1, dat$param$circ.BI.prob)
  dat$attr$circ[newIds[intersect(newH,newM)]] <- rbinom(length(intersect(newH,newM)), 1, dat$param$circ.H.prob)
  dat$attr$circ[newIds[intersect(newHI,newM)]] <- rbinom(length(intersect(newHI,newM)), 1, dat$param$circ.HI.prob)
  dat$attr$circ[newIds[intersect(newW,newM)]] <- rbinom(length(intersect(newW,newM)), 1, dat$param$circ.W.prob)
  
  dat$attr$circ[newIds[newF]] <- 0
  
  # Role
  dat$attr$role.class[newIds[newF]] <- "R"
  dat$attr$role.class[newIds[newM]] <- "I"
  
  if (length(new.msm) > 0) {

  
  dat$attr$role.class[newIds[intersect(newB,new.msm)]] <- sample(c("I", "R", "V"),
                                              length(intersect(newB,new.msm)), replace = TRUE,
                                              prob = dat$param$role.B.msm.prob)
  dat$attr$role.class[newIds[intersect(newBI,new.msm)]] <- sample(c("I", "R", "V"),
                                              length(intersect(newBI,new.msm)), replace = TRUE,
                                              prob = dat$param$role.BI.msm.prob)
  dat$attr$role.class[newIds[intersect(newH,new.msm)]] <- sample(c("I", "R", "V"),
                                              length(intersect(newH,new.msm)), replace = TRUE,
                                              prob = dat$param$role.H.msm.prob)
  dat$attr$role.class[newIds[intersect(newHI,new.msm)]] <- sample(c("I", "R", "V"),
                                              length(intersect(newHI,new.msm)), replace = TRUE,
                                               prob = dat$param$role.HI.msm.prob)
  dat$attr$role.class[newIds[intersect(newW,new.msm)]] <- sample(c("I", "R", "V"),
                                              length(intersect(newW,new.msm)), replace = TRUE,
                                              prob = dat$param$role.W.msm.prob)
  
  }
  
  if (length(new.msmf) > 0) {
  dat$attr$role.class[newIds[intersect(newB,new.msmf)]] <- sample(c("I", "R", "V"),
                                                             length(intersect(newB,new.msmf)), replace = TRUE,
                                                             prob = dat$param$role.B.msmf.prob)
  dat$attr$role.class[newIds[intersect(newBI,new.msmf)]] <- sample(c("I", "R", "V"),
                                                              length(intersect(newBI,new.msmf)), replace = TRUE,
                                                              prob = dat$param$role.BI.msmf.prob)
  dat$attr$role.class[newIds[intersect(newH,new.msmf)]] <- sample(c("I", "R", "V"),
                                                             length(intersect(newH,new.msmf)), replace = TRUE,
                                                             prob = dat$param$role.H.msmf.prob)
  dat$attr$role.class[newIds[intersect(newHI,new.msmf)]] <- sample(c("I", "R", "V"),
                                                              length(intersect(newHI,new.msmf)), replace = TRUE,
                                                              prob = dat$param$role.HI.msmf.prob)
  dat$attr$role.class[newIds[intersect(newW,new.msmf)]] <- sample(c("I", "R", "V"),
                                                             length(intersect(newW,new.msmf)), replace = TRUE,
                                                             prob = dat$param$role.W.msmf.prob)
  }

  ins.quot <- rep(NA,nBirths)
 
  ins.quot[dat$attr$role.class[newIds] == "I"]  <- 1
  ins.quot[dat$attr$role.class[newIds] == "R"]  <- 0
 
  ins.quot.v <- runif(sum(dat$attr$role.class[newIds] == "V",na.rm = TRUE))
  x <- which(dat$attr$role.class[newIds] == "V")
  ins.quot[x]  <-ins.quot.v
                                 
  dat$attr$ins.quot[newIds] <- ins.quot


  # CCR5
  ccr5.B.f.prob <- dat$param$ccr5.B.f.prob
  ccr5.BI.f.prob <- dat$param$ccr5.BI.f.prob
  ccr5.H.f.prob <- dat$param$ccr5.H.f.prob
  ccr5.HI.f.prob <- dat$param$ccr5.HI.f.prob
  ccr5.W.f.prob <- dat$param$ccr5.W.f.prob
  
  
  
  dat$attr$ccr5[newIds[intersect(newB,newF)]] <- sample(c("WW", "DW", "DD"),
                                        length(intersect(newB,newF)), replace = TRUE,
                                        prob = c(1 - sum(ccr5.B.f.prob),
                                                 ccr5.B.f.prob[2], ccr5.B.f.prob[1]))
  dat$attr$ccr5[newIds[intersect(newBI,newF)]] <- sample(c("WW", "DW", "DD"),
                                        length(intersect(newBI,newF)), replace = TRUE,
                                        prob = c(1 - sum(ccr5.BI.f.prob),
                                                 ccr5.BI.f.prob[2], ccr5.BI.f.prob[1]))
  dat$attr$ccr5[newIds[intersect(newH,newF)]] <- sample(c("WW", "DW", "DD"),
                                        length(intersect(newH,newF)), replace = TRUE,
                                        prob = c(1 - sum(ccr5.H.f.prob),
                                                 ccr5.H.f.prob[2], ccr5.H.f.prob[1]))
  dat$attr$ccr5[newIds[intersect(newHI,newF)]] <- sample(c("WW", "DW", "DD"),
                                        length(intersect(newHI,newF)), replace = TRUE,
                                        prob = c(1 - sum(ccr5.HI.f.prob),
                                                 ccr5.HI.f.prob[2], ccr5.HI.f.prob[1]))
  dat$attr$ccr5[newIds[intersect(newW,newF)]] <- sample(c("WW", "DW", "DD"),
                                        length(intersect(newW,newF)), replace = TRUE,
                                        prob = c(1 - sum(ccr5.W.f.prob),
                                                 ccr5.W.f.prob[2], ccr5.W.f.prob[1]))

  ccr5.B.m.prob <- dat$param$ccr5.B.m.prob
  ccr5.BI.m.prob <- dat$param$ccr5.BI.m.prob
  ccr5.H.m.prob <- dat$param$ccr5.H.m.prob
  ccr5.HI.m.prob <- dat$param$ccr5.HI.m.prob
  ccr5.W.m.prob <- dat$param$ccr5.W.m.prob
  
  
  
  dat$attr$ccr5[newIds[intersect(newB,newM)]] <- sample(c("WW", "DW", "DD"),
                                              length(intersect(newB,newM)), replace = TRUE,
                                              prob = c(1 - sum(ccr5.B.m.prob),
                                                       ccr5.B.m.prob[2], ccr5.B.m.prob[1]))
  dat$attr$ccr5[newIds[intersect(newBI,newM)]] <- sample(c("WW", "DW", "DD"),
                                               length(intersect(newBI,newM)), replace = TRUE,
                                               prob = c(1 - sum(ccr5.BI.m.prob),
                                                        ccr5.BI.m.prob[2], ccr5.BI.m.prob[1]))
  dat$attr$ccr5[newIds[intersect(newH,newM)]] <- sample(c("WW", "DW", "DD"),
                                              length(intersect(newH,newM)), replace = TRUE,
                                              prob = c(1 - sum(ccr5.H.m.prob),
                                                       ccr5.H.m.prob[2], ccr5.H.m.prob[1]))
  dat$attr$ccr5[newIds[intersect(newHI,newM)]] <- sample(c("WW", "DW", "DD"),
                                               length(intersect(newHI,newM)), replace = TRUE,
                                               prob = c(1 - sum(ccr5.HI.m.prob),
                                                        ccr5.HI.m.prob[2], ccr5.HI.m.prob[1]))
  dat$attr$ccr5[newIds[intersect(newW,newM)]] <- sample(c("WW", "DW", "DD"),
                                              length(intersect(newW,newM)), replace = TRUE,
                                              prob = c(1 - sum(ccr5.W.m.prob),
                                                       ccr5.W.m.prob[2], ccr5.W.m.prob[1]))

  # Degree
  dat$attr$deg.cohab[newIds] <- 0
  dat$attr$deg.pers[newIds] <- 0

  dat$attr$deg.cohab.c[newIds] <- 0
  dat$attr$deg.pers.c[newIds] <- 0

  # One-off risk group
  dat$attr$riskg[newIds] <- sample(1:5, nBirths, TRUE)

  # UAI group
  p1.msm <- dat$param$cond.pers.always.prob.msm
  p2.msm <- dat$param$cond.inst.always.prob.msm
  rho.msm <- dat$param$cond.always.prob.corr.msm
  uai.always.msm <- bindata::rmvbin(nBirths, c(p1.msm, p2.msm), bincorr = (1 - rho.msm) * diag(2) + rho.msm)
  dat$attr$cond.always.pers.msm[newIds] <- uai.always.msm[, 1]
  dat$attr$cond.always.inst.msm[newIds] <- uai.always.msm[, 2]
  
  # UVI group
  p1.het <- dat$param$cond.pers.always.prob.het
  p2.het <- dat$param$cond.inst.always.prob.het
  rho.het <- dat$param$cond.always.prob.corr.het
  uvi.always.het <- bindata::rmvbin(nBirths, c(p1.het, p2.het), bincorr = (1 - rho.het) * diag(2) + rho.het)
  dat$attr$cond.always.pers.het[newIds] <- uvi.always.het[, 1]
  dat$attr$cond.always.inst.het[newIds] <- uvi.always.het[, 2]

  # PrEP
  dat$attr$prepStat[newIds] <- 0
  
  

  
  
  return(dat)
}


