
#' @title Births Module for up to 5 race groups heterosexuals and MSM for preserving demographic distribution.
#'
#' @description Module function for births or entries into the sexually active
#'              population.
#'
#' @inheritParams aging_shamp
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
births_KTM <- function(dat, at){

  if(at > 1){
    
  ## Variables


  ## Process

  nBirths.gen <- max(dat$epi$dth.gen[at],0) 
  nBirths.dis <- max(dat$epi$dth.dis[at],0)
  nBirths.age <- max(dat$epi$dth.age[at],0)

  ##create population growth
  if (dat$param$p.growth == TRUE){
    
    if (at == 2){
      
      nBirths.gen <- nBirths.gen + dat$param$p.growth.size
    }
  }  


  nBirths <- nBirths.age + nBirths.gen +  nBirths.dis
  

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


#setBirthAttr_shamp <- function(dat, at, nBirths.gen, nBirths.age, nBirths.dis) {

setBirthAttr_shamp <- function(dat, at, nBirths.gen, nBirths.age, nBirths.dis) {
    
  # Parameters
  percent.male <- dat$param$percent.male


  
  ##Not replcing disease death at this time
  nBirths <- nBirths.gen + nBirths.age + nBirths.dis
  

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
  #dat$attr$arv.BI.pos[newIds] <- rep(0, nBirths)
  #dat$attr$arv.HI.pos[newIds] <- rep(0, nBirths)
  dat$attr$sex[newIds] <- rep(0, nBirths)
  dat$attr$race[newIds] <- rep(0, nBirths)
  dat$attr$age[newIds] <- rep(0, nBirths)
  dat$attr$sqrt.age[newIds] <- rep(0, nBirths)
  dat$attr$agesq[newIds] <- rep(0, nBirths)
  dat$attr$age.group[newIds] <- rep(0, nBirths)
  dat$attr$evertest[newIds] <- rep(0, nBirths)
  dat$attr$infected.gen[newIds]<-rep(NA, nBirths)
  
  dat$attr$age.inf[newIds]<-rep(NA, nBirths)
  dat$attr$age.diag[newIds]<-rep(NA, nBirths)
  dat$attr$age.group[newIds]<-rep(NA, nBirths)
  dat$attr$sqrt.age.adj[newIds]<-rep(NA, nBirths)
  
  dat$attr$age.inf[newIds]<-rep(NA, nBirths)
  dat$attr$age.diag[newIds]<-rep(NA, nBirths)
  

  dat$attr$arrival.time[newIds] <- rep(at, nBirths)


  sex<-rbinom(length(nBirths),1,dat$param$percent.male)
  sex<-ifelse(sex==1,"M","F")
  
  race<-rep("B",nBirths)
  age<-rep(dat$param$birth.age,nBirths)
  age.group<-rep(1,nBirths)

#asighn attributes to dat.
  dat$attr$sex[newIds] <- sex
  dat$attr$race[newIds] <- race
  dat$attr$age[newIds] <- age
  dat$attr$sqrt.age <- sqrt(dat$attr$age) 
  dat$attr$agesq <- dat$attr$age^2 


  dat$attr$sqrt.age.adj<-ifelse(dat$attr$sex=="M",dat$attr$sqrt.age,
                                ifelse(dat$attr$sex=="F",sqrt(dat$attr$age-5),dat$attr$sqrt.age))
  
  
  
  #AGEGROUP
  dat$attr$age.group<-ifelse(dat$att$age < 20 ,1,
                             ifelse(dat$attr$age >= 20 & dat$attr$age < 25,2,
                                    ifelse(dat$attr$age >= 25 & dat$attr$age < 30,3,
                                           ifelse(dat$attr$age >= 30 & dat$attr$age < 35,4,
                                                  ifelse(dat$attr$age >= 35,5,dat$attr$age.group)))))
  
  #SEX.AGE.GROUP.
  dat$attr$sex.age.group<-ifelse(dat$attr$sex == "M" & dat$attr$age < 20 ,"M.1",
                        ifelse(dat$attr$sex == "M" & age >= 20 & dat$attr$age < 25,"M.2",
                               ifelse(dat$attr$sex == "M" & dat$attr$age >= 25 & dat$attr$age < 30,"M.3",
                                      ifelse(dat$attr$sex == "M" & dat$attr$age >= 30 & dat$attr$age < 35,"M.4",
                                             ifelse(dat$attr$sex == "M" & dat$attr$age >= 35,"M.5",
                                                    ifelse(dat$attr$sex == "F" & dat$attr$age < 20 ,"M.1",
                                                           ifelse(dat$attr$sex == "F" & dat$attr$age >= 20 & dat$attr$age < 25,"M.2",
                                                                  ifelse(dat$attr$sex == "F" & dat$attr$age >= 25 & dat$attr$age < 30,"M.3",
                                                                         ifelse(dat$attr$sex == "F" & dat$attr$age >= 30 & dat$attr$age < 35,"M.4",
                                                                                ifelse(dat$attr$sex == "F" & dat$attr$age >= 35,"M.5",
                                                                                       dat$attr$sex.age.group))))))))))
  
                      
  
  newF<-which(sex=="F")
  newM<-which(sex=="M")

  
  newB<-which(race=="B")
  newBI<-which(race=="BI")
  newH<-which(race=="H")
  newHI<-which(race=="HI")
  newW<-which(race=="W")
  
  newB.m<-which(race=="B" & sex=="M")
  newBI.m<-which(race=="BI" & sex=="M")
  newH.m<-which(race=="H" & sex=="M")
  newHI.m<-which(race=="HI" & sex=="M")
  newW.m<-which(race=="W" & sex=="M")
  
  sex.ident<-rep(NA,nBirths)
  sex.ident[newF]<-"f"
  
  sex.ident[newB.m]<-sample(c("msm","msmf","msf"),length(newB.m),
                                    prob=c(dat$param$msm.frac,dat$param$msmf.frac.B,(1-(dat$param$msm.frac+dat$param$msmf.frac.B))),
                                           replace=TRUE)
  
  sex.ident[newBI.m]<-sample(c("msm","msmf","msf"),length(newBI.m),
                          prob=c(dat$param$msm.frac,dat$param$msmf.frac.BI,(1-(dat$param$msm.frac+dat$param$msmf.frac.BI))),
                          replace=TRUE)
  
  sex.ident[newH.m]<-sample(c("msm","msmf","msf"),length(newH.m),
                          prob=c(dat$param$msm.frac,dat$param$msmf.frac.H,(1-(dat$param$msm.frac+dat$param$msmf.frac.H))),
                          replace=TRUE)
  
  sex.ident[newHI.m]<-sample(c("msm","msmf","msf"),length(newHI.m),
                          prob=c(dat$param$msm.frac,dat$param$msmf.frac.HI,(1-(dat$param$msm.frac+dat$param$msmf.frac.HI))),
                          replace=TRUE)
  
  sex.ident[newW.m]<-sample(c("msm","msmf","msf"),length(newW.m),
                          prob=c(dat$param$msm.frac,dat$param$msmf.frac.W,(1-(dat$param$msm.frac+dat$param$msmf.frac.W))),
                          replace=TRUE)
  
  dat$attr$sex.ident[newIds] <- sex.ident
  dat$attr$msmf<-ifelse(dat$attr$sex.ident == "msmf",1,0)
  new.msm<-which(sex.ident=="msm")
  new.msmf<-which(sex.ident=="msmf")
  new.msf<-which(sex.ident=="msf")
 

 
  #Status
  dat$attr$status[newIds] <- rep(0, nBirths)
 
  #Set status for Black and Hispanic imigrant that are entering at an age greater than 18 based on prevelance parameters pos.entry.BI pos.entry.HI

  #BI.imigrating <- which(dat$attr$race[newIds] == "BI" & dat$attr$age[newIds] == dat$param$birth.age)
  #BI.status <- rbinom(length(BI.imigrating),1,prob=dat$param$pos.entry.BI)
  #dat$attr$status[newIds][BI.imigrating] <- BI.status
  #dat$attr$arv.BI.pos[newIds][BI.imigrating] <- BI.status * at
    
  #HI.imigrating <- which(dat$attr$race[newIds] == "HI" & dat$attr$age[newIds] == dat$param$birth.age)
  

 
  # dat$attr$inst.ai.class[newIds] <- sample(1:dat$param$num.inst.ai.classes,
  #                                         nBirths, replace = TRUE)
 
   #### females f. 
  
  if (length(newF) > 0) {
  dat$attr$tt.traj[newIds[intersect(newB,newF)]] <- sample(c(1, 2, 3, 4),
                                           length(intersect(newB,newF)), replace = TRUE,
                                           prob = dat$param$tt.traj.B.f.prob)
  }
  
  #### Males msf. 
  if (length(new.msf) > 0) {
  dat$attr$tt.traj[newIds[intersect(newB,new.msf)]] <- sample(c(1, 2, 3, 4),
                                                  length(intersect(newB,new.msf)), replace = TRUE,
                                                 prob = dat$param$tt.traj.B.msf.prob)
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
  
  dat$attr$cohab.lt[newIds] <- 0
  dat$attr$pers.lt[newIds] <- 0
  dat$attr$onetime.lt[newIds] <- 0
  



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
  
  
  #Kenya TM intervention
  dat$attr$partner.serv[newIds] <- 0
  dat$attr$partner.serv.time[newIds] <- 0
 
  return(dat)
}


