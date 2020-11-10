
#' @title HIV Testing Module for up to 5 race groups heterosexuals and MSM.
#'
#' @description Module function for HIV diagnostic testing of infected persons.
#'
#' @inheritParams aging_shamp
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
test_KTM <- function(dat, at) {

  ## Variables

  # Attributes
  diag.status <- dat$attr$diag.status
  inf.time <- dat$attr$inf.time
  sex <- dat$attr$sex
  age.group <- dat$attr$age.group
  status <- dat$attr$status
  inf.time <- dat$attr$inf.time
  stage <- dat$attr$stage
  evertest <- dat$attr$evertest
  age <- dat$attr$age
  lnt <- dat$attr$last.neg.test
  partner.serv.part <- dat$attr$partner.serv.part
  partner.serv.part.time <- dat$attr$partner.serv.part.time

  prepStat <- dat$attr$prepStat
  prep.elig.time <- dat$attr$prep.elig.time 
  prep.tst.int <- dat$param$prep.tst.int


  # Parameters
  test.prob.nevtest.m <- dat$param$test.prob.nevtest.m
  test.prob.nevtest.f <- dat$param$test.prob.nevtest.f
  test.prob.tested.m <- dat$param$test.prob.tested.m
  test.prob.tested.f <- dat$param$test.prob.tested.f
  
  twind.int.ab <- dat$param$test.window.int.ab
  twind.int.rna <- dat$param$test.window.int.rna
  tsincelntst <- at - dat$attr$last.neg.test
  
   
  seek.hc.AHI.prob <- dat$param$ seek.hc.AHI.prob
  sym.seek.prob <- dat$param$sym.seek.prob
  sym.test.prob.bl <- dat$param$sym.test.prob.b
  sym.test.prob.tm <- dat$param$sym.test.prob.tm
  
  partner.test.prob.bl <- dat$param$partner.test.prob.bl
  partner.test.prob.tm <- dat$param$partner.test.prob.tm
  PS.time <- dat$param$PS.time
  

  intervention_TM <- dat$param$intervention_TM
  
  ##Place holders
  tst.neg.ab <- NULL
  tst.neg.rna <- NULL
  tst.pos.ab <- NULL
  tst.pos.rna <- NULL
  tst.ab <- NULL
  tst.rna <- NULL
  prep.elig <- NULL
  clinic.test <- NULL
  tst.pos.clinic.test <- NULL
  
## Process

#females
  #First time testers
  #Age group 1
    elig.f.age1.nt <- which(sex == "F"  & age.group==1 & evertest==0 & partner.serv.part==0 &
                    (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age1.nt <- rep(test.prob.nevtest.f[1], length(elig.f.age1.nt))
    
    tst.f.age1.nt <- elig.f.age1.nt[rbinom(length(elig.f.age1.nt), 1, rates.f.age1.nt) == 1]
    
    #Age group 2
    elig.f.age2.nt <- which(sex == "F"  & age.group==2 & evertest==0 & partner.serv.part==0 &
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age2.nt <- rep(test.prob.nevtest.f[2], length(elig.f.age2.nt))
    
    tst.f.age2.nt <- elig.f.age2.nt[rbinom(length(elig.f.age2.nt), 1, rates.f.age2.nt) == 1]
    
    #Age group 3
    elig.f.age3.nt <- which(sex == "F"  & age.group==3 & evertest==0 & partner.serv.part==0 &
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age3.nt <- rep(test.prob.nevtest.f[3], length(elig.f.age3.nt))
    
    tst.f.age3.nt <- elig.f.age3.nt[rbinom(length(elig.f.age3.nt), 1, rates.f.age3.nt) == 1]
    
    #Repeat testers
    #Age group 1
    elig.f.age1.t <- which(sex == "F"  & age.group==1 & evertest==1 & partner.serv.part==0 &
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age1.t <- rep(test.prob.tested.f[1], length(elig.f.age1.t))
    
    tst.f.age1.t <- elig.f.age1.t[rbinom(length(elig.f.age1.t), 1, rates.f.age1.t) == 1]
    
    #Age group 2
    elig.f.age2.t <- which(sex == "F"  & age.group==2 & evertest==1 & partner.serv.part==0 &
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age2.t <- rep(test.prob.tested.f[2], length(elig.f.age2.t))
    
    tst.f.age2.t <- elig.f.age2.t[rbinom(length(elig.f.age2.t), 1, rates.f.age2.t) == 1]
    
    #Age group 3
    elig.f.age3.t <- which(sex == "F"  & age.group==3 & evertest==1 & partner.serv.part==0 &
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age3.t <- rep(test.prob.tested.f[3], length(elig.f.age3.t))
    
    tst.f.age3.t <- elig.f.age3.t[rbinom(length(elig.f.age3.t), 1, rates.f.age3.t) == 1]
    
    #Age group 4+
    elig.f.age4.t <- which(sex == "F"  & (age.group == 4 | age.group == 5) & evertest==1 & partner.serv.part==0 &
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age4.t <- rep(test.prob.tested.f[4], length(elig.f.age4.t))
    
    tst.f.age4.t <- elig.f.age4.t[rbinom(length(elig.f.age4.t), 1, rates.f.age4.t) == 1]
    

#males

    elig.m.age1.nt <- which(sex == "M"  & age.group==1 & evertest==0 & partner.serv.part==0 & 
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age1.nt <- rep(test.prob.nevtest.m[1], length(elig.m.age1.nt))
    
    tst.m.age1.nt <- elig.m.age1.nt[rbinom(length(elig.m.age1.nt), 1, rates.m.age1.nt) == 1]
    
    #Age group 2
    elig.m.age2.nt <- which(sex == "M"  & age.group==2 & evertest==0 & partner.serv.part==0 &
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age2.nt <- rep(test.prob.nevtest.m[2], length(elig.m.age2.nt))
    
    tst.m.age2.nt <- elig.m.age2.nt[rbinom(length(elig.m.age2.nt), 1, rates.m.age2.nt) == 1]
    
    #Age group 3
    elig.m.age3.nt <- which(sex == "M"  & age.group==3 & evertest==0 & partner.serv.part==0 &
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age3.nt <- rep(test.prob.nevtest.m[3], length(elig.m.age3.nt))
    
    tst.m.age3.nt <- elig.m.age3.nt[rbinom(length(elig.m.age3.nt), 1, rates.m.age3.nt) == 1]
    
    #Repeat testers
    #Age group 1
    elig.m.age1.t <- which(sex == "M"  & age.group==1 & evertest==1 & partner.serv.part==0 &
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age1.t <- rep(test.prob.tested.m[1], length(elig.m.age1.t))
    
    tst.m.age1.t <- elig.m.age1.t[rbinom(length(elig.m.age1.t), 1, rates.m.age1.t) == 1]
    
    #Age group 2
    elig.m.age2.t <- which(sex == "M"  & age.group==2 & evertest==1 & partner.serv.part==0 &
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age2.t <- rep(test.prob.tested.m[2], length(elig.m.age2.t))
    
    tst.m.age2.t <- elig.m.age2.t[rbinom(length(elig.m.age2.t), 1, rates.m.age2.t) == 1]
    
    #Age group 3
    elig.m.age3.t <- which(sex == "M"  & age.group==3 & evertest==1 & partner.serv.part==0 &
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age3.t <- rep(test.prob.tested.m[3], length(elig.m.age3.t))
    
    tst.m.age3.t <- elig.m.age3.t[rbinom(length(elig.m.age3.t), 1, rates.m.age3.t) == 1]
    
    #Age group 4
    elig.m.age4.t <- which(sex == "M"  & (age.group == 4 | age.group == 5) & evertest==1 & partner.serv.part==0 &
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age4.t <- rep(test.prob.tested.m[4], length(elig.m.age4.t))
    
    tst.m.age4.t <- elig.m.age4.t[rbinom(length(elig.m.age4.t), 1, rates.m.age4.t) == 1]
    
    

    tst.background <- c(tst.f.age1.nt,tst.f.age2.nt,tst.f.age3.nt,
                   tst.f.age1.t,tst.f.age2.t,tst.f.age3.t,tst.f.age4.t,
                   tst.m.age1.nt,tst.m.age2.nt,tst.m.age3.nt,
                   tst.m.age1.t,tst.m.age2.t,tst.m.age3.t,tst.m.age4.t)


    if(intervention_TM=="NONE"){


    ## Symptom based testing

    #Those with HIV
    sym_HIV <- which(status == 1 & inf.time >= (at - 3) & diag.status ==0  & partner.serv.part==0 & age < 40)
    selected <- rbinom(length(sym_HIV),1,seek.hc.AHI.prob*sym.test.prob.bl)
    sym_HIV <- sym_HIV[selected == 1]
    

    #Those with other illness
    sym<-which((diag.status == 0 | is.na(diag.status) == TRUE)  & age < 40)
    selected <- rbinom(length(sym),1,sym.seek.prob*sym.test.prob.bl)
    sym <- sym[selected == 1]
    
    #Count total number of background tests.  AHI symptom and other symptom tests will need to be removed from the background 
    #to be counted separately
    ntest<-length(tst.background)
    
    #remove individuals selected for testing as both background and AHI symptoms
    n.keep <- ntest-length(sym_HIV)-length(sym)
    
    tst.background <- tst.background[!(tst.background %in% sym_HIV)]
    tst.background <- tst.background[!(tst.background %in% sym)]

    if(length(tst.background) >= n.keep){
    tst.background <- sample(tst.background,n.keep,replace=FALSE)}
    
    if(length(tst.background) < n.keep){
      dat$epi$undertest[at] <- n.keep - length(tst.background)}
      else
        {dat$epi$undertest[at]<-0}
    
    dat$epi$n.tests.bg[at] <- max(0,length(tst.background))
   
     #PARTNER SERVICES 

    #Partner services testing
    tst.ps <- which(partner.serv.part==1)
    selected <- rbinom(length(tst.ps),1,partner.test.prob.bl)
    tst.ps <- tst.ps[selected == 1]
    
    #Select off partners that are already diagnosed
    exclude <- which(diag.status==1)
    PS.prior.diag <- intersect(tst.ps,exclude)
    tst.ps <- setdiff(tst.ps,exclude)
    
    ##TEST tst.background, tst.PREP, partner services -  Antibody tests
    ##sym_HIV, tst.sym.m, tst.sym.f.   -  RNA tests
    
    tst.ab <- c(tst.background, sym_HIV, sym, tst.ps)
    in.clinic <- c(sym_HIV, sym, tst.ps)
    
    tst.pos.ab <- tst.ab[status[tst.ab] == 1 & inf.time[tst.ab] <= at - twind.int.ab ]
    tst.neg.ab <- setdiff(tst.ab, tst.pos.ab)
  
    
    ##PUT IN INDICATOR FOR NEW INDEX FOR PARTNER SERVICES
    ##ZERO OUT THE INDICATOR AFTER THE OS MODULE IS RUN
    PS.index<-intersect(in.clinic,tst.pos.ab)
    
    ##Select all symptom based for background partner services
    dat$attr$PS.index.prev[PS.index]<-1
    
    
    ##Remove partner that have been tested from list of those being tracked for testing and stop the clock
    dat$attr$partner.serv.part[tst.ps] <- 0
    dat$attr$partner.serv.part.time[tst.ps] <- 0
    
    ##Remove followed-up and already diagnosed from continued follow-up that have been tested from list of those being tracked for testing and stop the clock
    dat$attr$partner.serv.part[PS.prior.diag] <- 0
    dat$attr$partner.serv.part.time[PS.prior.diag] <- 0
    
    ##Advance the clock on partners still being followed
    dat$attr$partner.serv.part.time <- ifelse(dat$attr$partner.serv.part.time >=1, dat$attr$partner.serv.part.time + 1, 
                                              ifelse(dat$attr$partner.serv.part.time > PS.time, 0, dat$attr$partner.serv.part.time))
    
    ##Remove those that have been lost to follow-up from the partner services list
    
    
    dat$attr$partner.serv.part <- ifelse(dat$attr$partner.serv.part.time == 0, 0, dat$attr$partner.serv.part)
    
    ##TEST tst.background, sym_HIV, sym   -  all Antibody tests
    
    tst.ab <- c(tst.background, sym_HIV, sym)
    
    tst.pos.ab <- tst.ab[status[tst.ab] == 1 & inf.time[tst.ab] <= at - twind.int.ab ]
    tst.neg.ab <- setdiff(tst.ab, tst.pos.ab)
    
    }
    

    if(intervention_TM=="TMP"){
      #THE SYMPT TESTS THAT WOULD BE DONE ANTWAY ARE BACKGROUND BUT WE NEED TO ADD NEW TESTS FOR BETTER PITC 
      ##WHO HAS HIV SYMPTOMS AND SELECT THEM FOR TESTING
      
      #Those with HIV
      sym_HIV <- which(status == 1 & inf.time >= (at - 3) & diag.status == 0  & partner.serv.part==0  & age < 40)
      selected <- rbinom(length(sym_HIV),1,seek.hc.AHI.prob*sym.test.prob.tm)
      sym_HIV.bg.count<-sum(rbinom(length(sym_HIV),1,seek.hc.AHI.prob*sym.test.prob.bl))
      sym_HIV <- sym_HIV[selected == 1]
      
      
      #Those with other illness
      sym<-which((diag.status == 0 | is.na(diag.status) == TRUE)  & age < 40)
      selected <- rbinom(length(sym),1,sym.seek.prob*sym.test.prob.tm)
      sym.bg.count<- sum(rbinom(length(sym),1,sym.seek.prob*sym.test.prob.bl))
      sym <- sym[selected == 1]
      
      #remove individuals selected for testing as both background and none AHI symptoms from list for background
      ntest<-length(tst.background)
      n.keep <- ntest - sym.bg.count - sym_HIV.bg.count
      
      tst.background <- tst.background[!(tst.background %in% sym_HIV)]
      tst.background <- tst.background[!(tst.background %in% sym)]
      

      if(length(tst.background) >= n.keep){
        tst.background <- sample(tst.background,n.keep,replace=FALSE)}
      
      if(length(tst.background) < n.keep){
        dat$epi$undertest[at] <- n.keep - length(tst.background)}
        else
          {dat$epi$undertest[at]<-0}
      
      dat$epi$n.tests.bg[at] <- max(0,length(tst.background))

      
      #PARTNER SERVICES
 
      
      #Partner services testing
      tst.ps <- which(partner.serv.part==1)
      selected <- rbinom(length(tst.ps),1,partner.test.prob.tm)
      tst.ps <- tst.ps[selected == 1]
      
      #Select off partners that are already diagnosed
      exclude <- which(diag.status==1)
      PS.prior.diag <- intersect(tst.ps,exclude)
      tst.ps <- setdiff(tst.ps,exclude)
      
      #PREP
      
      tst.prep <- which(prepStat == 1 & lnt == at- prep.tst.int)
      
      
      ##TEST tst.background, tst.PREP, partner services -  Antibody tests
      ##sym_HIV, tst.sym.m, tst.sym.f.   -  RNA tests
      
      tst.ab <- c(tst.background, tst.prep)
      tst.rna <- c(sym_HIV, sym, tst.ps)
      
      tst.pos.ab <- tst.ab[status[tst.ab] == 1 & inf.time[tst.ab] <= at - twind.int.ab ]
      tst.pos.rna <- tst.rna[status[tst.rna] == 1 & inf.time[tst.rna] <= at - twind.int.rna ]
      tst.neg.ab <- setdiff(tst.ab, tst.pos.ab)
      tst.neg.rna <- setdiff(tst.rna, tst.pos.rna)
      
      ##If they are testing because of partner services and are negative they are eligible for PrEP
      prep.elig <- intersect(tst.ps,tst.neg.rna)
      prep.elig.time[prep.elig] <- at
      
      ##PUT IN INDICATOR FOR NEW INDEX FOR PARTNER SERVICES
      ##ZERO OUT THE INDICATOR AFTER THE OS MODULE IS RUN
      acute <- which(inf.time[tst.pos.rna] >= at - twind.int.ab )
      prev <- which(inf.time[tst.pos.rna] < at - twind.int.ab )

      
      dat$attr$PS.index.acute[tst.pos.rna[acute]]<-1
      dat$attr$PS.index.prev[tst.pos.rna[prev]]<-1
      
      
      ##Remove partner that have been tested from list of those being tracked for testing and stop the clock
      dat$attr$partner.serv.part[tst.ps] <- 0
      dat$attr$partner.serv.part.time[tst.ps] <- 0
      
      ##Remove followed-up and already diagnosed from continued follow-up that have been tested from list of those being tracked for testing and stop the clock
      dat$attr$partner.serv.part[PS.prior.diag] <- 0
      dat$attr$partner.serv.part.time[PS.prior.diag] <- 0
      
      ##Advance the clock on partners still being followed
      dat$attr$partner.serv.part.time <- ifelse(dat$attr$partner.serv.part.time >=1, dat$attr$partner.serv.part.time + 1, 
                                                ifelse(dat$attr$partner.serv.part.time > PS.time, 0, dat$attr$partner.serv.part.time))
    
      ##Remove those that have been lost to follow-up from the partner services list

      
      dat$attr$partner.serv.part <- ifelse(dat$attr$partner.serv.part.time == 0, 0, dat$attr$partner.serv.part)
                                                
    }
  
    if(intervention_TM=="TESTS"){
      #THE SYMPT TESTS THAT WOULD BE DONE ANTWAY ARE BACKGROUND BUT WE NEED TO ADD NEW TESTS FOR BETTER PITC 
      ##WHO HAS HIV SYMPTOMS AND SELECT THEM FOR TESTING
      
      #Those with HIV
      sym_HIV <- which(status == 1 & inf.time >= (at - 3) & diag.status == 0  & partner.serv.part==0  & age < 40)
      selected <- rbinom(length(sym_HIV),1,seek.hc.AHI.prob*sym.test.prob.tm)
      sym_HIV.bg.count<-sum(rbinom(length(sym_HIV),1,seek.hc.AHI.prob*sym.test.prob.bl))
      sym_HIV <- sym_HIV[selected == 1]
      
      
      #Those with other illness
      sym<-which((diag.status == 0 | is.na(diag.status) == TRUE)  & age < 40)
      selected <- rbinom(length(sym),1,sym.seek.prob*sym.test.prob.tm)
      sym.bg.count<- sum(rbinom(length(sym),1,sym.seek.prob*sym.test.prob.bl))
      sym <- sym[selected == 1]
      
      #remove individuals selected for testing as both background and none AHI symptoms from list for background
      ntest<-length(tst.background)
      n.keep <- ntest - sym.bg.count - sym_HIV.bg.count
      
      tst.background <- tst.background[!(tst.background %in% sym_HIV)]
      tst.background <- tst.background[!(tst.background %in% sym)]
      
      
      if(length(tst.background) >= n.keep){
        tst.background <- sample(tst.background,n.keep,replace=FALSE)}
      
      if(length(tst.background) < n.keep){
        dat$epi$undertest[at] <- n.keep - length(tst.background)}
      else
      {dat$epi$undertest[at]<-0}
      
      dat$epi$n.tests.bg[at] <- max(0,length(tst.background))
      
      
      #PARTNER SERVICES
      
      
      #Partner services testing
      tst.ps <- which(partner.serv.part==1)
      selected <- rbinom(length(tst.ps),1,partner.test.prob.tm)
      tst.ps <- tst.ps[selected == 1]
      
      #Select off partners that are already diagnosed
      exclude <- which(diag.status==1)
      PS.prior.diag <- intersect(tst.ps,exclude)
      tst.ps <- setdiff(tst.ps,exclude)
      
      #PREP
      
      tst.prep <- which(prepStat == 1 & lnt == at- prep.tst.int)
      
      
      ##TEST tst.background, tst.PREP, partner services -  Antibody tests
      ##sym_HIV, tst.sym.m, tst.sym.f.   -  RNA tests
      
      tst.ab <- c(tst.background, tst.prep,sym_HIV, sym, tst.ps)
      clinic.test <-c(sym_HIV, sym, tst.ps)
      tst.rna <- NULL
      
      tst.pos.ab <- tst.ab[status[tst.ab] == 1 & inf.time[tst.ab] <= at - twind.int.ab ]
      tst.pos.clinic.test <- clinic.test[status[clinic.test] == 1 & inf.time[clinic.test] <= at - twind.int.ab ]
      tst.neg.ab <- setdiff(tst.ab, tst.pos.ab)
  
      ##If they are testing because of partner services and are negative they are eligible for PrEP
      prep.elig <- intersect(tst.ps,tst.neg.ab)
      prep.elig.time[prep.elig] <- at
      
      ##PUT IN INDICATOR FOR NEW INDEX FOR PARTNER SERVICES
      ##ZERO OUT THE INDICATOR AFTER THE OS MODULE IS RUN
      acute <- which(inf.time[tst.pos.clinic.test] >= at - twind.int.ab )
      prev <- which(inf.time[tst.pos.clinic.test] < at - twind.int.ab )
      
      
      dat$attr$PS.index.acute[tst.pos.clinic.test[acute]]<-1
      dat$attr$PS.index.prev[tst.pos.clinic.test[prev]]<-1
      
      
      ##Remove partner that have been tested from list of those being tracked for testing and stop the clock
      dat$attr$partner.serv.part[tst.ps] <- 0
      dat$attr$partner.serv.part.time[tst.ps] <- 0
      
      ##Remove followed-up and already diagnosed from continued follow-up that have been tested from list of those being tracked for testing and stop the clock
      dat$attr$partner.serv.part[PS.prior.diag] <- 0
      dat$attr$partner.serv.part.time[PS.prior.diag] <- 0
      
      ##Advance the clock on partners still being followed
      dat$attr$partner.serv.part.time <- ifelse(dat$attr$partner.serv.part.time >=1, dat$attr$partner.serv.part.time + 1, 
                                                ifelse(dat$attr$partner.serv.part.time > PS.time, 0, dat$attr$partner.serv.part.time))
      
      ##Remove those that have been lost to follow-up from the partner services list
      
      
      dat$attr$partner.serv.part <- ifelse(dat$attr$partner.serv.part.time == 0, 0, dat$attr$partner.serv.part)
      
    }
    

  # Attributes
  dat$attr$last.neg.test[tst.neg.ab] <- at
  dat$attr$last.neg.test[tst.neg.rna] <- at
  dat$attr$diag.status[tst.pos.ab] <- 1 
  dat$attr$diag.status[tst.pos.rna]  <- 1
  dat$attr$diag.time[tst.pos.ab] <- at
  dat$attr$diag.time[tst.pos.rna] <- at
  dat$attr$evertest[tst.ab] <- 1
  dat$attr$evertest[tst.rna] <- 1
  
  ##PS
  PS.pos<-which(dat$attr$diag.status[tst.ps]==1)
  PS.pos <- tst.ps[PS.pos]
  PS.neg<-setdiff(tst.ps,PS.pos)
  
  
  dat$attr$PS.diag.pos.time[PS.pos] <- at
  dat$attr$PS.diag.neg[PS.neg] <- 1
  dat$attr$diag.status[PS.neg] <- 0
  
  ##
  dat$attr$age.diag[tst.pos.ab] <- dat$attr$age[tst.pos.ab]
  dat$attr$age.diag[tst.pos.rna] <- dat$attr$age[tst.pos.rna]
  
  # Prevalence

  
  dat$epi$partners.found[at] <- max(0,length(tst.ps)+length(PS.prior.diag))
  dat$epi$partners.positive[at] <- max(0,sum(dat$attr$diag.status[tst.ps])) 
  dat$epi$partners.negative[at] <- max(0,length(tst.ps)-sum(dat$attr$diag.status[tst.ps]))
  dat$epi$PS.prior.diag[at] <- max(0, length(PS.prior.diag))
 
  dat$epi$n.tests[at] <- sum(length(tst.neg.ab),length(tst.neg.rna),length(tst.pos.ab),length(tst.pos.rna))
  dat$epi$n.tests.ab[at] <- sum(length(tst.neg.ab),length(tst.pos.ab))
  dat$epi$n.tests.rna[at] <- sum(length(tst.neg.rna),length(tst.pos.rna))
  dat$epi$n.tests.tst.ps[at] <- max(0,length(tst.ps))
  
  acute <- which(inf.time[tst.pos.rna] >= at - twind.int.ab )
  acute <- max(length(acute),0)
  dat$epi$diag.prevalent[at] <- sum(length(tst.pos.ab),length(tst.pos.rna)) - acute
  dat$epi$diag.acute[at] <- acute
  
  dat$epi$incid.diag[at] <- sum(length(tst.pos.ab),length(tst.pos.rna))
  
  age1<- dat$attr$age[tst.pos.ab]
  age2<- dat$attr$age[tst.pos.rna]
  dat$epi$mean.age.diag[at] <- max(0,mean(c(age1,age2)))
  
  dat$attr$prepElig[prep.elig]<-1
  
  
  return(dat)
}
