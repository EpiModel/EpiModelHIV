
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
  sex <- dat$attr$sex
  age.group <- dat$attr$age.group
  status <- dat$attr$status
  inf.time <- dat$attr$inf.time
  stage <- dat$attr$stage
  evertest <- dat$attr$evertest

  partner.serv <- dat$attr$partner.serv
  partner.serv.time <- dat$attr$partner.serv.time
  cel.temp <- dat$cel.temp
  cel.complete <- dat$cel.temp
  
  prepStat <- dat$attr$prepStat
  prep.tst.int <- dat$param$prep.tst.int

  # Parameters
  test.prob.nevtest.m <- dat$param$test.prob.nevtest.m
  test.prob.nevtest.f <- dat$param$test.prob.nevtest.f
  test.prob.tested.m <- dat$param$test.prob.tested.m
  test.prob.tested.f <- dat$param$test.prob.tested.f
  
  twind.int.ab <- dat$param$test.window.int.ab
  twind.int.rna <- dat$param$test.window.int.rna

  seek.hc.AHI.prob <- dat$param$seek.hc.AHI.prob
  sym.prob <- dat$param$sym.prob
  seek.hc.prob.m <- dat$param$seek.hc.prob.m
  seek.hc.prob.f <- dat$param$seek.hc.prob.f
  PITC.prob <- dat$param$PITC.prob
  PITC.TMP.prob <- dat$param$PITC.TMP.prob
  
  partner.test.prob <- dat$param$partner.test.prob
  
  

## Process

#females
  #First time testers
  #Age group 1
    elig.f.age1.nt <- which(sex == "F"  & age.group==1 & evertest==0 & partner.serv==00
                    (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age1.nt <- rep(test.prob.nevtest.f[1], length(elig.f.age1.nt))
    
    tst.f.age1.nt <- elig.f.age1.nt[rbinom(length(elig.f.age1.nt), 1, rates.f.age1.nt) == 1]
    
    #Age group 2
    elig.f.age2.nt <- which(sex == "F"  & age.group==2 & evertest==0 & partner.serv==0
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age2.nt <- rep(test.prob.nevtest.f[2], length(elig.f.age2.nt))
    
    tst.f.age2.nt <- elig.f.age2.nt[rbinom(length(elig.f.age2.nt), 1, rates.f.age2.nt) == 1]
    
    #Age group 3
    elig.f.age3.nt <- which(sex == "F"  & age.group==3 & evertest==0 & partner.serv==0
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age3.nt <- rep(test.prob.nevtest.f[3], length(elig.f.age3.nt))
    
    tst.f.age3.nt <- elig.f.age3.nt[rbinom(length(elig.f.age3.nt), 1, rates.f.age3.nt) == 1]
    
    #Repeat testers
    #Age group 1
    elig.f.age1.t <- which(sex == "F"  & age.group==1 & evertest==1 & partner.serv==0
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age1.t <- rep(test.prob.tested.f[1], length(elig.f.age1.t))
    
    tst.f.age1.t <- elig.f.age1.n[rbinom(length(elig.f.age1.n), 1, rates.f.age1.n) == 1]
    
    #Age group 2
    elig.f.age2.t <- which(sex == "F"  & age.group==2 & evertest==1 & partner.serv==0
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age2.t <- rep(test.prob.tested.f[2], length(elig.f.age2.t))
    
    tst.f.age2.t <- elig.f.age2.t[rbinom(length(elig.f.age2.t), 1, rates.f.age2.t) == 1]
    
    #Age group 3
    elig.f.age3.t <- which(sex == "F"  & age.group==3 & evertest==1 & partner.serv==0
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age3.t <- rep(test.prob.tested.f[3], length(elig.f.age3.t))
    
    tst.f.age3.t <- elig.f.age3.t[rbinom(length(elig.f.age3.t), 1, rates.f.age3.t) == 1]
    
    #Age group 4+
    elig.f.age4.t <- which(sex == "F"  & age.group>=4 & evertest==1 & partner.serv==0
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.f.age4.t <- rep(test.prob.tested.f[4], length(elig.f.age4.t))
    
    tst.f.age4.t <- elig.f.age4.t[rbinom(length(elig.f.age4.t), 1, rates.f.age4.t) == 1]
    

#males

    elig.m.age1.nt <- which(sex == "F"  & age.group==1 & evertest==0 & partner.serv==0
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age1.nt <- rep(test.prob.nevtest.f[1], length(elig.m.age1.nt))
    
    tst.m.age1.nt <- elig.m.age1.nt[rbinom(length(elig.m.age1.nt), 1, rates.m.age1.nt) == 1]
    
    #Age group 2
    elig.m.age2.nt <- which(sex == "F"  & age.group==2 & evertest==0 & partner.serv==0
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age2.nt <- rep(test.prob.nevtest.f[2], length(elig.m.age2.nt))
    
    tst.m.age2.nt <- elig.m.age2.nt[rbinom(length(elig.m.age2.nt), 1, rates.m.age2.nt) == 1]
    
    #Age group 3
    elig.m.age3.nt <- which(sex == "F"  & age.group==3 & evertest==0 & partner.serv==0
                            (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age3.nt <- rep(test.prob.nevtest.f[3], length(elig.m.age3.nt))
    
    tst.m.age3.nt <- elig.m.age3.nt[rbinom(length(elig.m.age3.nt), 1, rates.m.age3.nt) == 1]
    
    #Repeat testers
    #Age group 1
    elig.m.age1.t <- which(sex == "F"  & age.group==1 & evertest==1 & partner.serv==0
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age1.t <- rep(test.prob.tested.f[1], length(elig.m.age1.t))
    
    tst.m.age1.t <- elig.m.age1.n[rbinom(length(elig.m.age1.n), 1, rates.m.age1.n) == 1]
    
    #Age group 2
    elig.m.age2.t <- which(sex == "F"  & age.group==2 & evertest==1 & partner.serv==0
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age2.t <- rep(test.prob.tested.f[2], length(elig.m.age2.t))
    
    tst.m.age2.t <- elig.m.age2.t[rbinom(length(elig.m.age2.t), 1, rates.m.age2.t) == 1]
    
    #Age group 3
    elig.m.age3.t <- which(sex == "F"  & age.group==3 & evertest==1 & partner.serv==0
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age3.t <- rep(test.prob.tested.f[3], length(elig.m.age3.t))
    
    tst.m.age3.t <- elig.m.age3.t[rbinom(length(elig.m.age3.t), 1, rates.m.age3.t) == 1]
    
    #Age group 4
    elig.m.age4.t <- which(sex == "F"  & age.group>=4 & evertest==1 & partner.serv==0
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    rates.m.age4.t <- rep(test.prob.tested.f[4], length(elig.m.age4.t))
    
    tst.m.age4.t <- elig.m.age4.t[rbinom(length(elig.m.age4.t), 1, rates.m.age4.t) == 1]
    
    

    tst.background <- c(tst.f.age1.nt,tst.f.age2.nt,tst.f.age3.nt,
                   tst.f.age1.t,tst.f.age2.t,tst.f.age3.t,tst.f.age4.t,
                   tst.m.age1.nt,tst.m.age2.nt,tst.m.age3.nt,
                   tst.m.age1.t,tst.m.age2.t,tst.m.age3.t,tst.m.age4.t)


    #seek.hc.AHI.prob
    #sym.prob
    #seek.hc.prob.f
    #seek.hc.prob.m
    #PITC.prob

    ##TEST
    if(intervention_TM=="None"){
    #ALL TESTS ARE BACKGROUND SO SHOULD BE SUBTRACTED OUT OF BACKGROUND AND NO PARTNER SERVICES  
    ##WHO HAS HIV SYMPTOMS AND SELECT THEM FOR TESTING
      
    sym_HIV <- which(status == 1 & inf.time > (at - 3) & diag.status ==1  & partner.serv==0)
    select <- rbinom(length(sym_HIV),1,seek.hc.AHI.prob*PITC.prob)
    sym_HIV <- sym_HIV[select == 1]

    #Count total number of background tests.  AHI symptom and other symptom tests will need to be removed from the background 
    #to be counted seperately
    ntest<-length(tst.background)
    
    #remove individuals selected for testing as both background and AHI symptoms
    tst.background <- tst.background[!(tst.background %in% sym_HIV)]

    
    ##WHICH ARE TESTING DUE TO HAVING SYMPTOMS (None AHI)
    tst.sym.f <- which(sex == "F"  & partner.serv==0 & (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    select <- rbinom(length(tst.sym.f),1,(sym.prob * seek.hc.prob * PITC.prob))
    tst.sym.f <- tst.sym.f[select ==1] 
    
    tst.sym.m <- which(sex == "M"  & partner.serv==0 & (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    select <- rbinom(length(tst.sym.m),1,(sym.prob * seek.hc.prob * PITC.prob))
    tst.sym.m <- tst.sym.m[select ==1]
    
    #remove individuals selected for testing as both background and none AHI symptoms from list for background
    tst.background <- tst.background[!(tst.background %in% tst.sym.f)]
    tst.background <- tst.background[!(tst.background %in% tst.sym.m)]
    
    n.keep <- ntest-length(sym_HIV)-length(tst.sym.m)-length(tst.sym.f)
    tst.background <- sample(tst.background,n.keep,replace=FALSE)
    
    ##TEST tst.background, sym_HIV, tst.sym.m, tst.sym.f   -  all Antibody tests
    
    tst.ab <- c(tst.background, sym_HIV, tst.sym.m, tst.sym.f)
    
    tst.pos.ab <- tst.ab[status[tst.ab] == 1 & inf.time[tst.ab] <= at - twind.int.ab ]
    tst.neg.ab <- setdiff(tst.ab, tst.pos.ab)
    
    }
    
    if(intervention_TM=="OPT"){
      #THE SYMPT TESTS THAT WOULD BE DONE ANTWAY ARE BACKGROUND BUT WE NEED TO ADD NEW TESTS FOR BETTER PITC 
      ##WHO HAS HIV SYMPTOMS AND SELECT THEM FOR TESTING
      
      sym_HIV <- which(status == 1 & inf.time > (at - 3) & diag.status ==1  & partner.serv==0)
      select1 <- rbinom(length(sym_HIV),1,seek.hc.AHI.prob*PITC.prob)
      select2 <- rbinom(length(sym_HIV),1,seek.hc.AHI.prob*PITC.TMP.prob)
      sym_HIV.bg <- sum(select1)
      sym_HIV <- sym_HIV[select2 == 1]
      
      #Count total number of background tests.  AHI symptom and other symptom tests will need to be removed from the background 
      #to be counted seperately
      ntest<-length(tst.background)
      
      #remove individuals selected for testing as both background and AHI symptoms
      tst.background <- tst.background[!(tst.background %in% sym_HIV)]
      
      
      ##WHICH ARE TESTING DUE TO HAVING SYMPTOMS (None AHI)
      tst.sym.f <- which(sex == "F"  & partner.serv==0 & (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
      select1 <- rbinom(length(tst.sym.f),1,(sym.prob * seek.hc.prob * PITC.prob))
      select2 <- rbinom(length(tst.sym.f),1,(sym.prob * seek.hc.prob * PITC.TMP.prob))
      tst.sym.f <- tst.sym.f[select2 ==1] 
      tst.sym.f.bg <- sum(select1)
      
      tst.sym.m <- which(sex == "M"  & partner.serv==0 & (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
      select1 <- rbinom(length(tst.sym.m),1,(sym.prob * seek.hc.prob * PITC.prob))
      select2 <- rbinom(length(tst.sym.m),1,(sym.prob * seek.hc.prob * PITC.TMP.prob))
      tst.sym.m <- tst.sym.m[select ==2]
      tst.sym.m.bg <- sum(select1)
      
      #remove individuals selected for testing as both background and none AHI symptoms from list for background
      tst.background <- tst.background[!(tst.background %in% tst.sym.f)]
      tst.background <- tst.background[!(tst.background %in% tst.sym.m)]
      
      n.keep <- ntest - sym_HIV.bg - tst.sym.f.bg - tst.sym.f.bg
      tst.background <- sample(tst.background,n.keep,replace=FALSE)
      
      ##TEST tst.background, sym_HIV, tst.sym.m, tst.sym.f.   -  all Antibody tests
      
      tst.ab <- c(tst.background, sym_HIV, tst.sym.m, tst.sym.f)
      
      tst.pos.ab <- tst.ab[status[tst.ab] == 1 & inf.time[tst.ab] <= at - twind.int.ab ]
      tst.neg.ab <- setdiff(tst.ab, tst.pos.ab)
    }
    
    if(intervention_TM=="TMP"){
      #THE SYMPT TESTS THAT WOULD BE DONE ANTWAY ARE BACKGROUND BUT WE NEED TO ADD NEW TESTS FOR BETTER PITC 
      ##WHO HAS HIV SYMPTOMS AND SELECT THEM FOR TESTING
      
      sym_HIV <- which(status == 1 & inf.time > (at - 3) & diag.status ==1  & partner.serv==0)
      select1 <- rbinom(length(sym_HIV),1,seek.hc.AHI.prob*PITC.prob)
      select2 <- rbinom(length(sym_HIV),1,seek.hc.AHI.prob*PITC.TMP.prob)
      sym_HIV.bg <- sum(select1)
      sym_HIV <- sym_HIV[select2 == 1]
      
      #Count total number of background tests.  AHI symptom and other symptom tests will need to be removed from the background 
      #to be counted seperately
      ntest<-length(tst.background)
      
      #remove individuals selected for testing as both background and AHI symptoms
      tst.background <- tst.background[!(tst.background %in% sym_HIV)]
      
      
      ##WHICH ARE TESTING DUE TO HAVING SYMPTOMS (None AHI)
      tst.sym.f <- which(sex == "F"  & partner.serv==0 & (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
      select1 <- rbinom(length(tst.sym.f),1,(sym.prob * seek.hc.prob * PITC.prob))
      select2 <- rbinom(length(tst.sym.f),1,(sym.prob * seek.hc.prob * PITC.TMP.prob))
      tst.sym.f <- tst.sym.f[select2 ==1] 
      tst.sym.f.bg <- sum(select1)
      
      tst.sym.m <- which(sex == "M"  & partner.serv==0 & (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
      select1 <- rbinom(length(tst.sym.m),1,(sym.prob * seek.hc.prob * PITC.prob))
      select2 <- rbinom(length(tst.sym.m),1,(sym.prob * seek.hc.prob * PITC.TMP.prob))
      tst.sym.m <- tst.sym.m[select ==2]
      tst.sym.m.bg <- sum(select1)
      
      #remove individuals selected for testing as both background and none AHI symptoms from list for background
      tst.background <- tst.background[!(tst.background %in% tst.sym.f)]
      tst.background <- tst.background[!(tst.background %in% tst.sym.m)]
      
      n.keep <- ntest - sym_HIV.bg - tst.sym.f.bg - tst.sym.f.bg
      tst.background <- sample(tst.background,n.keep,replace=FALSE)

      
      #PARTNER SERVICES AND PREP
      
      # PrEP testing
      tst.prep <- which((diag.status == 0 | is.na(diag.status)) &
                          prepStat == 1  &
                          tsincelntst >= prep.tst.int)
      

      ##TEST tst.background, tst.PREP, partner services -  Antibody tests
      ##sym_HIV, tst.sym.m, tst.sym.f.   -  RNA tests
      
      tst.ab <- c(tst.background, tst.prep, tst.ps)
      tst.rna <- c(sym_HIV, tst.sym.m, tst.sym.f)
      
      tst.pos.ab <- tst.ab[status[tst.ab] == 1 & inf.time[tst.ab] <= at - twind.int.ab ]
      tst.pos.rna <- tst.rna[status[tst.rna] == 1 & inf.time[tst.rna] <= at - twind.int.rna ]
      tst.neg.ab <- setdiff(tst.ab, tst.pos.ab)
      tst.neg.rna <- setdiff(tst.rna, tst.pos.rna)
      
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
  
  dat$attr$age.diag[tst.pos.ab] <- dat$attr$age[tst.pos.ab]
  dat$attr$age.diag[tst.pos.rna] <- dat$attr$age[tst.pos.rna]
  
  # Prevalence
  dat$epi$n.tests <- sum(length(tst.neg.ab),length(tst.neg.rna),length(tst.pos.ab),length(tst.pos.rna))
  dat$epi$n.tests.ab <- sum(length(tst.neg.ab),length(tst.pos.ab))
  dat$epi$n.tests.rna <- sum(length(tst.neg.rna),length(tst.pos.rna))

  acute <- which(inf.time[tst.pos.rna] <= at - twind.int.ab )
  actute <- max(length(acute),0)
  dat$epi$diag.prevalent <- sum(length(tst.neg.ab),length(tst.neg.rna),length(tst.pos.ab),length(tst.pos.rna)) - acute
  dat$epi$diag.acute <- acute
  
  return(dat)
}

