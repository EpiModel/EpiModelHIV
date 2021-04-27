
#' @title HIV Testing Module for hometesting project:Testing interventions
#'
#' @description Module function for HIV diagnostic testing of infected persons using clinic and hometests following 
#' four testing types (None testers (can test under supplimentations interventions), opportunistic only testers
#' regular testers and risk based testers. All tester use opportunistic tests and there replacement or 
#' supplementation effects opportunistic tests for all tester types.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This testing module supports two testing interventions, replacement of exsisting test with home based tests
#' and supplimentation of existing tests with homebased tests.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{last.neg.test},
#' \code{diag.status} and \code{diag.time} attributes.
#'
#' @keywords module
#' @export
#'



testing_interventions <- function(dat, at) {
  

  ## Attributes
  active <- dat$attr$active
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  tt.traj <- dat$attr$tt.traj
  status <- dat$attr$status
  inf.time <- dat$attr$inf.time
  ai.class <- dat$attr$ai.class
  tt <- dat$attr$tt
  
  
  #Parameters
  twind.int <- dat$param$test.window.int
  twind.hometest.int <- dat$param$hometest.window.int
  max.risk.test.time <- dat$param$max.risk.test.time   
  
  #Testing protocals
  Opportunity.replace <- dat$param$Opportunity.replace
  Opportunity.supp <- dat$param$Opportunity.supp
  Regular.replace <- dat$param$Regular.replace
  Regular.supp <- dat$param$Regular.supp
  Regular.supp2 <- dat$param$Regular.supp2
  Risk.replace <- dat$param$Risk.replace
  Risk.supp <- dat$param$Risk.supp
  Never.test.supp <- dat$param$Never.test.supp
  HT.start <- dat$param$HT.start
  
  hometest.rep.opp.frac <- dat$param$hometest.rep.opp.frac
  hometest.rep.reg.frac <- dat$param$hometest.rep.reg.frac
  hometest.rep.risk.frac <- dat$param$hometest.rep.risk.frac
  
  #Probabilities
  mean.test.B.opp.prob <- dat$param$mean.test.B.opp.prob 
  mean.test.W.opp.prob <- dat$param$mean.test.W.opp.prob 
  mean.test.B.opp.reg.prob <- dat$param$mean.test.B.opp.reg.prob
  mean.test.W.opp.reg.prob<- dat$param$mean.test.W.opp.reg.prob
  mean.test.B.opp.risk.prob <- dat$param$mean.test.B.opp.risk.prob 
  mean.test.W.opp.risk.prob <-dat$param$ mean.test.W.opp.risk.prob
  
  supp.opp.home.test.prob <- dat$param$supp.opp.home.test.prob
  supp.reg.home.test.prob <- dat$param$supp.reg.home.test.prob
  supp.risk.home.test.prob <- dat$param$supp.risk.home.test.prob
  supp.nevertest.home.test.prob <- dat$param$supp.nevertest.home.test.prob
  
  #Regular tester testing intervals
  mean.test.B.low.int <- dat$param$mean.test.B.low.int
  mean.test.B.high.int <- dat$param$mean.test.B.high.int
  mean.test.W.low.int <- dat$param$mean.test.W.low.int
  mean.test.W.high.int <- dat$param$mean.test.W.high.int
  
  #Regular tester testing intervals when supplimented with extra annual test
  mean.test.B.low.int.supp2 <- dat$param$mean.test.B.low.int.supp2
  mean.test.B.high.int.supp2 <- dat$param$mean.test.B.high.int.supp2
  mean.test.W.low.int.supp2 <- dat$param$mean.test.W.low.int.supp2
  mean.test.W.high.int.supp2 <- dat$param$mean.test.W.high.int.supp2

  
  #Tester type specific opportunistic testing interval
  mean.test.opp.NO.W.int  <- dat$param$mean.test.opp.NO.W.int
  mean.test.opp.NO.B.int <- dat$param$mean.test.opp.NO.B.int
  mean.test.opp.ReT.W.int <- dat$param$mean.test.opp.ReT.W.int
  mean.test.opp.ReT.B.int <- dat$param$mean.test.opp.ReT.B.int
  mean.test.opp.Risk.W.int <- dat$param$mean.test.opp.Risk.W.int
  mean.test.opp.Risk.B.int  <- dat$param$mean.test.opp.Risk.B.int
 
  #Risk type specific testing interval
  mean.test.risk.CAI.nonmain.B.int <- dat$param$mean.test.risk.CAI.nonmain.B.int
  mean.test.risk.CAI.nonmain.W.int <- dat$param$mean.test.risk.CAI.nonmain.W.int
  mean.test.risk.AI.known.sd.B.int <- dat$param$mean.test.risk.AI.known.sd.B.int
  mean.test.risk.AI.known.sd.W.int <- dat$param$mean.test.risk.AI.known.sd.W.int
  mean.test.risk.newmain.B.int <- dat$param$mean.test.risk.newmain.B.int
  mean.test.risk.newmain.W.int <- dat$param$mean.test.risk.newmain.W.int
  
  #Risk based testing probabilities
  mean.test.risk.CAI.nonmain.B.prob <- dat$param$mean.test.risk.CAI.nonmain.B.prob
  mean.test.risk.CAI.nonmain.W.prob <- dat$param$mean.test.risk.CAI.nonmain.W.prob
  mean.test.risk.AI.known.sd.B.prob <- dat$param$mean.test.risk.AI.known.sd.B.prob
  mean.test.risk.AI.known.sd.W.prob <- dat$param$mean.test.risk.AI.known.sd.W.prob
  mean.test.risk.newmain.B.prob <- dat$param$mean.test.risk.newmain.B.prob
  mean.test.risk.newmain.W.prob <- dat$param$mean.test.risk.newmain.W.prob
  
  ##PrEP testing 
  prepStat <- dat$attr$prepStat 

  
#####  tests
#  if(at > 1){
    #Opportunity.replace <- TRUE
    #Opportunity.supp <- TRUE
    #Regular.replace <- TRUE
    #Regular.supp <- TRUE
    #Risk.replace <- TRUE
    #Risk.supp <- TRUE
    #Never.test.supp <- TRUE
#  }
  
  
  ## Process

  #Clear the vectors with NULL values
  #NULL VALUES
    tst.all <- tst.B <- tst.W <- tst.B.real <- tst.W.real<- tst.all.real <- 
    tst.neg <- tst.neg.clinic <- tst.neg.home <- tst.neg.home.s <- tst.neg.not.supp <- 
    tst.pos <- tst.pos.B <- tst.pos.W <- tst.pos.clinic <- tst.pos.home <- tst.pos.home.s <-
    tst.pos.B.home <- tst.pos.B.clinic <- tst.pos.B.home.s <-
    tst.pos.W.home <- tst.pos.W.clinic <- tst.pos.W.home.s <-
    NO.supp.B.real <- NO.supp.W.real <- ReT.supp.B.real <- ReT.supp.W.real <-
    Risk.supp.B.real <- Risk.supp.W.real <- NN.supp.B.real <- NN.supp.W.real <-
    clinictestIDs.B_Opp <- clinictestIDs.W_Opp <- hometestIDs.B_Opp <- hometestIDs.W_Opp <- tst.pos.home.comb <- NULL 
  

  # OPPORTUNISTIC ONLY TESTERS "NO" 
  # Frequency of testing is based on time since last opportnistic test being greater than 
  # a fixed interval for opportunistic tests (mean.test.opp.NO.R.int)
  
  #Get time since last test
  tsincelntst <- at - dat$attr$last.neg.test.opp
  tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]
  
  #Select those eligible to test based on time since last test
  # if the time since their last test is larger than the mean test interval, then there is a chance they test.
  tst.B <- which(active == 1 & race == "B" & tt == "NO" &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= mean.test.opp.NO.B.int
                   & prepStat == 0)  
  tst.W <- which(active == 1 & race == "W" & tt== "NO" &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= mean.test.opp.NO.W.int
                   & prepStat == 0)
  
  # Select a sample of those eligible to take a test to really test based on the testing rate mean.test.B.opp.rate
  tst.B.real <- tst.B[runif(length(tst.B)) < mean.test.B.opp.prob]
  tst.W.real <- tst.W[runif(length(tst.W)) < mean.test.W.opp.prob]

  #The clinic based opportunistic testers
  clinictestIDs.B_Opp <- tst.B.real
  clinictestIDs.W_Opp <- tst.W.real
  
    
  ## REPLACEMENT: Opportunistic test replacement with home tests
  ## If tests are being replaced the testers that are selected are split into two groups (clinic tests and home tests)
  ## The fraction of the tests that are home tests is determined by hometest.rep.frac.

  ##Do the lengths of the home and clinic sum to the original clinic
  length(clinictestIDs.B_Opp)
  length(clinictestIDs.W_Opp)
  
  if (Opportunity.replace == TRUE & at >= HT.start){
      hometestIDs.B_Opp <- tst.B.real[runif(length(tst.B.real)) < hometest.rep.opp.frac]
      clinictestIDs.B_Opp <- setdiff(tst.B.real,hometestIDs.B_Opp)
  
  hometestIDs.W_Opp <- tst.W.real[runif(length(tst.W.real)) < hometest.rep.opp.frac]
  clinictestIDs.W_Opp <- setdiff(tst.W.real,hometestIDs.W_Opp)}
  
  ## SUPPLIMENT: Opportunistic only testers take additional home tests
  if (Opportunity.supp == TRUE & at >= HT.start){
    #Select those eligible for a supplimental test
    NO.supp.B <- which(active == 1 & race == "B" & tt == "NO" &
                     (diag.status == 0 | is.na(diag.status)) & prepStat == 0)  
    NO.supp.W <- which(active == 1 & race == "W" & tt== "NO" &
                     (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    
    ##Remove IDS already taking an opp test.
    NO.supp.B <- setdiff(NO.supp.B,tst.B.real)
    NO.supp.W <- setdiff(NO.supp.W,tst.W.real)
    
    #Select those that will take a supplimental home test based on probability supp.home.test.prob
    NO.supp.B.real <- NO.supp.B[runif(length(NO.supp.B)) < supp.opp.home.test.prob]
    NO.supp.W.real <- NO.supp.W[runif(length(NO.supp.W)) < supp.opp.home.test.prob]    
  }
  

  ##Take the tests
  #clinic test
  
  tst.pos.B.clinic <- clinictestIDs.B_Opp[status[clinictestIDs.B_Opp] == 1 & inf.time[clinictestIDs.B_Opp] <= at - twind.int]
  tst.neg.B.clinic <- setdiff(clinictestIDs.B_Opp, tst.pos.B.clinic)
  
  tst.pos.W.clinic <- clinictestIDs.W_Opp[status[clinictestIDs.W_Opp] == 1 & inf.time[clinictestIDs.W_Opp] <= at - twind.int]
  tst.neg.W.clinic <- setdiff(clinictestIDs.W_Opp, tst.pos.W.clinic)
  
  tst.pos.clinic <- c(tst.pos.B.clinic, tst.pos.W.clinic)
  tst.neg.clinic <- c(tst.neg.B.clinic, tst.neg.W.clinic)
  
 
  ##Home tests: replacements for opportunistic tests
  if (Opportunity.replace == TRUE & at >= HT.start){
  #home test
  tst.pos.B.home <- hometestIDs.B_Opp[status[hometestIDs.B_Opp] == 1 & inf.time[hometestIDs.B_Opp] <= at - twind.hometest.int]
  tst.neg.B.home <- setdiff(hometestIDs.B_Opp, tst.pos.B.home)
  
  tst.pos.W.home <- hometestIDs.W_Opp[status[hometestIDs.W_Opp] == 1 & inf.time[hometestIDs.W_Opp] <= at - twind.hometest.int]
  tst.neg.W.home <- setdiff(hometestIDs.W_Opp, tst.pos.W.home)
  
  tst.pos.home <- c(tst.pos.B.home, tst.pos.W.home)
  tst.neg.home <- c(tst.neg.B.home, tst.neg.W.home)
  
  }
  
  ##Home tests: supplimental tests
  if (Opportunity.supp == TRUE & at >= HT.start){
    #home test
    tst.pos.B.home.s <- NO.supp.B.real[status[NO.supp.B.real] == 1 & inf.time[NO.supp.B.real] <= at - twind.hometest.int]
    tst.neg.B.home.s <- setdiff(NO.supp.B.real, tst.pos.B.home.s)
    
    tst.pos.W.home.s <- NO.supp.W.real[status[NO.supp.W.real] == 1 & inf.time[NO.supp.W.real] <= at - twind.hometest.int]
    tst.neg.W.home.s <- setdiff(NO.supp.W.real, tst.pos.W.home.s)
    
    tst.pos.home.s <- c(tst.pos.B.home.s, tst.pos.W.home.s)
    tst.neg.home.s <- c(tst.neg.B.home.s, tst.neg.W.home.s)
    
  }
  
  
  tst.neg <- c(tst.neg.clinic, tst.neg.home,tst.neg.home.s)
  tst.neg.not.supp <- c(tst.neg.clinic, tst.neg.home)
  tst.pos <- c(tst.pos.clinic, tst.pos.home,tst.pos.home.s)
  tst.pos.clinic <- c(tst.pos.clinic)
  tst.pos.home.comb <- c(tst.pos.home,tst.pos.home.s)
  tst.pos.B <- c(tst.pos.B.home, tst.pos.B.clinic,tst.pos.B.home.s)
  tst.pos.W <- c(tst.pos.W.home, tst.pos.W.clinic,tst.pos.W.home.s)
  
  tst.all.real <- c(tst.B.real, tst.W.real, NO.supp.B.real, NO.supp.W.real)
  tst.all.clinic <- c(tst.pos.clinic, tst.neg.clinic)
  tst.all.HT <- c(tst.pos.home, tst.neg.home, tst.pos.home.s, tst.neg.home.s)
  tst.all.opp.elig <- c(tst.B, tst.W)
  tst.all.opp.elig <- setdiff(tst.all.opp.elig,tst.pos)
  
  # Attributes
  dat$attr$diag.opp.clin.test[tst.pos.clinic] <- 1
  dat$attr$diag.opp.HT.test[tst.pos.home.comb] <- 1
  dat$attr$last.neg.test.opp[tst.all.opp.elig] <- at
  dat$attr$diag.status[tst.pos] <- 1
  dat$attr$diag.time[tst.pos] <- at
  

  ### Summary statistics
  dat$epi$opp.clin.test.num[at] <- max(0,length(tst.pos.clinic[!is.na(tst.pos.clinic)]) + length(tst.neg.clinic[!is.na(tst.neg.clinic)])) 
  dat$epi$opp.home.test.num[at] <- max(0,length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.neg.home[!is.na(tst.neg.home)]) + 
                                         length(tst.pos.home.s[!is.na(tst.pos.home.s)]) + length(tst.neg.home.s[!is.na(tst.neg.home.s)]))
  dat$epi$opp.clin.test.pos[at] <- max(0,length(tst.pos.clinic[!is.na(tst.pos.clinic)]))
  dat$epi$opp.home.test.pos[at] <- max(0,length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.pos.home.s[!is.na(tst.pos.home.s)]))

  ###########ADD THESE COUNTS
  dat$attr$Opp.Test.Num[tst.all.real] <- dat$attr$Opp.Test.Num[tst.all.real] + 1
  dat$attr$Opp.Clin.Test.Num[tst.all.clinic] <- dat$attr$Opp.Clin.Test.Num[tst.all.clinic] + 1
  #dat$attr$Opp.HT.Test.Num[tst.all.HT] <- dat$attr$Opp.HT.Test.Num[tst.all.HT] + 1
  
  
###########################################################################################################

##Regular testers
  
  #Clear the vectors with NULL values
  #NULL VALUES
  tst.all <- tst.B <- tst.W <- tst.B.real <- tst.W.real<- tst.all.real <- 
    tst.neg <- tst.neg.clinic <- tst.neg.home <- tst.neg.home.s <- tst.neg.not.supp <- 
    tst.pos <- tst.pos.B <- tst.pos.W <- tst.pos.clinic <- tst.pos.home <- tst.pos.home.s <-
    tst.pos.B.home <- tst.pos.B.clinic <- tst.pos.B.home.s <-
    tst.pos.W.home <- tst.pos.W.clinic <- tst.pos.W.home.s <-
    NO.supp.B.real <- NO.supp.W.real <- ReT.supp.B.real <- ReT.supp.W.real <-
    Risk.supp.B.real <- Risk.supp.W.real <- NN.supp.B.real <- NN.supp.W.real <-
    clinictestIDs.B_Opp <- clinictestIDs.W_Opp <- hometestIDs.B_Opp <- hometestIDs.W_Opp <- tst.pos.home.comb <- NULL 
  
  #Update current status
  diag.status <- dat$attr$diag.status
  
  ##Select those that will test.
  # Regular tests
  tsincelntst <- at - dat$attr$last.neg.test
  tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]
  
  tst.B.l <- which(active == 1 & race == "B" & ai.class == "Low_AI" & tt == "ReT" &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= mean.test.B.low.int & prepStat == 0)  # if the time since their last test is larger than the mean test interval, then there is a chance they test.
  tst.B.h <- which(active == 1 & race == "B" & ai.class == "High_AI" & tt == "ReT" &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= mean.test.B.high.int & prepStat == 0)  # if the time since their last test is larger than the mean test interval, then there is a chance they test.
  tst.B <- c(tst.B.l, tst.B.h)
  
  tst.W.l <- which(active == 1 & race == "W" & ai.class == "Low_AI" & tt== "ReT" &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= mean.test.W.low.int & prepStat == 0)
  tst.W.h <- which(active == 1 & race == "W" & ai.class == "High_AI" & tt== "ReT" &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= mean.test.W.high.int & prepStat == 0)
  tst.W <- c(tst.W.l, tst.W.h)
  
  tst.all <- c(tst.W, tst.B)
  
  #The clinic based regular testers
  clinictestIDs.B_Opp <- tst.B
  clinictestIDs.W_Opp <- tst.W
  
  
  ## REPLACEMENT: Regular test replacement with home tests
  ## If tests are being replaced the testers that are selected are split into two groups (clinic tests and home tests)
  ## The fraction of the tests that are home tests is determined by hometest.rep.reg.frac.
  
  if (Regular.replace == TRUE & at >= HT.start){
    hometestIDs.B_Opp <- tst.B[runif(length(tst.B)) < hometest.rep.reg.frac]
    clinictestIDs.B_Opp <- setdiff(tst.B,hometestIDs.B_Opp)
    
    hometestIDs.W_Opp <- tst.W[runif(length(tst.W)) < hometest.rep.reg.frac]
    clinictestIDs.W_Opp <- setdiff(tst.W,hometestIDs.W_Opp)}
  
  ## SUPPLIMENT: Regular testers may add additioanl tests at hazard supp.reg.home.test.prob
  if (Regular.supp == TRUE & at >= HT.start){
    #Select those eligible for a supplimental test
    ReT.supp.B <- which(active == 1 & race == "B" & tt == "ReT" &
                         (diag.status == 0 | is.na(diag.status)) & prepStat == 0)  
    ReT.supp.W <- which(active == 1 & race == "W" & tt== "ReT" &
                         (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    
    ##Remove IDS already taking an opp test.
    ReT.supp.B <- setdiff(ReT.supp.B,tst.B)
    ReT.supp.W <- setdiff(ReT.supp.W,tst.W)
    
    #Select those that will take a supplimental home test based on probability supp.home.test.prob
    ReT.supp.B.real <- ReT.supp.B[runif(length(ReT.supp.B)) < supp.reg.home.test.prob]
    ReT.supp.W.real <- ReT.supp.W[runif(length(ReT.supp.W)) < supp.reg.home.test.prob]    
  }
  
  ##Supplement with one more test per year by decreasing the regular testing interval (mean.test....supp2)
  ##and using a home test for the additional test with prob hometest.rep.in.reg.supp2
  if (Regular.supp2 == TRUE  & at >= HT.start){
    #NULL VALUES
    tst.all <- tst.B <- tst.W <- tst.B.real <- tst.W.real<- tst.all.real <- 
      tst.neg <- tst.neg.clinic <- tst.neg.home <- tst.neg.home.s <- tst.neg.not.supp <- 
      tst.pos <- tst.pos.B <- tst.pos.W <- tst.pos.clinic <- tst.pos.home <- tst.pos.home.s <-
      tst.pos.B.home <- tst.pos.B.clinic <- tst.pos.B.home.s <-
      tst.pos.W.home <- tst.pos.W.clinic <- tst.pos.W.home.s <-
      NO.supp.B.real <- NO.supp.W.real <- ReT.supp.B.real <- ReT.supp.W.real <-
      Risk.supp.B.real <- Risk.supp.W.real <- NN.supp.B.real <- NN.supp.W.real <-
      clinictestIDs.B_Opp <- clinictestIDs.W_Opp <- hometestIDs.B_Opp <- hometestIDs.W_Opp <- tst.pos.home.comb <- NULL 
    
    
    tsincelntst <- at - dat$attr$last.neg.test
    tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]
    
    tst.B.l <- which(active == 1 & race == "B" & ai.class == "Low_AI" & tt == "ReT" &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= mean.test.B.low.int.supp2 & prepStat == 0)  # if the time since their last test is larger than the mean test interval, then there is a chance they test.
    tst.B.h <- which(active == 1 & race == "B" & ai.class == "High_AI" & tt == "ReT" &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= mean.test.B.high.int.supp2 & prepStat == 0)  # if the time since their last test is larger than the mean test interval, then there is a chance they test.

    tst.W.l <- which(active == 1 & race == "W" & ai.class == "Low_AI" & tt== "ReT" &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= mean.test.W.low.int.supp2 & prepStat == 0)
    tst.W.h <- which(active == 1 & race == "W" & ai.class == "High_AI" & tt== "ReT" &
                       (diag.status == 0 | is.na(diag.status)) &
                       tsincelntst >= mean.test.W.high.int.supp2 & prepStat == 0)

    
    
    ## Replace the increased fractions of tests with home tests

      tst.B.l.home <- tst.B.l[runif(length(tst.B.l)) < 1/(52/mean.test.B.low.int.supp2)]
      tst.B.l.clin <- setdiff(tst.B.l,tst.B.l.home)
      
      tst.B.h.home <- tst.B.h[runif(length(tst.B.h)) < 1/(52/mean.test.B.high.int.supp2)]
      tst.B.h.clin <- setdiff(tst.B.h,tst.B.h.home)
      
      tst.W.l.home <- tst.W.l[runif(length(tst.W.l)) < 1/(52/mean.test.W.low.int.supp2)]
      tst.W.l.clin <- setdiff(tst.W.l,tst.W.l.home)
      
      tst.W.h.home <- tst.W.h[runif(length(tst.W.h)) < 1/(52/mean.test.W.high.int.supp2)]
      tst.W.h.clin <- setdiff(tst.W.h,tst.W.h.home)
      
      
      hometestIDs.B_Opp <- c(tst.B.l.home, tst.B.h.home)
      clinictestIDs.B_Opp <- c(tst.B.l.clin, tst.B.h.clin)
  
      hometestIDs.W_Opp <- c(tst.W.l.home, tst.W.h.home)
      clinictestIDs.W_Opp <- c(tst.W.l.clin, tst.W.h.clin)
    
  }
  
 #TAKE the assigned tests  
  
  #clinic test
  
  tst.pos.B.clinic <- clinictestIDs.B_Opp[status[clinictestIDs.B_Opp] == 1 & inf.time[clinictestIDs.B_Opp] <= at - twind.int]
  tst.neg.B.clinic <- setdiff(clinictestIDs.B_Opp, tst.pos.B.clinic)
  
  tst.pos.W.clinic <- clinictestIDs.W_Opp[status[clinictestIDs.W_Opp] == 1 & inf.time[clinictestIDs.W_Opp] <= at - twind.int]
  tst.neg.W.clinic <- setdiff(clinictestIDs.W_Opp, tst.pos.W.clinic)
  
  tst.pos.clinic <- c(tst.pos.B.clinic, tst.pos.W.clinic)
  tst.neg.clinic <- c(tst.neg.B.clinic, tst.neg.W.clinic)
  
  
  if(Regular.replace == TRUE & at >= HT.start){
    #home test
    tst.pos.B.home <- hometestIDs.B_Opp[status[hometestIDs.B_Opp] == 1 & inf.time[hometestIDs.B_Opp] <= at - twind.hometest.int]
    tst.neg.B.home <- setdiff(hometestIDs.B_Opp, tst.pos.B.home)
    
    tst.pos.W.home <- hometestIDs.W_Opp[status[hometestIDs.W_Opp] == 1 & inf.time[hometestIDs.W_Opp] <= at - twind.hometest.int]
    tst.neg.W.home <- setdiff(hometestIDs.W_Opp, tst.pos.W.home)
    
    tst.pos.home <- c(tst.pos.B.home, tst.pos.W.home)
    tst.neg.home <- c(tst.neg.B.home, tst.neg.W.home)
  }
  
  if(Regular.supp == TRUE & at >= HT.start){
    
    tst.pos.B.home.s <- ReT.supp.B.real[status[ReT.supp.B.real] == 1 & inf.time[ReT.supp.B.real] <= at - twind.hometest.int]
    tst.neg.B.home.s <- setdiff(ReT.supp.B.real, tst.pos.B.home.s)
    
    tst.pos.W.home.s <- ReT.supp.W.real[status[ReT.supp.W.real] == 1 & inf.time[ReT.supp.W.real] <= at - twind.hometest.int]
    tst.neg.W.home.s <- setdiff(ReT.supp.W.real, tst.pos.W.home.s)
    
    tst.pos.home.s <- c(tst.pos.B.home.s, tst.pos.W.home.s)
    tst.neg.home.s <- c(tst.neg.B.home.s, tst.neg.W.home.s)
  }
 
  if(Regular.supp2 == TRUE & at >= HT.start){
    #home test
    tst.pos.B.home.s <- hometestIDs.B_Opp[status[hometestIDs.B_Opp] == 1 & inf.time[hometestIDs.B_Opp] <= at - twind.hometest.int]
    tst.neg.B.home.s <- setdiff(hometestIDs.B_Opp, tst.pos.B.home.s)
    
    tst.pos.W.home.s <- hometestIDs.W_Opp[status[hometestIDs.W_Opp] == 1 & inf.time[hometestIDs.W_Opp] <= at - twind.hometest.int]
    tst.neg.W.home.s <- setdiff(hometestIDs.W_Opp, tst.pos.W.home.s)
    
    tst.pos.home.s <- c(tst.pos.B.home.s, tst.pos.W.home.s)
    tst.neg.home.s <- c(tst.neg.B.home.s, tst.neg.W.home.s)
  }
  
  tst.neg <- c(tst.neg.clinic, tst.neg.home,tst.neg.home.s)
  tst.neg.not.supp <- c(tst.neg.clinic, tst.neg.home)
  tst.pos <- c(tst.pos.clinic, tst.pos.home,tst.pos.home.s)
  tst.pos.clinic <- c(tst.pos.clinic)
  tst.pos.home.comb <- c(tst.pos.home,tst.pos.home.s)
  tst.pos.B <- c(tst.pos.B.home, tst.pos.B.clinic,tst.pos.B.home.s)
  tst.pos.W <- c(tst.pos.W.home, tst.pos.W.clinic,tst.pos.W.home.s)
  
  tst.all.real <- c(tst.B, tst.W, ReT.supp.B.real, ReT.supp.W.real)
  tst.all.ReT.elig <- c(tst.B, tst.W)
  
  tst.all.clinic <- c(tst.pos.clinic, tst.neg.clinic)
  tst.all.HT <- c(tst.pos.home, tst.neg.home, tst.pos.home.s, tst.neg.home.s)

  ## Output
  
  # Attributes
  dat$attr$diag.reg.clin.test[tst.pos.clinic] <- 1
  dat$attr$diag.reg.HT.test[tst.pos.home.comb] <- 1

  
  dat$attr$last.neg.test[tst.neg.not.supp] <- at
  if(Regular.supp2 == TRUE & at >= HT.start){
    dat$attr$last.neg.test[tst.neg] <- at
  }
  
  dat$attr$last.neg.test.in.clinic[tst.all.clinic] <- at
  
  dat$attr$diag.status[tst.pos] <- 1
  dat$attr$diag.time[tst.pos] <- at
  
  dat$attr$Reg.Test.Num[tst.all.real] <- dat$attr$Reg.Test.Num[tst.all.real ]  + 1
  dat$attr$Reg.Clin.Test.Num[tst.all.clinic] <- dat$attr$Reg.Clin.Test.Num[tst.all.clinic]  + 1
  dat$attr$Reg.HT.Test.Num[tst.all.HT] <- dat$attr$Reg.HT.Test.Num[tst.all.HT]  + 1
  
  
  ### Summary statistics
  dat$epi$reg.clin.test.num[at] <- max(0,length(tst.pos.clinic[!is.na(tst.pos.clinic)]) + length(tst.neg.clinic[!is.na(tst.neg.clinic)]))
  dat$epi$reg.home.test.num[at] <- max(0,length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.neg.home[!is.na(tst.neg.home)]) + 
                                         length(tst.pos.home.s[!is.na(tst.pos.home.s)]) + length(tst.neg.home.s[!is.na(tst.neg.home.s)]))
  dat$epi$reg.clin.test.pos[at] <- max(0,length(tst.pos.clinic[!is.na(tst.pos.clinic)]))
  dat$epi$reg.home.test.pos[at] <- max(0,length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.pos.home.s[!is.na(tst.pos.home.s)]))

############################################################################  
##Opportunuistic tests for regular testers
  
  #Clear the vectors with NULL values
  #NULL VALUES
    tst.all <- tst.B <- tst.W <- tst.B.real <- tst.W.real<- tst.all.real <- 
    tst.neg <- tst.neg.clinic <- tst.neg.home <- tst.neg.home.s <- tst.neg.not.supp <- 
    tst.pos <- tst.pos.B <- tst.pos.W <- tst.pos.clinic <- tst.pos.home <- tst.pos.home.s <-
    tst.pos.B.home <- tst.pos.B.clinic <- tst.pos.B.home.s <-
    tst.pos.W.home <- tst.pos.W.clinic <- tst.pos.W.home.s <-
    NO.supp.B.real <- NO.supp.W.real <- ReT.supp.B.real <- ReT.supp.W.real <-
    Risk.supp.B.real <- Risk.supp.W.real <- NN.supp.B.real <- NN.supp.W.real <-
    clinictestIDs.B_Opp <- clinictestIDs.W_Opp <- hometestIDs.B_Opp <- hometestIDs.W_Opp <- tst.pos.home.comb <- NULL 
  
    #Update current diagnosis status
    diag.status <- dat$attr$diag.status
  
  # Regular testers only "ReT" 
  # Frequency of testing is based on time since last opportnistic test being greater than 
  # a fixed interval for opportunistic tests for regular testers (mean.test.opp.NO.R.int)
  
  #Get time since last test
  tsincelntst <- at - dat$attr$last.neg.test.opp
  tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]
  
  #Select those eligible to test based on time since last test
  # if the time since their last test is larger than the mean test interval, then there is a chance they test.
  tst.B <- which(active == 1 & race == "B" & tt == "ReT" &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= mean.test.opp.ReT.B.int
                 & prepStat == 0)  
  tst.W <- which(active == 1 & race == "W" & tt== "ReT" &
                   (diag.status == 0 | is.na(diag.status)) &
                   tsincelntst >= mean.test.opp.ReT.W.int
                 & prepStat == 0)
  
  # Select a sample of those eligible to take a test to really test based on the testing rate mean.test.B.opp.rate
  tst.B.real <- tst.B[runif(length(tst.B)) < mean.test.B.opp.reg.prob]
  tst.W.real <- tst.W[runif(length(tst.W)) < mean.test.W.opp.reg.prob]
  
  #The clinic based opportunistic testers
  clinictestIDs.B_Opp <- tst.B.real
  clinictestIDs.W_Opp <- tst.W.real
  
  
  ## REPLACEMENT: Opportunistic test replacement with home tests
  ## If tests are being replaced the testers that are selected are split into two groups (clinic tests and home tests)
  ## The fraction of the tests that are home tests is determined by hometest.rep.frac.
  
  if (Opportunity.replace == TRUE  & at >= HT.start){
    hometestIDs.B_Opp <- tst.B.real[runif(length(tst.B.real)) < hometest.rep.opp.frac]
    clinictestIDs.B_Opp <- setdiff(tst.B.real,hometestIDs.B_Opp)
    
    hometestIDs.W_Opp <- tst.W.real[runif(length(tst.W.real)) < hometest.rep.opp.frac]
    clinictestIDs.W_Opp <- setdiff(tst.W.real,hometestIDs.W_Opp)}
  
  ## SUPPLIMENT: Opportunistic only testers take additional home tests
  if (Opportunity.supp == TRUE & at >= HT.start){
    #Select those eligible for a supplimental test
    ReT.supp.B <- which(active == 1 & race == "B" & tt == "ReT" &
                         (diag.status == 0 | is.na(diag.status)) & prepStat == 0)  
    ReT.supp.W <- which(active == 1 & race == "W" & tt== "ReT" &
                         (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
    
    ##Remove IDS already taking an opp test.
    ReT.supp.B <- setdiff(ReT.supp.B,tst.B.real)
    ReT.supp.W <- setdiff(ReT.supp.W,tst.W.real)
    
    #Select those that will take a supplimental home test based on probability supp.home.test.prob
    ReT.supp.B.real <- ReT.supp.B[runif(length(ReT.supp.B)) < supp.opp.home.test.prob]
    ReT.supp.W.real <- ReT.supp.W[runif(length(ReT.supp.W)) < supp.opp.home.test.prob]    
  }
  
  
  ##Take the tests
  #clinic test
  
  tst.pos.B.clinic <- clinictestIDs.B_Opp[status[clinictestIDs.B_Opp] == 1 & inf.time[clinictestIDs.B_Opp] <= at - twind.int]
  tst.neg.B.clinic <- setdiff(clinictestIDs.B_Opp, tst.pos.B.clinic)
  
  tst.pos.W.clinic <- clinictestIDs.W_Opp[status[clinictestIDs.W_Opp] == 1 & inf.time[clinictestIDs.W_Opp] <= at - twind.int]
  tst.neg.W.clinic <- setdiff(clinictestIDs.W_Opp, tst.pos.W.clinic)
  
  tst.pos.clinic <- c(tst.pos.B.clinic, tst.pos.W.clinic)
  tst.neg.clinic <- c(tst.neg.B.clinic, tst.neg.W.clinic)
  
  
  ##Home tests: replacements for opportunistic tests
  if (Opportunity.replace==TRUE & at >= HT.start){
    #home test
    tst.pos.B.home <- hometestIDs.B_Opp[status[hometestIDs.B_Opp] == 1 & inf.time[hometestIDs.B_Opp] <= at - twind.hometest.int]
    tst.neg.B.home <- setdiff(hometestIDs.B_Opp, tst.pos.B.home)
    
    tst.pos.W.home <- hometestIDs.W_Opp[status[hometestIDs.W_Opp] == 1 & inf.time[hometestIDs.W_Opp] <= at - twind.hometest.int]
    tst.neg.W.home <- setdiff(hometestIDs.W_Opp, tst.pos.W.home)
    
    tst.pos.home <- c(tst.pos.B.home, tst.pos.W.home)
    tst.neg.home <- c(tst.neg.B.home, tst.neg.W.home)
    
  }
  
  ##Home tests: supplimental tests
  if (Opportunity.supp == TRUE & at >= HT.start){
    #home test
    tst.pos.B.home.s <- ReT.supp.B.real[status[ReT.supp.B.real] == 1 & inf.time[ReT.supp.B.real] <= at - twind.hometest.int]
    tst.neg.B.home.s <- setdiff(ReT.supp.B.real, tst.pos.B.home.s)
    
    tst.pos.W.home.s <- ReT.supp.W.real[status[ReT.supp.W.real] == 1 & inf.time[ReT.supp.W.real] <= at - twind.hometest.int]
    tst.neg.W.home.s <- setdiff(ReT.supp.W.real, tst.pos.W.home.s)
    
    tst.pos.home.s <- c(tst.pos.B.home.s, tst.pos.W.home.s)
    tst.neg.home.s <- c(tst.neg.B.home.s, tst.neg.W.home.s)
    
  }
  
  
  tst.neg <- c(tst.neg.clinic, tst.neg.home,tst.neg.home.s)
  tst.neg.not.supp <- c(tst.neg.clinic, tst.neg.home)
  tst.pos <- c(tst.pos.clinic, tst.pos.home,tst.pos.home.s)
  tst.pos.clinic <- c(tst.pos.clinic)
  tst.pos.home.comb <- c(tst.pos.home,tst.pos.home.s)
  tst.pos.B <- c(tst.pos.B.home, tst.pos.B.clinic,tst.pos.B.home.s)
  tst.pos.W <- c(tst.pos.W.home, tst.pos.W.clinic,tst.pos.W.home.s)
  
  tst.all.real <- c(tst.B.real, tst.W.real, ReT.supp.B.real, ReT.supp.W.real)
  tst.all.opp.elig <- c(tst.B, tst.W)
  tst.all.opp.elig <- setdiff(tst.all.opp.elig,tst.pos)
  tst.all.clinic <- c(tst.pos.clinic, tst.neg.clinic)
  tst.all.HT <- c(tst.pos.home, tst.neg.home, tst.pos.home.s, tst.neg.home.s)
  
  # Attributes
  dat$attr$diag.opp.clin.test[tst.pos.clinic] <- 1
  dat$attr$diag.opp.HT.test[tst.pos.home.comb] <- 1
  dat$attr$last.neg.test.opp[tst.all.opp.elig] <- at
  dat$attr$diag.status[tst.pos] <- 1
  dat$attr$diag.time[tst.pos] <- at
  
  dat$attr$Opp.Test.Num[tst.all.real] <- dat$attr$Opp.Test.Num[tst.all.real] + 1
  dat$attr$Opp.Clin.Test.Num[tst.all.clinic] <- dat$attr$Opp.Clin.Test.Num[tst.all.clinic] + 1
  #dat$attr$Opp.HT.Test.Num[tst.all.HT] <- dat$attr$Opp.HT.Test.Num[tst.all.HT] + 1

  
  ### Summary statistics
  dat$epi$opp.clin.test.num[at] <- dat$epi$opp.clin.test.num[at] + length(tst.pos.clinic[!is.na(tst.pos.clinic)]) + length(tst.neg.clinic[!is.na(tst.neg.clinic)])
  dat$epi$opp.home.test.num[at] <- dat$epi$opp.home.test.num[at] + length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.pos.home[!is.na(tst.pos.home)]) + 
                                    length(tst.pos.home.s[!is.na(tst.pos.home.s)]) + length(tst.neg.home.s[!is.na(tst.neg.home.s)])
  dat$epi$opp.clin.test.pos[at] <- dat$epi$opp.clin.test.pos[at] + length(tst.pos.clinic[!is.na(tst.pos.clinic)])
  dat$epi$opp.home.test.pos[at] <- dat$epi$opp.home.test.pos[at] + length(tst.pos.home.comb[!is.na(tst.pos.home.comb)])
  
  
  ##########################################################################################
  ##########################################################################################
  #Clear the vectors with NULL values
  #NULL VALUES
  tst.all <- tst.B <- tst.W <- tst.B.real <- tst.W.real<- tst.all.real <- 
    tst.neg <- tst.neg.clinic <- tst.neg.home <- tst.neg.home.s <- tst.neg.not.supp <- 
    tst.pos <- tst.pos.B <- tst.pos.W <- tst.pos.clinic <- tst.pos.home <- tst.pos.home.s <-
    tst.pos.B.home <- tst.pos.B.clinic <- tst.pos.B.home.s <-
    tst.pos.W.home <- tst.pos.W.clinic <- tst.pos.W.home.s <-
    NO.supp.B.real <- NO.supp.W.real <- ReT.supp.B.real <- ReT.supp.W.real <-
    Risk.supp.B.real <- Risk.supp.W.real <- NN.supp.B.real <- NN.supp.W.real <-
    clinictestIDs.B_Opp <- clinictestIDs.W_Opp <- hometestIDs.B_Opp <- hometestIDs.W_Opp <- tst.pos.home.comb <- NULL 
 
  ############################################################################
   #### For Risk testers, test bases on risk status
###########################################################
  
 # CAI in non-main partnership
  

   
   tst.B <- which(active == 1 & race == "B" & dat$attr$ind.uai.nonmain == 1 &
                    (diag.status == 0 | is.na(diag.status)) &
                    tt == "Risk"& prepStat == 0 & dat$attr$ind.uai.nonmain.clock >= mean.test.risk.CAI.nonmain.B.int
                  & ((dat$attr$Risk.Test.Num < ceiling((at-dat$attr$arrival.time)/52) * max.risk.test.time) | is.na(dat$attr$Risk.Test.Num)))
   
   tst.W <- which(active == 1 & race == "W" & dat$attr$ind.uai.nonmain == 1 &
                    (diag.status == 0 | is.na(diag.status)) &
                    tt == "Risk"& prepStat == 0 & dat$attr$ind.uai.nonmain.clock >= mean.test.risk.CAI.nonmain.W.int
                  & ((dat$attr$Risk.Test.Num < ceiling((at-dat$attr$arrival.time)/52) * max.risk.test.time) | is.na(dat$attr$Risk.Test.Num))) 

   tst.all<-c(tst.B, tst.W)
   
   # Select a sample of those eligible to take a test to really test based on the testing rate mean.test.B.opp.rate
   tst.B.real <- tst.B[runif(length(tst.B)) < mean.test.risk.CAI.nonmain.B.prob]
   tst.W.real <- tst.W[runif(length(tst.W)) < mean.test.risk.CAI.nonmain.W.prob]
   
   #The clinic based opportunistic testers
   clinictestIDs.B_Risk <- tst.B.real
   clinictestIDs.W_Risk <- tst.W.real

   
   ## REPLACEMENT: Risk based test replacement with home tests
   ## If tests are being replaced the testers that are selected are split into two groups (clinic tests and home tests)
   ## The fraction of the tests that are home tests is determined by hometest.rep.risk.frac.
   
   if (Risk.replace==TRUE & at >= HT.start){
     hometestIDs.B_Risk <- tst.B.real[runif(length(tst.B.real)) < hometest.rep.risk.frac]
     clinictestIDs.B_Risk <- setdiff(tst.B.real,hometestIDs.B_Risk)
     
     hometestIDs.W_Risk <- tst.W.real[runif(length(tst.W.real)) < hometest.rep.risk.frac]
     clinictestIDs.W_Risk <- setdiff(tst.W.real,hometestIDs.W_Risk)}
   
   ## SUPPLIMENT: Risk testers may add additional probability of testing with a home test
   if (Risk.supp==TRUE & at >= HT.start){
     #Select those take a supplimental test
     
     Risk.supp.B <- tst.B[runif(length(tst.B)) < supp.risk.home.test.prob]
     Risk.supp.W <- tst.W[runif(length(tst.W)) < supp.risk.home.test.prob]
     
     ##Remove IDS already taking an opp test.
     Risk.supp.B.real <- setdiff(Risk.supp.B,tst.B.real)
     Risk.supp.W.real <- setdiff(Risk.supp.W,tst.W.real)
     
   
   }
   
   
   
   #TAKE the assigned tests  
   
   #clinic test
   
   tst.pos.B.clinic <- clinictestIDs.B_Risk[status[clinictestIDs.B_Risk] == 1 & inf.time[clinictestIDs.B_Risk] <= at - twind.int]
   tst.neg.B.clinic <- setdiff(clinictestIDs.B_Risk, tst.pos.B.clinic)
   
   tst.pos.W.clinic <- clinictestIDs.W_Risk[status[clinictestIDs.W_Risk] == 1 & inf.time[clinictestIDs.W_Risk] <= at - twind.int]
   tst.neg.W.clinic <- setdiff(clinictestIDs.W_Risk, tst.pos.W.clinic)
   
   tst.pos.clinic <- c(tst.pos.B.clinic, tst.pos.W.clinic)
   tst.neg.clinic <- c(tst.neg.B.clinic, tst.neg.W.clinic)
   
   
   if(Risk.replace == TRUE & at >= HT.start){
     #home test
     tst.pos.B.home <- hometestIDs.B_Risk[status[hometestIDs.B_Risk] == 1 & inf.time[hometestIDs.B_Risk] <= at - twind.hometest.int]
     tst.neg.B.home <- setdiff(hometestIDs.B_Risk, tst.pos.B.home)
     
     tst.pos.W.home <- hometestIDs.W_Risk[status[hometestIDs.W_Risk] == 1 & inf.time[hometestIDs.W_Risk] <= at - twind.hometest.int]
     tst.neg.W.home <- setdiff(hometestIDs.W_Risk, tst.pos.W.home)
     
     tst.pos.home <- c(tst.pos.B.home, tst.pos.W.home)
     tst.neg.home <- c(tst.neg.B.home, tst.neg.W.home)
   }
   
   if(Risk.supp == TRUE & at >= HT.start){
     
     tst.pos.B.home.s <- Risk.supp.B.real[status[Risk.supp.B.real] == 1 & inf.time[Risk.supp.B.real] <= at - twind.hometest.int]
     tst.neg.B.home.s <- setdiff(Risk.supp.B.real, tst.pos.B.home.s)
     
     tst.pos.W.home.s <- Risk.supp.W.real[status[Risk.supp.W.real] == 1 & inf.time[Risk.supp.W.real] <= at - twind.hometest.int]
     tst.neg.W.home.s <- setdiff(Risk.supp.W.real, tst.pos.W.home.s)
     
     tst.pos.home.s <- c(tst.pos.B.home.s, tst.pos.W.home.s)
     tst.neg.home.s <- c(tst.neg.B.home.s, tst.neg.W.home.s)
   }
   
   tst.neg <- c(tst.neg.clinic, tst.neg.home,tst.neg.home.s)
   tst.neg.not.supp <- c(tst.neg.clinic, tst.neg.home)
   tst.pos <- c(tst.pos.clinic, tst.pos.home,tst.pos.home.s)
   tst.pos.clinic <- c(tst.pos.clinic)
   tst.pos.home.comb <- c(tst.pos.home,tst.pos.home.s)
   tst.pos.B <- c(tst.pos.B.home, tst.pos.B.clinic,tst.pos.B.home.s)
   tst.pos.W <- c(tst.pos.W.home, tst.pos.W.clinic,tst.pos.W.home.s)
   
   tst.all.real <- c(tst.B.real, tst.W.real, Risk.supp.B.real, Risk.supp.W.real)
   tst.all.Risk.elig <- c(tst.B, tst.W)
   
   tst.all.clinic <- c(tst.pos.clinic, tst.neg.clinic)
   tst.all.HT <- c(tst.pos.home, tst.neg.home, tst.pos.home.s, tst.neg.home.s)
   
   
   ## Output
   
   # Attributes
   dat$attr$diag.risk.clin.test[tst.pos.clinic] <- 1
   dat$attr$diag.risk.HT.test[tst.pos.home.comb] <- 1
   
   dat$attr$diag.status[tst.pos] <- 1
   dat$attr$diag.time[tst.pos] <- at
   dat$attr$last.neg.test.in.clinic[tst.all.clinic] <- at
   
   

   dat$attr$Risk.Test.Num[tst.all.real] <- max(0,dat$attr$Risk.Test.Num[tst.all.real] + 1)
   dat$attr$Risk.Clin.Test.Num[tst.all.clinic] <- max(0,dat$attr$Risk.Clin.Test.Num[tst.all.clinic] + 1)
   dat$attr$Risk.HT.Test.Num[tst.all.HT ] <- max(0,dat$attr$Risk.HT.Test.Num[tst.all.HT] + 1)

   
   
   ### Summary statistics
   dat$epi$risk.clin.test.num[at] <- max(0,length(tst.pos.clinic[!is.na(tst.pos.clinic)]) + length(tst.neg.clinic[!is.na(tst.neg.clinic)]))
   dat$epi$risk.home.test.num[at] <- max(0,length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.neg.home[!is.na(tst.neg.home)]) + 
                                           length(tst.pos.home.s[!is.na(tst.pos.home.s)]) + length(tst.neg.home.s[!is.na(tst.neg.home.s)]))
   dat$epi$risk.clin.test.pos[at] <- max(0,length(tst.pos.clinic[!is.na(tst.pos.clinic)]))
   dat$epi$risk.home.test.pos[at] <- max(0,length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.pos.home.s[!is.na(tst.pos.home.s)]))
   
   
   # Reset Risk based testing attrinutes.
   #All for those that tested, specific risk catagory for those that could have tested but did not.
   
   dat$attr$ind.uai.nonmain[tst.all] <- NA
   dat$attr$ind.uai.nonmain.clock[tst.all] <- NA
   
   dat$attr$ind.uai.known.sd[tst.all.real] <- NA
   dat$attr$ind.uai.known.sd.clock[tst.all.real] <- NA 
   
   dat$attr$ind.newmain[tst.all.real] <- NA
   dat$attr$ind.newmain.clock[tst.all.real] <- NA
   
   
   #Clear the vectors with NULL values
   #NULL VALUES
   tst.all <- tst.B <- tst.W <- tst.B.real <- tst.W.real<- tst.all.real <- 
     tst.neg <- tst.neg.clinic <- tst.neg.home <- tst.neg.home.s <- tst.neg.not.supp <- 
     tst.pos <- tst.pos.B <- tst.pos.W <- tst.pos.clinic <- tst.pos.home <- tst.pos.home.s <-
     tst.pos.B.home <- tst.pos.B.clinic <- tst.pos.B.home.s <-
     tst.pos.W.home <- tst.pos.W.clinic <- tst.pos.W.home.s <-
     NO.supp.B.real <- NO.supp.W.real <- ReT.supp.B.real <- ReT.supp.W.real <-
     Risk.supp.B.real <- Risk.supp.W.real <- NN.supp.B.real <- NN.supp.W.real <-
     clinictestIDs.B_Opp <- clinictestIDs.W_Opp <- hometestIDs.B_Opp <- hometestIDs.W_Opp <- tst.pos.home.comb <- NULL

  ########################################################## 
  ### CAI within known serodiscordant partnership
 
   tst.B <- which(active == 1 & race == "B" & dat$attr$ind.uai.known.sd == 1 &
                    (diag.status == 0 | is.na(diag.status)) &
                    tt == "Risk"& prepStat == 0 & dat$attr$ind.uai.known.sd.clock >= mean.test.risk.AI.known.sd.B.int
                  & ((dat$attr$Risk.Test.Num < ceiling((at-dat$attr$arrival.time)/52) * max.risk.test.time) | is.na(dat$attr$Risk.Test.Num)))
   tst.W <- which(active == 1 & race == "W" & dat$attr$ind.uai.known.sd == 1 &
                    (diag.status == 0 | is.na(diag.status)) &
                    tt == "Risk"& prepStat == 0 & dat$attr$ind.uai.known.sd.clock >= mean.test.risk.AI.known.sd.W.int
                  & ((dat$attr$Risk.Test.Num < ceiling((at-dat$attr$arrival.time)/52) * max.risk.test.time) | is.na(dat$attr$Risk.Test.Num)))
   
   tst.all<-c(tst.B, tst.W)
   
   # Select a sample of those eligible to take a test to really test based on the testing rate mean.test.B.opp.rate
   tst.B.real <- tst.B[runif(length(tst.B)) < mean.test.risk.AI.known.sd.B.prob]
   tst.W.real <- tst.W[runif(length(tst.W)) < mean.test.risk.AI.known.sd.W.prob]
   
   #The clinic based opportunistic testers
   clinictestIDs.B_Risk <- tst.B.real
   clinictestIDs.W_Risk <- tst.W.real
   
   
   ## REPLACEMENT: Risk based test replacement with home tests
   ## If tests are being replaced the testers that are selected are split into two groups (clinic tests and home tests)
   ## The fraction of the tests that are home tests is determined by hometest.rep.risk.frac.
   
   if (Risk.replace==TRUE & at >= HT.start){
     hometestIDs.B_Risk <- tst.B.real[runif(length(tst.B.real)) < hometest.rep.risk.frac]
     clinictestIDs.B_Risk <- setdiff(tst.B.real,hometestIDs.B_Risk)
     
     hometestIDs.W_Risk <- tst.W.real[runif(length(tst.W.real)) < hometest.rep.risk.frac]
     clinictestIDs.W_Risk <- setdiff(tst.W.real,hometestIDs.W_Risk)}
   
   ## SUPPLIMENT: Risk testers may add additional probability of testing with a home test
   if (Risk.supp==TRUE & at >= HT.start){
     #Select those take a supplimental test
     
     Risk.supp.B <- tst.B[runif(length(tst.B)) < supp.risk.home.test.prob]
     Risk.supp.W <- tst.W[runif(length(tst.W)) < supp.risk.home.test.prob]
     
     ##Remove IDS already taking an opp test.
     Risk.supp.B.real <- setdiff(Risk.supp.B,tst.B.real)
     Risk.supp.W.real <- setdiff(Risk.supp.W,tst.W.real)
     
     
   }
   
   
   
   #TAKE the assigned tests  
   
   #clinic test
   
   tst.pos.B.clinic <- clinictestIDs.B_Risk[status[clinictestIDs.B_Risk] == 1 & inf.time[clinictestIDs.B_Risk] <= at - twind.int]
   tst.neg.B.clinic <- setdiff(clinictestIDs.B_Risk, tst.pos.B.clinic)
   
   tst.pos.W.clinic <- clinictestIDs.W_Risk[status[clinictestIDs.W_Risk] == 1 & inf.time[clinictestIDs.W_Risk] <= at - twind.int]
   tst.neg.W.clinic <- setdiff(clinictestIDs.W_Risk, tst.pos.W.clinic)
   
   tst.pos.clinic <- c(tst.pos.B.clinic, tst.pos.W.clinic)
   tst.neg.clinic <- c(tst.neg.B.clinic, tst.neg.W.clinic)
   
   
   if(Risk.replace == TRUE & at >= HT.start){
     #home test
     tst.pos.B.home <- hometestIDs.B_Risk[status[hometestIDs.B_Risk] == 1 & inf.time[hometestIDs.B_Risk] <= at - twind.hometest.int]
     tst.neg.B.home <- setdiff(hometestIDs.B_Risk, tst.pos.B.home)
     
     tst.pos.W.home <- hometestIDs.W_Risk[status[hometestIDs.W_Risk] == 1 & inf.time[hometestIDs.W_Risk] <= at - twind.hometest.int]
     tst.neg.W.home <- setdiff(hometestIDs.W_Risk, tst.pos.W.home)
     
     tst.pos.home <- c(tst.pos.B.home, tst.pos.W.home)
     tst.neg.home <- c(tst.neg.B.home, tst.neg.W.home)
   }
   
   if(Risk.supp == TRUE & at >= HT.start){
     
     tst.pos.B.home.s <- Risk.supp.B.real[status[Risk.supp.B.real] == 1 & inf.time[Risk.supp.B.real] <= at - twind.hometest.int]
     tst.neg.B.home.s <- setdiff(Risk.supp.B.real, tst.pos.B.home.s)
     
     tst.pos.W.home.s <- Risk.supp.W.real[status[Risk.supp.W.real] == 1 & inf.time[Risk.supp.W.real] <= at - twind.hometest.int]
     tst.neg.W.home.s <- setdiff(Risk.supp.W.real, tst.pos.W.home.s)
     
     tst.pos.home.s <- c(tst.pos.B.home.s, tst.pos.W.home.s)
     tst.neg.home.s <- c(tst.neg.B.home.s, tst.neg.W.home.s)
   }
   
   tst.neg <- c(tst.neg.clinic, tst.neg.home,tst.neg.home.s)
   tst.neg.not.supp <- c(tst.neg.clinic, tst.neg.home)
   tst.pos <- c(tst.pos.clinic, tst.pos.home,tst.pos.home.s)
   tst.pos.clinic <- c(tst.pos.clinic)
   tst.pos.home.comb <- c(tst.pos.home,tst.pos.home.s)
   tst.pos.B <- c(tst.pos.B.home, tst.pos.B.clinic,tst.pos.B.home.s)
   tst.pos.W <- c(tst.pos.W.home, tst.pos.W.clinic,tst.pos.W.home.s)
   
   tst.all.real <- c(tst.B.real, tst.W.real, Risk.supp.B.real, Risk.supp.W.real)
   tst.all.Risk.elig <- c(tst.B, tst.W)
   tst.all.clinic <- c(tst.pos.clinic, tst.neg.clinic)
   tst.all.HT <- c(tst.pos.home, tst.neg.home, tst.pos.home.s, tst.neg.home.s)
   
   ## Output
   
   # Attributes
   dat$attr$diag.risk.clin.test[tst.pos.clinic] <- 1
   dat$attr$diag.risk.HT.test[tst.pos.home.comb] <- 1
   
   dat$attr$diag.status[tst.pos] <- 1
   dat$attr$diag.time[tst.pos] <- at
   dat$attr$last.neg.test.in.clinic[tst.all.clinic] <- at
   
   
   dat$attr$Risk.Test.Num[tst.all.real] <- max(0,dat$attr$Risk.Test.Num[tst.all.real] + 1)
   dat$attr$Risk.Clin.Test.Num[tst.all.clinic] <- max(0,dat$attr$Risk.Clin.Test.Num[tst.all.clinic] + 1)
   dat$attr$Risk.HT.Test.Num[tst.all.HT ] <- max(0,dat$attr$Risk.HT.Test.Num[tst.all.HT] + 1)
   
   
   ### Summary statistics
   dat$epi$risk.clin.test.num[at] <- max(0, dat$epi$risk.clin.test.num[at] + (length(tst.pos.clinic[!is.na(tst.pos.clinic)]) + length(tst.neg.clinic[!is.na(tst.neg.clinic)])))
   dat$epi$risk.home.test.num[at] <- max(0, dat$epi$risk.home.test.num[at] + (length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.neg.home[!is.na(tst.neg.home)]) + 
                                     length(tst.pos.home.s[!is.na(tst.pos.home.s)]) + length(tst.neg.home.s[!is.na(tst.neg.home.s)])))
   dat$epi$risk.clin.test.pos[at] <- max(0, dat$epi$risk.clin.test.pos[at] + (length(tst.pos.clinic[!is.na(tst.pos.clinic)])))
   dat$epi$risk.home.test.pos[at] <- max(0, dat$epi$risk.home.test.pos[at] + (length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.pos.home.s[!is.na(tst.pos.home.s)])))
   
   # Reset Risk based testing attrinutes.
   #All for those that tested, specific risk catagory for those that could have tested but did not.
   
   dat$attr$ind.uai.nonmain[tst.all.real] <- NA
   dat$attr$ind.uai.nonmain.clock[tst.all.real] <- NA
   
   dat$attr$ind.uai.known.sd[tst.all] <- NA
   dat$attr$ind.uai.known.sd.clock[tst.all] <- NA 
   
   dat$attr$ind.newmain[tst.all.real] <- NA
   dat$attr$ind.newmain.clock[tst.all.real] <- NA
   
   #Clear the vectors with NULL values
   #NULL VALUES
   tst.all <- tst.B <- tst.W <- tst.B.real <- tst.W.real<- tst.all.real <- 
     tst.neg <- tst.neg.clinic <- tst.neg.home <- tst.neg.home.s <- tst.neg.not.supp <- 
     tst.pos <- tst.pos.B <- tst.pos.W <- tst.pos.clinic <- tst.pos.home <- tst.pos.home.s <-
     tst.pos.B.home <- tst.pos.B.clinic <- tst.pos.B.home.s <-
     tst.pos.W.home <- tst.pos.W.clinic <- tst.pos.W.home.s <-
     NO.supp.B.real <- NO.supp.W.real <- ReT.supp.B.real <- ReT.supp.W.real <-
     Risk.supp.B.real <- Risk.supp.W.real <- NN.supp.B.real <- NN.supp.W.real <-
     clinictestIDs.B_Opp <- clinictestIDs.W_Opp <- hometestIDs.B_Opp <- hometestIDs.W_Opp <- tst.pos.home.comb <- NULL
   
  ######################################################################### 
  ## Acquisition of new main partner 

   tst.B <- which(active == 1 & race == "B" & dat$attr$ind.newmain.sd == 1 &
                    (diag.status == 0 | is.na(diag.status)) &
                    tt == "Risk"& prepStat == 0 & dat$attr$ind.newmain.clock >= mean.test.risk.newmain.B.int
                  & ((dat$attr$Risk.Test.Num < ceiling((at-dat$attr$arrival.time)/52) * max.risk.test.time) | is.na(dat$attr$Risk.Test.Num))) 
   tst.W <- which(active == 1 & race == "W" & dat$attr$ind.newmain == 1 &
                    (diag.status == 0 | is.na(diag.status)) &
                    tt == "Risk"& prepStat == 0 & dat$attr$ind.newmain.clock >= mean.test.risk.newmain.W.int
                  & ((dat$attr$Risk.Test.Num < ceiling((at-dat$attr$arrival.time)/52) * max.risk.test.time) | is.na(dat$attr$Risk.Test.Num)))  
   
   tst.all<-c(tst.B, tst.W)
   
   # Select a sample of those eligible to take a test to really test based on the testing rate mean.test.B.opp.rate
   tst.B.real <- tst.B[runif(length(tst.B)) < mean.test.risk.newmain.B.prob]
   tst.W.real <- tst.W[runif(length(tst.W)) < mean.test.risk.newmain.W.prob]
   
   #The clinic based opportunistic testers
   clinictestIDs.B_Risk <- tst.B.real
   clinictestIDs.W_Risk <- tst.W.real
   
   
   ## REPLACEMENT: Risk based test replacement with home tests
   ## If tests are being replaced the testers that are selected are split into two groups (clinic tests and home tests)
   ## The fraction of the tests that are home tests is determined by hometest.rep.risk.frac.
   
   if (Risk.replace==TRUE & at >= HT.start){
     hometestIDs.B_Risk <- tst.B.real[runif(length(tst.B.real)) < hometest.rep.risk.frac]
     clinictestIDs.B_Risk <- setdiff(tst.B.real,hometestIDs.B_Risk)
     
     hometestIDs.W_Risk <- tst.W.real[runif(length(tst.W.real)) < hometest.rep.risk.frac]
     clinictestIDs.W_Risk <- setdiff(tst.W.real,hometestIDs.W_Risk)}
   
   ## SUPPLIMENT: Risk testers may add additional probability of testing with a home test
   if (Risk.supp==TRUE & at >= HT.start){
     #Select those take a supplimental test
     
     Risk.supp.B <- tst.B[runif(length(tst.B)) < supp.risk.home.test.prob]
     Risk.supp.W <- tst.W[runif(length(tst.W)) < supp.risk.home.test.prob]
     
     ##Remove IDS already taking an opp test.
     Risk.supp.B.real <- setdiff(Risk.supp.B,tst.B.real)
     Risk.supp.W.real <- setdiff(Risk.supp.W,tst.W.real)
     
     
   }
   
   
   
   #TAKE the assigned tests  
   
   #clinic test
   
   tst.pos.B.clinic <- clinictestIDs.B_Risk[status[clinictestIDs.B_Risk] == 1 & inf.time[clinictestIDs.B_Risk] <= at - twind.int]
   tst.neg.B.clinic <- setdiff(clinictestIDs.B_Risk, tst.pos.B.clinic)
   
   tst.pos.W.clinic <- clinictestIDs.W_Risk[status[clinictestIDs.W_Risk] == 1 & inf.time[clinictestIDs.W_Risk] <= at - twind.int]
   tst.neg.W.clinic <- setdiff(clinictestIDs.W_Risk, tst.pos.W.clinic)
   
   tst.pos.clinic <- c(tst.pos.B.clinic, tst.pos.W.clinic)
   tst.neg.clinic <- c(tst.neg.B.clinic, tst.neg.W.clinic)
   
   
   if(Risk.replace == TRUE & at >= HT.start){
     #home test
     tst.pos.B.home <- hometestIDs.B_Risk[status[hometestIDs.B_Risk] == 1 & inf.time[hometestIDs.B_Risk] <= at - twind.hometest.int]
     tst.neg.B.home <- setdiff(hometestIDs.B_Risk, tst.pos.B.home)
     
     tst.pos.W.home <- hometestIDs.W_Risk[status[hometestIDs.W_Risk] == 1 & inf.time[hometestIDs.W_Risk] <= at - twind.hometest.int]
     tst.neg.W.home <- setdiff(hometestIDs.W_Risk, tst.pos.W.home)
     
     tst.pos.home <- c(tst.pos.B.home, tst.pos.W.home)
     tst.neg.home <- c(tst.neg.B.home, tst.neg.W.home)
   }
   
   if(Risk.supp == TRUE & at >= HT.start){
     
     tst.pos.B.home.s <- Risk.supp.B.real[status[Risk.supp.B.real] == 1 & inf.time[Risk.supp.B.real] <= at - twind.hometest.int]
     tst.neg.B.home.s <- setdiff(Risk.supp.B.real, tst.pos.B.home.s)
     
     tst.pos.W.home.s <- Risk.supp.W.real[status[Risk.supp.W.real] == 1 & inf.time[Risk.supp.W.real] <= at - twind.hometest.int]
     tst.neg.W.home.s <- setdiff(Risk.supp.W.real, tst.pos.W.home.s)
     
     tst.pos.home.s <- c(tst.pos.B.home.s, tst.pos.W.home.s)
     tst.neg.home.s <- c(tst.neg.B.home.s, tst.neg.W.home.s)
   }
   
   tst.neg <- c(tst.neg.clinic, tst.neg.home,tst.neg.home.s)
   tst.neg.not.supp <- c(tst.neg.clinic, tst.neg.home)
   tst.pos <- c(tst.pos.clinic, tst.pos.home,tst.pos.home.s)
   tst.pos.clinic <- c(tst.pos.clinic)
   tst.pos.home.comb <- c(tst.pos.home,tst.pos.home.s)
   tst.pos.B <- c(tst.pos.B.home, tst.pos.B.clinic,tst.pos.B.home.s)
   tst.pos.W <- c(tst.pos.W.home, tst.pos.W.clinic,tst.pos.W.home.s)
   
   tst.all.real <- c(tst.B.real, tst.W.real, Risk.supp.B.real, Risk.supp.W.real)
   tst.all.Risk.elig <- c(tst.B, tst.W)
   tst.all.clinic <- c(tst.pos.clinic, tst.neg.clinic)
   tst.all.HT <- c(tst.pos.home, tst.neg.home, tst.pos.home.s, tst.neg.home.s)
   
   
   ## Output
   
   # Attributes
   dat$attr$diag.risk.clin.test[tst.pos.clinic] <- 1
   dat$attr$diag.risk.HT.test[tst.pos.home.comb] <- 1
   
   dat$attr$diag.status[tst.pos] <- 1
   dat$attr$diag.time[tst.pos] <- at
   dat$attr$last.neg.test.in.clinic[tst.all.clinic] <- at
   
   dat$attr$Risk.Test.Num[tst.all.real] <- max(0,dat$attr$Risk.Test.Num[tst.all.real] + 1)
   dat$attr$Risk.Clin.Test.Num[tst.all.clinic] <- max(0,dat$attr$Risk.Clin.Test.Num[tst.all.clinic] + 1)
   dat$attr$Risk.HT.Test.Num[tst.all.HT ] <- max(0,dat$attr$Risk.HT.Test.Num[tst.all.HT] + 1)
   
   
   ### Summary statistics
   dat$epi$risk.clin.test.num[at] <- max(0, dat$epi$risk.clin.test.num[at] + (length(tst.pos.clinic[!is.na(tst.pos.clinic)]) + length(tst.neg.clinic[!is.na(tst.neg.clinic)])))
   dat$epi$risk.home.test.num[at] <- max(0, dat$epi$risk.home.test.num[at] + (length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.neg.home[!is.na(tst.neg.home)]) + 
                                     length(tst.pos.home.s[!is.na(tst.pos.home.s)]) + length(tst.neg.home.s[!is.na(tst.neg.home.s)])))
   dat$epi$risk.clin.test.pos[at] <- max(0, dat$epi$risk.clin.test.pos[at] + (length(tst.pos.clinic[!is.na(tst.pos.clinic)])))
   dat$epi$risk.home.test.pos[at] <- max(0, dat$epi$risk.home.test.pos[at] + (length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.pos.home.s[!is.na(tst.pos.home.s)])))
   
   # Reset Risk based testing attrinutes.
   #All for those that tested, specific risk catagory for those that could have tested but did not.
   
   dat$attr$ind.uai.nonmain[tst.all.real] <- NA
   dat$attr$ind.uai.nonmain.clock[tst.all.real] <- NA
   
   dat$attr$ind.uai.known.sd[tst.all.real] <- NA
   dat$attr$ind.uai.known.sd.clock[tst.all.real] <- NA 
   
   dat$attr$ind.newmain[tst.all] <- NA
   dat$attr$ind.newmain.clock[tst.all] <- NA
   
   ##############################################################################################################################
   ##############################################################################################################################
   ######  Opportunistic tests for Risk based testers.
   
   ##Opportunuistic tests for regular testers
   
   #Clear the vectors with NULL values
   #NULL VALUES
   tst.all <- tst.B <- tst.W <- tst.B.real <- tst.W.real<- tst.all.real <- 
     tst.neg <- tst.neg.clinic <- tst.neg.home <- tst.neg.home.s <- tst.neg.not.supp <- 
     tst.pos <- tst.pos.B <- tst.pos.W <- tst.pos.clinic <- tst.pos.home <- tst.pos.home.s <-
     tst.pos.B.home <- tst.pos.B.clinic <- tst.pos.B.home.s <-
     tst.pos.W.home <- tst.pos.W.clinic <- tst.pos.W.home.s <-
     NO.supp.B.real <- NO.supp.W.real <- ReT.supp.B.real <- ReT.supp.W.real <-
     Risk.supp.B.real <- Risk.supp.W.real <- NN.supp.B.real <- NN.supp.W.real <-
     clinictestIDs.B_Opp <- clinictestIDs.W_Opp <- hometestIDs.B_Opp <- hometestIDs.W_Opp <- tst.pos.home.comb <- NULL 

   
   # Risk testers only "Risk" 
   # Frequency of testing is based on time since last opportnistic test being greater than 
   # a fixed interval for opportunistic tests for risk testers (mean.test.opp.Risk.R.int)
   
   #Get time since last test
   tsincelntst <- at - dat$attr$last.neg.test.opp
   tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]
   
   #Select those eligible to test based on time since last test
   # if the time since their last test is larger than the mean test interval, then there is a chance they test.
   tst.B <- which(active == 1 & race == "B" & tt == "Risk" &
                    (diag.status == 0 | is.na(diag.status)) &
                    tsincelntst >= mean.test.opp.Risk.B.int
                  & prepStat == 0)  
   tst.W <- which(active == 1 & race == "W" & tt== "Risk" &
                    (diag.status == 0 | is.na(diag.status)) &
                    tsincelntst >= mean.test.opp.Risk.W.int
                  & prepStat == 0)
   
   # Select a sample of those eligible to take a test to really test based on the testing rate mean.test.B.opp.rate
   tst.B.real <- tst.B[runif(length(tst.B)) < mean.test.B.opp.reg.prob]
   tst.W.real <- tst.W[runif(length(tst.W)) < mean.test.W.opp.reg.prob]
   
   #The clinic based opportunistic testers
   clinictestIDs.B_Opp <- tst.B.real
   clinictestIDs.W_Opp <- tst.W.real
   
   
   ## REPLACEMENT: Opportunistic test replacement with home tests
   ## If tests are being replaced the testers that are selected are split into two groups (clinic tests and home tests)
   ## The fraction of the tests that are home tests is determined by hometest.rep.frac.
   
   if (Opportunity.replace == TRUE & at >= HT.start){
     hometestIDs.B_Opp <- tst.B.real[runif(length(tst.B.real)) < hometest.rep.opp.frac]
     clinictestIDs.B_Opp <- setdiff(tst.B.real,hometestIDs.B_Opp)
     
     hometestIDs.W_Opp <- tst.W.real[runif(length(tst.W.real)) < hometest.rep.opp.frac]
     clinictestIDs.W_Opp <- setdiff(tst.W.real,hometestIDs.W_Opp)}
   
   ## SUPPLIMENT: Opportunistic only testers take additional home tests
   if (Opportunity.supp == TRUE & at >= HT.start){
     #Select those eligible for a supplimental test
     Risk.supp.B <- which(active == 1 & race == "B" & tt == "Risk" &
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)  
     Risk.supp.W <- which(active == 1 & race == "W" & tt== "Risk" &
                           (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
     
     ##Remove IDS already taking an opp test.
     Risk.supp.B <- setdiff(Risk.supp.B,tst.B.real)
     Risk.supp.W <- setdiff(Risk.supp.W,tst.W.real)
     
     #Select those that will take a supplimental home test based on probability supp.home.test.prob
     Risk.supp.B.real <- Risk.supp.B[runif(length(Risk.supp.B)) < supp.opp.home.test.prob]
     Risk.supp.W.real <- Risk.supp.W[runif(length(Risk.supp.W)) < supp.opp.home.test.prob]    
   }
   
   
   ##Take the tests
   #clinic test
   
   tst.pos.B.clinic <- clinictestIDs.B_Opp[status[clinictestIDs.B_Opp] == 1 & inf.time[clinictestIDs.B_Opp] <= at - twind.int]
   tst.neg.B.clinic <- setdiff(clinictestIDs.B_Opp, tst.pos.B.clinic)
   
   tst.pos.W.clinic <- clinictestIDs.W_Opp[status[clinictestIDs.W_Opp] == 1 & inf.time[clinictestIDs.W_Opp] <= at - twind.int]
   tst.neg.W.clinic <- setdiff(clinictestIDs.W_Opp, tst.pos.W.clinic)
   
   tst.pos.clinic <- c(tst.pos.B.clinic, tst.pos.W.clinic)
   tst.neg.clinic <- c(tst.neg.B.clinic, tst.neg.W.clinic)
   
   
   ##Home tests: replacements for opportunistic tests
   if (Opportunity.replace == TRUE & at >= HT.start){
     #home test
     tst.pos.B.home <- hometestIDs.B_Opp[status[hometestIDs.B_Opp] == 1 & inf.time[hometestIDs.B_Opp] <= at - twind.hometest.int]
     tst.neg.B.home <- setdiff(hometestIDs.B_Opp, tst.pos.B.home)
     
     tst.pos.W.home <- hometestIDs.W_Opp[status[hometestIDs.W_Opp] == 1 & inf.time[hometestIDs.W_Opp] <= at - twind.hometest.int]
     tst.neg.W.home <- setdiff(hometestIDs.W_Opp, tst.pos.W.home)
     
     tst.pos.home <- c(tst.pos.B.home, tst.pos.W.home)
     tst.neg.home <- c(tst.neg.B.home, tst.neg.W.home)
     
   }
   
   ##Home tests: supplimental tests
   if (Opportunity.supp == TRUE & at >= HT.start){
     #home test
     tst.pos.B.home.s <- Risk.supp.B.real[status[Risk.supp.B.real] == 1 & inf.time[Risk.supp.B.real] <= at - twind.hometest.int]
     tst.neg.B.home.s <- setdiff(Risk.supp.B.real, tst.pos.B.home.s)
     
     tst.pos.W.home.s <- Risk.supp.W.real[status[Risk.supp.W.real] == 1 & inf.time[Risk.supp.W.real] <= at - twind.hometest.int]
     tst.neg.W.home.s <- setdiff(Risk.supp.W.real, tst.pos.W.home.s)
     
     tst.pos.home.s <- c(tst.pos.B.home.s, tst.pos.W.home.s)
     tst.neg.home.s <- c(tst.neg.B.home.s, tst.neg.W.home.s)
     
   }
   
   
   tst.neg <- c(tst.neg.clinic, tst.neg.home,tst.neg.home.s)
   tst.neg.not.supp <- c(tst.neg.clinic, tst.neg.home)
   tst.pos <- c(tst.pos.clinic, tst.pos.home,tst.pos.home.s)
   tst.pos.clinic <- c(tst.pos.clinic)
   tst.pos.home.comb <- c(tst.pos.home,tst.pos.home.s)
   tst.pos.B <- c(tst.pos.B.home, tst.pos.B.clinic,tst.pos.B.home.s)
   tst.pos.W <- c(tst.pos.W.home, tst.pos.W.clinic,tst.pos.W.home.s)
   
   tst.all.real <- c(tst.B.real, tst.W.real, Risk.supp.B.real, Risk.supp.W.real)
   tst.all.opp.elig <- c(tst.B, tst.W)
   tst.all.opp.elig <- setdiff(tst.all.opp.elig,tst.pos)
   tst.all.clinic <- c(tst.pos.clinic, tst.neg.clinic)
   tst.all.HT <- c(tst.pos.home, tst.neg.home, tst.pos.home.s, tst.neg.home.s)
   
   
   # Attributes
   dat$attr$diag.opp.clin.test[tst.pos.clinic] <- 1
   dat$attr$diag.opp.HT.test[tst.pos.home.comb] <- 1
   dat$attr$last.neg.test.opp[tst.all.opp.elig] <- at
   dat$attr$diag.status[tst.pos] <- 1
   dat$attr$diag.time[tst.pos] <- at
   
   dat$attr$Opp.Test.Num[tst.all.real] <- dat$attr$Opp.Test.Num[tst.all.real] + 1
   dat$attr$Opp.Clin.Test.Num[tst.all.clinic] <- dat$attr$Opp.Clin.Test.Num[tst.all.clinic] + 1
   #dat$attr$Opp.HT.Test.Num[tst.all.HT] <- dat$attr$Opp.HT.Test.Num[tst.all.HT] + 1

   
   ### Summary statistics
   dat$epi$opp.clin.test.num[at] <- dat$epi$opp.clin.test.num[at] + length(tst.pos.clinic[!is.na(tst.pos.clinic)]) + length(tst.neg.clinic[!is.na(tst.neg.clinic)])
   dat$epi$opp.home.test.num[at] <- dat$epi$opp.home.test.num[at] + length(tst.pos.home[!is.na(tst.pos.home)]) + length(tst.neg.home[!is.na(tst.neg.home)]) + 
                                    length(tst.pos.home.s[!is.na(tst.pos.home.s)]) + length(tst.neg.home.s[!is.na(tst.neg.home.s)])
   dat$epi$opp.clin.test.pos[at] <- dat$epi$opp.clin.test.pos[at] + length(tst.pos.clinic[!is.na(tst.pos.clinic)])
   dat$epi$opp.home.test.pos[at] <- dat$epi$opp.home.test.pos[at] + length(tst.pos.home.comb[!is.na(tst.pos.home.comb)])
   
   
   ##############################################################
   #########  NEVER TESTERS
   
  
   if (Never.test.supp == TRUE & at >= HT.start){
     
     
     #Clear the vectors with NULL values
     #NULL VALUES
     tst.all <- tst.B <- tst.W <- tst.B.real <- tst.W.real<- tst.all.real <- 
       tst.neg <- tst.neg.clinic <- tst.neg.home <- tst.neg.home.s <- tst.neg.not.supp <- 
       tst.pos <- tst.pos.B <- tst.pos.W <- tst.pos.clinic <- tst.pos.home <- tst.pos.home.s <-
       tst.pos.B.home <- tst.pos.B.clinic <- tst.pos.B.home.s <-
       tst.pos.W.home <- tst.pos.W.clinic <- tst.pos.W.home.s <-
       NO.supp.B.real <- NO.supp.W.real <- ReT.supp.B.real <- ReT.supp.W.real <-
       Risk.supp.B.real <- Risk.supp.W.real <- NN.supp.B.real <- NN.supp.W.real <-
       clinictestIDs.B_Opp <- clinictestIDs.W_Opp <- hometestIDs.B_Opp <- hometestIDs.W_Opp <- tst.pos.home.comb <- NULL 
     
     
     # Neveer Testers TESTERS "NN" 
     # Frequency of testing is based on time since last opportnistic test being greater than 
     # a fixed interval for opportunistic tests (mean.test.opp.NO.R.int)
     

     #Select those eligible to test based on time since last test
     # if the time since their last test is larger than the mean test interval, then there is a chance they test.
     tst.B <- which(active == 1 & race == "B" & tt == "NN" &
                      (diag.status == 0 | is.na(diag.status)) & prepStat == 0)  
     tst.W <- which(active == 1 & race == "W" & tt== "NN" &
                      (diag.status == 0 | is.na(diag.status)) & prepStat == 0)
     
     # Select a sample of those eligible to take a test to really test based on the testing rate mean.test.B.opp.rate
     NN.supp.B.real <- tst.B[runif(length(tst.B)) < supp.nevertest.home.test.prob]
     NN.supp.W.real <- tst.W[runif(length(tst.W)) < supp.nevertest.home.test.prob]
     

       #Take the home test
       tst.pos.B.home.s <- NN.supp.B.real[status[NN.supp.B.real] == 1 & inf.time[NN.supp.B.real] <= at - twind.hometest.int]
       tst.neg.B.home.s <- setdiff(NN.supp.B.real, tst.pos.B.home.s)
       
       tst.pos.W.home.s <- NN.supp.W.real[status[NN.supp.W.real] == 1 & inf.time[NN.supp.W.real] <= at - twind.hometest.int]
       tst.neg.W.home.s <- setdiff(NN.supp.W.real, tst.pos.W.home.s)
       
       tst.pos.home.s <- c(tst.pos.B.home.s, tst.pos.W.home.s)
       tst.neg.home.s <- c(tst.neg.B.home.s, tst.neg.W.home.s)
       
   
     
     tst.neg <- c(tst.neg.home.s)
     tst.pos <- c(tst.pos.home.s)
     tst.pos.home <- c(tst.pos.home.s)
     tst.pos.B <- c(tst.pos.B.home.s)
     tst.pos.W <- c(tst.pos.W.home.s)
     
     tst.all.real <- c(NN.supp.B.real, NN.supp.W.real)
     tst.all.HT <- c(tst.pos.home.s, tst.neg.home.s)

     # Attributes
     dat$attr$diag.NN.HT.test[tst.pos.home] <- 1
     dat$attr$diag.status[tst.pos] <- 1
     dat$attr$diag.time[tst.pos] <- at
     
     
     ### Summary statistics
     dat$epi$NN.home.test.num[at] <- max(0,length(tst.pos.home.s[!is.na(tst.pos.home.s)]) + length(tst.neg.home.s[!is.na(tst.neg.home.s)]))
     dat$epi$NN.home.test.pos[at] <- max(0,length(tst.pos.home.s[!is.na(tst.pos.home.s)]))
     
     ###########ADD THESE COUNTS
     dat$attr$NN.HT.Test.Num[tst.all.HT ] <- dat$attr$NN.HT.Test.Num[tst.all.HT] + 1
}
     
   
  return(dat)
}




