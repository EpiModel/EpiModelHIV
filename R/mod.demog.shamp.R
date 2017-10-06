
#' @title Demography Check Module
#'
#' @description Module updated demographic clasifications that change over time.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
demogupdate_shamp <- function(dat, at) {

  #AGECAT
  dat$attr$agecat<-ifelse(dat$att$age < 26 ,"18-25",
                           ifelse(dat$attr$age > 25 & dat$attr$age < 36,"26-35",
                                  ifelse(dat$attr$age > 35 & dat$attr$age < 46,"36-45",
                                         ifelse(dat$attr$age > 45,"46-60",dat$attr$agecat))))
  

  #Create a new factor for age by sex group by pers.c.
  
  dat$attr$race.sex.pers<-rep(NA,length(dat$attr$age)) 
  dat$attr$race.sex.pers<-ifelse(dat$attr$race=="B" & dat$attr$sex=="M" & dat$attr$deg.pers.c==0,"Ref.0",
                                  ifelse(dat$attr$race=="BI" & dat$attr$sex=="M" & dat$attr$deg.pers.c==0,"Ref.0",
                                  ifelse(dat$attr$race=="H" & dat$attr$sex=="M" & dat$attr$deg.pers.c==0,"Ref.0",
                                  ifelse(dat$attr$race=="HI" & dat$attr$sex=="M" & dat$attr$deg.pers.c==0,"Ref.0",
                                  ifelse(dat$attr$race=="W" & dat$attr$sex=="M" & dat$attr$deg.pers.c==0,"Ref.0",
                                  ifelse(dat$attr$race=="B" & dat$attr$sex=="F" & dat$attr$deg.pers.c==0,"Ref.0",
                                  ifelse(dat$attr$race=="BI" & dat$attr$sex=="F" & dat$attr$deg.pers.c==0,"Ref.0",
                                  ifelse(dat$attr$race=="H" & dat$attr$sex=="F" & dat$attr$deg.pers.c==0,"Ref.0",
                                  ifelse(dat$attr$race=="HI" & dat$attr$sex=="F" & dat$attr$deg.pers.c==0,"Ref.0",
                                  ifelse(dat$attr$race=="W" & dat$attr$sex=="F" & dat$attr$deg.pers.c==0,"Ref.0",
                                  ifelse(dat$attr$race=="B" & dat$attr$sex=="M" & dat$attr$deg.pers.c==1,"B.M.p1",
                                  ifelse(dat$attr$race=="BI" & dat$attr$sex=="M" & dat$attr$deg.pers.c==1,"BI.M.p1",
                                  ifelse(dat$attr$race=="H" & dat$attr$sex=="M" & dat$attr$deg.pers.c==1,"H.M.p1",
                                  ifelse(dat$attr$race=="HI" & dat$attr$sex=="M" & dat$attr$deg.pers.c==1,"HI.M.p1",
                                  ifelse(dat$attr$race=="W" & dat$attr$sex=="M" & dat$attr$deg.pers.c==1,"W.M.p1",
                                  ifelse(dat$attr$race=="B" & dat$attr$sex=="F" & dat$attr$deg.pers.c==1,"B.F.p1",
                                  ifelse(dat$attr$race=="BI" & dat$attr$sex=="F" & dat$attr$deg.pers.c==1,"BI.F.p1",
                                  ifelse(dat$attr$race=="H" & dat$attr$sex=="F" & dat$attr$deg.pers.c==1,"H.F.p1",
                                  ifelse(dat$attr$race=="HI" & dat$attr$sex=="F" & dat$attr$deg.pers.c==1,"HI.F.p1",
                                  ifelse(dat$attr$race=="W" & dat$attr$sex=="F" & dat$attr$deg.pers.c==1,"W.F.p1",
                                  dat$attr$race.sex.pers))))))))))))))))))))
  
  
  
  #Create a new factor for age by sex group by main.c 
  dat$attr$race.sex.cohab<-rep(NA,length(dat$attr$age))
  dat$attr$race.sex.cohab<-ifelse(dat$attr$race=="B" & dat$attr$sex=="M" & dat$attr$deg.cohab.c==0,"Ref.0",
                                   ifelse(dat$attr$race=="BI" & dat$attr$sex=="M" & dat$attr$deg.cohab.c==0,"Ref.0",
                                   ifelse(dat$attr$race=="H" & dat$attr$sex=="M" & dat$attr$deg.cohab.c==0,"Ref.0",
                                   ifelse(dat$attr$race=="HI" & dat$attr$sex=="M" & dat$attr$deg.cohab.c==0,"Ref.0",
                                   ifelse(dat$attr$race=="W" & dat$attr$sex=="M" & dat$attr$deg.cohab.c==0,"Ref.0",
                                   ifelse(dat$attr$race=="B" & dat$attr$sex=="F" & dat$attr$deg.cohab.c==0,"Ref.0",
                                   ifelse(dat$attr$race=="BI" & dat$attr$sex=="F" & dat$attr$deg.cohab.c==0,"Ref.0",
                                   ifelse(dat$attr$race=="H" & dat$attr$sex=="F" & dat$attr$deg.cohab.c==0,"Ref.0",
                                   ifelse(dat$attr$race=="HI" & dat$attr$sex=="F" & dat$attr$deg.cohab.c==0,"Ref.0",
                                   ifelse(dat$attr$race=="W" & dat$attr$sex=="F" & dat$attr$deg.cohab.c==0,"Ref.0",
                                   ifelse(dat$attr$race=="B" & dat$attr$sex=="M" & dat$attr$deg.cohab.c==1,"B.M.p1",
                                   ifelse(dat$attr$race=="BI" & dat$attr$sex=="M" & dat$attr$deg.cohab.c==1,"BI.M.p1",
                                   ifelse(dat$attr$race=="H" & dat$attr$sex=="M" & dat$attr$deg.cohab.c==1,"H.M.p1",
                                   ifelse(dat$attr$race=="HI" & dat$attr$sex=="M" & dat$attr$deg.cohab.c==1,"HI.M.p1",
                                   ifelse(dat$attr$race=="W" & dat$attr$sex=="M" & dat$attr$deg.cohab.c==1,"W.M.p1",
                                   ifelse(dat$attr$race=="B" & dat$attr$sex=="F" & dat$attr$deg.cohab.c==1,"B.F.p1",
                                   ifelse(dat$attr$race=="BI" & dat$attr$sex=="F" & dat$attr$deg.cohab.c==1,"BI.F.p1",
                                   ifelse(dat$attr$race=="H" & dat$attr$sex=="F" & dat$attr$deg.cohab.c==1,"H.F.p1",
                                   ifelse(dat$attr$race=="HI" & dat$attr$sex=="F" & dat$attr$deg.cohab.c==1,"HI.F.p1",
                                   ifelse(dat$attr$race=="W" & dat$attr$sex=="F" & dat$attr$deg.cohab.c==1,"W.F.p1",
                                   dat$attr$race.sex.cohab))))))))))))))))))))
  
  
  
  #  Recalculate the demog.cat for new demo.
  
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
  
  
  
  
  return(dat)
}
