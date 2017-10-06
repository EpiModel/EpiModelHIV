
# SHAMP  -----------------------------------------------------------------

#' @title Import and Check Data
#'
#' @description Imports the data for simulations using ergm.ego and SHAMP function.  Ego and alter data are checked to verify that all 
#'              of the required fields are present.  Egos and alters with missing data are dropped. 
#'
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{data.params}, which can be passed to
#' EpiModel function \code{params}
#' EpiModel function \code{setup}
#' .
#'
#' @keywords shamp
#'
#' @export
input_shamp <- function(data, data.params, immigration, msm.msmf) {
 

data.params$rates_c <-NULL 
data.params$durs_c <-NULL
data.params$rates_p <-NULL 
data.params$durs_p <-NULL
data.params$msm.frac <-NULL
data.params$msmf.frac <-NULL
data.params$age.adj<-NULL
  
##Check to make sure there is data for egos and all three alters in the expected data.frames.
names(data) 
x<- names(data)  
  if (any(x=="egos")==FALSE){
    stop("data list must include observed egos data.frame 'egos'",
         call. = FALSE)}

if (any(x=="altersCohab")==FALSE){
  stop("data must include alters  data.frame 'altersCohab'",
       call. = FALSE)}

if (any(x=="altersPers")==FALSE){
  stop("data must include alters  data.frame 'altersPers'",
       call. = FALSE)}

if (any(x=="altersOT")==FALSE){
  stop("data must include alters  data.frame 'altersOT'",
       call. = FALSE)}

##Check to make sure we have all of the expected attributes for egos.

fields<-c("weight","ego","sex","age","sqrt.age","sex.ident","immigrant","role.class","race","deg.cohab", "deg.pers")
fields.egos<-names(data$egos)

for (i in 1:length(fields)){
  if(any(fields.egos==fields[i]) == FALSE){
  stop("missing required field in egos",
       call. = FALSE)}}


  ##Check to make sure we have all of the expected attributes for altersCohab, altersPers, and alterOT
  
  fields<-c("ego","len")

  fields.altersCohab<-names(data$altersCohab)
  for (i in 1:length(fields)){
    if(any(fields.altersCohab==fields[i]) == FALSE){
      stop("missing required field in altersCohab",
           call. = FALSE)}}
    
  fields.altersPers<-names(data$altersPers)
  for (i in 1:length(fields)){
    if(any(fields.altersPers==fields[i]) == FALSE){
      stop("missing required field in altersPers",
           call. = FALSE)}}
    
  fields.altersOT<-names(data$altersOT)
  for (i in 1:length(fields)){
    if(any(fields.altersOT==fields[i]) == FALSE){
      stop("missing required field in altersOT'",
           call. = FALSE)}}
    
  
##Drop egos and alters with missing data.
  
##egos with missing data.
  
  fields<-c("weight","ego","sex","age","sqrt.age","sex.ident","immigrant","role.class","race","deg.cohab", "deg.pers")
  
  ids.list<-NULL
  for(i in 1:length(fields)){
  temp<-is.na(data$egos[i])
  ids<-which(temp==TRUE)
  ego.ids<-data$egos$ego[ids]
  ids.list<-c(ids.list,ego.ids)}
  
  if(length(ego.ids>0)){warning("Egos and associated alters deleted: missing required field",call. = FALSE)
    print(length(ids.list))}

  data$egos<-data$egos[data$egos$ego %in% ids.list == FALSE,]
  data$altersCohab<-data$altersCohab[data$altersCohab$ego %in% ids.list == FALSE,]
  data$altersPers<-data$altersPers[data$altersPers$ego %in% ids.list == FALSE,]
  data$altersOT<-data$altersOT[data$altersOT$ego %in% ids.list == FALSE,]

##Alters with missing data.
##Main
fields<-names(data$altersCohab)

alter.list<-NULL
for(i in 1:length(fields)){
  temp<-is.na(data$altersCohab[i])
  ids<-which(temp==TRUE)
  alter.ids<-data$altersCohab$ego[ids]
  alter.list<-c(alter.list,alter.ids)}

if(length(alter.list>0)){warning("Cohab alters deleted: missing values",call. = FALSE)
    print(length(alter.list))}
    data$altersCohab<-data$altersCohab[data$altersCohab$ego %in% alter.list == FALSE,]


##Casual.

fields<-names(data$altersPers)

alter.list<-NULL
for(i in 1:length(fields)){
  temp<-is.na(data$altersPers[i])
  ids<-which(temp==TRUE)
  alter.ids<-data$altersPers$ego[ids]
  alter.list<-c(alter.list,alter.ids)}

if(length(alter.list>0)){warning("Casual alters deleted: missing values",call. = FALSE)
   print(length(alter.list))}
    data$altersPers<-data$altersPers[data$altersPers$ego %in% alter.list == FALSE,]



##one time.

##If not imputing race of missing alters move BI and HI assignment for egos to after the missing data checks.
   
##If using immigration set Black and Hispanic immigrants to BI and HI.
    if(immigration==TRUE){
      data$egos$race<-ifelse(data$egos$race=="B" & data$egos$immigrant=="Yes","BI",
                             ifelse(data$egos$race=="H" & data$egos$immigrant=="Yes","HI",
                                    data$egos$race))
      
      data$altersCohab$race<-ifelse(data$altersCohab$race=="B" & data$altersCohab$immigrant=="Yes","BI",
                                    ifelse(data$altersCohab$race=="H" & data$altersCohab$immigrant=="Yes","HI",
                                           data$altersCohab$race))
      
      data$altersPers$race<-ifelse(data$altersPers$race=="B" & data$altersPers$immigrant=="Yes","BI",
                                   ifelse(data$altersPers$race=="H" & data$altersPers$immigrant=="Yes","HI",
                                          data$altersPers$race))
      
      data$altersOT$race<-ifelse(data$altersOT$race=="B" & data$altersOT$immigrant=="Yes","BI",
                                 ifelse(data$altersOT$race=="H" & data$altersOT$immigrant=="Yes","HI",
                                        data$altersOT$race))
    }
    
    
    {warning("OT alters race imputed)",call. = FALSE)}

mis.r <- which(is.na(data$altersOT$race)==TRUE)
for (i in 1:length(mis.r)){
  ego<- data$altersOT$ego[mis.r[i]]
  ego.race<-data$egos$race[data$egos$ego==ego]
  if (ego.race=="B") {imp.race<-sample(c("H","W","B"),1,prob=c(0.073,	0.152, 	0.775))}
  if (ego.race=="BI") {imp.race<-sample(c("H","W","B","BI"),1,prob=c(0.073,	0.152,0.3875,0.3875))}
  if (ego.race=="H") {imp.race<-sample(c("H","W","B"),1,prob=c(0.714,	0.230,	0.056))}
  if (ego.race=="HI") {imp.race<-sample(c("H","HI", "W","B"),1,prob=c(0.357,0.357,	0.230,	0.056))}
  if (ego.race=="W") {imp.race<-sample(c("H","W","B"),1,prob=c(0.097,	0.888,	0.016))}
  data$altersOT$race[mis.r[i]]<-imp.race
  
}

fields<-names(data$altersOT)

alter.list<-NULL

for(i in 1:length(fields)){
  temp<-is.na(data$altersOT[i])
  ids<-which(temp==TRUE)
  alter.ids<-data$altersOT$ego[ids]
  alter.list<-c(alter.list,alter.ids)}


if(length(alter.list>0)){warning("OT alters deleted: missing values",call. = FALSE)
  print(length(alter.list))}


    data$altersOT<-data$altersOT[data$altersOT$ego %in% alter.list == FALSE,]
    
##Data checks completed############################################

###  Calculate the paramters needed for estimation that are not line data.
    
#Age, sqrt.age and age.adj.
 
    data$altersCohab$age <- ifelse(data$altersCohab$age > 60, 60, data$altersCohab$age)
    data$altersCohab$age <- ifelse(data$altersCohab$age < 18, 18, data$altersCohab$age)
    data$altersCohab$sqrt.age <- sqrt(data$altersCohab$age)
    
    data$altersPers$age <- ifelse(data$altersPers$age > 60, 60, data$altersPers$age)
    data$altersPers$age <- ifelse(data$altersPers$age < 18, 18, data$altersPers$age)
    data$altersPers$sqrt.age <- sqrt(data$altersPers$age)
    
    data$altersOT$age <- ifelse(data$altersOT$age > 60, 60, data$altersOT$age)
    data$altersOT$age <- ifelse(data$altersOT$age < 18, 18, data$altersOT$age)
    data$altersOT$sqrt.age <- sqrt(data$altersOT$age)
    
    ids.me<-sample(data$egos$ego[data$egos$sex=="M"],100000,replace=TRUE,prob=data$egos$weight[data$egos$sex=="M"])
    x<-which((ids.me %in% data$altersCohab$ego) == TRUE)
    ids.me<-ids.me[x]      
    
    sqrt.age.me<-rep(NA,length(ids.me))
    sqrt.age.fa<-rep(NA,length(ids.me))
  
    for(i in 1:length(ids.me)){
      sqrt.age.me[i] <- data$egos$sqrt.age[data$egos$ego==ids.me[i]]
      sqrt.age.fa[i] <- sample(c(data$altersCohab$sqrt.age[data$altersCohab$ego == ids.me[i]],
                                 data$altersCohab$sqrt.age[data$altersCohab$ego == ids.me[i]]),size=1)
    }
  
    

    ids.fe<-sample(data$egos$ego[data$egos$sex=="F"],100000,replace=TRUE,prob=data$egos$weight[data$egos$sex=="F"])
    x<-which((ids.fe %in% data$altersCohab$ego) == TRUE)
    ids.fe<-ids.fe[x]      
    
    sqrt.age.fe<-rep(NA,length(ids.fe))
    sqrt.age.ma<-rep(NA,length(ids.fe))
    
    for(i in 1:length(ids.fe)){
      sqrt.age.fe[i] <- data$egos$sqrt.age[data$egos$ego==ids.fe[i]]
      sqrt.age.ma[i] <- sample(c(data$altersCohab$sqrt.age[data$altersCohab$ego == ids.fe[i]],
                                 data$altersCohab$sqrt.age[data$altersCohab$ego == ids.fe[i]]),size=1)
    }
    
    

    agediff.1<-mean(sqrt.age.me)-mean(sqrt.age.fa)
    agediff.2<-mean(sqrt.age.ma)-mean(sqrt.age.fe)
    age.adj <- mean(agediff.1,agediff.2)
    data.params$age.adj<-age.adj
    
    
    
    #sqrt.age.adj
    data$egos$sqrt.age.adj<-ifelse(data$egos$sex=="M",data$egos$sqrt.age,
                                  ifelse(data$egos$sex=="F",data$egos$sqrt.age + age.adj, data$egos$sqrt.age))
    data$altersCohab$sqrt.age.adj<-ifelse(data$altersCohab$sex=="M",data$altersCohab$sqrt.age,
                                   ifelse(data$altersCohab$sex=="F",data$altersCohab$sqrt.age  + age.adj, data$altersCohab$sqrt.age))
    data$altersPers$sqrt.age.adj<-ifelse(data$altersPers$sex=="M",data$altersPers$sqrt.age,
                                   ifelse(data$altersPers$sex=="F",data$altersPers$sqrt.age  + age.adj, data$altersPers$sqrt.age))
    data$altersOT$sqrt.age.adj<-ifelse(data$altersOT$sex=="M",data$altersOT$sqrt.age,
                                   ifelse(data$altersOT$sex=="F",data$altersOT$sqrt.age  + age.adj, data$altersOT$sqrt.age))
  
  #Limit Cohab and persistent to a 2 by 2 set (Cohab 0:1, Pers 0:1)  
  ##Create deg.pers.c and deg.cohab.c for capping main at 1 and pers at 1 for model fit.
  data$egos$deg.pers<-as.numeric(data$egos$deg.pers)
  data$egos$deg.pers.c<-as.numeric(data$egos$deg.pers)
  data$egos$deg.pers.c<-ifelse(data$egos$deg.pers.c > 0, 1,data$egos$deg.pers.c)
  data$egos$deg.cohab<-as.numeric(data$egos$deg.cohab)
  data$egos$deg.cohab.c<-data$egos$deg.cohab
  data$egos$deg.cohab.c<-ifelse(data$egos$deg.cohab.c > 0, 1,data$egos$deg.cohab.c)
    
  
  # Mean durations
  data.params$durs_c <- (min(data$altersCohab$len)+(2*(median(data$altersCohab$len)))+max(data$altersCohab$len))/4
  data.params$durs_c <- (data.params$durs_c/12)*365
  data.params$durs_p <- (min(data$altersPers$len)+(2*(median(data$altersPers$len)))+max(data$altersPers$len))/4
  data.params$durs_p <- (data.params$durs_p/12)*365 
  

  
  ##If using msm.msmf set the population proportions to those in the data for msm and msmf.
  if(msm.msmf==TRUE){

    data.params$msm.frac<-max(0,sum(data$egos$sex.ident=="msm")/sum(data$egos$sex=="M"))
    data.params$msmf.frac<-max(0,sum(data$egos$sex.ident=="msmf")/sum(data$egos$sex=="M"))
  }
  


  ##Make vector with age, race and sex specific proportions for demographic consistency in birth and deaths.
  
  data$egos$demog.cat<-rep(NA,length(data$egos$ego))
  
  sex.groups<-sort(unique(data$egos$sex))
  for (i in 1:(length(sex.groups))){
    data$egos$demog.cat<-ifelse(data$egos$sex==sex.groups[i],i*1000,data$egos$demog.cat)      
  }
  
  race.groups<-sort(unique(data$egos$race))
  for (i in 1:(length(race.groups))){
    data$egos$demog.cat<-ifelse(data$egos$race==race.groups[i],data$egos$demog.cat+(i*100),data$egos$demog.cat)      
  }
  
  
  age.groups<-sort(unique(data$egos$age))
  for (i in 1:(length(age.groups))){
    data$egos$demog.cat<-ifelse (data$egos$age==age.groups[i],data$egos$demog.cat+(age.groups[i]),data$egos$demog.cat)      
  }
  

  temp<-table(data$egos$demog.cat)
  data.params$demog.list<-as.numeric(names(temp))
  data.params$demog.dist<-as.numeric(temp)
  data.params$demog.dist<-data.params$demog.dist/(sum(data.params$demog.dist))
  data.params$sex.groups<-sex.groups
  data.params$race.groups<-race.groups
  data.params$age.groups<-age.groups
  
  ##Weight the demog.dist vector
  if(any(data$egos$weight>0)){
    for(i in 1:(length(data.params$demog.list))){
      ids<-which(data$egos$demog.cat==data.params$demog.list[i])
      weight<-mean(data$egos$weight[ids])
      data.params$demog.dist[i]<-data.params$demog.dist[i]*weight}
  } 
  
  data.params$demog.dist<-data.params$demog.dist/mean(data$egos$weight)
  #class(data.params) <- "data.params"
  
  #Create an age catagory
  data$egos$agecat<-rep(NA,length(data$egos$age))
  data$egos$agecat<-ifelse(data$egos$age < 26 ,"18-25",
                           ifelse(data$egos$age > 25 & data$egos$age < 36,"26-35",
                           ifelse(data$egos$age > 35 & data$egos$age < 46,"36-45",
                           ifelse(data$egos$age > 45,"46-60",data$egos$agecat))))
  

  #Create a new factor for age by sex group.
  data$egos$race.sex<-rep(NA,length(data$egos$ego))
  data$egos$race.sex<-ifelse(data$egos$race=="B" & data$egos$sex=="M","B.M",
                             ifelse(data$egos$race=="BI" & data$egos$sex=="M","BI.M",
                             ifelse(data$egos$race=="H" & data$egos$sex=="M","H.M",
                             ifelse(data$egos$race=="HI" & data$egos$sex=="M","HI.M",
                             ifelse(data$egos$race=="W" & data$egos$sex=="M","W.M",
                             ifelse(data$egos$race=="B" & data$egos$sex=="F","B.F",
                             ifelse(data$egos$race=="BI" & data$egos$sex=="F","BI.F",
                             ifelse(data$egos$race=="H" & data$egos$sex=="F","H.F",
                             ifelse(data$egos$race=="HI" & data$egos$sex=="F","HI.F",
                             ifelse(data$egos$race=="W" & data$egos$sex=="F","W.F",
                                    data$egos$race.sex))))))))))

  data$altersCohab$race.sex<-rep(NA,length(data$altersCohab$race))
  data$altersCohab$race.sex<-ifelse(data$altersCohab$race=="B" & data$altersCohab$sex=="M","B.M",
                             ifelse(data$altersCohab$race=="BI" & data$altersCohab$sex=="M","BI.M",
                             ifelse(data$altersCohab$race=="H" & data$altersCohab$sex=="M","H.M",
                             ifelse(data$altersCohab$race=="HI" & data$altersCohab$sex=="M","HI.M",
                             ifelse(data$altersCohab$race=="W" & data$altersCohab$sex=="M","W.M",
                             ifelse(data$altersCohab$race=="B" & data$altersCohab$sex=="F","B.F",
                             ifelse(data$altersCohab$race=="BI" & data$altersCohab$sex=="F","BI.F",
                             ifelse(data$altersCohab$race=="H" & data$altersCohab$sex=="F","H.F",
                             ifelse(data$altersCohab$race=="HI" & data$altersCohab$sex=="F","HI.F",
                             ifelse(data$altersCohab$race=="W" & data$altersCohab$sex=="F","W.F",
                                    data$altersCohab$race.sex))))))))))
  
  data$altersPers$race.sex<-rep(NA,length(data$altersPers$race))  
  data$altersPers$race.sex<-ifelse(data$altersPers$race=="B" & data$altersPers$sex=="M","B.M",
                              ifelse(data$altersPers$race=="BI" & data$altersPers$sex=="M","BI.M",
                              ifelse(data$altersPers$race=="H" & data$altersPers$sex=="M","H.M",
                              ifelse(data$altersPers$race=="HI" & data$altersPers$sex=="M","HI.M",
                              ifelse(data$altersPers$race=="W" & data$altersPers$sex=="M","W.M",
                              ifelse(data$altersPers$race=="B" & data$altersPers$sex=="F","B.F",
                              ifelse(data$altersPers$race=="BI" & data$altersPers$sex=="F","BI.F",
                              ifelse(data$altersPers$race=="H" & data$altersPers$sex=="F","H.F",
                              ifelse(data$altersPers$race=="HI" & data$altersPers$sex=="F","HI.F",
                              ifelse(data$altersPers$race=="W" & data$altersPers$sex=="F","W.F",
                                     data$altersPers$race.sex))))))))))
  
  data$altersOT$race.sex<-rep(NA,length(data$altersOT$race))   
  data$altersOT$race.sex<-ifelse(data$altersOT$race=="B" & data$altersOT$sex=="M","B.M",
                          ifelse(data$altersOT$race=="BI" & data$altersOT$sex=="M","BI.M",
                          ifelse(data$altersOT$race=="H" & data$altersOT$sex=="M","H.M",
                          ifelse(data$altersOT$race=="HI" & data$altersOT$sex=="M","HI.M",
                          ifelse(data$altersOT$race=="W" & data$altersOT$sex=="M","W.M",
                          ifelse(data$altersOT$race=="B" & data$altersOT$sex=="F","B.F",
                          ifelse(data$altersOT$race=="BI" & data$altersOT$sex=="F","BI.F",
                          ifelse(data$altersOT$race=="H" & data$altersOT$sex=="F","H.F",
                          ifelse(data$altersOT$race=="HI" & data$altersOT$sex=="F","HI.F",
                          ifelse(data$altersOT$race=="W" & data$altersOT$sex=="F","W.F",
                                 data$altersOT$race.sex))))))))))
  
  

  #Create a new factor for age by sex group by pers.c.
  
data$egos$race.sex.pers<-rep(NA,length(data$egos$ego)) 
data$egos$race.sex.pers<-ifelse(data$egos$race=="B" & data$egos$sex=="M" & data$egos$deg.pers.c==0,"Ref.0",
ifelse(data$egos$race=="BI" & data$egos$sex=="M" & data$egos$deg.pers.c==0,"Ref.0",
ifelse(data$egos$race=="H" & data$egos$sex=="M" & data$egos$deg.pers.c==0,"Ref.0",
ifelse(data$egos$race=="HI" & data$egos$sex=="M" & data$egos$deg.pers.c==0,"Ref.0",
ifelse(data$egos$race=="W" & data$egos$sex=="M" & data$egos$deg.pers.c==0,"Ref.0",
ifelse(data$egos$race=="B" & data$egos$sex=="F" & data$egos$deg.pers.c==0,"Ref.0",
ifelse(data$egos$race=="BI" & data$egos$sex=="F" & data$egos$deg.pers.c==0,"Ref.0",
ifelse(data$egos$race=="H" & data$egos$sex=="F" & data$egos$deg.pers.c==0,"Ref.0",
ifelse(data$egos$race=="HI" & data$egos$sex=="F" & data$egos$deg.pers.c==0,"Ref.0",
ifelse(data$egos$race=="W" & data$egos$sex=="F" & data$egos$deg.pers.c==0,"Ref.0",
ifelse(data$egos$race=="B" & data$egos$sex=="M" & data$egos$deg.pers.c==1,"B.M.p1",
ifelse(data$egos$race=="BI" & data$egos$sex=="M" & data$egos$deg.pers.c==1,"BI.M.p1",
ifelse(data$egos$race=="H" & data$egos$sex=="M" & data$egos$deg.pers.c==1,"H.M.p1",
ifelse(data$egos$race=="HI" & data$egos$sex=="M" & data$egos$deg.pers.c==1,"HI.M.p1",
ifelse(data$egos$race=="W" & data$egos$sex=="M" & data$egos$deg.pers.c==1,"W.M.p1",
ifelse(data$egos$race=="B" & data$egos$sex=="F" & data$egos$deg.pers.c==1,"B.F.p1",
ifelse(data$egos$race=="BI" & data$egos$sex=="F" & data$egos$deg.pers.c==1,"BI.F.p1",
ifelse(data$egos$race=="H" & data$egos$sex=="F" & data$egos$deg.pers.c==1,"H.F.p1",
ifelse(data$egos$race=="HI" & data$egos$sex=="F" & data$egos$deg.pers.c==1,"HI.F.p1",
ifelse(data$egos$race=="W" & data$egos$sex=="F" & data$egos$deg.pers.c==1,"W.F.p1",
       data$egos$race.sex.pers))))))))))))))))))))
  
  

  #Create a new factor for age by sex group by main.c 
data$egos$race.sex.cohab<-rep(NA,length(data$egos$ego))
data$egos$race.sex.cohab<-ifelse(data$egos$race=="B" & data$egos$sex=="M" & data$egos$deg.cohab.c==0,"Ref.0",
ifelse(data$egos$race=="BI" & data$egos$sex=="M" & data$egos$deg.cohab.c==0,"Ref.0",
ifelse(data$egos$race=="H" & data$egos$sex=="M" & data$egos$deg.cohab.c==0,"Ref.0",
ifelse(data$egos$race=="HI" & data$egos$sex=="M" & data$egos$deg.cohab.c==0,"Ref.0",
ifelse(data$egos$race=="W" & data$egos$sex=="M" & data$egos$deg.cohab.c==0,"Ref.0",
ifelse(data$egos$race=="B" & data$egos$sex=="F" & data$egos$deg.cohab.c==0,"Ref.0",
ifelse(data$egos$race=="BI" & data$egos$sex=="F" & data$egos$deg.cohab.c==0,"Ref.0",
ifelse(data$egos$race=="H" & data$egos$sex=="F" & data$egos$deg.cohab.c==0,"Ref.0",
ifelse(data$egos$race=="HI" & data$egos$sex=="F" & data$egos$deg.cohab.c==0,"Ref.0",
ifelse(data$egos$race=="W" & data$egos$sex=="F" & data$egos$deg.cohab.c==0,"Ref.0",
ifelse(data$egos$race=="B" & data$egos$sex=="M" & data$egos$deg.cohab.c==1,"B.M.p1",
ifelse(data$egos$race=="BI" & data$egos$sex=="M" & data$egos$deg.cohab.c==1,"BI.M.p1",
ifelse(data$egos$race=="H" & data$egos$sex=="M" & data$egos$deg.cohab.c==1,"H.M.p1",
ifelse(data$egos$race=="HI" & data$egos$sex=="M" & data$egos$deg.cohab.c==1,"HI.M.p1",
ifelse(data$egos$race=="W" & data$egos$sex=="M" & data$egos$deg.cohab.c==1,"W.M.p1",
ifelse(data$egos$race=="B" & data$egos$sex=="F" & data$egos$deg.cohab.c==1,"B.F.p1",
ifelse(data$egos$race=="BI" & data$egos$sex=="F" & data$egos$deg.cohab.c==1,"BI.F.p1",
ifelse(data$egos$race=="H" & data$egos$sex=="F" & data$egos$deg.cohab.c==1,"H.F.p1",
ifelse(data$egos$race=="HI" & data$egos$sex=="F" & data$egos$deg.cohab.c==1,"HI.F.p1",
ifelse(data$egos$race=="W" & data$egos$sex=="F" & data$egos$deg.cohab.c==1,"W.F.p1",
       data$egos$race.sex.cohab))))))))))))))))))))


  #Create a new factor for age by sex group degree (1).
 
  return(list(data.params,data))

}


