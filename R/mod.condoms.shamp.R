
#' @title Condom Use Module for 5 race groups and het/msm.
#'
#' @description Module function stochastically simulates potential condom use
#'              for each act on the discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' For each act on the discordant edgelist, condom use is stochastically simulated
#' based on the partnership type, racial combination of the dyad and whether the
#' relationship is male-male (ai) or hetersexual (vi). Other
#' modifiers for the probability of condom use in that pair are diagnosis of
#' disease, disclosure of status, and full or partial HIV viral suppression
#' given HIV anti-retroviral therapy.
#'
#' @return
#' Updates the discordant edgelist with a \code{uai} variable indicating whether
#' condoms were used in that act for male-male dyads and a \code{uvi} variable 
#' indicating whether condoms were used in that act for heterosexual dyads.
#'
#' @keywords module SHAMP msm het vi ai
#' @export
#'
condoms_shamp <- function(dat, at) {

  for (type in c("main", "pers", "inst")) {

    ## Variables ##

    # Attributes
    uid <- dat$attr$uid
    diag.status <- dat$attr$diag.status
    race <- dat$attr$race
    sex.ident <- dat$attr$sex.ident
    sex <- dat$attr$sex
    status<-dat$attr$status

    # Parameters

    if (type == "main") {
      cond.B.prob.msm <- dat$param$cond.main.B.prob.msm
      cond.BI.prob.msm <- dat$param$cond.main.BI.prob.msm
      cond.H.prob.msm <- dat$param$cond.main.H.prob.msm
      cond.HI.prob.msm <- dat$param$cond.main.HI.prob.msm
      cond.W.prob.msm <- dat$param$cond.main.W.prob.msm
      
      cond.B.prob.het <- dat$param$cond.main.B.prob.het
      cond.BI.prob.het <- dat$param$cond.main.BI.prob.het
      cond.H.prob.het <- dat$param$cond.main.H.prob.het
      cond.HI.prob.het <- dat$param$cond.main.HI.prob.het
      cond.W.prob.het <- dat$param$cond.main.W.prob.het
      
      diag.beta.msm <- dat$param$cond.diag.main.beta.msm
      discl.beta.msm <- dat$param$cond.discl.main.beta.msm
      diag.beta.het <- dat$param$cond.diag.main.beta.het
      discl.beta.het <- dat$param$cond.discl.main.beta.het
      cond.always.msm <- NULL
      cond.always.het <- NULL
      ptype <- 1
    }
    if (type == "pers") {
      cond.B.prob.msm <- dat$param$cond.pers.B.prob.msm
      cond.BI.prob.msm <- dat$param$cond.pers.BI.prob.msm
      cond.H.prob.msm <- dat$param$cond.pers.H.prob.msm
      cond.HI.prob.msm <- dat$param$cond.pers.HI.prob.msm
      cond.W.prob.msm <- dat$param$cond.pers.W.prob.msm
      cond.always.msm <- dat$attr$cond.always.pers.msm
      
      cond.B.prob.het <- dat$param$cond.pers.B.prob.het
      cond.BI.prob.het <- dat$param$cond.pers.BI.prob.het
      cond.H.prob.het <- dat$param$cond.pers.H.prob.het
      cond.HI.prob.het <- dat$param$cond.pers.HI.prob.het
      cond.W.prob.het <- dat$param$cond.pers.W.prob.het
      cond.always.het <- dat$attr$cond.always.pers.het
      
      diag.beta.msm <- dat$param$cond.diag.pers.beta.msm
      discl.beta.msm <- dat$param$cond.discl.pers.beta.msm
      diag.beta.het <- dat$param$cond.diag.pers.beta.het
      discl.beta.het <- dat$param$cond.discl.pers.beta.het
     
      ptype <- 2
    }
    
    if (type == "inst") {
      cond.B.prob.msm <- dat$param$cond.inst.B.prob.msm
      cond.BI.prob.msm <- dat$param$cond.inst.BI.prob.msm
      cond.H.prob.msm <- dat$param$cond.inst.H.prob.msm
      cond.HI.prob.msm <- dat$param$cond.inst.HI.prob.msm
      cond.W.prob.msm <- dat$param$cond.inst.W.prob.msm
      cond.always.msm <- dat$attr$cond.always.inst.msm
      
      cond.B.prob.het <- dat$param$cond.inst.B.prob.het
      cond.BI.prob.het <- dat$param$cond.inst.BI.prob.het
      cond.H.prob.het <- dat$param$cond.inst.H.prob.het
      cond.HI.prob.het <- dat$param$cond.inst.HI.prob.het
      cond.W.prob.het <- dat$param$cond.inst.W.prob.het
      cond.always.het <- dat$attr$cond.always.inst.het
      
      diag.beta.msm <- dat$param$cond.diag.inst.beta.msm
      discl.beta.msm <- dat$param$cond.discl.inst.beta.msm
      diag.beta.het <- dat$param$cond.diag.inst.beta.het
      discl.beta.het <- dat$param$cond.discl.inst.beta.het

      ptype <- 3
    }

    el <- dat$temp$el

    
      ## Handle ai and vi seperately

    elt.vi <- el[el[, "ptype"] == ptype & el[, "vi"] > 0, ]
    elt.ai <- el[el[, "ptype"] == ptype & el[, "ai"] > 0, ]
    
  ##VI.


    if(nrow(elt.vi) > 0){

    ## Process ##

    # Base condom probs
    cond.prob <- rep(NA, nrow(elt.vi))
    race.p1 <- race[elt.vi[, 1]]
    race.p2 <- race[elt.vi[, 2]]
    sex.p1 <- sex[elt.vi[,1]]
    sex.p2 <- sex[elt.vi[,2]]
    

    num.het <- (sex.p1 == "F") + (sex.p2 == "F")
   
    cond.prob[race.p1=="B" & race.p2=="B" & num.het==1] <- cond.B.prob.het
    cond.prob[race.p1=="B" & race.p2=="BI" & num.het==1] <- mean(c(cond.B.prob.het,cond.BI.prob.het))
    cond.prob[race.p1=="B" & race.p2=="H" & num.het==1] <- mean(c(cond.B.prob.het,cond.H.prob.het))
    cond.prob[race.p1=="B" & race.p2=="HI" & num.het==1] <- mean(c(cond.B.prob.het,cond.HI.prob.het))
    cond.prob[race.p1=="B" & race.p2=="W" & num.het==1] <- mean(c(cond.B.prob.het,cond.W.prob.het))
    cond.prob[race.p1=="BI" & race.p2=="B" & num.het==1] <- mean(c(cond.B.prob.het,cond.BI.prob.het))
    cond.prob[race.p1=="H" & race.p2=="B" & num.het==1] <- mean(c(cond.B.prob.het,cond.H.prob.het))
    cond.prob[race.p1=="HI" & race.p2=="B" & num.het==1] <- mean(c(cond.B.prob.het,cond.HI.prob.het))
    cond.prob[race.p1=="W" & race.p2=="B" & num.het==1] <- mean(c(cond.B.prob.het,cond.W.prob.het))   
    
    cond.prob[race.p1=="BI" & race.p2=="BI" & num.het==1] <- cond.BI.prob.het
    cond.prob[race.p1=="BI" & race.p2=="H" & num.het==1] <- mean(c(cond.BI.prob.het,cond.H.prob.het))
    cond.prob[race.p1=="BI" & race.p2=="HI" & num.het==1] <- mean(c(cond.BI.prob.het,cond.HI.prob.het))
    cond.prob[race.p1=="BI" & race.p2=="W" & num.het==1] <- mean(c(cond.BI.prob.het,cond.W.prob.het))
    cond.prob[race.p1=="H" & race.p2=="BI" & num.het==1] <- mean(c(cond.BI.prob.het,cond.H.prob.het))
    cond.prob[race.p1=="HI" & race.p2=="BI" & num.het==1] <- mean(c(cond.BI.prob.het,cond.HI.prob.het))
    cond.prob[race.p1=="W" & race.p2=="BI" & num.het==1] <- mean(c(cond.BI.prob.het,cond.W.prob.het)) 

    cond.prob[race.p1=="H" & race.p2=="H" & num.het==1] <- cond.H.prob.het
    cond.prob[race.p1=="H" & race.p2=="HI" & num.het==1] <- mean(c(cond.H.prob.het,cond.HI.prob.het))
    cond.prob[race.p1=="H" & race.p2=="W" & num.het==1] <- mean(c(cond.H.prob.het,cond.W.prob.het))
    cond.prob[race.p1=="HI" & race.p2=="H" & num.het==1] <- mean(c(cond.H.prob.het,cond.HI.prob.het))
    cond.prob[race.p1=="W" & race.p2=="H" & num.het==1] <- mean(c(cond.H.prob.het,cond.W.prob.het))
    
    cond.prob[race.p1=="HI" & race.p2=="HI" & num.het==1] <- cond.HI.prob.het
    cond.prob[race.p1=="HI" & race.p2=="W" & num.het==1] <- mean(c(cond.HI.prob.het,cond.W.prob.het))
    cond.prob[race.p1=="W" & race.p2=="HI" & num.het==1] <- mean(c(cond.HI.prob.het,cond.W.prob.het))
    
    cond.prob[race.p1=="W" & race.p2=="W" & num.het==1] <- cond.HI.prob.het
    
   
    # Transform to UI logit

    uvi.prob <- 1 - cond.prob
    uvi.logodds <- log(uvi.prob / (1 - uvi.prob))


    # Diagnosis modifier
    pos.diag <- diag.status[elt.vi[, 1]]
    isDx <- which(pos.diag == 1)
    uvi.logodds[isDx] <- uvi.logodds[isDx] + diag.beta.het

    # Disclosure modifier
    isDiscord <- which((elt.vi[, "st1"] - elt.vi[, "st2"]) == 1)
    delt <- matrix(elt.vi[isDiscord, ],ncol=7)
    discl.list <- dat$temp$discl.list
    disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
    delt.cdl <- uid[delt[, 1]] * 1e7 + uid[delt[, 2]]
    discl.disc <- (delt.cdl %in% disclose.cdl)

    discl <- rep(NA, nrow(elt.vi))
    discl[isDiscord] <- discl.disc

    isDisc <- which(discl == 1)
    uvi.logodds[isDisc] <- uvi.logodds[isDisc] + discl.beta.het
    
    # Back transform to prob
    old.uvi.prob <- uvi.prob
    uvi.prob <- exp(uvi.logodds) / (1 + exp(uvi.logodds))
    
    uvi.prob[is.na(uvi.prob) & old.uvi.prob == 0] <- 0
    uvi.prob[is.na(uvi.prob) & old.uvi.prob == 1] <- 1
    
    # UVI group
    if (type %in% c("pers", "inst")) {
      ca1 <- cond.always.het[elt.vi[, 1]]
      ca2 <- cond.always.het[elt.vi[, 2]]
      uvi.prob <- ifelse(ca1 == 1 | ca2 == 1, 0, uvi.prob)
      if (type == "pers") {
        dat$epi$cprob.always.pers.het[at] <- mean(uvi.prob == 0)
      } else {
        dat$epi$cprob.always.inst.het[at] <- mean(uvi.prob == 0)
      }
    }
    
    vi.vec <- elt.vi[, "vi"]
    pos <- rep(elt.vi[, "p1"], vi.vec)
    neg <- rep(elt.vi[, "p2"], vi.vec)
    ptype <- rep(elt.vi[, "ptype"], vi.vec)
    
    uvi.prob.peract <- rep(uvi.prob, vi.vec)
    uvi <- rbinom(length(pos), 1, uvi.prob.peract)
    


    if (type == "main") {
      pid <- rep(1:length(vi.vec), vi.vec)
      uai<-rep(0,length(uvi))
      al <- cbind(pos, neg, ptype, uai, uvi, pid)
    } else {
      pid <- rep(max(al[, "pid"]) + (1:length(vi.vec)), vi.vec)
      uai<-rep(0,length(uvi))
      tmp.al <- cbind(pos, neg, ptype, uai, uvi, pid)
      al <- rbind(al, tmp.al)
    }
    max.pid.vi<-max(pid)
    
    }
    
### AI
  
    if(nrow(elt.ai) > 0){
    
    ## Process ##
    
    # Base condom probs
    cond.prob <- rep(NA, nrow(elt.ai))
    race.p1 <- race[elt.ai[, 1]]
    race.p2 <- race[elt.ai[, 2]]
    sex.p1 <- sex[elt.ai[,1]]
    sex.p2 <- sex[elt.ai[,2]]
    
    
    num.het <- (sex.p1 == "F") + (sex.p2 == "F")
    
    cond.prob[race.p1=="B" & race.p2=="B" & num.het==0] <- cond.B.prob.msm
    cond.prob[race.p1=="B" & race.p2=="BI" & num.het==0] <- mean(c(cond.B.prob.msm,cond.BI.prob.msm))
    cond.prob[race.p1=="B" & race.p2=="H" & num.het==0] <- mean(c(cond.B.prob.msm,cond.H.prob.msm))
    cond.prob[race.p1=="B" & race.p2=="HI" & num.het==0] <- mean(c(cond.B.prob.msm,cond.HI.prob.msm))
    cond.prob[race.p1=="B" & race.p2=="W" & num.het==0] <- mean(c(cond.B.prob.msm,cond.W.prob.msm))
    cond.prob[race.p1=="BI" & race.p2=="B" & num.het==0] <- mean(c(cond.B.prob.msm,cond.BI.prob.msm))
    cond.prob[race.p1=="H" & race.p2=="B" & num.het==0] <- mean(c(cond.B.prob.msm,cond.H.prob.msm))
    cond.prob[race.p1=="HI" & race.p2=="B" & num.het==0] <- mean(c(cond.B.prob.msm,cond.HI.prob.msm))
    cond.prob[race.p1=="W" & race.p2=="B" & num.het==0] <- mean(c(cond.B.prob.msm,cond.W.prob.msm))   
    
    cond.prob[race.p1=="BI" & race.p2=="BI" & num.het==0] <- cond.BI.prob.msm
    cond.prob[race.p1=="BI" & race.p2=="H" & num.het==0] <- mean(c(cond.BI.prob.msm,cond.H.prob.msm))
    cond.prob[race.p1=="BI" & race.p2=="HI" & num.het==0] <- mean(c(cond.BI.prob.msm,cond.HI.prob.msm))
    cond.prob[race.p1=="BI" & race.p2=="W" & num.het==0] <- mean(c(cond.BI.prob.msm,cond.W.prob.msm))
    cond.prob[race.p1=="H" & race.p2=="BI" & num.het==0] <- mean(c(cond.BI.prob.msm,cond.H.prob.msm))
    cond.prob[race.p1=="HI" & race.p2=="BI" & num.het==0] <- mean(c(cond.BI.prob.msm,cond.HI.prob.msm))
    cond.prob[race.p1=="W" & race.p2=="BI" & num.het==0] <- mean(c(cond.BI.prob.msm,cond.W.prob.msm)) 
    
    cond.prob[race.p1=="H" & race.p2=="H" & num.het==0] <- cond.H.prob.msm
    cond.prob[race.p1=="H" & race.p2=="HI" & num.het==0] <- mean(c(cond.H.prob.msm,cond.HI.prob.msm))
    cond.prob[race.p1=="H" & race.p2=="W" & num.het==0] <- mean(c(cond.H.prob.msm,cond.W.prob.msm))
    cond.prob[race.p1=="HI" & race.p2=="H" & num.het==0] <- mean(c(cond.H.prob.msm,cond.HI.prob.msm))
    cond.prob[race.p1=="W" & race.p2=="H" & num.het==0] <- mean(c(cond.H.prob.msm,cond.W.prob.msm))
    
    cond.prob[race.p1=="HI" & race.p2=="HI" & num.het==0] <- cond.HI.prob.msm
    cond.prob[race.p1=="HI" & race.p2=="W" & num.het==0] <- mean(c(cond.HI.prob.msm,cond.W.prob.msm))
    cond.prob[race.p1=="W" & race.p2=="HI" & num.het==0] <- mean(c(cond.HI.prob.msm,cond.W.prob.msm))
    
    cond.prob[race.p1=="W" & race.p2=="W" & num.het==0] <- cond.W.prob.msm
    
    
    # Transform to UI logit
    uai.prob <- 1 - cond.prob
    uai.logodds <- log(uai.prob / (1 - uai.prob))
    

    # Diagnosis modifier
    pos.diag <- diag.status[elt.ai[, 1]]
    isDx <- which(pos.diag == 1)
    uai.logodds[isDx] <- uai.logodds[isDx] + diag.beta.msm

    # Disclosure modifier
    isDiscord <- which((elt.ai[, "st1"] - elt.ai[, "st2"]) == 1)
    delt <- matrix(elt.vi[isDiscord, ],ncol=7)
    discl.list <- dat$temp$discl.list
    disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
    delt.cdl <- uid[delt[, 1]] * 1e7 + uid[delt[, 2]]
    discl.disc <- (delt.cdl %in% disclose.cdl)
    
    discl <- rep(NA, nrow(elt.ai))
    discl[isDiscord] <- discl.disc
    
    isDisc <- which(discl == 1)
    uai.logodds[isDisc] <- uai.logodds[isDisc] + discl.beta.msm

    # Back transform to prob
    old.uai.prob <- uai.prob
    uai.prob <- exp(uai.logodds) / (1 + exp(uai.logodds))
    
    uai.prob[is.na(uai.prob) & old.uai.prob == 0] <- 0
    uai.prob[is.na(uai.prob) & old.uai.prob == 1] <- 1
    
    # UAI group
    if (type %in% c("pers", "inst")) {
      ca1 <- cond.always.msm[elt.ai[, 1]]
      ca2 <- cond.always.msm[elt.ai[, 2]]
      uai.prob <- ifelse(ca1 == 1 | ca2 == 1, 0, uai.prob)
      if (type == "pers") {
        dat$epi$cprob.always.pers.msm[at] <- mean(uai.prob == 0)
      } else {
        dat$epi$cprob.always.inst.msm[at] <- mean(uai.prob == 0)
      }
    }
    

    ai.vec <- elt.ai[, "ai"]
    pos <- rep(elt.ai[, "p1"], ai.vec)
    neg <- rep(elt.ai[, "p2"], ai.vec)
    ptype <- rep(elt.ai[, "ptype"], ai.vec)
    
    uai.prob.peract <- rep(uai.prob, ai.vec)
    uai <- rbinom(length(pos), 1, uai.prob.peract)
    

    if (type == "main") {
      pid <- rep(1:length(ai.vec), ai.vec)
      if(max.pid.vi > 1){pid<-pid+max.pid.vi}
      uvi <- rep(0,length(uai))
      al <- cbind(pos, neg, ptype, uai, uvi, pid)
    } else {
      pid <- rep(max(al[, "pid"]) + (1:length(ai.vec)), ai.vec)
      if(max.pid.vi > 1){pid<-pid+max.pid.vi}
      uvi <- rep(0,length(uai))
      tmp.al <- cbind(pos, neg, ptype, uai, uvi, pid)
      al <- rbind(al, tmp.al)
    }
    }

  } # end ptype loop

  dat$temp$al <- al

  return(dat)
}
