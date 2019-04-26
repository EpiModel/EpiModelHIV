
#' @title Condom Use Module
#'
#' @description Module function stochastically simulates potential condom use
#'              for each act on the discordant edgelist.
#'
#' @inheritParams aging_camplc
#'
#' @details
#' For each act on the discordant edgelist, condom use is stochastically simulated
#' based on the partnership type and racial combination of the dyad. Other
#' modifiers for the probability of condom use in that pair are diagnosis of
#' disease, disclosure of status, and full or partial HIV viral suppression
#' given HIV anti-retroviral therapy.
#'
#' @return
#' Updates the discordant edgelist with a \code{uai} variable indicating whether
#' condoms were used in that act.
#'
#' @keywords module msm
#' @export
#'
condoms_msm <- function(dat, at) {

  for (type in c("main", "pers", "asmm", "inst")) {

    ## Variables ##

    # Attributes
    uid <- dat$attr$uid
    diag.status <- dat$attr$diag.status
    race <- dat$attr$race
    age <- dat$attr$age
    cond.int.active <- dat$attr$cond.int.active

    # Parameters
    cond.rr.BB <- dat$param$cond.rr.BB
    cond.rr.BW <- dat$param$cond.rr.BW
    cond.rr.WW <- dat$param$cond.rr.WW
    cond.asmm.by.age <- dat$param$cond.asmm.by.age
    cond.edu <- dat$param$cond.edu
    cond.post.edu.rr <- dat$param$cond.post.edu.rr

    if (type == "main") {
      cond.BB.prob <- dat$param$cond.main.BB.prob
      cond.BW.prob <- dat$param$cond.main.BW.prob
      cond.WW.prob <- dat$param$cond.main.WW.prob
      diag.beta <- dat$param$cond.diag.main.beta
      discl.beta <- dat$param$cond.discl.main.beta
      cond.always <- NULL
      ptype <- 1
    }
    if (type == "pers") {
      cond.BB.prob <- dat$param$cond.pers.BB.prob
      cond.BW.prob <- dat$param$cond.pers.BW.prob
      cond.WW.prob <- dat$param$cond.pers.WW.prob
      cond.post.edu.prob <- dat$param$cond.post.edu.prob
      diag.beta <- dat$param$cond.diag.pers.beta
      discl.beta <- dat$param$cond.discl.pers.beta
      cond.always <- dat$attr$cond.always.pers
      ptype <- 2
    }
    if (type == "asmm") {
      cond.BB.prob <- dat$param$cond.asmm.BB.prob
      cond.BW.prob <- dat$param$cond.asmm.BW.prob
      cond.WW.prob <- dat$param$cond.asmm.WW.prob
      cond.13.15 <- dat$param$cond.asmm.13.15
      cond.16.17 <- dat$param$cond.asmm.16.17
      cond.18 <- dat$param$cond.asmm.18
      cond.post.edu.prob <- dat$param$cond.post.edu.prob
      diag.beta <- dat$param$cond.diag.asmm.beta
      discl.beta <- dat$param$cond.discl.asmm.beta
      cond.always <- dat$attr$cond.always.asmm
      ptype <- 3
    }
    if (type == "inst") {
      cond.BB.prob <- dat$param$cond.inst.BB.prob
      cond.BW.prob <- dat$param$cond.inst.BW.prob
      cond.WW.prob <- dat$param$cond.inst.WW.prob
      cond.post.edu.prob <- dat$param$cond.post.edu.prob
      diag.beta <- dat$param$cond.diag.inst.beta
      discl.beta <- dat$param$cond.discl.inst.beta
      cond.always <- dat$attr$cond.always.inst
      ptype <- 4
    }

    el <- dat$temp$el
    elt <- el[el[, "ptype"] == ptype, ]

    ## Process ##


    # Base condom probs
    race.p1 <- race[elt[, 1]]
    race.p2 <- race[elt[, 2]]
    num.B <- (race.p1 == "B") + (race.p2 == "B")
    cond.prob <- (num.B == 2) * (cond.BB.prob * cond.rr.BB) +
      (num.B == 1) * (cond.BW.prob * cond.rr.BW) +
      (num.B == 0) * (cond.WW.prob * cond.rr.WW)
    

    age.p1 <- pmin(floor(age[elt[, 1]]), floor(age[elt[, 2]]))
    age.p2 <- pmax(floor(age[elt[, 1]]), floor(age[elt[, 2]]))
    
    age.p1 <- ifelse(age.p1 < 16,1,
                     ifelse(age.p1 == 16 | age.p1 == 17,2,
                            ifelse(age.p1 == 18,3,99)))
    
    age.p2 <- ifelse(age.p2 < 16,1,
                     ifelse(age.p2 == 16 | age.p2 == 17,2,
                            ifelse(age.p2 == 18,3,99)))
    
    # Change the probabilities for relationships with an asmm if using age specific probabilities  
    if (cond.asmm.by.age == TRUE & ptype == 3){
      
    cond.prob <- ifelse(age.p1 == 1 & age.p2 == 1, cond.13.15,
           ifelse(age.p1 == 1 & age.p2 == 2, (cond.13.15 + cond.16.17)/2,
                  ifelse(age.p1 == 1 & age.p2 == 3, (cond.13.15 + cond.18)/2,
                         ifelse(age.p1 == 1 & age.p2 == 99, cond.13.15,
                                ifelse(age.p1 == 2 & age.p2 == 2, cond.16.17,
                                       ifelse(age.p1 == 2 & age.p2 == 3, (cond.16.17 + cond.16.17)/2,
                                              ifelse(age.p1 == 2 & age.p2 == 99, cond.16.17,
                                                     ifelse(age.p1 == 3 & (age.p2 == 3 | age.p2 == 99), cond.18 , cond.prob))))))))

    }
    
    
    if (cond.asmm.by.age == TRUE & cond.edu == TRUE & (ptype == 2 | ptype == 3 | ptype == 4)){
      
      p1 <- pmin(elt[, 1],elt[, 2])
      p2 <- pmax(elt[, 1],elt[, 2])
      p1.c <- cond.int.active[p1]
      p2.c <- cond.int.active[p2]
      
      # Update condom probs for education effect
      
      cond.prob <-        ifelse(age.p1 == 1 & p1.c == 1 & age.p2 == 1 & p2.c ==0, (cond.13.15 * cond.post.edu.rr + cond.13.15)/2,
                          ifelse(age.p1 == 1 & p1.c == 0 & age.p2 == 1 & p2.c ==1, (cond.13.15 * cond.post.edu.rr + cond.13.15)/2,
                          ifelse(age.p1 == 1 & p1.c == 1 & age.p2 == 1 & p2.c ==1, (cond.13.15 * cond.post.edu.rr),
                          ifelse(age.p1 == 1 & p1.c == 1 & age.p2 == 2 & p2.c == 0, (cond.13.15 * cond.post.edu.rr + cond.16.17)/2,
                          ifelse(age.p1 == 1 & p1.c == 0 & age.p2 == 2 & p2.c == 1, (cond.13.15 + cond.16.17 * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 1 & p1.c == 1 & age.p2 == 2 & p2.c == 1, (cond.13.15 * cond.post.edu.rr + cond.16.17 * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 1 & p1.c == 1 & age.p2 == 3 & p2.c == 0, (cond.13.15 * cond.post.edu.rr + cond.18)/2,
                          ifelse(age.p1 == 1 & p1.c == 0 & age.p2 == 3 & p2.c == 1, (cond.13.15 + cond.18 * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 1 & p1.c == 1 & age.p2 == 3 & p2.c == 1, (cond.13.15 * cond.post.edu.rr + cond.18 * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 1 & p1.c == 1 & age.p2 == 99 & p2.c == 0, (cond.13.15 * cond.post.edu.rr + cond.WW.prob)/2,
                          ifelse(age.p1 == 1 & p1.c == 0 & age.p2 == 99 & p2.c == 1, (cond.13.15 + cond.WW.prob * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 1 & p1.c == 1 & age.p2 == 99 & p2.c == 1, (cond.13.15 * cond.post.edu.rr + cond.WW.prob * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 2 & p1.c == 1 & age.p2 == 2 & p2.c == 0, (cond.16.17 * cond.post.edu.rr + cond.16.17)/2,
                          ifelse(age.p1 == 2 & p1.c == 0 & age.p2 == 2 & p2.c == 1, (cond.16.17 + cond.16.17 * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 2 & p1.c == 1 & age.p2 == 2 & p2.c == 1, (cond.16.17 * cond.post.edu.rr),
                          ifelse(age.p1 == 2 & p1.c == 1 & age.p2 == 3 & p2.c == 0, (cond.16.17 * cond.post.edu.rr + cond.18)/2,
                          ifelse(age.p1 == 2 & p1.c == 0 & age.p2 == 3 & p2.c == 1, (cond.16.17 + cond.18 * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 2 & p1.c == 1 & age.p2 == 3 & p2.c == 1, (cond.16.17 * cond.post.edu.rr + cond.18 * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 2 & p1.c == 1 & age.p2 == 99 & p2.c == 0, (cond.16.17 * cond.post.edu.rr + cond.WW.prob)/2,
                          ifelse(age.p1 == 2 & p1.c == 0 & age.p2 == 99 & p2.c == 1, (cond.16.17 + cond.WW.prob * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 2 & p1.c == 1 & age.p2 == 99 & p2.c == 1, (cond.16.17 * cond.post.edu.rr + cond.WW.prob * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 3 & p1.c == 1 & age.p2 == 3 & p2.c == 0, (cond.18 * cond.post.edu.rr + cond.18)/2,
                          ifelse(age.p1 == 3 & p1.c == 0 & age.p2 == 3 & p2.c == 1, (cond.18 + cond.18 * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 3 & p1.c == 1 & age.p2 == 3 & p2.c == 1, (cond.18 * cond.post.edu.rr + cond.18 * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 3 & p1.c == 1 & age.p2 == 99 & p2.c == 0, (cond.18 * cond.post.edu.rr + cond.WW.prob)/2,
                          ifelse(age.p1 == 3 & p1.c == 0 & age.p2 == 99 & p2.c == 1, (cond.18 + cond.WW.prob * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 3 & p1.c == 1 & age.p2 == 99 & p2.c == 1, (cond.18 * cond.post.edu.rr + cond.WW.prob * cond.post.edu.rr)/2,
                          ifelse(age.p1 == 99 & p1.c == 1 & age.p2 == 99 & p2.c == 0, (cond.WW.prob * cond.post.edu.rr + cond.WW.prob)/2, 
                          ifelse(age.p1 == 99 & p1.c == 0 & age.p2 == 99 & p2.c == 1, (cond.WW.prob + cond.WW.prob * cond.post.edu.rr)/2, 
                          ifelse(age.p1 == 99 & p1.c == 1 & age.p2 == 99 & p2.c == 1, (cond.WW.prob * cond.post.edu.rr),cond.prob))))))))))))))))))))))))))))))
      
    }
    
    
    
    
 
    # Transform to UAI logit
    uai.prob <- 1 - cond.prob
    uai.logodds <- log(uai.prob / (1 - uai.prob))

    # Diagnosis modifier
    pos.diag <- diag.status[elt[, 1]]
    isDx <- which(pos.diag == 1)
    uai.logodds[isDx] <- uai.logodds[isDx] + diag.beta

    # Disclosure modifier
    isDiscord <- which((elt[, "st1"] - elt[, "st2"]) == 1)
    
    if (length(isDiscord) > 0){
    delt <- elt[isDiscord,, drop=FALSE]
    discl.list <- dat$temp$discl.list
    disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
    delt.cdl <- uid[delt[, 1]] * 1e7 + uid[delt[, 2]]
    discl.disc <- (delt.cdl %in% disclose.cdl)

    discl <- rep(NA, nrow(elt))
    discl[isDiscord] <- discl.disc

    isDisc <- which(discl == 1)
    uai.logodds[isDisc] <- uai.logodds[isDisc] + discl.beta
    }

    # Back transform to prob
    old.uai.prob <- uai.prob
    uai.prob <- exp(uai.logodds) / (1 + exp(uai.logodds))

    uai.prob[is.na(uai.prob) & old.uai.prob == 0] <- 0
    uai.prob[is.na(uai.prob) & old.uai.prob == 1] <- 1

    # UAI group
    if (type %in% c("pers", "inst")) {
      ca1 <- cond.always[elt[, 1]]
      ca2 <- cond.always[elt[, 2]]
      uai.prob <- ifelse(ca1 == 1 | ca2 == 1, 0, uai.prob)
      if (type == "pers") {
        dat$epi$cprob.always.pers[at] <- mean(uai.prob == 0)
      } else {
        dat$epi$cprob.always.inst[at] <- mean(uai.prob == 0)
      }
    }

    # PrEP Status (risk compensation)
    rcomp.prob <- dat$param$rcomp.prob
    rcomp.adh.groups <- dat$param$rcomp.adh.groups
    if (rcomp.prob > 0) {

      prepStat <- dat$attr$prepStat
      prepClass <- dat$attr$prepClass

      idsRC <- which((prepStat[elt[, 1]] == 1 & prepClass[elt[, 1]] %in% rcomp.adh.groups) |
                       (prepStat[elt[, 2]] == 1 & prepClass[elt[, 2]] %in% rcomp.adh.groups))

      if (dat$param$rcomp.main.only == TRUE & ptype > 1) {
        idsRC <- NULL
      }
      if (dat$param$rcomp.discl.only == TRUE) {
        idsRC <- intersect(idsRC, isDisc)
      }
      uai.prob[idsRC] <- 1 - (1 - uai.prob[idsRC]) * (1 - rcomp.prob)
    }

    ai.vec <- elt[, "ai"]
    pos <- rep(elt[, "p1"], ai.vec)
    neg <- rep(elt[, "p2"], ai.vec)
    ptype <- rep(elt[, "ptype"], ai.vec)

    uai.prob.peract <- rep(uai.prob, ai.vec)
    uai <- rbinom(length(pos), 1, uai.prob.peract)

    if (type == "main") {
      pid <- rep(1:length(ai.vec), ai.vec)
      al <- cbind(pos, neg, ptype, uai, pid)
    } else {
      pid <- rep(max(al[, "pid"]) + (1:length(ai.vec)), ai.vec)
      tmp.al <- cbind(pos, neg, ptype, uai, pid)
      al <- rbind(al, tmp.al)
    }

  } # end ptype loop

  dat$temp$al <- al

  if (at == 2) {
    dat$epi$ai.events <- rep(NA, 2)
    dat$epi$uai.events <- rep(NA, 2)
  }
  dat$epi$ai.events[at] <- nrow(al)
  dat$epi$uai.events[at] <- sum(al[, "uai"])

  return(dat)
}
