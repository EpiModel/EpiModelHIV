
#' @title Death Module for up to 5 race groups heterosexuals and MSM.
#'
#' @description Module function for simulting both general and disease-related
#'              deaths among population members.
#'
#' @inheritParams aging_shamp
#'
#' @details
#' Deaths are divided into two categories: general deaths, for which demographic
#' data on age-specific mortality rates applies; and disease-related diseases,
#' for which the rate of death is a function of progression to end-stage AIDS.
#' Which nodes have died is determined stochastically for general deaths using
#' draws from a binomial distribution, and deterministically for disease-related
#' deaths after nodes have reach a maximum viral load value set in the
#' \code{vl.fatal} parameter.
#'
#' @return
#' This function returns the updated \code{dat} object accounting for deaths.
#' The deaths are deactivated from the main and casual networks, as those are in
#' \code{networkDynamic} class objects; dead nodes are not deleted from the
#' instant network until the \code{\link{simnet_shamp}} module for bookkeeping
#' purposes.
#'
#' @keywords module shamp msm het
#' @export
#'
deaths_KTM <- function(dat, at) {

if(dat$param$death_stats==TRUE) {
  ##Tracking disease related stats at death.
  death.stats <- dat$death.stats
  diag.time <- dat$attr$diag.time
  inf.time <- dat$attr$inf.time
  cum.time.on.tx <- dat$attr$cum.time.on.tx
  stage4.onset.time <- dat$attr$stage4.onset.time
  age.inf <- dat$attr$age.inf
  age.diag <- dat$attr$age.diag
}
  
  ## General deaths
  age <- floor(dat$attr$age)
  age.exact <- dat$attr$age
  race <- dat$attr$race
  sex <- dat$attr$sex
  poi <- dat$attr$poi
  deg.cohab.c  <- dat$attr$deg.cohab.c
  deg.pers.c  <- dat$attr$deg.pers.c
  diag.status <- dat$attr$diag.status
  tx.status <- dat$attr$tx.status

  alive.B.f <- which(race == "B" & sex == "F")
  age.B.f <- age[alive.B.f]
  death.B.f.prob <- dat$param$asmr.B.f[age.B.f]
  deaths.B.f <- alive.B.f[rbinom(length(death.B.f.prob), 1, death.B.f.prob) == 1]


  alive.B.m <- which(race == "B" & sex == "M")
  age.B.m <- age[alive.B.m]
  death.B.m.prob <- dat$param$asmr.B.m[age.B.m]
  deaths.B.m <- alive.B.m[rbinom(length(death.B.m.prob), 1, death.B.m.prob) == 1]
  
  #Those over 40 with no ties or ties with other > 40
  remove<-which(poi==1 & deg.cohab.c==0 & deg.pers.c==0)
  
  if(dat$param$kill.poi.rels == TRUE){
    kill<-NULL
    for(i in 1:length(dat$temp.kill)){
      j<-which(dat$attr$uid==dat$temp.kill[i])
      kill<-c(kill,j)
      
    }
    remove<-c(remove,kill)
    dat$temp.kill<-NULL
  }
    

  
  dth.gen <- c(deaths.B.f, deaths.B.m)
  #Select general deaths under 40 for replacement
  dth.gen.rep <- which(age[dth.gen] < 40)
  dth.gen.rep <- length(dth.gen.rep)

  ## Disease deaths
  dth.dis <- which(dat$attr$stage == 4 &
                   dat$attr$vl >= dat$param$vl.fatal)

  #Select dis deaths under 40 for replacement
  dth.dis.rep <- which(age[dth.dis] < 40)
    
  dth.all <- NULL
  dth.all <- unique(c(dth.gen, dth.dis,remove))
  
  ##Those that will age out this time step.
  dth.age <- which(age >= dat$param$exit.age)
  dth.age.pos <- which((age >= dat$param$exit.age) & diag.status ==1)
  dth.age.pos.ontx <- which((age >= dat$param$exit.age) & diag.status == 1  & tx.status == 1)

  
  if (length(dth.all) > 0) {
    dat$attr$active[dth.all] <- 0
    for (i in 1:3) {
      dat$el[[i]] <- tergmLite::delete_vertices(dat$el[[i]], dth.all)
    }
    dat$attr <- deleteAttr(dat$attr, dth.all)
    if (unique(sapply(dat$attr, length)) != attributes(dat$el[[1]])$n) {
      stop("mismatch between el and attr length in death mod")
    }
  }



  ## Summary Output
  dat$epi$dth.gen[at] <- max(0,dth.gen.rep)
  dat$epi$dth.dis[at] <- max(0,length(dth.dis.rep))
  dat$epi$dth.age[at] <-max(0,length(dth.age))
  dat$epi$dth.remove[at] <-max(0,length(remove))
  
  dat$epi$dth.dis.poi[at] <- max(0,sum(age[dth.dis] < 40, na.rm = TRUE), na.rm = TRUE)
  age.poi <- age.exact[dth.dis] 
  poi <- which(age.poi < 40)
  age.poi <- age.poi[poi]
  dat$epi$dth.dis.age.poi[at] <-max(0,mean(age.poi, na.rm = TRUE))
  
  
  dat$epi$dth.age.pos[at] <- max(0,length(dth.age.pos))
  dat$epi$dth.age.pos.ontx[at] <- max(0,length(dth.age.pos.ontx)) 
    
  if(dat$param$death_stats==TRUE){
  
  if (at == 1){
    death.stats$age <- age[dth.dis.rep]
    death.stats$race <- race[dth.dis.rep]
    death.stats$age.inf <- age.inf[dth.dis.rep]
    death.stats$age.diag <- age.diag[dth.dis.rep]
    death.stats$diag.status <- diag.status[dth.dis.rep]
    death.stats$diag.time <- diag.time[dth.dis.rep]
    death.stats$duration.diagnosed <- at - diag.time[dth.dis.rep]
    death.stats$inf.time <- inf.time[dth.dis.rep]
    death.stats$duration.positive <- at - inf.time[dth.dis.rep]
    death.stats$on.treatment <- tx.status[dth.dis.rep]
    death.stats$duration.on.treatment <- cum.time.on.tx[dth.dis.rep]
    death.stats$stage4.onset.time <- stage4.onset.time[dth.dis.rep]
    death.stats$diag.to.stage4 <- stage4.onset.time[dth.dis.rep] - diag.time[dth.dis.rep]  
    death.stats$stage4.to.death <- at - stage4.onset.time[dth.dis.rep]
    
    death.stats$age.out$race <- race[dth.age.pos]
    death.stats$age.out$age.inf <- age.inf[dth.age.pos]
    death.stats$age.out$age.diag <- age.diag[dth.age.pos]
    death.stats$age.out$diag.status <- diag.status[dth.age.pos]
    death.stats$age.out$diag.time <- diag.time[dth.age.pos]
    death.stats$age.out$duration.diagnosed <- at - diag.time[dth.age.pos]
    death.stats$age.out$inf.time <- inf.time[dth.age.pos]
    death.stats$age.out$duration.positive <- at - inf.time[dth.age.pos]
    death.stats$age.out$on.treatment <- tx.status[dth.age.pos]
    death.stats$age.out$duration.on.treatment <- cum.time.on.tx[dth.age.pos]
    death.stats$age.out$stage4.onset.time <- stage4.onset.time[dth.age.pos]
    death.stats$age.out$diag.to.stage4 <- stage4.onset.time[dth.age.pos] - diag.time[dth.age.pos]  
    death.stats$age.out$stage4.to.ageout <- at -stage4.onset.time[dth.age.pos]  
    
    
  }
  
  if (at > 1){

    
    death.stats$age <- c(death.stats$age, age[dth.dis.rep])
    death.stats$race <- c(death.stats$race, race[dth.dis.rep])
    death.stats$diag.status <- c(death.stats$diag.status, diag.status[dth.dis.rep])
    death.stats$diag.time <- c(death.stats$diag.time, diag.time[dth.dis.rep])
    death.stats$duration.diagnosed <- c(death.stats$duration.diagnosed, at-diag.time[dth.dis.rep])
    death.stats$inf.time <- c(death.stats$inf.time, inf.time[dth.dis.rep])
    death.stats$duration.positive <- c(death.stats$duration.positive, at - inf.time[dth.dis.rep])
    death.stats$on.treatment <- c(death.stats$on.treatment, tx.status[dth.dis.rep])
    death.stats$duration.on.treatment <- c(death.stats$duration.on.treatment, cum.time.on.tx[dth.dis.rep])
    death.stats$stage4.onset.time <- c(death.stats$stage4.onset.time, stage4.onset.time[dth.dis.rep])
    death.stats$diag.to.stage4 <- c(death.stats$diag.to.stage4, stage4.onset.time[dth.dis.rep] - diag.time[dth.dis.rep])  
    death.stats$stage4.to.death <- c(death.stats$stage4.to.death, at - stage4.onset.time[dth.dis.rep])
    
    death.stats$age.out$race <- c(death.stats$age.out$race, race[dth.age.pos])
    death.stats$age.out$diag.status <- c(death.stats$age.out$diag.status, diag.status[dth.age.pos])
    death.stats$age.out$diag.time <- c(death.stats$age.out$diag.time, diag.time[dth.age.pos])
    death.stats$age.out$duration.diagnosed <- c(death.stats$age.out$duration.diagnosed, at - diag.time[dth.age.pos])
    death.stats$age.out$inf.time <- c(death.stats$age.out$inf.time, inf.time[dth.age.pos])
    death.stats$age.out$duration.positive <- c(death.stats$age.out$duration.positive, at - inf.time[dth.age.pos])
    death.stats$age.out$on.treatment <- c(death.stats$age.out$on.treatment, tx.status[dth.age.pos])
    death.stats$age.out$duration.on.treatment <- c(death.stats$age.out$duration.on.treatment, cum.time.on.tx[dth.age.pos])
    death.stats$age.out$stage4.onset.time <- c(death.stats$age.out$stage4.onset.time, stage4.onset.time[dth.age.pos])
    death.stats$age.out$diag.to.stage4 <- c(death.stats$age.out$diag.to.stage4, stage4.onset.time[dth.age.pos] - diag.time[dth.age.pos])  
    death.stats$age.out$stage4.to.ageout <- c(death.stats$age.out$stage4.to.ageout, at - stage4.onset.time[dth.age.pos]  )
    
  }
  
  dat$death.stats <- death.stats
  
  }
  
  return(dat)
}

