
#' @title MDD assignment sub-module
#'
#' @description Function to assign MDD (major depression disorder) - the overall sum of an
#' individual's lifetime MDD risk as a function of race, stigma experience and HIV
#' status
#'
#' @inheritParams aging_msm
#'
#' @export
#'



mddassign_msm <- function(dat, at) {

  ## Input
  # Core Attributes
  active            <- get_attr(dat, "active")
  status            <- get_attr(dat, "status")
  diag.status       <- get_attr(dat, "diag.status")
  diag.time         <- get_attr(dat, "diag.time")
  race              <- get_attr(dat, "race")
  sexStigmaClass    <- get_attr(dat, "sexStigmaClass")
  mdd.status        <- get_attr(dat, "mdd.status")
  mdd.diag.status   <- get_attr(dat, "mdd.diag.status")
  mdd.tx.status     <- get_attr(dat, "mdd.tx.status")

  # Parameters
  mdd.prob.base           <- get_param(dat, "mdd.prob.base")
  mddcoef.stigma    <- get_param(dat, "mddcoef.stigma")
  mddcoef.race      <- get_param(dat, "mddcoef.race")
  mddcoef.hiv       <- get_param(dat, "mddcoef.hiv")



  #A. assign mdd for any missing mdd
  #---------------------------------------------------
  #get all who need mdd assignment
  idsEligMDDAssign <- which(
    active == 1 &
      is.na(mdd.status)
    )

  #get mdd.prob by race, stigma class and hiv status
  if (length(idsEligMDDAssign) > 0) {

    mdd.prob <- rep(mdd.prob.base, length(idsEligMDDAssign))
    tp <- qlogis(mdd.prob)

    tp <- tp + log(mddcoef.stigma[sexStigmaClass[idsEligMDDAssign]])
    tp <- tp + log(mddcoef.race[race[idsEligMDDAssign]])
    tp <- tp + log(mddcoef.hiv[ status[idsEligMDDAssign] + 1])
    mdd.prob <- plogis(tp)

    #determine those with mdd==1 and update dat attributes
    mdd_status_new <- as.numeric(runif(length(idsEligMDDAssign)) < mdd.prob)
    idsMDDAssigned_new <- idsEligMDDAssign[mdd_status_new == 1]

    #update dat attributes
    dat <- set_attr(dat, "mdd.status", mdd_status_new, posit_ids = idsEligMDDAssign)
    dat <- set_attr(dat, "mdd.diag.status", 0, posit_ids = idsMDDAssigned_new)
    dat <- set_attr(dat, "mdd.tx.status", 0, posit_ids = idsMDDAssigned_new)
  }



  #B. mdd transitions (0 --> 1) among newly diagnosed
  #---------------------------------------------------
  #get all who may need mdd re-assignment
  newdiag <- which(diag.status == 1 & diag.time == at - 1)

  idsEligMDDtransition <- which(
    active == 1 &
      diag.status == 1 &
      diag.time == at - 1 &
      get_attr(dat, "mdd.status") == 0)

  #get mdd.prob by race, stigma class and diagnosed hiv status
  if (length(idsEligMDDtransition) > 0) {

    mdd.prob <- tp <- rep(0, length(idsEligMDDtransition))

    tp <- rep(-2.23852622321257000, length(idsEligMDDtransition))
    tp <- tp + log(mddcoef.stigma[sexStigmaClass[idsEligMDDtransition]])
    tp <- tp + log(mddcoef.race[race[idsEligMDDtransition]])
    tp <- tp + log(mddcoef.hiv[status[idsEligMDDtransition] + 1])
    mdd.prob <- plogis(tp)

    #determine those with mdd==1 reassignments
    mdd_status_reassign <- as.numeric(runif(length(idsEligMDDtransition)) < mdd.prob)
    idsMDDAssigned_transitions <- idsEligMDDtransition[mdd_status_reassign == 1]

    #update dat attributes
    dat <- set_attr(dat, "mdd.status", mdd_status_reassign, posit_ids = idsEligMDDtransition)
    dat <- set_attr(dat, "mdd.diag.status", 0, posit_ids = idsMDDAssigned_transitions)
    dat <- set_attr(dat, "mdd.tx.status", 0, posit_ids = idsMDDAssigned_transitions)

    #epi trackers
    dat <- set_epi(dat, "mdd.change.num", at,  length(mdd_status_reassign))
    dat <- set_epi(dat, "mdd.change.prp", at,  length(which(mdd_status_reassign == 1))
                   / length(mdd_status_reassign))
  }


  return(dat)
}






#' @title MDE determination sub-module
#'
#' @description Function to determine and track major depression episodes (MDE) among
#' those with MDD. Transitions set as a function of HIV status where appropriate.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'



mde_msm <- function(dat, at){

  ## Input
  # attributes
  active            <- get_attr(dat, "active")
  age               <- get_attr(dat, "age")
  race              <- get_attr(dat, "race")
  status            <- get_attr(dat, "status")
  diag.status       <- get_attr(dat, "diag.status")
  diag.time         <- get_attr(dat, "diag.time")
  mdd.status        <- get_attr(dat, "mdd.status")
  mdd.diag.status   <- get_attr(dat, "mdd.diag.status")
  mdd.tx.status     <- get_attr(dat, "mdd.tx.status")
  mdd.lasttx.cumlt  <- get_attr(dat, "mdd.lasttx.cumlt")
  mdeCurr           <- get_attr(dat, "mdeCurr")
  mdeCurr.starttime <- get_attr(dat, "mdeCurr.starttime")
  mde.last.endtime  <- get_attr(dat, "mde.last.endtime")
  mde.counter       <- get_attr(dat, "mde.counter")
  mde.symseverity   <- get_attr(dat, "mde.symseverity")
  mde.roleimpair    <- get_attr(dat, "mde.roleimpair")
  prepElig          <- get_attr(dat, "prepElig")


  # parameters
  mde.start.prob    <- get_param(dat, "mde.start.prob")
  mde.recurr.int   <- get_param(dat, "mde.recurr.int")
  mde.symsevgrp.dist<- get_param(dat, "mde.symsevgrp.dist")
  mde.sevimp.symgrp.dist  <- get_param(dat, "mde.sevimp.symgrp.dist")





  ## MDE starts
  #-------------------------------------------------
  tsinceLastMDEended <- at - mde.last.endtime

  #type = first mde (Elig: active, has mdd, never mde)
  mde.start.prob.wk <- c(1 - ((1-mde.start.prob[1])^(1/52)),
                         1 - ((1-mde.start.prob[2])^(1/52)))

  idsEligMDEStart_new <- which(
    active == 1 &
      mdd.status == 1 &
      (is.na(mde.counter) | mde.counter == 0))

  prob.start.new <- mde.start.prob.wk[status [idsEligMDEStart_new] + 1]
  mde_started_new <- as.numeric(runif(length(idsEligMDEStart_new)) < prob.start.new)
  idsMDEStarted_new <- idsEligMDEStart_new[mde_started_new == 1]


  #type = recurrence ff an immediate past episode with no tx or inadequate tx
  idsEligMDEStart_noadeqtx <- which(
    active == 1 &
      mdd.status == 1 &
      mdeCurr == 0 &
      mde.counter >= 1 &
      tsinceLastMDEended >= 9 &
      mdd.tx.status == 0 &
      mdd.lasttx.cumlt < 12)

  prob.recurr.noadeqtx <- 1 - 2^(-1 / mde.recurr.int[1])
  mde_recurr_noadeqtx <- as.numeric(runif(length(idsEligMDEStart_noadeqtx))
                                       < prob.recurr.noadeqtx)
  idsMDEStarted_noadeqtx <- idsEligMDEStart_noadeqtx[mde_recurr_noadeqtx == 1]


  #type = recurrence ff an immediate past epi that was adequately tx
  idsEligMDEStart_adeqtx <- which(
    active == 1 &
      mdd.status == 1 &
      mdeCurr == 0 &
      mde.counter >= 1 &
      tsinceLastMDEended >= 9 &
      mdd.tx.status == 0 &
      mdd.lasttx.cumlt >= 12)

  prob.recurr.adeqtx <- 1 - 2^(-1 / mde.recurr.int[2])
  mde_recurr_adeqtx <- as.numeric(runif(length(idsEligMDEStart_adeqtx))
                                  < prob.recurr.adeqtx)
  idsMDEStarted_adeqtx <- idsEligMDEStart_adeqtx[mde_recurr_adeqtx == 1]


  #all mde starts
  idsMDEStarted <- c(idsMDEStarted_new,
                     idsMDEStarted_noadeqtx,
                     idsMDEStarted_adeqtx)


  #update attributes
  if (length(idsMDEStarted) > 0) {
    mdeCurr[idsMDEStarted]             <- 1
    mdeCurr.starttime[idsMDEStarted]   <- at
    mde.counter[idsMDEStarted] <- ifelse(is.na(mde.counter[idsMDEStarted]), 1,
                                       mde.counter[idsMDEStarted] + 1)

    #assign mde symptom severity and role impairment
    mde.symseverity[idsMDEStarted] <- sample(
      x = 1:4,
      size = length(idsMDEStarted),
      replace = T,
      prob = mde.symsevgrp.dist
    )

    prob.severeimp <- mde.sevimp.symgrp.dist[mde.symseverity[idsMDEStarted]]
    has_severeimp <- as.numeric(runif(length(idsMDEStarted)) < prob.severeimp)
    idsMDEStarted_severeimp <- idsMDEStarted[has_severeimp == 1]

    mde.roleimpair[idsMDEStarted_severeimp] <- 1
  }


  #track new starts vs recurrences
  mde_newstarts_age      <- age[idsMDEStarted_new]
  mde_recurrences       <- c(idsMDEStarted_noadeqtx, idsMDEStarted_adeqtx)
  mde_recurr_remissint <- tsinceLastMDEended[mde_recurrences]



  # update dat attributes
  #-------------------------------------------------
  dat <- set_attr(dat, "mdeCurr", mdeCurr)
  dat <- set_attr(dat, "mdeCurr.starttime", mdeCurr.starttime)
  dat <- set_attr(dat, "mde.counter", mde.counter)
  dat <- set_attr(dat, "mde.last.endtime", mde.last.endtime)
  dat <- set_attr(dat, "mde.symseverity", mde.symseverity)
  dat <- set_attr(dat, "mde.roleimpair", mde.roleimpair)



  # epi trackers
  #-------------------------------------------------
  dat <- set_epi(dat, "mde.ever.numall", at, length(which(mde.counter >= 1)))
  dat <- set_epi(dat, "mde.ever.prpall", at, length(which(mde.counter >= 1)) /
                   length(which(active == 1)))
  dat <- set_epi(dat, "mde.ever.prpmdd", at, length(which(mde.counter >= 1)) /
                   length(which(mdd.status == 1)))
  dat <- set_epi(dat, "mde.ever.prpmddhiv1", at,
                 length(which(mde.counter >= 1 & status == 1)) /
                   length(which(mdd.status == 1 & status == 1)))
  dat <- set_epi(dat, "mde.ever.prpmddhiv0", at,
                 length(which(mde.counter >= 1 & status == 0)) /
                   length(which(mdd.status == 1 & status == 0)))

  dat <- set_epi(dat, "mde.epinum.min", at, suppressWarnings(min(mde.counter, na.rm = T)))
  dat <- set_epi(dat, "mde.epinum.median", at, median(mde.counter, na.rm = T))
  dat <- set_epi(dat, "mde.epinum.max", at, suppressWarnings(max(mde.counter, na.rm = T)))

  dat <- set_epi(dat, "mde.starts.prpnew", at, length(idsMDEStarted_new)/
                   length(idsMDEStarted))

  dat <- set_epi(dat, "mde.start.numall", at, length(idsMDEStarted))
  dat <- set_epi(dat, "mde.start.numhiv0", at,
                 length(intersect(idsMDEStarted, which(status == 0))))
  dat <- set_epi(dat, "mde.start.numhiv1", at,
                 length(intersect(idsMDEStarted, which(status == 1))))

  dat <- set_epi(dat, "mde.start.numsym1", at,
                 length(intersect(idsMDEStarted, which(mde.symseverity == 1))))
  dat <- set_epi(dat, "mde.start.numsym2", at,
                 length(intersect(idsMDEStarted, which(mde.symseverity == 2))))
  dat <- set_epi(dat, "mde.start.numsym3", at,
                 length(intersect(idsMDEStarted, which(mde.symseverity == 3))))
  dat <- set_epi(dat, "mde.start.numsym4", at,
                 length(intersect(idsMDEStarted, which(mde.symseverity == 4))))

  dat <- set_epi(dat, "mde.start.new.numall", at, length(idsMDEStarted_new))
  dat <- set_epi(dat, "mde.start.new.numhiv0", at,
                 length(intersect(idsMDEStarted_new, which(status == 0))))
  dat <- set_epi(dat, "mde.start.new.numhiv1", at,
                 length(intersect(idsMDEStarted_new, which(status == 1))))

  dat <- set_epi(dat, "mde.start.new.age.min", at,
                 suppressWarnings(min(mde_newstarts_age, na.rm = T)))
  dat <- set_epi(dat, "mde.start.new.age.median", at,
                 median(mde_newstarts_age, na.rm = T))
  dat <- set_epi(dat, "mde.start.new.age.max", at,
                 suppressWarnings(max(mde_newstarts_age, na.rm = T)))
  dat <- set_epi(dat, "mde.start.new.age.less30", at,
                 suppressWarnings(length(which(mde_newstarts_age < 30)))/
                   length(mde_newstarts_age))

  dat <- set_epi(dat, "mde.start.recurr.numall", at, length(mde_recurrences))
  dat <- set_epi(dat, "mde.start.recurr.numhiv0", at,
                 length(intersect(mde_recurrences, which(status == 0))))
  dat <- set_epi(dat, "mde.start.recurr.numhiv1", at,
                 length(intersect(mde_recurrences, which(status == 1))))

  dat <- set_epi(dat, "mde.start.recurr.dur.min", at,
                 suppressWarnings(min(mde_recurr_remissint, na.rm = T)))
  dat <- set_epi(dat, "mde.start.recurr.dur.median", at,
                 median(mde_recurr_remissint, na.rm = T))
  dat <- set_epi(dat, "mde.start.recurr.dur.max", at,
                 suppressWarnings(max(mde_recurr_remissint, na.rm = T)))

  dat <- set_epi(dat, "mde.active.numall", at, length(which(mdeCurr == 1)))
  dat <- set_epi(dat, "mde.active.prpmdd", at, length(which(mdeCurr == 1)) /
                   length(which(mdd.status == 1)))
  dat <- set_epi(dat, "mde.active.prpmddhiv0", at,
                 length(which(mdeCurr == 1 & diag.status == 0)) /
                   length(which(mdd.status == 1 & diag.status == 0)))
  dat <- set_epi(dat, "mde.active.prpmddhiv1", at,
                 length(which(mdeCurr == 1 & diag.status == 1)) /
                   length(which(mdd.status == 1 & diag.status == 1)))



return(dat)
}








#' @title MDD Care
#'
#' @description Function to determine and track major depression disorder (MDD) diagnosis.
#' and treatment among MSM with MDD/MDE. Transitions set as a function of HIV status where
#' appropriate.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'



mddcare_msm <- function(dat, at){

  ## Input
  # attributes
  active            <- get_attr(dat, "active")
  age               <- get_attr(dat, "age")
  race              <- get_attr(dat, "race")
  status            <- get_attr(dat, "status")
  diag.status       <- get_attr(dat, "diag.status")
  diag.time         <- get_attr(dat, "diag.time")
  last.neg.test <- get_attr(dat, "last.neg.test")

  mdd.status        <- get_attr(dat, "mdd.status")
  mdd.diag.status   <- get_attr(dat, "mdd.diag.status")
  mdd.diag.time     <- get_attr(dat, "mdd.diag.time")
  mdd.scrn.counter  <- get_attr(dat, "mdd.scrn.counter")
  mdd.lastscrn.time <- get_attr(dat, "mdd.lastscrn.time")
  mdd.tx.status     <- get_attr(dat, "mdd.tx.status")
  mdd.lasttx.cumlt  <- get_attr(dat, "mdd.lasttx.cumlt")
  mdd.txepi.counter <- get_attr(dat, "mdd.txepi.counter")
  mdd.txinit.time   <- get_attr(dat, "mdd.txinit.time")
  mdd.currtx.tottime<- get_attr(dat, "mdd.currtx.tottime")
  mdd.currtx.stage  <- get_attr(dat, "mdd.currtx.stage")
  mdd.txintv.status <- get_attr(dat, "mdd.txintv.status")
  mdd.txintvinit.time<- get_attr(dat, "mdd.txintvinit.time")
  mdd.txintv.last.endtime<- get_attr(dat, "mdd.txintv.last.endtime")

  mdeCurr           <- get_attr(dat, "mdeCurr")
  mdeCurr.starttime <- get_attr(dat, "mdeCurr.starttime")
  mde.last.endtime  <- get_attr(dat, "mde.last.endtime")
  mde.counter       <- get_attr(dat, "mde.counter")
  mde.symseverity   <- get_attr(dat, "mde.symseverity")
  mde.roleimpair    <- get_attr(dat, "mde.roleimpair")
  prepElig          <- get_attr(dat, "prepElig")
  prepStat          <- get_attr(dat, "prepStat")
  prepStartTime     <- get_attr(dat, "prepStartTime")



  # parameters
  mdd.diag.gen.prob        <- get_param(dat, "mdd.diag.gen.prob")
  mddscrnuptk.pstat.prob   <- get_param(dat, "mddscrnuptk.pstat.prob")
  mddscrnuptk.pind.prob    <- get_param(dat, "mddscrnuptk.pind.prob")
  mdd.txinit.prob.reg      <- get_param(dat, "mdd.txinit.prob.reg")
  mdd.txinit.prob.prep     <- get_param(dat, "mdd.txinit.prob.prep")
  mde.spontres.int         <- get_param(dat, "mde.spontres.int")
  mde.txremiss.prob        <- get_param(dat, "mde.txremiss.prob")
  mdd.txltfu.prob          <- get_param(dat, "mdd.txltfu.prob")
  minadeq.txdur            <- get_param(dat, "minadeq.txdur")




  #if(at >= get_param(dat, "mhefx.start") & length(which(last.neg.test==at) > 10)) browser()

  #MDD Screening
  #--------------------------------------------------------
  tsinceCurrMDEStarted <- at - mdeCurr.starttime

  #path A: screening via general route (standard - self/clinician referral)
  #presenting for depression screening in this route is expected to be dependent on
  #symptom severity and presence of any role impairment in affected person.
  #Also more likely that individuals going through a first episode will seek care unless
  #its a particularly severe episode or there is some impairment to their social, relational
  #or occupational roles
  mdd.diag.gen.probwk <- 1 - ((1-mdd.diag.gen.prob)^(1/52))

  idsEligMDDDiag_gen <- which(
    active == 1 &
      mdeCurr == 1 &
      tsinceCurrMDEStarted >= 2 &
      (mde.symseverity > 1 | (mde.symseverity == 1 & mde.roleimpair == 1)) &
      (is.na(mdd.tx.status) | mdd.tx.status == 0))

  diag.prob.gen <- mdd.diag.gen.probwk[status [idsEligMDDDiag_gen] + 1]
  mdd_diag_gen <- as.numeric(runif(length(idsEligMDDDiag_gen)) < diag.prob.gen)
  idsMDDDiag_gen <- idsEligMDDDiag_gen[mdd_diag_gen == 1]



  #path B: screening within context of PrEP [intervention part 1]
  #not referral-based, so expecting screening to not be based on
  #symptom severity, role impairment or being already in mdd care (if routine screening)

  #B1: routine mental health screening implemented for PrEP starters only
  #set these epi trackers to 0 to prevent NAs in early return cases
  dat <- set_epi(dat, "mdd.scrndat.numpstat", at, 0)
  dat <- set_epi(dat, "mdd.diagat.numpstat", at, 0)

  mddscrnuptk.pstat.probwk <- mddscrnuptk.pstat.prob #1 - ((1-mddscrnuptk.pstat.prob)^(1/52))

  idsEligMDDDiag_pstat <- which(
    active == 1 &
      prepStat == 1 &
      prepStartTime == at)
  idsEligMDDDiag_pstat <- setdiff(idsEligMDDDiag_pstat, idsMDDDiag_gen)

  idsMDDDiag_pstat <- numeric()

  if(at >= get_param(dat, "mh.scrninterv.start")) {
    scrn.prob.prep <- mddscrnuptk.pstat.probwk
    mdd_scrn_prep <- as.numeric(runif(length(idsEligMDDDiag_pstat)) < scrn.prob.prep)

    scrnd.pstat <- idsEligMDDDiag_pstat[mdd_scrn_prep == 1]
    dat <- set_epi(dat, "mdd.scrndat.numpstat", at, length(scrnd.pstat))

    idsMDDDiag_pstat <- scrnd.pstat[which(
      (mdeCurr[scrnd.pstat] == 1 & tsinceCurrMDEStarted[scrnd.pstat] >= 2) |
        (mdeCurr[scrnd.pstat] == 0 & at - mde.last.endtime[scrnd.pstat] < 52))]
  }


  #B2: routine mental health screening implemented at time of HIV testing (those who test neg)
  # set epi trackers to 0 to prevent NAs in early return cases
  dat <- set_epi(dat, "mdd.scrndat.numpind", at, 0)
  dat <- set_epi(dat, "mdd.diagat.numpind", at, 0)

  mddscrnuptk.pind.probwk <- mddscrnuptk.pind.prob #1 - ((1-mddscrnuptk.pind.prob)^(1/52))

  idsEligMDDDiag_pind <- which(
    active == 1 &
      last.neg.test == at)
  idsEligMDDDiag_pind <- setdiff(idsEligMDDDiag_pind, idsMDDDiag_gen)

  idsMDDDiag_pind <- numeric()

  if(at >= get_param(dat, "mh.scrninterv.start")) {
    scrn.prob.prep <- mddscrnuptk.pind.probwk
    mdd_scrn_prep <- as.numeric(runif(length(idsEligMDDDiag_pind)) < scrn.prob.prep)

    scrnd.pind <- idsEligMDDDiag_pind[mdd_scrn_prep == 1]
    dat <- set_epi(dat, "mdd.scrndat.numpind", at, length(scrnd.pind))

    idsMDDDiag_pind <- scrnd.pind[which(
      (mdeCurr[scrnd.pind] == 1 & tsinceCurrMDEStarted[scrnd.pind] >= 2) |
        (mdeCurr[scrnd.pind] == 0 & at - mde.last.endtime[scrnd.pind] < 52))]
  }



  #update diagnosis-related attributes
  idsMDDDiag <- c(idsMDDDiag_gen, idsMDDDiag_pstat, idsMDDDiag_pind)
  idsMDDDiag_first <- idsMDDDiag[which(mdd.scrn.counter[idsMDDDiag] == 0)]

  if (length(idsMDDDiag) > 0) {
    mdd.scrn.counter[idsMDDDiag] <- mdd.scrn.counter[idsMDDDiag] + 1
    mdd.lastscrn.time[idsMDDDiag]    <- at

    mdd.diag.status[idsMDDDiag_first]      <- 1
    mdd.diag.time[idsMDDDiag_first]        <- at
  }




  #MDD treatment initiation
  #--------------------------------------------------------
  tsincelastMDDScrn <- at - mdd.lastscrn.time
  mdd.txinit.prob.regwk <- 1 - ((1-mdd.txinit.prob.reg)^(1/52))
  mdd.txinit.prob.prepwk <- 1 - ((1-mdd.txinit.prob.prep)^(1/52))


  #treatment 1: [Intervention part 2] linkage to mental health care with added sexual
  #stigma mitigating component. Per sec. lit. interventions that help are those that aim
  #at increasing resilience and equipping affected persons with adaptive coping strategies
  diag_prep <- c(idsMDDDiag_pstat, idsMDDDiag_pind)
  idsEligMDDtxInit_txintv <- intersect(diag_prep, which(
    # (mde.counter > 1 |
    #    (mde.counter == 1 & mde.symseverity > 1) |
    #    (mde.counter == 1 & mde.symseverity == 1 & mde.roleimpair == 1)) &
      (is.na(mdd.tx.status) | mdd.tx.status == 0)
  ))

  idsMDDtxInit_txintv <- numeric()
  if(at >= get_param(dat, "mh.txinterv.start")) {
    prob.txinit_txintv <- mdd.txinit.prob.prepwk[race [idsEligMDDtxInit_txintv]]
    mdd_txinit_txintv <- as.numeric(runif(length(idsEligMDDtxInit_txintv)) < prob.txinit_txintv)
    idsMDDtxInit_txintv <- idsEligMDDtxInit_txintv[mdd_txinit_txintv == 1]
  }



  #treatment 2: [standard mdd care] linkage to regular depression care (no stigma component)
  idsEligMDDtxInit_txreg <- which(
    active == 1 &
      mdeCurr == 1 &
      tsinceCurrMDEStarted >= 2 &
      (mde.counter > 1 |
         (mde.counter == 1 & mde.symseverity > 1) |
         (mde.counter == 1 & mde.symseverity == 1 & mdd.lasttx.cumlt > 0) |
         (mde.counter == 1 & mde.symseverity == 1 & mde.roleimpair == 1)) &
      mdd.diag.status == 1 &
      (is.na(mdd.tx.status) | mdd.tx.status == 0)
  )
  idsEligMDDtxInit_txreg <- setdiff(idsEligMDDtxInit_txreg, idsMDDtxInit_txintv)

  prob.txInit_txreg <- mdd.txinit.prob.regwk[race [idsEligMDDtxInit_txreg]]
  mdd_txInit_txreg <- as.numeric(runif(length(idsEligMDDtxInit_txreg)) < prob.txInit_txreg)
  idsMDDtxInit_txreg <- idsEligMDDtxInit_txreg[mdd_txInit_txreg == 1]



  idsMDDtxInit <- c(idsMDDtxInit_txreg, idsMDDtxInit_txintv)
  ######
  #if(length(intersect(idsMDDtxInit, which(race==2))) > 0) browser()


  #update tx-related attributes
  if (length(idsMDDtxInit) > 0) {
    mdd.tx.status[idsMDDtxInit]      <- 1
    mdd.txinit.time[idsMDDtxInit]    <- at

    if(length(idsMDDtxInit_txintv) > 0) {
      mdd.txintv.status[idsMDDtxInit_txintv]       <- 1
      mdd.txintvinit.time[idsMDDtxInit_txintv]     <- at
      mdd.txintv.last.endtime[idsMDDtxInit_txintv] <- NA
    }

    mdd.txepi.counter[idsMDDtxInit] <- mdd.txepi.counter[idsMDDtxInit] + 1

    #reset total time on tx for last episode to zero
    mdd.lasttx.cumlt[idsMDDtxInit]  <- 0
  }

  #changes to tx stage and mde.severity/mde.roleimpair with continued tx
  ontx <- which(mdd.tx.status == 1)

  mdd.currtx.tottime[ontx] <- ifelse(mdd.txinit.time[ontx] == at, 0,
    mdd.currtx.tottime[ontx] + 1)

  mdd.currtx.stage[which(mdd.tx.status == 1 & mdd.currtx.tottime <= 4 * 4)] <- 1
  mdd.currtx.stage[which(mdd.tx.status == 1 &
                           mdd.currtx.tottime > 4 * 4 &
                           mdd.currtx.tottime < 9 * 4)] <- 2
  mdd.currtx.stage[which(mdd.tx.status == 1 & mdd.currtx.tottime >= 9 * 4)]<- 3

  stgtransit.acu.conti <- which(mdd.tx.status == 1 & mdd.currtx.tottime == (4*4) + 1)
  mde.roleimpair[stgtransit.acu.conti] <- 0
  mde.symseverity[stgtransit.acu.conti] <- 1

  # stgtransit.conti.mainte <- which(mdd.tx.status == 1 & mdd.currtx.tottime == (9*4) + 1)
  # mde.symseverity[stgtransit.conti.mainte] <- 1




  #MDE remission (spont. res, tx remission & deaths)
  #-------------------------------------------------
  tsincelastMDDScrn <- at - mdd.lastscrn.time

  #mde stop type 1: spontaneous resolution (depends on mde.symseverity)
  #time to spontaneous resolution varies by symptom severity
  idsEligMDEStop_spontres <- which(
    active == 1 &
      mdeCurr == 1 &
      tsinceCurrMDEStarted >= 2 &
      mdd.tx.status == 0
  )

  mde.spontres.prob <- 1 - 2^(-1 / mde.spontres.int[mde.symseverity[idsEligMDEStop_spontres]])
  mde_stopped_spontres <- as.numeric(runif(length(idsEligMDEStop_spontres)) < mde.spontres.prob)
  idsMDEStopped_spontres <- idsEligMDEStop_spontres[mde_stopped_spontres == 1]


  #mde stop type 2a: remission following receipt of minimally adequate treatment
  #remission among those on tx set to occur while in any phase of tx (but has to occur
  #after individual has received minimally adequate tx i.e. >= 8wks of tx). However, prob of
  #remission is turned up if in continuation or maintenance phase of tx.
  idsEligMDEStop_txremiss <- which(
    active == 1 &
      mdeCurr == 1 &
      tsinceCurrMDEStarted >= 2 &
      mdd.tx.status == 1 &
      mdd.currtx.tottime >= minadeq.txdur
  )

  prob.remiss <- 1 - ((1 - mde.txremiss.prob)^(1/52))
  # prob.remiss.base <- 1 - ((1 - mde.txremiss.prob)^(1/52))
  # prob.remiss.elig <- rep(prob.remiss.base, length(idsEligMDEStop_txremiss))
  # #turn up prob of remission if in continuation or maintenance mdd tx phase
  # txphase.elig <- mdd.currtx.stage[idsEligMDEStop_txremiss]
  # mult <- c(1,1.5,2)
  # mult.elig <- mult[txphase.elig]
  # prob.remiss <- plogis(qlogis(prob.remiss.elig) + log(mult.elig))

  mde_stopped_remiss <- as.numeric(runif(length(idsEligMDEStop_txremiss)) < prob.remiss)
  idsMDEStopped_txremiss <- idsEligMDEStop_txremiss[mde_stopped_remiss == 1]


  #mde stop type 2b: remission following completed tx (stage 3 folks only)
  #setting a hard stop to mde after individual has received up to 3 years of treatment.
  idsMDEStopped_endmainte <- which(
    active == 1 &
      mdeCurr == 1 &
      mdd.tx.status == 1 &
      mdd.currtx.stage == 3 &
      mdd.currtx.tottime >= 3*52
  )

  dat <- set_epi(dat, "mde.stop.txremis.prpactiveontx", at,
                 length(c(idsMDEStopped_txremiss, idsMDEStopped_endmainte)) /
                   length(which(active == 1 & mdeCurr == 1 & tsinceCurrMDEStarted >= 2 &
                                        mdd.tx.status == 1)))


  #stop type 3: mde stop due to death (suicide or other)
  idsMDEStopped_dth <- which(
    active == 0 &
      mdeCurr == 1)




  #all mde stops
  idsMDEStopped <- c(idsMDEStopped_spontres,
                     idsMDEStopped_txremiss,
                     idsMDEStopped_endmainte,
                     idsMDEStopped_dth)


  mdeCurr.dur <- tsinceCurrMDEStarted[idsMDEStopped]


  #update attributes
  if (length(idsMDEStopped) > 0) {
    mdeCurr[idsMDEStopped]             <- 0
    mde.last.endtime[idsMDEStopped]    <- at
    mdeCurr.starttime[idsMDEStopped]   <- NA

    #reset symptom severity/role impairment = 0
    mde.symseverity[idsMDEStopped]     <- 0
    mde.roleimpair[idsMDEStopped]      <- 0
  }




  #MDD treatment end
  #-------------------------------------------------
  #tx discontinuation: MDE remission after starting treatment (any tx stage)
  idsMDDtxStopped_remiss <- which(
    active == 1 &
      mdd.tx.status == 1 &
      mde.symseverity == 0 &
      mde.roleimpair == 0 &
      mde.last.endtime == at
  )


  #tx discontinuation: end of maintenance phase
  idsMDDtxStopped_endmainte <- idsMDEStopped_endmainte


  #tx discontinuation: premature drop out from tx (w/o mde remission, any tx stage)
  mdd.txltfu.probwk <- 1 - ((1-mdd.txltfu.prob)^(1/52))

  idsEligMDDtxStop_ltfu <- which(
    active == 1 &
      mdeCurr == 1 &
      mdd.tx.status == 1
  )

  prob.ltfu <- mdd.txltfu.probwk
  mdd_tx_ltfu <- as.numeric(runif(length(idsEligMDDtxStop_ltfu)) < prob.ltfu)
  idsMDDtxStopped_ltfu <- idsEligMDDtxStop_ltfu[mdd_tx_ltfu == 1]


  #tx discontinuation: death
  idsMDDtxStopped_dth <- which(!active & mdd.tx.status == 1)


  #all tx stops
  idsMDDtxStopped <- c(idsMDDtxStopped_remiss,
                       idsMDDtxStopped_endmainte,
                       idsMDDtxStopped_ltfu,
                       idsMDDtxStopped_dth)

  cuml.timeontx <- at - mdd.txinit.time[idsMDDtxStopped]

  #if(length(idsMDDtxStopped) > 0) browser()

  #update attributes
  if(length(idsMDDtxStopped) > 0){
    mdd.tx.status[idsMDDtxStopped]          <- 0
    mdd.lasttx.cumlt[idsMDDtxStopped]       <- cuml.timeontx

    mdd.txintv.status[idsMDDtxStopped]      <- 0
    mdd.txintvinit.time[idsMDDtxStopped]    <- NA
    mdd.txintv.last.endtime[idsMDDtxStopped]<- at

    mdd.txinit.time[idsMDDtxStopped]        <- 0
    mdd.currtx.stage[idsMDDtxStopped]       <- 0
    mdd.currtx.tottime[idsMDDtxStopped]     <- 0

  }





  # update dat attributes
  #-------------------------------------------------
  #diagnosis
  dat <- set_attr(dat, "mdd.diag.status", mdd.diag.status)
  dat <- set_attr(dat, "mdd.diag.time", mdd.diag.time)
  dat <- set_attr(dat, "mdd.scrn.counter", mdd.scrn.counter)
  dat <- set_attr(dat, "mdd.lastscrn.time", mdd.lastscrn.time)
  #tx init and stops
  dat <- set_attr(dat, "mdd.tx.status", mdd.tx.status)
  dat <- set_attr(dat, "mdd.txepi.counter", mdd.txepi.counter)
  dat <- set_attr(dat, "mdd.txinit.time", mdd.txinit.time)
  dat <- set_attr(dat, "mdd.lasttx.cumlt", mdd.lasttx.cumlt)
  dat <- set_attr(dat, "mdd.currtx.tottime", mdd.currtx.tottime)
  dat <- set_attr(dat, "mdd.currtx.stage", mdd.currtx.stage)
  dat <- set_attr(dat, "mdd.txintv.status", mdd.txintv.status)
  dat <- set_attr(dat, "mdd.txintvinit.time", mdd.txintvinit.time)
  dat <- set_attr(dat, "mdd.txintv.last.endtime", mdd.txintv.last.endtime)
  #mde end
  dat <- set_attr(dat, "mdeCurr", mdeCurr)
  dat <- set_attr(dat, "mdeCurr.starttime", mdeCurr.starttime)
  dat <- set_attr(dat, "mde.last.endtime", mde.last.endtime)
  dat <- set_attr(dat, "mde.symseverity", mde.symseverity)
  dat <- set_attr(dat, "mde.roleimpair", mde.roleimpair)



  #epi trackers
  #-------------------------------------------------
  #diagnosis
  dat <- set_epi(dat, "mdd.eligdiagat.numgen", at, length(idsEligMDDDiag_gen))
  dat <- set_epi(dat, "mdd.diagat.numgen", at, length(idsMDDDiag_gen))
  dat <- set_epi(dat, "mdd.diagat.prpeliggen", at, length(idsMDDDiag_gen) /
                   length(idsEligMDDDiag_gen))

  dat <- set_epi(dat, "mdd.PrEPstartsat", at,
                 length(which(prepStat == 1 & prepStartTime==at)))
  dat <- set_epi(dat, "mdd.eligscrnat.numpstat", at, length(idsEligMDDDiag_pstat))
  dat <- set_epi(dat, "mdd.diagat.numpstat", at, length(idsMDDDiag_pstat))
  dat <- set_epi(dat, "mdd.diagat.prpeligpstat", at, length(idsMDDDiag_pstat) /
                   length(idsEligMDDDiag_pstat))

  dat <- set_epi(dat, "mdd.PrEPeligat", at,length(which(prepElig == 1)))
  dat <- set_epi(dat, "mdd.eligscrnat.numpind", at, length(idsEligMDDDiag_pind))
  dat <- set_epi(dat, "mdd.diagat.numpind", at, length(idsMDDDiag_pind))
  dat <- set_epi(dat, "mdd.diagat.prpeligpind", at, length(idsMDDDiag_pind) /
                   length(idsEligMDDDiag_pind))

  dat <- set_epi(dat, "mdd.diagat.numall", at, length(idsMDDDiag))
  dat <- set_epi(dat, "mdd.diagat.numhiv0", at,
                 length(intersect(idsMDDDiag, which(status == 0))))
  dat <- set_epi(dat, "mdd.diagat.numhiv1", at,
                 length(intersect(idsMDDDiag, which(status == 1))))
  dat <- set_epi(dat, "mdd.diagat.numfirst", at, length(idsMDDDiag_first))

  dat <- set_epi(dat, "mdd.diagever.numall", at, length(which(mdd.diag.status == 1)))
  dat <- set_epi(dat, "mdd.diagever.numhiv0", at,
                 length(which(mdd.diag.status == 1 & status == 0)))
  dat <- set_epi(dat, "mdd.diagever.numhiv1", at,
                 length(which(mdd.diag.status == 1 & status == 0)))

  dat <- set_epi(dat, "mdd.diagever.prpmddall", at, length(which(mdd.diag.status == 1)) /
                   length(which(mdd.status == 1)))
  dat <- set_epi(dat, "mdd.diagever.prpmddhiv0", at,
                 length(which(mdd.diag.status == 1 & mdd.status == 1 & status == 0)) /
                   length(which(mdd.status == 1 & status == 0)))
  dat <- set_epi(dat, "mdd.diagever.prpmddhiv1", at,
                 length(which(mdd.diag.status == 1 & mdd.status == 1 & status == 1)) /
                   length(which(mdd.status == 1 & status == 1)))

  dat <- set_epi(dat, "mdd.diagever.numscrn.med", at,
                 median(mdd.scrn.counter[which(mdd.diag.status==1)], na.rm = T))
  dat <- set_epi(dat, "mdd.diagever.numscrn.min", at,
                 suppressWarnings(min(mdd.scrn.counter[which(mdd.diag.status==1)], na.rm = T)))
  dat <- set_epi(dat, "mdd.diagever.numscrn.max", at,
                 suppressWarnings(max(mdd.scrn.counter[which(mdd.diag.status==1)], na.rm = T)))


  #tx initiation
  dat <- set_epi(dat, "mdd.eligtxstart.numtxreg", at, length(idsEligMDDtxInit_txreg))
  dat <- set_epi(dat, "mdd.txstart.numtxreg", at, length(idsMDDtxInit_txreg))
  dat <- set_epi(dat, "mdd.txstart.prpeligtxreg", at, length(idsMDDtxInit_txreg) /
                   length(idsEligMDDtxInit_txreg))

  dat <- set_epi(dat, "mdd.eligtxstart.numtxintv", at, length(idsEligMDDtxInit_txintv))
  dat <- set_epi(dat, "mdd.txstart.numtxintv", at, length(idsMDDtxInit_txintv))
  dat <- set_epi(dat, "mdd.txstart.prpeligtxintv", at, length(idsMDDtxInit_txintv) /
                   length(idsEligMDDtxInit_txintv))

  dat <- set_epi(dat, "mdd.txstart.numall", at, length(idsMDDtxInit))

  dat <- set_epi(dat, "mdd.txcurr.numall", at,
                 length(which(mdd.tx.status == 1 &
                                (mde.roleimpair == 1 | mde.symseverity > 1))) )
  dat <- set_epi(dat, "mdd.txelig.numall", at,
                 length(which(mdeCurr == 1 & tsinceCurrMDEStarted >=2 &
                                (mde.roleimpair == 1 | mde.symseverity > 1) &
                                mdd.diag.status == 1)))

  dat <- set_epi(dat, "mdd.txcurr.prpactiveeverdiag", at,
                 length(which(mdd.tx.status == 1 &
                          (mde.roleimpair == 1 | mde.symseverity > 1))) /
                   length(which(mdeCurr == 1 & tsinceCurrMDEStarted >=2 &
                                  (mde.roleimpair == 1 | mde.symseverity > 1) &
                                     mdd.diag.status == 1)))

  dat <- set_epi(dat, "mdd.txcurr.numrace1", at,
                 length(which(mdd.tx.status == 1 &
                                (mde.roleimpair == 1 | mde.symseverity > 1) &
                                race == 1)))
  dat <- set_epi(dat, "mdd.txcurr.numrace2", at,
                 length(which(mdd.tx.status == 1 &
                                (mde.roleimpair == 1 | mde.symseverity > 1) &
                                race == 2)))
  dat <- set_epi(dat, "mdd.txcurr.numrace3", at,
                 length(which(mdd.tx.status == 1 &
                                (mde.roleimpair == 1 | mde.symseverity > 1) &
                                race == 3)))

  dat <- set_epi(dat, "mdd.txelig.numrace1", at,
                 length(which(mdeCurr == 1 & tsinceCurrMDEStarted >=2 &
                                (mde.roleimpair == 1 | mde.symseverity > 1) &
                                mdd.diag.status == 1 & race == 1)))
  dat <- set_epi(dat, "mdd.txelig.numrace2", at,
                 length(which(mdeCurr == 1 & tsinceCurrMDEStarted >=2 &
                                (mde.roleimpair == 1 | mde.symseverity > 1) &
                                mdd.diag.status == 1 & race == 1)))
  dat <- set_epi(dat, "mdd.txelig.numrace3", at,
                 length(which(mdeCurr == 1 & tsinceCurrMDEStarted >=2 &
                                (mde.roleimpair == 1 | mde.symseverity > 1) &
                                mdd.diag.status == 1 & race == 1)))

  dat <- set_epi(dat, "mdd.txcurr.prpactiveeverdiag.race1", at,
                 length(which(mdd.tx.status == 1 &
                                (mde.roleimpair == 1 | mde.symseverity > 1) &
                                race == 1)) /
                   length(which(mdeCurr == 1 & tsinceCurrMDEStarted >=2 &
                                  (mde.roleimpair == 1 | mde.symseverity > 1) &
                                  mdd.diag.status == 1 & race == 1)))
  dat <- set_epi(dat, "mdd.txcurr.prpactiveeverdiag.race2", at,
                 length(which(mdd.tx.status == 1 &
                                (mde.roleimpair == 1 | mde.symseverity > 1) &
                                race == 2)) /
                   length(which(mdeCurr == 1 & tsinceCurrMDEStarted >=2 &
                                  (mde.roleimpair == 1 | mde.symseverity > 1) &
                                  mdd.diag.status == 1 & race == 2)))
  dat <- set_epi(dat, "mdd.txcurr.prpactiveeverdiag.race3", at,
                 length(which(mdd.tx.status == 1 &
                                (mde.roleimpair == 1 | mde.symseverity > 1) &
                                race == 3)) /
                   length(which(mdeCurr == 1 & tsinceCurrMDEStarted >=2 &
                                  (mde.roleimpair == 1 | mde.symseverity > 1) &
                                  mdd.diag.status == 1 & race == 3)))

  dat <- set_epi(dat, "mdd.txinittime.min", at,
                suppressWarnings(min(mdd.txinit.time[idsMDDtxInit], na.rm = T)))
  dat <- set_epi(dat, "mdd.txinittime.med", at,
                 median(mdd.txinit.time[idsMDDtxInit], na.rm = T))
  dat <- set_epi(dat, "mdd.txinittime.max", at,
                 suppressWarnings(max(mdd.txinit.time[idsMDDtxInit], na.rm = T)))

  dat <- set_epi(dat, "mdd.mhcare.numpcvdneed", at,
                 length(which(
                   mdeCurr == 1 & tsinceCurrMDEStarted >= 2 &
                     (mde.symseverity >= 2 |
                        (mde.symseverity == 1 & mde.roleimpair == 1)))))
  dat <- set_epi(dat, "mdd.mhcare.numrcvanytx", at,
                 length(which(
                   mdeCurr == 1 & tsinceCurrMDEStarted >= 2 &
                     (mde.symseverity >= 2 |
                        (mde.symseverity == 1 & mde.roleimpair == 1)) &
                     mdd.tx.status == 1)))
  dat <- set_epi(dat, "mdd.mhcare.numrcvminadeqtx", at,
                 length(which(
                   mdeCurr == 1 & tsinceCurrMDEStarted >= 2 &
                     (mde.symseverity >= 2 |
                        (mde.symseverity == 1 & mde.roleimpair == 1))  &
                     mdd.tx.status == 1 & mdd.currtx.tottime >= minadeq.txdur)))

  dat <- set_epi(dat, "mdd.mhcare.prppcvd_rcvanytx", at,
                 length(which(
                   mdeCurr == 1 & tsinceCurrMDEStarted >= 2 &
                     (mde.symseverity >= 2 |
                        (mde.symseverity == 1 & mde.roleimpair == 1)) &
                     mdd.tx.status == 1))/
                   length(which(
                     mdeCurr == 1 & tsinceCurrMDEStarted >= 2 &
                       (mde.symseverity >= 2 |
                          (mde.symseverity == 1 & mde.roleimpair == 1)))))
  dat <- set_epi(dat, "mdd.mhcare.prppcvd_rcvminadeqtx", at,
                 length(which(
                   mdeCurr == 1 & tsinceCurrMDEStarted >= 2 &
                     (mde.symseverity >= 2 |
                        (mde.symseverity == 1 & mde.roleimpair == 1))  &
                     mdd.tx.status == 1 & mdd.currtx.tottime >= minadeq.txdur)) /
                   length(which(
                     mdeCurr == 1 & tsinceCurrMDEStarted >= 2 &
                       (mde.symseverity >= 2 |
                          (mde.symseverity == 1 & mde.roleimpair == 1)))))

  #mde stops
  dat <- set_epi(dat, "mde.stop.numall", at, length(idsMDEStopped))
  dat <- set_epi(dat, "mde.stop.numspontres", at, length(idsMDEStopped_spontres))
  dat <- set_epi(dat, "mde.stop.numtxremiss", at, length(idsMDEStopped_txremiss))
  dat <- set_epi(dat, "mde.stop.numendmainte", at, length(idsMDEStopped_endmainte))
  dat <- set_epi(dat, "mde.stop.numdth", at, length(idsMDEStopped_dth))

  dat <- set_epi(dat, "mde.stop.spontres.prpactiveuntx", at,
                 length(idsMDEStopped_spontres) /
                   length(which(active==1 & mdeCurr == 1 & tsinceCurrMDEStarted >= 2 &
                                  mdd.tx.status == 0)))

  dat <- set_epi(dat, "mde.epidur.min", at,
                 suppressWarnings(min(mdeCurr.dur, na.rm = T)))
  dat <- set_epi(dat, "mde.epidur.median", at,
                 median(mdeCurr.dur, na.rm = T))
  dat <- set_epi(dat, "mde.epidur.max", at,
                 suppressWarnings(max(mdeCurr.dur, na.rm = T)))

  #mde tx stops
  dat <- set_epi(dat, "mdd.txstop.numall", at, length(idsMDDtxStopped))
  dat <- set_epi(dat, "mdd.txstop.numremiss", at, length(idsMDDtxStopped_remiss))
  dat <- set_epi(dat, "mdd.txstop.numendmainte", at, length(idsMDDtxStopped_endmainte))
  dat <- set_epi(dat, "mdd.txstop.numltfu", at, length(idsMDDtxStopped_ltfu))
  dat <- set_epi(dat, "mdd.txstop.numdth", at, length(idsMDDtxStopped_dth))


  return(dat)
}






#' @title MDD-related suicide attempts
#'
#' @description Function to determine and track suicide attempts related to major
#' depression disorder among network nodes (untreated/inadequately treated mdd ONLY)
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'



mddsuitry_msm <- function(dat, at){

  ## Input
  # attributes
  active            <- get_attr(dat, "active")
  age               <- get_attr(dat, "age")
  race              <- get_attr(dat, "race")
  status            <- get_attr(dat, "status")
  mdd.status        <- get_attr(dat, "mdd.status")
  mdeCurr           <- get_attr(dat, "mdeCurr")
  mdeCurr.starttime <- get_attr(dat, "mdeCurr.starttime")
  mde.symseverity   <- get_attr(dat, "mde.symseverity")

  mdd.tx.status     <- get_attr(dat, "mdd.tx.status")
  mdd.lasttx.cumlt  <- get_attr(dat, "mdd.lasttx.cumlt")
  mdd.currtx.tottime<- get_attr(dat, "mdd.currtx.tottime")

  mdd.suitry.status <- get_attr(dat, "mdd.suitry.status")
  mdd.suitry.counter<- get_attr(dat, "mdd.suitry.counter")
  mdd.suitry.time   <- get_attr(dat, "mdd.suitry.time")


  # parameters
  mdd.suitry.prob       <- get_param(dat, "mdd.suitry.prob")



  #attempts
  tsinceMDEstarted <- at - mdeCurr.starttime

  idsEligMDDSuiTry <- which(
    active == 1 & (
      (mdeCurr == 1 & tsinceMDEstarted >= 2 & mde.symseverity > 2 & mdd.tx.status == 0) |
        (mdeCurr == 1 & mdd.tx.status == 1 & mdd.currtx.tottime < 8))
  )

  prob.suitry <- 1 - ((1-mdd.suitry.prob)^(1/52))
  mdd_suitried <- as.numeric(runif(length(idsEligMDDSuiTry)) < prob.suitry)
  idsMDDSuiTried <- idsEligMDDSuiTry[mdd_suitried == 1]


  #update attributes
  if(length(idsMDDSuiTried) > 0) {
    mdd.suitry.status[idsMDDSuiTried]   <- 1
    mdd.suitry.counter[idsMDDSuiTried] <- ifelse(is.na(mdd.suitry.counter[idsMDDSuiTried]),
                                                 1, mdd.suitry.counter[idsMDDSuiTried] + 1)
    mdd.suitry.time[idsMDDSuiTried] <- at
  }

  dat <- set_attr(dat, "mdd.suitry.status", mdd.suitry.status)
  dat <- set_attr(dat, "mdd.suitry.counter", mdd.suitry.counter)
  dat <- set_attr(dat, "mdd.suitry.time", mdd.suitry.time)



  #epi tracker
  dat <- set_epi(dat, "mdd.eligsuitry.numall", at, length(idsEligMDDSuiTry))
  dat <- set_epi(dat, "mdd.suitry.numall", at, length(idsMDDSuiTried))
  dat <- set_epi(dat, "mdd.suitry.numhiv0", at, length(intersect(idsMDDSuiTried, which(status == 0))))
  dat <- set_epi(dat, "mdd.suitry.numhiv1", at, length(intersect(idsMDDSuiTried, which(status == 1))))

  dat <- set_epi(dat, "mdd.suitry.numsev1", at,
                 length(intersect(idsMDDSuiTried, which(mde.symseverity == 1))))
  dat <- set_epi(dat, "mdd.suitry.numsev2", at,
                 length(intersect(idsMDDSuiTried, which(mde.symseverity == 2))))
  dat <- set_epi(dat, "mdd.suitry.numsev3", at,
                 length(intersect(idsMDDSuiTried, which(mde.symseverity == 3))))
  dat <- set_epi(dat, "mdd.suitry.numsev4", at,
                 length(intersect(idsMDDSuiTried, which(mde.symseverity == 4))))

  dat <- set_epi(dat, "mdd.suitry.nummin", at,
                 suppressWarnings(min(mdd.suitry.counter[idsMDDSuiTried], na.rm = T)))
  dat <- set_epi(dat, "mdd.suitry.nummed", at,
                 median(mdd.suitry.counter[idsMDDSuiTried], na.rm = T))
  dat <- set_epi(dat, "mdd.suitry.nummax", at,
                 suppressWarnings(max(mdd.suitry.counter[idsMDDSuiTried], na.rm = T)))

  return(dat)
}
