
#' @title Prevalence Calculations within Time Steps
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Summary statistic calculations are of two broad forms: prevalence and
#' incidence. This function establishes the summary statistic vectors for both
#' prevalence and incidence at time 1, and then calculates the prevalence
#' statistics for times 2 onward. Incidence statistics (e.g., number of new
#' infections or deaths) are calculated within the modules as they depend on
#' vectors that are not stored external to the module.
#'
#' @return
#' This function returns the \code{dat} object with an updated summary of current
#' attributes stored in \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
prevalence_msm <- function(dat, at) {

  ## Variables

  # Attributes

  active <- dat$attr$active
  race <- dat$attr$race
  status <- dat$attr$status
  prepStat <- dat$attr$prepStat
  prepElig <- dat$attr$prepElig
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT
  syphilis <- dat$attr$syphilis
  inf.role <- dat$attr$inf.role
  inf.time <- dat$attr$inf.time
  rGC.infTime <- dat$attr$rGC.infTime
  uGC.infTime <- dat$attr$uGC.infTime
  rCT.infTime <- dat$attr$rCT.infTime
  uCT.infTime <- dat$attr$uCT.infTime
  syph.infTime <- dat$attr$syph.infTime
  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  rCT.sympt <- dat$attr$rCT.sympt
  uCT.sympt <- dat$attr$uCT.sympt
  stage.syph <- dat$attr$stage.syph
  stage.prim.sympt <- dat$attr$stage.prim.sympt
  stage.seco.sympt <- dat$attr$stage.seco.sympt
  stage.earlat.sympt <- dat$attr$stage.earlat.sympt
  stage.latelat.sympt <- dat$attr$stage.latelat.sympt
  stage.latelatelat.sympt <- dat$attr$stage.latelatelat.sympt
  stage.tert.sympt <- dat$attr$stage.tert.sympt

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  if (at == 1) {
    dat$epi$num <- rNA
    dat$epi$num.B <- rNA
    dat$epi$num.W <- rNA
    dat$epi$s.num <- rNA
    dat$epi$i.num <- rNA
    dat$epi$i.num.B <- rNA
    dat$epi$i.num.W <- rNA
    dat$epi$i.prev <- rNA
    dat$epi$i.prev.B <- rNA
    dat$epi$i.prev.W <- rNA
    dat$epi$incid <- rNA
    dat$epi$ir100 <- rNA

    dat$epi$prepCurr <- rNA
    dat$epi$prepCov <- rNA
    dat$epi$prepElig <- rNA
    dat$epi$prepStart <- rNA
    dat$epi$i.num.prep0 <- rNA
    dat$epi$i.num.prep1 <- rNA
    
    dat$epi$time.hivneg <- rep(0, nsteps)
    dat$epi$time.on.prep <- rep(0, nsteps)
    dat$epi$time.off.prep <- rep(0, nsteps)
    dat$epi$stage.time.ar.ndx <- rep(0, nsteps)
    dat$epi$stage.time.ar.dx <- rep(0, nsteps)
    dat$epi$stage.time.ar.art <- rep(0, nsteps)
    dat$epi$stage.time.af.ndx <- rep(0, nsteps)
    dat$epi$stage.time.af.dx <- rep(0, nsteps)
    dat$epi$stage.time.af.art <- rep(0, nsteps)
    dat$epi$stage.time.chronic.ndx <- rep(0, nsteps)
    dat$epi$stage.time.chronic.dx <- rep(0, nsteps)
    dat$epi$stage.time.chronic.art <- rep(0, nsteps)
    dat$epi$stage.time.aids.ndx <- rep(0, nsteps)
    dat$epi$stage.time.aids.dx <- rep(0, nsteps)
    dat$epi$stage.time.aids.art <- rep(0, nsteps)
    
    dat$epi$hivtests.prep <- rep(0, nsteps)
    dat$epi$hivtests.nprep <- rep(0, nsteps)
    dat$epi$totalhivtests.prep <- rep(0, nsteps)
    dat$epi$totalhivtests <- rep(0, nsteps)
    
    dat$epi$totalrGCsympttests <- rep(0, nsteps)
    dat$epi$totaluGCsympttests <- rep(0, nsteps)
    dat$epi$totalGCsympttests <- rep(0, nsteps)
    dat$epi$totalrCTsympttests <- rep(0, nsteps)
    dat$epi$totaluCTsympttests <- rep(0, nsteps)
    dat$epi$totalCTsympttests <- rep(0, nsteps)
    dat$epi$totalsyphsympttests <- rep(0, nsteps)
    dat$epi$totalstisympttests <- rep(0, nsteps)
    
    dat$epi$totalrGCasympttests <- rep(0, nsteps)
    dat$epi$totaluGCasympttests <- rep(0, nsteps)
    dat$epi$totalGCasympttests <- rep(0, nsteps)
    dat$epi$totalrCTasympttests <- rep(0, nsteps)
    dat$epi$totaluCTasympttests <- rep(0, nsteps)
    dat$epi$totalCTasympttests <- rep(0, nsteps)
    dat$epi$totalsyphasympttests <- rep(0, nsteps)
    dat$epi$totalstiasympttests <- rep(0, nsteps)
    
    dat$epi$totalrGCasympttests.prep <- rep(0, nsteps)
    dat$epi$totaluGCasympttests.prep <- rep(0, nsteps)
    dat$epi$totalGCasympttests.prep <- rep(0, nsteps)
    dat$epi$totalrCTasympttests.prep <- rep(0, nsteps)
    dat$epi$totaluCTasympttests.prep <- rep(0, nsteps)
    dat$epi$totalCTasympttests.prep <- rep(0, nsteps)
    dat$epi$totalsyphasympttests.prep <- rep(0, nsteps)
    dat$epi$totalstiasympttests.prep <- rep(0, nsteps)

    dat$epi$prev.rgc <- rNA
    dat$epi$prev.ugc <- rNA
    dat$epi$prev.gc <- rNA
    dat$epi$prev.gc.sympt <- rNA
    dat$epi$prev.gc.dual <- rNA

    dat$epi$prev.rct <- rNA
    dat$epi$prev.uct <- rNA
    dat$epi$prev.ct <- rNA
    dat$epi$prev.ct.sympt <- rNA
    dat$epi$prev.ct.dual <- rNA

    dat$epi$prev.rgcct <- rNA
    dat$epi$prev.ugcct <- rNA
    
    dat$epi$prev.syph <- rNA
    dat$epi$prev.stage.prim <- rNA
    dat$epi$prev.stage.seco <- rNA
    dat$epi$prev.stage.earlat <- rNA
    dat$epi$prev.stage.latelat <- rNA
    dat$epi$prev.stage.latelatelat <- rNA
    dat$epi$prev.stage.tert <- rNA
    dat$epi$prev.earlysyph <- rNA
    dat$epi$prev.latesyph <- rNA
    
    #HIV Coinfection
    dat$epi$prev.syph.hivneg <- rNA
    dat$epi$prev.syph.hivpos <- rNA
    
    dat$epi$prev.gc.hivneg <- rNA
    dat$epi$prev.gc.hivpos <- rNA
    
    dat$epi$prev.ct.hivneg <- rNA
    dat$epi$prev.ct.hivpos <- rNA

    dat$epi$incid.rgc <- rNA
    dat$epi$incid.ugc <- rNA
    dat$epi$incid.gc <- rNA
    dat$epi$incid.rct <- rNA
    dat$epi$incid.uct <- rNA
    dat$epi$incid.ct <- rNA
    dat$epi$incid.syph <- rNA
    
    dat$epi$ir100.rgc <- rNA
    dat$epi$ir100.ugc <- rNA
    dat$epi$ir100.gc <- rNA
    dat$epi$ir100.rct <- rNA
    dat$epi$ir100.uct <- rNA
    dat$epi$ir100.ct <- rNA
    dat$epi$ir100.syph <- rNA
    
    dat$epi$ir100.sti <- rNA
    dat$epi$incid.gcct.prep <- rNA
    dat$epi$incid.syph.prep <- rNA
    
    #PAF
    dat$epi$sti_paf <- rNA
    dat$epi$sti_u_paf <- rNA
    dat$epi$sti_r_paf <- rNA
    dat$epi$sti_syph_paf <- rNA
    dat$epi$sti_u_sympt_paf <- rNA
    dat$epi$sti_u_asympt_paf <- rNA
    dat$epi$sti_r_sympt_paf <- rNA
    dat$epi$sti_r_asympt_paf <- rNA
    dat$epi$sti_syph_sympt_paf <- rNA
    dat$epi$sti_syph_asympt_paf <- rNA
    dat$epi$sti_hiv_sum <- rNA
    dat$epi$sti_u_hiv_sum <- rNA
    dat$epi$sti_r_hiv_sum <- rNA
    dat$epi$sti_syph_hiv_sum <- rNA
    dat$epi$hiv_sum <- rNA
    dat$epi$sti_u_sympt_hiv_sum <- rNA
    dat$epi$sti_u_asympt_hiv_sum <- rNA
    dat$epi$sti_r_sympt_hiv_sum <- rNA
    dat$epi$sti_r_asympt_hiv_sum <- rNA
    dat$epi$sti_syph_sympt_hiv_sum <- rNA
    dat$epi$sti_syph_asympt_hiv_sum <- rNA

    dat$epi$recov.rgc <- rNA
    dat$epi$recov.ugc <- rNA
    dat$epi$recov.rct <- rNA
    dat$epi$recov.uct <- rNA
    dat$epi$recov.prim.syph <- rNA
    dat$epi$recov.seco.syph <- rNA
    dat$epi$recov.earlat.syph <- rNA
    dat$epi$recov.latelat.syph <- rNA
    dat$epi$recov.latelatelat.syph <- rNA
    dat$epi$recov.tert.syph <- rNA
    dat$epi$recov.earlysyph <- rNA
    dat$epi$recov.latesyph <- rNA
    dat$epi$recov.syphilis <- rNA
    
    dat$epi$trans.main <- rNA
    dat$epi$trans.casl <- rNA
    dat$epi$trans.inst <- rNA

    dat$epi$txGC <- rNA
    dat$epi$txCT <- rNA
    dat$epi$txsyph <- rNA
    
    dat$epi$stiactiveind <- rNA  
    dat$epi$recentpartners <- rNA
    dat$epi$recentSTI <- rNA
    dat$epi$newpartner <- rNA
    dat$epi$concurrpart <- rNA
    dat$epi$partnersti <- rNA
    dat$epi$uai.nmain <- rNA
    dat$epi$uai.any <- rNA
    
    dat$epi$eptpartelig <- rNA
    dat$epi$eptprovided <- rNA
    dat$epi$eptTx <- rNA
    dat$epi$eptprop_provided <- rNA
    dat$epi$eptprop_tx <- rNA
    
    ##########
  }

  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$num.B[at] <- sum(race == "B", na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == "W", na.rm = TRUE)
  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$i.num.B[at] <- sum(status == 1 & race == "B", na.rm = TRUE)
  dat$epi$i.num.W[at] <- sum(status == 1 & race == "W", na.rm = TRUE)
  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.B[at] <- dat$epi$i.num.B[at] / dat$epi$num.B[at]
  dat$epi$i.prev.W[at] <- dat$epi$i.num.W[at] / dat$epi$num.W[at]
  dat$epi$ir100[at] <- (dat$epi$incid[at] / sum(status == 0, na.rm = TRUE)) * 5200

  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  dat$epi$i.num.prep0[at] <- sum((is.na(prepStat) | prepStat == 0) & status == 1, na.rm = TRUE)
  dat$epi$i.num.prep1[at] <- sum(prepStat == 1 & status == 1, na.rm = TRUE)
  dat$epi$i.prev.prep0[at] <- dat$epi$i.num.prep0[at] /
    sum((is.na(prepStat) | prepStat == 0), na.rm = TRUE)
  if (at == 1) {
    dat$epi$i.prev.prep1[1] <- 0
  } else {
    dat$epi$i.prev.prep1[at] <- dat$epi$i.num.prep1[at] / sum(prepStat == 1, na.rm = TRUE)
  }

  dat$epi$time.hivneg[at] <- sum(dat$attr$time.hivneg, na.rm = TRUE)
  dat$epi$time.on.prep[at] <- sum(dat$attr$time.on.prep, na.rm = TRUE)
  dat$epi$time.off.prep[at] <- sum(dat$attr$time.off.prep, na.rm = TRUE)
  dat$epi$stage.time.ar.ndx[at] <- sum(dat$attr$stage.time.ar.ndx, na.rm = TRUE)
  dat$epi$stage.time.ar.dx[at] <- sum(dat$attr$stage.time.ar.dx, na.rm = TRUE)
  dat$epi$stage.time.ar.art[at] <- sum(dat$attr$stage.time.ar.art, na.rm = TRUE)
  dat$epi$stage.time.af.ndx[at] <- sum(dat$attr$stage.time.af.ndx, na.rm = TRUE)
  dat$epi$stage.time.af.dx[at] <- sum(dat$attr$stage.time.af.dx, na.rm = TRUE)
  dat$epi$stage.time.af.art[at] <- sum(dat$attr$stage.time.af.art, na.rm = TRUE)
  dat$epi$stage.time.chronic.ndx[at] <- sum(dat$attr$stage.time.chronic.ndx, na.rm = TRUE)
  dat$epi$stage.time.chronic.dx[at] <- sum(dat$attr$stage.time.chronic.dx, na.rm = TRUE)
  dat$epi$stage.time.chronic.art[at] <- sum(dat$attr$stage.time.chronic.art, na.rm = TRUE) 
  dat$epi$stage.time.aids.ndx[at] <- sum(dat$attr$stage.time.aids.ndx, na.rm = TRUE)
  dat$epi$stage.time.aids.dx[at] <- sum(dat$attr$stage.time.aids.dx, na.rm = TRUE)
  dat$epi$stage.time.aids.art[at] <- sum(dat$attr$stage.time.aids.art, na.rm = TRUE)
  dat$epi$stisympttests[at] <- sum(dat$epi$syphsympttests[at], dat$epi$CTsympttests[at], dat$epi$GCsympttests[at], na.rm = TRUE)
  dat$epi$stiasympttests[at] <- sum(dat$epi$syphasympttests[at], dat$epi$CTasympttests[at], dat$epi$GCasympttests[at], na.rm = TRUE)
  dat$epi$stiasympttests.prep[at] <- sum(dat$epi$syphasympttests.prep[at], dat$epi$CTasympttests.prep[at], dat$epi$GCasympttests.prep[at], na.rm = TRUE)
  dat$epi$stiasympttests.pos[at] <- sum(dat$epi$syphasympttests.pos[at], dat$epi$CTasympttests.pos[at], dat$epi$GCasympttests.pos[at], na.rm = TRUE)
  dat$epi$stiasympttests.pos.prep[at] <- sum(dat$epi$syphasympttests.pos.prep[at], dat$epi$CTasympttests.pos.prep[at], dat$epi$GCasympttests.pos.prep[at], na.rm = TRUE)
  
  # Number of tests: Total resets at PrEP or STI testing start time
  if (at == 2 | at == dat$param$prep.start | at == dat$param$stitest.start) {
      
      dat$epi$totalhivtests[at] <- sum(dat$epi$hivtests.prep[at], dat$epi$hivtests.nprep[at])
      dat$epi$totalhivtests.prep[at] <- dat$epi$hivtests.prep[at]
      
      dat$epi$totalrGCsympttests[at] <- dat$epi$rGCsympttests[at]
      dat$epi$totaluGCsympttests[at] <- dat$epi$uGCsympttests[at]
      dat$epi$totalGCsympttests[at] <- dat$epi$GCsympttests[at]
      dat$epi$totalrCTsympttests[at] <- dat$epi$rCTsympttests[at]
      dat$epi$totaluCTsympttests[at] <- dat$epi$uCTsympttests[at]
      dat$epi$totalCTsympttests[at] <- dat$epi$CTsympttests[at]
      dat$epi$totalsyphsympttests[at] <- dat$epi$syphsympttests[at]
      dat$epi$totalstisympttests[at] <- sum(dat$epi$syphsympttests[at], dat$epi$CTsympttests[at], dat$epi$GCsympttests[at], na.rm = TRUE)
      
      dat$epi$totalrGCasympttests[at] <- dat$epi$rGCasympttests[at]
      dat$epi$totaluGCasympttests[at] <- dat$epi$uGCasympttests[at]
      dat$epi$totalGCasympttests[at] <- dat$epi$GCasympttests[at]
      dat$epi$totalrGCasympttests.pos[at] <- dat$epi$rGCasympttests.pos[at]
      dat$epi$totaluGCasympttests.pos[at] <- dat$epi$uGCasympttests.pos[at]
      dat$epi$totalGCasympttests.pos[at] <- dat$epi$GCasympttests.pos[at]
      
      dat$epi$totalrCTasympttests[at] <- dat$epi$rCTasympttests[at]
      dat$epi$totaluCTasympttests[at] <- dat$epi$uCTasympttests[at]
      dat$epi$totalCTasympttests[at] <- dat$epi$CTsympttests[at]
      dat$epi$totalrCTasympttests.pos[at] <- dat$epi$rCTasympttests.pos[at]
      dat$epi$totaluCTasympttests.pos[at] <- dat$epi$uCTasympttests.pos[at]
      dat$epi$totalCTasympttests.pos[at] <- dat$epi$CTasympttests.pos[at]
      
      dat$epi$totalsyphasympttests[at] <- dat$epi$syphasympttests[at]
      dat$epi$totalsyphasympttests.pos[at] <- dat$epi$syphasympttests.pos[at]
      
      dat$epi$totalstiasympttests[at] <- sum(dat$epi$syphasympttests[at], dat$epi$CTasympttests[at], dat$epi$GCasympttests[at], na.rm = TRUE)
      dat$epi$totalstiasympttests.pos[at] <- sum(dat$epi$syphasympttests.pos[at], dat$epi$CTasympttests.pos[at], dat$epi$GCasympttests.pos[at], na.rm = TRUE)
      
      dat$epi$totalrGCasympttests.prep[at] <- dat$epi$rGCasympttests.prep[at]
      dat$epi$totaluGCasympttests.prep[at] <- dat$epi$uGCasympttests.prep[at]
      dat$epi$totalGCasympttests.prep[at] <- dat$epi$GCasympttests.prep[at]
      dat$epi$totalrGCasympttests.pos.prep[at] <- dat$epi$rGCasympttests.pos.prep[at]
      dat$epi$totaluGCasympttests.pos.prep[at] <- dat$epi$uGCasympttests.pos.prep[at]
      dat$epi$totalGCasympttests.pos.prep[at] <- dat$epi$GCasympttests.pos.prep[at]
      
      dat$epi$totalrCTasympttests.prep[at] <- dat$epi$rCTasympttests.prep[at]
      dat$epi$totaluCTasympttests.prep[at] <- dat$epi$uCTasympttests.prep[at]
      dat$epi$totalCTasympttests.prep[at] <- dat$epi$CTasympttests.prep[at]
      dat$epi$totalrCTasympttests.pos.prep[at] <- dat$epi$rCTasympttests.pos.prep[at]
      dat$epi$totaluCTasympttests.pos.prep[at] <- dat$epi$uCTasympttests.pos.prep[at]
      dat$epi$totalCTasympttests.pos.prep[at] <- dat$epi$CTasympttests.pos.prep[at]
      
      dat$epi$totalsyphasympttests.prep[at] <- dat$epi$syphasympttests.prep[at]
      dat$epi$totalsyphasympttests.pos.prep[at] <- dat$epi$syphasympttests.pos.prep[at]
      
      dat$epi$totalstiasympttests.prep[at] <- sum(dat$epi$syphasympttests.prep[at], dat$epi$CTasympttests.prep[at], dat$epi$GCasympttests.prep[at], na.rm = TRUE)
      dat$epi$totalstiasympttests.pos.prep[at] <- sum(dat$epi$syphasympttests.pos.prep[at], dat$epi$CTasympttests.pos.prep[at], dat$epi$GCasympttests.pos.prep[at], na.rm = TRUE)
      
  }
  if ((at > 2 & at < dat$param$prep.start) | (at > dat$param$prep.start) |
      (at > 2 & at < dat$param$stitest.start) | (at > dat$param$stitest.start)) {
      
      dat$epi$totalhivtests[at] <- dat$epi$hivtests.nprep[at] + dat$epi$hivtests.prep[at] + dat$epi$totalhivtests[at - 1]
      dat$epi$totalhivtests.prep[at] <- dat$epi$hivtests.prep[at] + dat$epi$totalhivtests.prep[at - 1]
      
      dat$epi$totalrGCsympttests[at] <- dat$epi$rGCsympttests[at] + dat$epi$totalrGCsympttests[at - 1] 
      dat$epi$totaluGCsympttests[at] <- dat$epi$uGCsympttests[at] + dat$epi$totaluGCsympttests[at - 1]
      dat$epi$totalGCsympttests[at] <- dat$epi$GCsympttests[at] + dat$epi$totalGCsympttests[at - 1]
      dat$epi$totalrCTsympttests[at] <- dat$epi$rCTsympttests[at] + dat$epi$totalrCTsympttests[at - 1]
      dat$epi$totaluCTsympttests[at] <- dat$epi$uCTsympttests[at] + dat$epi$totaluCTsympttests[at - 1]
      dat$epi$totalCTsympttests[at] <- dat$epi$CTsympttests[at] + dat$epi$totalCTsympttests[at - 1]
      dat$epi$totalsyphsympttests[at] <- dat$epi$syphsympttests[at] + dat$epi$totalsyphsympttests[at - 1]
      dat$epi$totalstisympttests[at] <- dat$epi$stisympttests[at] + dat$epi$totalstisympttests[at - 1]
      
      dat$epi$totalrGCasympttests[at] <- dat$epi$rGCasympttests[at] + dat$epi$totalrGCasympttests[at - 1]
      dat$epi$totaluGCasympttests[at] <- dat$epi$uGCasympttests[at] + dat$epi$totaluGCasympttests[at - 1]
      dat$epi$totalGCasympttests[at] <- dat$epi$GCasympttests[at] + dat$epi$totalGCasympttests[at - 1]
      dat$epi$totalrCTasympttests[at] <- dat$epi$rCTasympttests[at] + dat$epi$totalrCTasympttests[at - 1]
      dat$epi$totaluCTasympttests[at] <- dat$epi$uCTasympttests[at] + dat$epi$totaluCTasympttests[at - 1]
      dat$epi$totalCTasympttests[at] <- dat$epi$CTasympttests[at] + dat$epi$totalCTasympttests[at - 1]
      dat$epi$totalsyphasympttests[at] <- dat$epi$syphasympttests[at] + dat$epi$totalsyphasympttests[at - 1]
      dat$epi$totalstiasympttests[at] <- dat$epi$stiasympttests[at] + dat$epi$totalstiasympttests[at - 1]
      
      dat$epi$totalrGCasympttests.prep[at] <- dat$epi$rGCasympttests.prep[at] + dat$epi$totalrGCasympttests.prep[at - 1]
      dat$epi$totaluGCasympttests.prep[at] <- dat$epi$uGCasympttests.prep[at] + dat$epi$totaluGCasympttests.prep[at - 1]
      dat$epi$totalGCasympttests.prep[at] <- dat$epi$GCasympttests.prep[at] + dat$epi$totalGCasympttests.prep[at - 1]
      dat$epi$totalrCTasympttests.prep[at] <- dat$epi$rCTasympttests.prep[at] + dat$epi$totalrCTasympttests.prep[at - 1] 
      dat$epi$totaluCTasympttests.prep[at] <- dat$epi$uCTasympttests.prep[at] + dat$epi$totaluCTasympttests.prep[at - 1]
      dat$epi$totalCTasympttests.prep[at] <- dat$epi$CTasympttests.prep[at] + dat$epi$totalCTasympttests.prep[at - 1]
      dat$epi$totalsyphasympttests.prep[at] <- dat$epi$syphasympttests.prep[at] + dat$epi$totalsyphasympttests.prep[at - 1]
      dat$epi$totalstiasympttests.prep[at] <- dat$epi$stiasympttests.prep[at] + dat$epi$totalstiasympttests.prep[at - 1]
      
  }
 
  # STI Prevalence
  
  dat$epi$prev.rgc[at] <- sum(rGC == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ugc[at] <- sum(uGC == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.gc[at] <- sum((rGC == 1 | uGC == 1), na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.gc.sympt[at] <- sum((rGC.sympt == 1 | uGC.sympt == 1)) / dat$epi$num[at]
  dat$epi$prev.gc.dual[at] <- sum((rGC == 1 & uGC == 1), na.rm = TRUE) / dat$epi$num[at]

  dat$epi$prev.rct[at] <- sum(rCT == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.uct[at] <- sum(uCT == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ct[at] <- sum((rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ct.sympt[at] <- sum((rCT.sympt == 1 | uCT.sympt == 1)) / dat$epi$num[at]
  dat$epi$prev.ct.dual[at] <- sum((rCT == 1 & uCT == 1), na.rm = TRUE) / dat$epi$num[at]

  dat$epi$prev.rgcct[at] <- sum(rGC == 1 | rCT == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$prev.ugcct[at] <- sum(uGC == 1 | uCT == 1, na.rm = TRUE) / dat$epi$num[at]
  
  dat$epi$prev.stage.incub[at] <- length(which(stage.syph == 1)) / length(which(syphilis == 1))
  dat$epi$prev.stage.prim[at] <- length(which(stage.syph == 2)) / length(which(syphilis == 1))
  dat$epi$prev.stage.incubprim[at] <- length(which(stage.syph == 1 | stage.syph == 2)) / length(which(syphilis == 1))
  dat$epi$prev.stage.seco[at] <- length(which(stage.syph == 3)) / length(which(syphilis == 1))
  dat$epi$prev.stage.earlat[at] <- length(which(stage.syph == 4)) / length(which(syphilis == 1))
  dat$epi$prev.stage.latelat[at] <- length(which(stage.syph == 5)) / length(which(syphilis == 1))
  dat$epi$prev.stage.latelatelat[at] <- length(which(stage.syph == 6)) / length(which(syphilis == 1))
  dat$epi$prev.stage.alllatelat[at] <- length(which(stage.syph %in% c(5,6))) / length(which(syphilis == 1))
  dat$epi$prev.stage.tert[at] <- length(which(stage.syph == 7)) / length(which(syphilis == 1))
  dat$epi$prev.earlysyph[at] <- length(which(stage.syph %in% c(1, 2, 3, 4))) / length(which(syphilis == 1))
  dat$epi$prev.latesyph[at] <- length(which(stage.syph %in% c(5, 6, 7))) / length(which(syphilis == 1))
  dat$epi$prev.syph[at] <- length(which(syphilis == 1)) / dat$epi$num[at]
  dat$epi$prev.primsecosyph[at] <- length(which(stage.syph %in% c(1, 2, 3))) / dat$epi$num[at]
  
  uGC.prev <- which(uGC == 1 & uGC.infTime < at)
  uCT.prev <- which(uCT == 1 & uCT.infTime < at)
  rGC.prev <- which(rGC == 1 & rGC.infTime < at)
  rCT.prev <- which(rCT == 1 & rCT.infTime < at)

  # Prevalence of HIV/STI overlap
  dat$epi$prev.primsecosyph.hivneg[at] <- length(intersect(which(status == 0), which(stage.syph %in% c(1, 2, 3)))) / dat$epi$s.num[at]
  dat$epi$prev.primsecosyph.hivpos[at] <- length(intersect(which(status == 1), which(stage.syph %in% c(1, 2, 3)))) / dat$epi$i.num[at]
  dat$epi$prev.syph.hivneg[at] <- length(intersect(which(status == 0), which(syphilis == 1))) / dat$epi$s.num[at]
  dat$epi$prev.syph.hivpos[at] <- length(intersect(which(status == 1), which(syphilis == 1))) / dat$epi$i.num[at]
  
  dat$epi$prev.gc.hivneg[at] <- length(intersect(which(status == 0), which((rGC == 1 | uGC == 1)))) / dat$epi$s.num[at]
  dat$epi$prev.gc.hivpos[at] <- length(intersect(which(status == 1), which((rGC == 1 | uGC == 1)))) / dat$epi$i.num[at]
  
  dat$epi$prev.ct.hivneg[at] <- length(intersect(which(status == 0), which((rCT == 1 | uCT == 1)))) / dat$epi$s.num[at]
  dat$epi$prev.ct.hivpos[at] <- length(intersect(which(status == 1), which((rCT == 1 | uCT == 1)))) / dat$epi$i.num[at]
  
  dat$epi$prev.hiv.primsecosyphpos[at] <- length(intersect(which(status == 1), which(stage.syph %in% c(1, 2, 3)))) /
                                    length(which(stage.syph %in% c(1, 2, 3)))
  dat$epi$prev.hiv.primsecosyphneg[at] <- length(intersect(which(status == 1), which(stage.syph %in% c(1, 2, 3)))) /
      length(which(stage.syph %in% c(1, 2, 3)))
  
  dat$epi$prev.hiv.syphpos[at] <- length(intersect(which(status == 1), which(syphilis == 1))) /
      length(which(syphilis == 1))
  dat$epi$prev.hiv.syphneg[at] <- length(intersect(which(status == 1), which(syphilis == 1))) /
      length(which(syphilis == 1))
  
  
  dat$epi$prev.hiv.gcpos[at] <- length(intersect(which(status == 1), which((rGC == 1 | uGC == 1)))) /
                                sum((rGC == 1 | uGC == 1), na.rm = TRUE)
  dat$epi$prev.hiv.gcneg[at] <- length(intersect(which(status == 1), which((rGC == 0 & uGC == 0)))) /
      sum((rGC == 0 & uGC == 0), na.rm = TRUE)
  
  dat$epi$prev.hiv.ctpos[at] <- length(intersect(which(status == 1), which((rCT == 1 | uCT == 1)))) / 
                                sum((rCT == 1 | uGC == 1), na.rm = TRUE)
  dat$epi$prev.hiv.ctneg[at] <- length(intersect(which(status == 1), which((rCT == 0 & uCT == 0)))) / 
      sum((rCT == 0 & uGC == 0), na.rm = TRUE)
  
  dat$epi$prev.rgc.hivpos[at] <- length(intersect(which(status == 1), which(rGC == 1))) / dat$epi$i.num[at]
  dat$epi$prev.ugc.hivpos[at] <- length(intersect(which(status == 1), which(uGC == 1))) / dat$epi$i.num[at]
  dat$epi$prev.rct.hivpos[at] <- length(intersect(which(status == 1), which(rCT == 1))) / dat$epi$i.num[at]
  dat$epi$prev.uct.hivpos[at] <- length(intersect(which(status == 1), which(uCT == 1))) / dat$epi$i.num[at]
  
  dat$epi$prev.rgc.hivneg[at] <- length(intersect(which(status == 1), which(rGC == 1))) / dat$epi$s.num[at]
  dat$epi$prev.ugc.hivneg[at] <- length(intersect(which(status == 1), which(uGC == 1))) / dat$epi$s.num[at]
  dat$epi$prev.rct.hivneg[at] <- length(intersect(which(status == 1), which(rCT == 1))) / dat$epi$s.num[at]
  dat$epi$prev.uct.hivneg[at] <- length(intersect(which(status == 1), which(uCT == 1))) / dat$epi$s.num[at]

  # Site-specific STI incidence rates        
  dat$epi$ir100.rgc[at] <- (dat$epi$incid.rgc[at] / sum(rGC == 0, na.rm = TRUE)) * 5200
  dat$epi$ir100.ugc[at] <- (dat$epi$incid.ugc[at] / sum(uGC == 0, na.rm = TRUE)) * 5200
  dat$epi$ir100.gc[at] <- (dat$epi$incid.gc[at] /
                             (sum(rGC == 0, na.rm = TRUE) +
                                sum(uGC == 0, na.rm = TRUE))) * 5200

  dat$epi$ir100.rct[at] <- (dat$epi$incid.rct[at] / sum(rCT == 0, na.rm = TRUE)) * 5200
  dat$epi$ir100.uct[at] <- (dat$epi$incid.uct[at] / sum(uCT == 0, na.rm = TRUE)) * 5200
  dat$epi$ir100.ct[at] <- (dat$epi$incid.ct[at] /
                             (sum(rCT == 0, na.rm = TRUE) +
                                sum(uCT == 0, na.rm = TRUE))) * 5200
  
  dat$epi$ir100.syph[at] <- (dat$epi$incid.syph[at] / sum(syphilis == 0 , na.rm = TRUE)) * 5200

  dat$epi$prev.sti[at] <- sum(rGC == 1 | uGC == 1 |
                                rCT == 1 | uCT == 1 | syphilis == 1 , na.rm = TRUE) / dat$epi$num[at]
  dat$epi$ir100.sti[at] <- ((dat$epi$incid.ct[at] + dat$epi$incid.gc[at] + dat$epi$incid.syph[at]) /
                              (sum(rGC == 0, na.rm = TRUE) +
                                 sum(uGC == 0, na.rm = TRUE) +
                                 sum(rCT == 0, na.rm = TRUE) +
                                 sum(uCT == 0, na.rm = TRUE) +
                                 sum(syphilis == 0, na.rm = TRUE))) * 5200

  dat$epi$ir100.sti.prep[at] <- (dat$epi$incid.gcct.prep[at] + dat$epi$incid.syph.prep[at] /
                                  (sum(rGC == 0 & prepStat == 1, na.rm = TRUE) +
                                   sum(uGC == 0 & prepStat == 1, na.rm = TRUE) +
                                   sum(rCT == 0 & prepStat == 1, na.rm = TRUE) +
                                   sum(uCT == 0 & prepStat == 1, na.rm = TRUE) +
                                       sum(syphilis == 0 & prepStat == 1, na.rm = TRUE))) * 5200

  # PAF
  syph.prev <- which(syphilis == 1 & syph.infTime < at)
  rGC.prev <- which(rGC == 1 & rGC.infTime < at)
  uGC.prev <- which(uGC == 1 & uGC.infTime < at)
  rCT.prev <- which(rCT == 1 & rCT.infTime < at)
  uCT.prev <- which(uCT == 1 & uCT.infTime < at)
  sti.prev <- unique(c(uGC.prev, uCT.prev, rGC.prev, rCT.prev, syph.prev))
  u.sti.prev <- unique(c(uGC.prev, uCT.prev))
  r.sti.prev <- unique(c(rGC.prev, rCT.prev))
  u.sti.sympt <- which(uGC.sympt == 1 | uCT.sympt == 1)
  u.sti.asympt <- setdiff(u.sti.prev, u.sti.sympt)
  r.sti.sympt <- which(rGC.sympt == 1 | rCT.sympt == 1)
  r.sti.asympt <- setdiff(r.sti.prev, r.sti.sympt)
  
  syph.sympt <- which(stage.prim.sympt == 1 | stage.seco.sympt == 1 | stage.earlat.sympt == 1 | stage.latelat.sympt == 1 | stage.latelatelat.sympt == 1 | stage.tert.sympt == 1)
  syph.asympt <- setdiff(syph.prev, syph.sympt)
  
  # PAF values
  if (at >= 2) {
      
      dat$epi$sti_paf[at] <- length(which(inf.time[sti.prev] == at)) / length(which(inf.time == at))
      dat$epi$sti_u_paf[at] <- length(which(inf.time[u.sti.prev] == at & inf.role[u.sti.prev] == 1)) / length(which(inf.time == at))
      dat$epi$sti_r_paf[at] <- length(which(inf.time[r.sti.prev] == at & inf.role[r.sti.prev] == 0)) / length(which(inf.time == at))
      dat$epi$sti_syph_paf[at] <- length(which(inf.time[syph.prev] == at)) / length(which(inf.time == at))
      dat$epi$sti_u_sympt_paf[at] <- length(which(inf.time[u.sti.sympt] == at)) / length(which(inf.time == at))
      dat$epi$sti_u_asympt_paf[at] <- length(which(inf.time[u.sti.asympt] == at)) / length(which(inf.time == at))
      dat$epi$sti_r_sympt_paf[at] <- length(which(inf.time[r.sti.sympt] == at)) / length(which(inf.time == at))
      dat$epi$sti_r_asympt_paf[at] <- length(which(inf.time[r.sti.asympt] == at)) / length(which(inf.time == at))
      dat$epi$sti_syph_sympt_paf[at] <- length(which(inf.time[syph.sympt] == at)) / length(which(inf.time == at))
      dat$epi$sti_syph_asympt_paf[at] <- length(which(inf.time[syph.asympt] == at)) / length(which(inf.time == at))
  
  # Sums at each time step (numerators from paf formulae)
      dat$epi$sti_hiv_sum[at] <- length(which(inf.time[sti.prev] == at))
      dat$epi$sti_u_hiv_sum[at] <- length(which(inf.time[u.sti.prev] == at & inf.role[u.sti.prev] == 1))
      dat$epi$sti_r_hiv_sum[at] <- length(which(inf.time[r.sti.prev] == at & inf.role[r.sti.prev] == 0))
      dat$epi$sti_syph_hiv_sum[at] <- length(which(inf.time[syph.prev] == at))
      dat$epi$hiv_sum[at] <- length(which(inf.time == at))
     
      dat$epi$sti_u_sympt_hiv_sum[at] <- length(which(inf.time[u.sti.sympt] == at))
      dat$epi$sti_u_asympt_hiv_sum[at] <- length(which(inf.time[u.sti.asympt] == at))
      dat$epi$sti_r_sympt_hiv_sum[at] <- length(which(inf.time[r.sti.sympt] == at))
      dat$epi$sti_r_asympt_hiv_sum[at] <- length(which(inf.time[r.sti.asympt] == at))
      dat$epi$sti_syph_sympt_hiv_sum[at] <- length(which(inf.time[syph.sympt] == at))
      dat$epi$sti_syph_asympt_hiv_sum[at] <- length(which(inf.time[syph.asympt] == at))
  }
  
  return(dat)
}


#' @title Prevalence Module
#'
#' @description Module function to calculate and store summary statistics for
#'              disease prevalence, demographics, and other epidemiological
#'              outcomes.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
prevalence_het <- function(dat, at) {

  status <- dat$attr$status
  male <- dat$attr$male
  age <- dat$attr$age

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  # Initialize vectors
  if (at == 1) {
    dat$epi$i.num <- rNA
    dat$epi$num <- rNA

    dat$epi$i.num.male <- rNA
    dat$epi$i.num.feml <- rNA
    dat$epi$i.prev.male <- rNA
    dat$epi$i.prev.feml <- rNA

    dat$epi$num.male <- rNA
    dat$epi$num.feml <- rNA
    dat$epi$meanAge <- rNA
    dat$epi$propMale <- rNA

    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA
  }

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$num[at] <- length(status)

  dat$epi$i.num.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.num.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE)
  dat$epi$i.prev.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE) /
    sum(male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE) /
    sum(male == 0, na.rm = TRUE)

  dat$epi$num.male[at] <- sum(male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(male == 0, na.rm = TRUE)
  dat$epi$meanAge[at] <- mean(age, na.rm = TRUE)
  dat$epi$propMale[at] <- mean(male, na.rm = TRUE)

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}
