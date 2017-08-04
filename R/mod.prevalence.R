
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
#' This function returns the \code{dat} object with an updated summary of
#' current attributes stored in \code{dat$epi}.
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
  syph.sympt <- dat$attr$syph.sympt

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
    dat$epi$deathage <- rNA

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
    dat$epi$stage.time.af.ndx <- rep(0, nsteps)
    dat$epi$stage.time.af.dx <- rep(0, nsteps)
    dat$epi$stage.time.early.chronic.ndx <- rep(0, nsteps)
    dat$epi$stage.time.early.chronic.dx.yrone <- rep(0, nsteps)
    dat$epi$stage.time.early.chronic.dx.yrstwotolate <- rep(0, nsteps)
    dat$epi$stage.time.early.chronic.art <- rep(0, nsteps)
    dat$epi$stage.time.late.chronic.ndx <- rep(0, nsteps)
    dat$epi$stage.time.late.chronic.dx <- rep(0, nsteps)
    dat$epi$stage.time.late.chronic.art <- rep(0, nsteps)
    dat$epi$stage.time.aids.ndx <- rep(0, nsteps)
    dat$epi$stage.time.aids.dx <- rep(0, nsteps)
    dat$epi$stage.time.aids.art <- rep(0, nsteps)

    dat$epi$hivtests.prep <- rep(0, nsteps)
    dat$epi$hivtests.nprep <- rep(0, nsteps)
    dat$epi$hivtests.pos <- rep(0, nsteps)

    dat$epi$rGCsympttests <- rep(0, nsteps)
    dat$epi$uGCsympttests <- rep(0, nsteps)
    dat$epi$rCTsympttests <- rep(0, nsteps)
    dat$epi$uCTsympttests <- rep(0, nsteps)
    dat$epi$syphsympttests <- rep(0, nsteps)

    dat$epi$rGCasympttests <- rep(0, nsteps)
    dat$epi$uGCasympttests <- rep(0, nsteps)
    dat$epi$GCasympttests <- rep(0, nsteps)
    dat$epi$rGCasympttests.pos <- rep(0, nsteps)
    dat$epi$uGCasympttests.pos <- rep(0, nsteps)
    dat$epi$GCasympttests.pos <- rep(0, nsteps)

    dat$epi$rCTasympttests <- rep(0, nsteps)
    dat$epi$uCTasympttests <- rep(0, nsteps)
    dat$epi$CTasympttests <- rep(0, nsteps)
    dat$epi$rCTasympttests.pos <- rep(0, nsteps)
    dat$epi$uCTasympttests.pos <- rep(0, nsteps)
    dat$epi$CTasympttests.pos <- rep(0, nsteps)

    dat$epi$syphasympttests <- rep(0, nsteps)
    dat$epi$syphasympttests.pos <- rep(0, nsteps)

    dat$epi$stiasympttests <- rep(0, nsteps)
    dat$epi$stiasympttests.pos <- rep(0, nsteps)

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

    #HIV/STI coinfection with conditional denominators
    dat$epi$prev.primsecosyph.hivneg <- rNA
    dat$epi$prev.primsecosyph.hivpos <- rNA
    dat$epi$prev.syph.hivneg <- rNA
    dat$epi$prev.syph.hivpos <- rNA
    dat$epi$prev.gc.hivneg <- rNA
    dat$epi$prev.gc.hivpos <- rNA
    dat$epi$prev.ct.hivneg <- rNA
    dat$epi$prev.ct.hivpos <- rNA
    dat$epi$prev.hiv.primsecosyphpos <- rNA
    dat$epi$prev.hiv.primsecosyphneg <- rNA
    dat$epi$prev.hiv.syphpos <- rNA
    dat$epi$prev.hiv.syphneg <- rNA
    dat$epi$prev.hiv.gcpos <- rNA
    dat$epi$prev.hiv.gcneg <- rNA
    dat$epi$prev.hiv.ctpos <- rNA
    dat$epi$prev.hiv.ctneg <- rNA
    dat$epi$prev.rgc.hivpos <- rNA
    dat$epi$prev.ugc.hivpos <- rNA
    dat$epi$prev.rct.hivpos <- rNA
    dat$epi$prev.uct.hivpos <- rNA
    dat$epi$prev.rgc.hivneg <- rNA
    dat$epi$prev.ugc.hivneg <- rNA
    dat$epi$prev.rct.hivneg <- rNA
    dat$epi$prev.uct.hivneg <- rNA

    #HIV/STI coinfection prevalence
    dat$epi$prev.rct.uct <- rNA
    dat$epi$prev.rgc.ugc <- rNA
    dat$epi$prev.rct.rgc <- rNA
    dat$epi$prev.uct.ugc <- rNA
    dat$epi$prev.gc.syph <- rNA
    dat$epi$prev.ct.syph <- rNA
    dat$epi$prev.gc.primsecosyph <- rNA
    dat$epi$prev.ct.primsecosyph <- rNA

    dat$epi$incid.rgc <- rNA
    dat$epi$incid.ugc <- rNA
    dat$epi$incid.gc <- rNA
    dat$epi$incid.rct <- rNA
    dat$epi$incid.uct <- rNA
    dat$epi$incid.ct <- rNA
    dat$epi$incid.syph <- rNA
    dat$epi$incid.sti <- rNA

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
    dat$epi$sum_GC <- rNA
    dat$epi$sum_CT <- rNA
    dat$epi$sum_syph <- rNA
    dat$epi$sum_urethral <- rNA
    dat$epi$sum_rectal <- rNA
    #2x2 for PAF
    #               HIV+
    #             STI+  STI-
    #HIV-   STI +  1    2
    #       STI -  3    4
    dat$epi$cell1_gc <- rNA
    dat$epi$cell2_gc <- rNA
    dat$epi$cell3_gc <- rNA
    dat$epi$cell4_gc <- rNA
    dat$epi$cell1_ct <- rNA
    dat$epi$cell2_ct <- rNA
    dat$epi$cell3_ct <- rNA
    dat$epi$cell4_ct <- rNA
    dat$epi$cell1_syph <- rNA
    dat$epi$cell2_syph <- rNA
    dat$epi$cell3_syph <- rNA
    dat$epi$cell4_syph <- rNA
    dat$epi$cell1_sti <- rNA
    dat$epi$cell2_sti <- rNA
    dat$epi$cell3_sti <- rNA
    dat$epi$cell4_sti <- rNA
    # dat$epi$sti_paf <- rNA
    # dat$epi$sti_u_paf <- rNA
    # dat$epi$sti_r_paf <- rNA
    # dat$epi$sti_syph_paf <- rNA
    # dat$epi$sti_u_sympt_paf <- rNA
    # dat$epi$sti_u_asympt_paf <- rNA
    # dat$epi$sti_r_sympt_paf <- rNA
    # dat$epi$sti_r_asympt_paf <- rNA
    # dat$epi$sti_syph_sympt_paf <- rNA
    # dat$epi$sti_syph_asympt_paf <- rNA
    # dat$epi$sti_hiv_sum <- rNA
    # dat$epi$sti_u_hiv_sum <- rNA
    # dat$epi$sti_r_hiv_sum <- rNA
    # dat$epi$sti_syph_hiv_sum <- rNA
    # dat$epi$hiv_sum <- rNA
    # dat$epi$sti_u_sympt_hiv_sum <- rNA
    # dat$epi$sti_u_asympt_hiv_sum <- rNA
    # dat$epi$sti_r_sympt_hiv_sum <- rNA
    # dat$epi$sti_r_asympt_hiv_sum <- rNA
    # dat$epi$sti_syph_sympt_hiv_sum <- rNA
    # dat$epi$sti_syph_asympt_hiv_sum <- rNA

    dat$epi$recov.rgc <- rNA
    dat$epi$recov.ugc <- rNA
    dat$epi$recov.rct <- rNA
    dat$epi$recov.uct <- rNA
    dat$epi$recov.earlysyph <- rNA
    dat$epi$recov.syphilis <- rNA

    dat$epi$trans.main <- rNA
    dat$epi$trans.pers <- rNA
    dat$epi$trans.inst <- rNA

    dat$epi$txGC <- rNA
    dat$epi$txCT <- rNA
    dat$epi$txsyph <- rNA
    dat$epi$txearlysyph <- rNA
    dat$epi$txlatesyph <- rNA
    dat$epi$txasympt <- rNA

    # dat$epi$stiactiveind <- rNA
    # dat$epi$recentpartners <- rNA
    # dat$epi$recentSTI <- rNA
    # dat$epi$newpartner <- rNA
    # dat$epi$concurrpart <- rNA
    # dat$epi$partnersti <- rNA
    # dat$epi$uai.nmain <- rNA
    # dat$epi$uai.any <- rNA

    dat$epi$eptCov <- rNA
    dat$epi$eptpartelig <- rNA
    dat$epi$eptpartprovided <- rNA
    dat$epi$eptpartuptake <- rNA
    dat$epi$eptTx <- rNA
    dat$epi$propindexeptElig <- rNA
    dat$epi$eptprop_provided <- rNA
    dat$epi$eptprop_tx <- rNA
    dat$epi$eptuninfectedprovided <- rNA
    dat$epi$eptuninfecteduptake <- rNA
    dat$epi$eptgcinfectsti <- rNA
    dat$epi$eptctinfectsti <- rNA
    dat$epi$eptgcinfecthiv <- rNA
    dat$epi$eptctinfecthiv <- rNA

    # Proportion with X partners
    dat$epi$zeropart <- rNA
    dat$epi$onepart <- rNA
    dat$epi$twopart <- rNA
    dat$epi$threepart <- rNA
    dat$epi$fourpart <- rNA
    dat$epi$fivepart <- rNA
    dat$epi$gtfivepart <- rNA
  }

  dat$epi$num[at] <- sum(race %in% c("B","W"), na.rm = TRUE)
  dat$epi$num.B[at] <- sum(race == "B", na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == "W", na.rm = TRUE)
  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$i.num.B[at] <- sum(status == 1 & race == "B", na.rm = TRUE)
  dat$epi$i.num.W[at] <- sum(status == 1 & race == "W", na.rm = TRUE)
  dat$epi$i.prev[at] <- ifelse(dat$epi$num[at] > 0, dat$epi$i.num[at] / dat$epi$num[at], 0)
  dat$epi$i.prev.B[at] <- ifelse(dat$epi$num[at] > 0, dat$epi$i.num.B[at] / dat$epi$num.B[at], 0)
  dat$epi$i.prev.W[at] <- ifelse(dat$epi$num[at] > 0, dat$epi$i.num.W[at] / dat$epi$num.W[at], 0)
  dat$epi$ir100[at] <- ifelse(sum(status == 0, na.rm = TRUE) > 0, (dat$epi$incid[at] /
                          sum(status == 0, na.rm = TRUE)) * 5200, 0)

  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  dat$epi$i.num.prep0[at] <- sum((is.na(prepStat) | prepStat == 0) & status == 1, na.rm = TRUE)
  dat$epi$i.num.prep1[at] <- sum(prepStat == 1 & status == 1, na.rm = TRUE)
  dat$epi$i.prev.prep0[at] <- dat$epi$i.num.prep0[at] / sum((is.na(prepStat) | prepStat == 0), na.rm = TRUE)

  if (at == 1) {
    dat$epi$i.prev.prep1[1] <- 0
  } else {
    dat$epi$i.prev.prep1[at] <- dat$epi$i.num.prep1[at] / sum(prepStat == 1, na.rm = TRUE)
  }

  dat$epi$time.on.prep[at] <- length(which(dat$attr$prepStat == 1))
  dat$epi$time.off.prep[at] <- length(which(dat$attr$prepStat == 0))
  dat$epi$stisympttests[at] <- sum(dat$epi$syphsympttests[at], dat$epi$CTsympttests[at],
                                   dat$epi$GCsympttests[at], na.rm = TRUE)
  dat$epi$stiasympttests[at] <- sum(dat$epi$syphasympttests[at], dat$epi$CTasympttests[at],
                                    dat$epi$GCasympttests[at], na.rm = TRUE)
  dat$epi$stiasympttests.prep[at] <- sum(dat$epi$syphasympttests.prep[at], dat$epi$CTasympttests.prep[at],
                                         dat$epi$GCasympttests.prep[at], na.rm = TRUE)
  dat$epi$stiasympttests.pos[at] <- sum(dat$epi$syphasympttests.pos[at], dat$epi$CTasympttests.pos[at],
                                        dat$epi$GCasympttests.pos[at], na.rm = TRUE)
  dat$epi$stiasympttests.pos.prep[at] <- sum(dat$epi$syphasympttests.pos.prep[at], dat$epi$CTasympttests.pos.prep[at],
                                           dat$epi$GCasympttests.pos.prep[at], na.rm = TRUE)

  # STI Prevalence

  dat$epi$prev.rgc[at] <- ifelse(dat$epi$num[at] > 0, sum(rGC == 1, na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.ugc[at] <- ifelse(dat$epi$num[at] > 0, sum(uGC == 1, na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.gc[at] <- ifelse(dat$epi$num[at] > 0, sum((rGC == 1 | uGC == 1), na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.gc.sympt[at] <- ifelse(dat$epi$num[at] > 0, sum((rGC.sympt == 1 | uGC.sympt == 1)) / dat$epi$num[at], 0)
  dat$epi$prev.gc.dual[at] <- ifelse(dat$epi$num[at] > 0, sum((rGC == 1 & uGC == 1), na.rm = TRUE) / dat$epi$num[at], 0)

  dat$epi$prev.rct[at] <- ifelse(dat$epi$num[at] > 0, sum(rCT == 1, na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.uct[at] <- ifelse(dat$epi$num[at] > 0, sum(uCT == 1, na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.ct[at] <- ifelse(dat$epi$num[at] > 0, sum((rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.ct.sympt[at] <- ifelse(dat$epi$num[at] > 0, sum((rCT.sympt == 1 | uCT.sympt == 1)) / dat$epi$num[at], 0)
  dat$epi$prev.ct.dual[at] <- ifelse(dat$epi$num[at] > 0, sum((rCT == 1 & uCT == 1), na.rm = TRUE) / dat$epi$num[at], 0)

  dat$epi$prev.rgcct[at] <- ifelse(dat$epi$num[at] > 0, sum(rGC == 1 | rCT == 1, na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.ugcct[at] <- ifelse(dat$epi$num[at] > 0, sum(uGC == 1 | uCT == 1, na.rm = TRUE) / dat$epi$num[at], 0)

  dat$epi$prev.stage.incub[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                         length(which(stage.syph == 1)) / length(which(syphilis == 1)), 0)
  dat$epi$prev.stage.prim[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                        length(which(stage.syph == 2)) / length(which(syphilis == 1)), 0)
  dat$epi$prev.stage.incubprim[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                             length(which(stage.syph == 1 | stage.syph == 2)) / length(which(syphilis == 1)), 0)
  dat$epi$prev.stage.seco[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                        length(which(stage.syph == 3)) / length(which(syphilis == 1)), 0)
  dat$epi$prev.stage.earlat[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                          length(which(stage.syph == 4)) / length(which(syphilis == 1)), 0)
  dat$epi$prev.stage.latelat[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                           length(which(stage.syph == 5)) / length(which(syphilis == 1)), 0)
  dat$epi$prev.stage.latelatelat[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                               length(which(stage.syph == 6)) / length(which(syphilis == 1)), 0)
  dat$epi$prev.stage.alllatelat[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                              length(which(stage.syph == 5 | stage.syph == 6)) / length(which(syphilis == 1)), 0)
  dat$epi$prev.stage.tert[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                        length(which(stage.syph == 7)) / length(which(syphilis == 1)), 0)
  dat$epi$prev.earlysyph[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                       length(which(stage.syph %in% c(1, 2, 3, 4))) / length(which(syphilis == 1)), 0)
  dat$epi$prev.latesyph[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                      length(which(stage.syph %in% c(5, 6, 7))) / length(which(syphilis == 1)), 0)
  dat$epi$prev.syph[at] <- ifelse(dat$epi$num[at] > 0, length(which(syphilis == 1)) / dat$epi$num[at], 0)
  dat$epi$prev.primsecosyph[at] <- ifelse(dat$epi$num[at] > 0, length(which(stage.syph %in% c(1, 2, 3))) / dat$epi$num[at], 0)

  # uGC.prev <- which(uGC == 1 & uGC.infTime < at)
  # uCT.prev <- which(uCT == 1 & uCT.infTime < at)
  # rGC.prev <- which(rGC == 1 & rGC.infTime < at)
  # rCT.prev <- which(rCT == 1 & rCT.infTime < at)

  # Prevalence of HIV/STI overlap (conditional denominators)
  dat$epi$prev.primsecosyph.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
                                                 length(intersect(which(status == 0), which(stage.syph %in% c(1, 2, 3)))) / dat$epi$s.num[at], 0)
  dat$epi$prev.primsecosyph.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
                                            length(intersect(which(status == 1), which(stage.syph %in% c(1, 2, 3)))) / dat$epi$i.num[at], 0)
  dat$epi$prev.syph.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
                                         length(intersect(which(status == 0), which(syphilis == 1))) / dat$epi$s.num[at], 0)
  dat$epi$prev.syph.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
                                         length(intersect(which(status == 1), which(syphilis == 1))) / dat$epi$i.num[at], 0)

  dat$epi$prev.gc.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
                                       length(intersect(which(status == 0), which((rGC == 1 | uGC == 1)))) / dat$epi$s.num[at], 0)
  dat$epi$prev.gc.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
                                       length(intersect(which(status == 1), which((rGC == 1 | uGC == 1)))) / dat$epi$i.num[at], 0)

  dat$epi$prev.ct.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
                                       length(intersect(which(status == 0),  which((rCT == 1 | uCT == 1)))) / dat$epi$s.num[at], 0)
  dat$epi$prev.ct.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
                                       length(intersect(which(status == 1), which((rCT == 1 | uCT == 1)))) / dat$epi$i.num[at], 0)

  dat$epi$prev.hiv.primsecosyphpos[at] <- ifelse(length(which(stage.syph %in% c(1, 2, 3))) > 0,
                                            length(intersect(which(status == 1), which(stage.syph %in% c(1, 2, 3)))) / length(which(stage.syph %in% c(1, 2, 3))), 0)
  dat$epi$prev.hiv.primsecosyphneg[at] <- ifelse(length(which(stage.syph %in% c(1, 2, 3))) > 0,
                                            length(intersect(which(status == 1), which(stage.syph %in% c(1, 2, 3)))) / length(which(stage.syph %in% c(1, 2, 3))), 0)

  dat$epi$prev.hiv.syphpos[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                         length(intersect(which(status == 1), which(syphilis == 1))) / length(which(syphilis == 1)), 0)
  dat$epi$prev.hiv.syphneg[at] <- ifelse(length(which(syphilis == 0)) > 0,
                                         length(intersect(which(status == 1), which(syphilis == 1))) / length(which(syphilis == 0)), 0)

  dat$epi$prev.hiv.gcpos[at] <- ifelse(sum((rGC == 1 | uGC == 1)) > 0,
                                       length(intersect(which(status == 1), which((rGC == 1 | uGC == 1)))) / sum((rGC == 1 | uGC == 1), na.rm = TRUE), 0)
  dat$epi$prev.hiv.gcneg[at] <- ifelse(sum((rGC == 0 & uGC == 0)) > 0,
                                       length(intersect(which(status == 1), which((rGC == 0 & uGC == 0)))) / sum((rGC == 0 & uGC == 0), na.rm = TRUE), 0)

  dat$epi$prev.hiv.ctpos[at] <- ifelse(sum((rCT == 1 | uCT == 1)) > 0,
                                       length(intersect(which(status == 1), which((rCT == 1 | uCT == 1)))) / sum((rCT == 1 | uCT == 1), na.rm = TRUE), 0)
  dat$epi$prev.hiv.ctneg[at] <- ifelse(sum((rCT == 0 & uCT == 0)) > 0,
                                       length(intersect(which(status == 1), which((rCT == 0 & uCT == 0)))) / sum((rCT == 0 & uCT == 0), na.rm = TRUE), 0)

  dat$epi$prev.rgc.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
                                        length(intersect(which(status == 1), which(rGC == 1))) / dat$epi$i.num[at], 0)
  dat$epi$prev.ugc.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
                                        length(intersect(which(status == 1),  which(uGC == 1))) / dat$epi$i.num[at], 0)
  dat$epi$prev.rct.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
                                        length(intersect(which(status == 1), which(rCT == 1))) / dat$epi$i.num[at], 0)
  dat$epi$prev.uct.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
                                        length(intersect(which(status == 1), which(uCT == 1))) / dat$epi$i.num[at], 0)

  dat$epi$prev.rgc.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
                                        length(intersect(which(status == 1), which(rGC == 1))) / dat$epi$s.num[at], 0)
  dat$epi$prev.ugc.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
                                        length(intersect(which(status == 1), which(uGC == 1))) / dat$epi$s.num[at], 0)
  dat$epi$prev.rct.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
                                        length(intersect(which(status == 1), which(rCT == 1))) / dat$epi$s.num[at], 0)
  dat$epi$prev.uct.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
                                        length(intersect(which(status == 1), which(uCT == 1))) / dat$epi$s.num[at], 0)

  # Co-infection prevalence
  dat$epi$prev.rct.uct[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rCT == 1), which(uCT == 1)))) / dat$epi$num[at], 0)
  dat$epi$prev.rgc.ugc[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rGC == 1), which(uGC == 1)))) / dat$epi$num[at], 0)
  dat$epi$prev.rct.rgc[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rGC == 1), which(rCT == 1))))  / dat$epi$num[at], 0)
  dat$epi$prev.uct.ugc[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(uCT == 1), which(uGC == 1)))) / dat$epi$num[at], 0)
  dat$epi$prev.gc.syph[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rGC == 1 | uGC == 1), which(syphilis == 1)))) / dat$epi$num[at], 0)
  dat$epi$prev.ct.syph[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rCT == 1 | uCT == 1), which(syphilis == 1)))) / dat$epi$num[at], 0)
  dat$epi$prev.gc.primsecosyph[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rGC == 1 | uGC == 1), which(stage.syph %in% c(1, 2))))) / dat$epi$num[at], 0)
  dat$epi$prev.ct.primsecosyph[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rCT == 1 | uCT == 1), which(stage.syph %in% c(1, 2))))) / dat$epi$num[at], 0)

  # Site-specific STI incidence rates
  dat$epi$ir100.rgc[at] <- ifelse(sum(rGC == 0, na.rm = TRUE) > 0, (dat$epi$incid.rgc[at] / sum(rGC == 0, na.rm = TRUE)) * 5200, 0)
  dat$epi$ir100.ugc[at] <- ifelse(sum(uGC == 0, na.rm = TRUE) > 0, (dat$epi$incid.ugc[at] / sum(uGC == 0, na.rm = TRUE)) * 5200, 0)
  dat$epi$ir100.gc[at] <- ifelse((sum(rGC == 0, na.rm = TRUE) + sum(uGC == 0, na.rm = TRUE)) > 0, (dat$epi$incid.gc[at] / (sum(rGC == 0, na.rm = TRUE) + sum(uGC == 0, na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.rct[at] <- ifelse(sum(rCT == 0, na.rm = TRUE) > 0, (dat$epi$incid.rct[at] / sum(rCT == 0, na.rm = TRUE)) * 5200, 0)
  dat$epi$ir100.uct[at] <- ifelse(sum(uCT == 0, na.rm = TRUE) > 0, (dat$epi$incid.uct[at] / sum(uCT == 0, na.rm = TRUE)) * 5200, 0)
  dat$epi$ir100.ct[at] <- ifelse((sum(rCT == 0, na.rm = TRUE) + sum(uCT == 0, na.rm = TRUE)) > 0, (dat$epi$incid.ct[at] / (sum(rCT == 0, na.rm = TRUE) + sum(uCT == 0, na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.syph[at] <- ifelse(sum(syphilis == 0, na.rm = TRUE) > 0, (dat$epi$incid.syph[at] / sum(syphilis == 0 , na.rm = TRUE)) * 5200, 0)

  dat$epi$prev.sti[at] <- ifelse(sum(rGC == 1 | uGC == 1 | rCT == 1 | uCT == 1 | syphilis == 1 , na.rm = TRUE) > 0,
                                 sum(rGC == 1 | uGC == 1 | rCT == 1 | uCT == 1 | syphilis == 1 , na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$ir100.sti[at] <- ifelse((sum(rGC == 0, na.rm = TRUE) + sum(uGC == 0, na.rm = TRUE) +
                                     sum(rCT == 0, na.rm = TRUE) + sum(uCT == 0, na.rm = TRUE) +
                                     sum(syphilis == 0, na.rm = TRUE)) > 0,
                                ((dat$epi$incid.ct[at] + dat$epi$incid.gc[at] + dat$epi$incid.syph[at]) /
                                (sum(rGC == 0, na.rm = TRUE) + sum(uGC == 0, na.rm = TRUE) +
                                 sum(rCT == 0, na.rm = TRUE) + sum(uCT == 0, na.rm = TRUE) +
                                 sum(syphilis == 0, na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.sti.prep[at] <- ifelse((sum(rGC == 0 & prepStat == 1, na.rm = TRUE) + sum(uGC == 0 & prepStat == 1, na.rm = TRUE) +
                                          sum(rCT == 0 & prepStat == 1, na.rm = TRUE) + sum(uCT == 0 & prepStat == 1, na.rm = TRUE) +
                                          sum(syphilis == 0 & prepStat == 1, na.rm = TRUE)) > 0,
                                  (dat$epi$incid.gcct.prep[at] + dat$epi$incid.syph.prep[at] /
                                  (sum(rGC == 0 & prepStat == 1, na.rm = TRUE) + sum(uGC == 0 & prepStat == 1, na.rm = TRUE) +
                                   sum(rCT == 0 & prepStat == 1, na.rm = TRUE) + sum(uCT == 0 & prepStat == 1, na.rm = TRUE) +
                                   sum(syphilis == 0 & prepStat == 1, na.rm = TRUE))) * 5200, 0)

  # PAF
  # syph.prev <- which(syphilis == 1 & syph.infTime < at)
  # rGC.prev <- which(rGC == 1 & rGC.infTime < at)
  # uGC.prev <- which(uGC == 1 & uGC.infTime < at)
  # rCT.prev <- which(rCT == 1 & rCT.infTime < at)
  # uCT.prev <- which(uCT == 1 & uCT.infTime < at)
  # sti.prev <- unique(c(uGC.prev, uCT.prev, rGC.prev, rCT.prev, syph.prev))
  # u.sti.prev <- unique(c(uGC.prev, uCT.prev))
  # r.sti.prev <- unique(c(rGC.prev, rCT.prev))
  # u.sti.sympt <- which(uGC.sympt == 1 | uCT.sympt == 1)
  # u.sti.asympt <- setdiff(u.sti.prev, u.sti.sympt)
  # r.sti.sympt <- which(rGC.sympt == 1 | rCT.sympt == 1)
  # r.sti.asympt <- setdiff(r.sti.prev, r.sti.sympt)
  #
  # syph.sympt <- which(syph.sympt == 1)
  # syph.asympt <- setdiff(syph.prev, syph.sympt)
  #
  # # PAF values
  # if (at >= 2) {
  #
  #     dat$epi$sti_paf[at] <- ifelse(length(which(inf.time == at)) > 0, length(which(inf.time[sti.prev] == at)) / length(which(inf.time == at)), 0)
  #     dat$epi$sti_u_paf[at] <- ifelse(length(which(inf.time == at)) > 0, length(which(inf.time[u.sti.prev] == at & inf.role[u.sti.prev] == 1)) / length(which(inf.time == at)), 0)
  #     dat$epi$sti_r_paf[at] <- ifelse(length(which(inf.time == at)) > 0, length(which(inf.time[r.sti.prev] == at & inf.role[r.sti.prev] == 0)) / length(which(inf.time == at)), 0)
  #     dat$epi$sti_syph_paf[at] <- ifelse(length(which(inf.time == at)) > 0, length(which(inf.time[syph.prev] == at)) / length(which(inf.time == at)), 0)
  #     dat$epi$sti_u_sympt_paf[at] <- ifelse(length(which(inf.time == at)) > 0, length(which(inf.time[u.sti.sympt] == at & inf.role[u.sti.sympt] == 1)) / length(which(inf.time == at)), 0)
  #     dat$epi$sti_u_asympt_paf[at] <- ifelse(length(which(inf.time == at)) > 0, length(which(inf.time[u.sti.asympt] == at & inf.role[u.sti.asympt] == 1)) / length(which(inf.time == at)), 0)
  #     dat$epi$sti_r_sympt_paf[at] <- ifelse(length(which(inf.time == at)) > 0, length(which(inf.time[r.sti.sympt] == at & inf.role[r.sti.sympt] == 0)) / length(which(inf.time == at)), 0)
  #     dat$epi$sti_r_asympt_paf[at] <- ifelse(length(which(inf.time == at)) > 0, length(which(inf.time[r.sti.asympt] == at & inf.role[r.sti.asympt] == 0)) / length(which(inf.time == at)), 0)
  #     dat$epi$sti_syph_sympt_paf[at] <- ifelse(length(which(inf.time == at)) > 0, length(which(inf.time[syph.sympt] == at)) / length(which(inf.time == at)), 0)
  #     dat$epi$sti_syph_asympt_paf[at] <- ifelse(length(which(inf.time == at)) > 0, length(which(inf.time[syph.asympt] == at)) / length(which(inf.time == at)), 0)
  #
  # # Sums at each time step (numerators from paf formulae)
  #     dat$epi$sti_hiv_sum[at] <- length(which(inf.time[sti.prev] == at))
  #     dat$epi$sti_u_hiv_sum[at] <- length(which(inf.time[u.sti.prev] == at & inf.role[u.sti.prev] == 1))
  #     dat$epi$sti_r_hiv_sum[at] <- length(which(inf.time[r.sti.prev] == at & inf.role[r.sti.prev] == 0))
  #     dat$epi$sti_syph_hiv_sum[at] <- length(which(inf.time[syph.prev] == at))
  #     dat$epi$hiv_sum[at] <- length(which(inf.time == at))
  #
  #     dat$epi$sti_u_sympt_hiv_sum[at] <- length(which(inf.time[u.sti.sympt] == at))
  #     dat$epi$sti_u_asympt_hiv_sum[at] <- length(which(inf.time[u.sti.asympt] == at))
  #     dat$epi$sti_r_sympt_hiv_sum[at] <- length(which(inf.time[r.sti.sympt] == at))
  #     dat$epi$sti_r_asympt_hiv_sum[at] <- length(which(inf.time[r.sti.asympt] == at))
  #     dat$epi$sti_syph_sympt_hiv_sum[at] <- length(which(inf.time[syph.sympt] == at))
  #     dat$epi$sti_syph_asympt_hiv_sum[at] <- length(which(inf.time[syph.asympt] == at))
  #}

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
  dat$epi$i.prev.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE) / sum(male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE) / sum(male == 0, na.rm = TRUE)

  dat$epi$num.male[at] <- sum(male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(male == 0, na.rm = TRUE)
  dat$epi$meanAge[at] <- mean(age, na.rm = TRUE)
  dat$epi$propMale[at] <- mean(male, na.rm = TRUE)

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) > (param$vl.acute.topeak + param$vl.acute.toset))
}
