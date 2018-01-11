
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
  race <- dat$attr$race
  status <- dat$attr$status
  prepStat <- dat$attr$prepStat
  prepElig <- dat$attr$prepElig
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT
  syphilis <- dat$attr$syphilis
  # rGC.sympt <- dat$attr$rGC.sympt
  # uGC.sympt <- dat$attr$uGC.sympt
  # rCT.sympt <- dat$attr$rCT.sympt
  # uCT.sympt <- dat$attr$uCT.sympt
  stage.syph <- dat$attr$stage.syph
  diag.status.syph <- dat$attr$diag.status.syph
  last.diag.time.syph <- dat$attr$last.diag.time.syph
  tslt.rgc <- dat$attr$time.since.last.test.rgc
  tslt.ugc <- dat$attr$time.since.last.test.ugc
  tslt.rct <- dat$attr$time.since.last.test.rct
  tslt.uct <- dat$attr$time.since.last.test.uct
  tslt.syph <- dat$attr$time.since.last.test.syph
  diag.status <- dat$attr$diag.status
  tt.traj.gc.hivpos <- dat$attr$tt.traj.gc.hivpos
  tt.traj.gc.hivneg <- dat$attr$tt.traj.gc.hivneg
  tt.traj.ct.hivneg <- dat$attr$tt.traj.ct.hivneg
  tt.traj.ct.hivpos <- dat$attr$tt.traj.ct.hivpos
  tt.traj.syph.hivpos <- dat$attr$tt.traj.syph.hivpos
  tt.traj.syph.hivneg <- dat$attr$tt.traj.syph.hivneg

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  if (at == 1) {

    # Population sizes and HIV incidence/prevalence
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

    # PrEP
    dat$epi$prepCurr <- rNA
    dat$epi$prepCov <- rNA
    dat$epi$prepElig <- rNA
    dat$epi$prepStart <- rNA
    dat$epi$i.num.prep0 <- rNA
    dat$epi$i.num.prep1 <- rNA

    #Time in health-related states
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

    # Number of HIV tests
    dat$epi$hivtests.prep <- rep(0, nsteps)
    dat$epi$hivtests.nprep <- rep(0, nsteps)
    dat$epi$hivtests.pos <- rep(0, nsteps)

    # Number of STI tests
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
    # dat$epi$rGCasympttests.hivpos <- rep(0, nsteps)
    # dat$epi$uGCasympttests.hivpos <- rep(0, nsteps)
    # dat$epi$GCasympttests.hivpos <- rep(0, nsteps)
    # dat$epi$rGCasympttests.pos.hivpos <- rep(0, nsteps)
    # dat$epi$uGCasympttests.pos.hivpos <- rep(0, nsteps)
    # dat$epi$GCasympttests.pos.hivpos <- rep(0, nsteps)
    # dat$epi$rGCasympttests.hivneg <- rep(0, nsteps)
    # dat$epi$uGCasympttests.hivneg <- rep(0, nsteps)
    # dat$epi$GCasympttests.hivneg <- rep(0, nsteps)
    # dat$epi$rGCasympttests.pos.hivneg <- rep(0, nsteps)
    # dat$epi$uGCasympttests.pos.hivneg <- rep(0, nsteps)
    # dat$epi$GCasympttests.pos.hivneg <- rep(0, nsteps)
    # dat$epi$rGCasympttests.prep <- rep(0, nsteps)
    # dat$epi$uGCasympttests.prep <- rep(0, nsteps)
    # dat$epi$GCasympttests.prep <- rep(0, nsteps)
    # dat$epi$rGCasympttests.pos.prep <- rep(0, nsteps)
    # dat$epi$uGCasympttests.pos.prep <- rep(0, nsteps)
    # dat$epi$GCasympttests.pos.prep <- rep(0, nsteps)

    dat$epi$rCTasympttests <- rep(0, nsteps)
    dat$epi$uCTasympttests <- rep(0, nsteps)
    dat$epi$CTasympttests <- rep(0, nsteps)
    dat$epi$rCTasympttests.pos <- rep(0, nsteps)
    dat$epi$uCTasympttests.pos <- rep(0, nsteps)
    dat$epi$CTasympttests.pos <- rep(0, nsteps)
    # dat$epi$rCTasympttests.hivpos <- rep(0, nsteps)
    # dat$epi$uCTasympttests.hivpos <- rep(0, nsteps)
    # dat$epi$CTasympttests.hivpos <- rep(0, nsteps)
    # dat$epi$rCTasympttests.pos.hivpos <- rep(0, nsteps)
    # dat$epi$uCTasympttests.pos.hivpos <- rep(0, nsteps)
    # dat$epi$CTasympttests.pos.hivpos <- rep(0, nsteps)
    # dat$epi$rCTasympttests.hivneg <- rep(0, nsteps)
    # dat$epi$uCTasympttests.hivneg <- rep(0, nsteps)
    # dat$epi$CTasympttests.hivneg <- rep(0, nsteps)
    # dat$epi$rCTasympttests.pos.hivneg <- rep(0, nsteps)
    # dat$epi$uCTasympttests.pos.hivneg <- rep(0, nsteps)
    # dat$epi$CTasympttests.pos.hivneg <- rep(0, nsteps)
    # dat$epi$rCTasympttests.prep <- rep(0, nsteps)
    # dat$epi$uCTasympttests.prep <- rep(0, nsteps)
    # dat$epi$CTasympttests.prep <- rep(0, nsteps)
    # dat$epi$rCTasympttests.pos.prep <- rep(0, nsteps)
    # dat$epi$uCTasympttests.pos.prep <- rep(0, nsteps)
    # dat$epi$CTasympttests.pos.prep <- rep(0, nsteps)

    dat$epi$syphasympttests <- rep(0, nsteps)
    dat$epi$syphasympttests.pos <- rep(0, nsteps)
    dat$epi$syphearlyasympttests.pos <- rep(0, nsteps)
    dat$epi$syphlateasympttests.pos <- rep(0, nsteps)
    # dat$epi$syphasympttests.hivpos <- rep(0, nsteps)
    # dat$epi$syphasympttests.pos.hivpos <- rep(0, nsteps)
    # dat$epi$syphearlyasympttests.pos.hivpos <- rep(0, nsteps)
    # dat$epi$syphlateasympttests.pos.hivpos <- rep(0, nsteps)
    # dat$epi$syphasympttests.hivneg <- rep(0, nsteps)
    # dat$epi$syphasympttests.pos.hivneg <- rep(0, nsteps)
    # dat$epi$syphearlyasympttests.pos.hivneg <- rep(0, nsteps)
    # dat$epi$syphlateasympttests.pos.hivneg <- rep(0, nsteps)
    # dat$epi$syphasympttests.prep <- rep(0, nsteps)
    # dat$epi$syphasympttests.pos.prep <- rep(0, nsteps)
    # dat$epi$syphearlyasympttests.pos.prep <- rep(0, nsteps)
    # dat$epi$syphlateasympttests.pos.prep <- rep(0, nsteps)

    dat$epi$stiasympttests <- rep(0, nsteps)
    dat$epi$stiasympttests.pos <- rep(0, nsteps)
    # dat$epi$stiasympttests.hivneg <- rep(0, nsteps)
    # dat$epi$stiasympttests.pos.hivneg <- rep(0, nsteps)
    # dat$epi$stiasympttests.hivpos <- rep(0, nsteps)
    # dat$epi$stiasympttests.pos.hivpos <- rep(0, nsteps)
    # dat$epi$stisympttests <- rep(0, nsteps)
    # dat$epi$stiasympttests.prep <- rep(0, nsteps)
    # dat$epi$stiasympttests.pos.prep <- rep(0, nsteps)

    #STI Testing due to HIV/Symptomatic STI
    dat$epi$rGC_hivdxtime <- rep(0, nsteps)
    dat$epi$uGC_hivdxtime <- rep(0, nsteps)
    dat$epi$rCT_hivdxtime <- rep(0, nsteps)
    dat$epi$uCT_hivdxtime <- rep(0, nsteps)
    dat$epi$syph_hivdxtime <- rep(0, nsteps)

    dat$epi$rGC_pos_hivdxtime <- rep(0, nsteps)
    dat$epi$uGC_pos_hivdxtime <- rep(0, nsteps)
    dat$epi$rCT_pos_hivdxtime <- rep(0, nsteps)
    dat$epi$uCT_pos_hivdxtime <- rep(0, nsteps)
    dat$epi$syph_pos_hivdxtime <- rep(0, nsteps)
    dat$epi$syph_earlypos_hivdxtime <- rep(0, nsteps)
    dat$epi$syph_latepos_hivdxtime <- rep(0, nsteps)

    dat$epi$rGC_symptstidxtime <- rep(0, nsteps)
    dat$epi$uGC_symptstidxtime <- rep(0, nsteps)
    dat$epi$rCT_symptstidxtime <- rep(0, nsteps)
    dat$epi$uCT_symptstidxtime <- rep(0, nsteps)
    dat$epi$syph_symptstidxtime <- rep(0, nsteps)

    dat$epi$rGC_pos_symptstidxtime <- rep(0, nsteps)
    dat$epi$uGC_pos_symptstidxtime <- rep(0, nsteps)
    dat$epi$rCT_pos_symptstidxtime <- rep(0, nsteps)
    dat$epi$uCT_pos_symptstidxtime <- rep(0, nsteps)
    dat$epi$syph_pos_symptstidxtime <- rep(0, nsteps)
    dat$epi$syph_earlypos_symptstidxtime <- rep(0, nsteps)
    dat$epi$syph_latepos_symptstidxtime <- rep(0, nsteps)

    # STI prevalence and coinfection prevalence
    dat$epi$prev.rgc <- rNA
    dat$epi$prev.ugc <- rNA
    dat$epi$prev.gc <- rNA
    # dat$epi$prev.gc.sympt <- rNA
    # dat$epi$prev.gc.dual <- rNA

    dat$epi$prev.rct <- rNA
    dat$epi$prev.uct <- rNA
    dat$epi$prev.ct <- rNA
    dat$epi$prev.gcct <- rNA
    # dat$epi$prev.ct.sympt <- rNA
    # dat$epi$prev.ct.dual <- rNA

    dat$epi$prev.rgcct <- rNA
    dat$epi$prev.ugcct <- rNA

    dat$epi$prev.syph <- rNA
    dat$epi$prev.stage.prim <- rNA
    dat$epi$prev.stage.seco <- rNA
    dat$epi$prev.stage.earlat <- rNA
    dat$epi$prev.stage.latelat <- rNA
    dat$epi$prev.stage.tert <- rNA
    dat$epi$prev.earlysyph <- rNA
    dat$epi$prev.latesyph <- rNA
    dat$epi$prev.primsecosyph <- rNA
    dat$epi$num.newearlydiagsyph <- rNA
    dat$epi$num.newlatediagsyph <- rNA
    dat$epi$early.late.syphratio <- rNA
    dat$epi$early.late.diagsyphratio <- rNA
    dat$epi$prev.dxhiv.dxipssyph <- rNA
    dat$epi$prev.dxhiv.atdxipssyph <- rNA

    # STI only by HIV serostatus
    dat$epi$prev.rgc.hivneg.only <- rNA
    dat$epi$prev.ugc.hivneg.only <- rNA
    # dat$epi$prev.gc.hivneg.only <- rNA
    dat$epi$prev.rct.hivneg.only <- rNA
    dat$epi$prev.uct.hivneg.only <- rNA
    # dat$epi$prev.ct.hivneg.only <- rNA
    dat$epi$prev.primsecosyph.hivneg.only <- rNA
    dat$epi$prev.rgc.hivpos.only <- rNA
    dat$epi$prev.ugc.hivpos.only <- rNA
    # dat$epi$prev.gc.hivpos.only <- rNA
    dat$epi$prev.rct.hivpos.only <- rNA
    dat$epi$prev.uct.hivpos.only <- rNA
    # dat$epi$prev.ct.hivpos.only <- rNA
    dat$epi$prev.primsecosyph.hivpos.only <- rNA

    # Multi STI
    dat$epi$prev.hivposmultsti <- rNA
    dat$epi$prev.hivnegmultsti <- rNA

    #HIV/STI coinfection with conditional HIV serostatus denominators
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
    # dat$epi$prev.rgc.hivpos <- rNA
    # dat$epi$prev.ugc.hivpos <- rNA
    # dat$epi$prev.rct.hivpos <- rNA
    # dat$epi$prev.uct.hivpos <- rNA
    # dat$epi$prev.rgc.hivneg <- rNA
    # dat$epi$prev.ugc.hivneg <- rNA
    # dat$epi$prev.rct.hivneg <- rNA
    # dat$epi$prev.uct.hivneg <- rNA

    #HIV/STI coinfection prevalence
    # dat$epi$prev.rct.uct <- rNA
    # dat$epi$prev.rgc.ugc <- rNA
    # dat$epi$prev.rct.rgc <- rNA
    # dat$epi$prev.uct.ugc <- rNA
    # dat$epi$prev.gc.syph <- rNA
    # dat$epi$prev.ct.syph <- rNA
    # dat$epi$prev.gc.primsecosyph <- rNA
    # dat$epi$prev.ct.primsecosyph <- rNA

    # STI incidence
    dat$epi$incid.rgc <- rNA
    dat$epi$incid.ugc <- rNA
    dat$epi$incid.gc <- rNA
    dat$epi$incid.gc.g1 <- rNA
    dat$epi$incid.gc.g2 <- rNA
    dat$epi$incid.rct <- rNA
    dat$epi$incid.uct <- rNA
    dat$epi$incid.ct <- rNA
    dat$epi$incid.ct.g1 <- rNA
    dat$epi$incid.ct.g2 <- rNA
    dat$epi$incid.syph <- rNA
    dat$epi$incid.syph.g1 <- rNA
    dat$epi$incid.syph.g2 <- rNA
    dat$epi$incid.sti <- rNA
    dat$epi$incid.sti.g1 <- rNA
    dat$epi$incid.sti.g2 <- rNA

    dat$epi$incid.gc.hivneg <- rNA
    dat$epi$incid.gc.hivpos <- rNA
    dat$epi$incid.ct.hivneg <- rNA
    dat$epi$incid.ct.hivpos <- rNA
    dat$epi$incid.syph.hivneg <- rNA
    dat$epi$incid.syph.hivpos <- rNA

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
    dat$epi$ir100.gcct <- rNA
    dat$epi$ir100.gcct.tttraj1 <- rNA
    dat$epi$ir100.gcct.tttraj2 <- rNA
    dat$epi$incid.gcct <- rNA
    dat$epi$incid.gcct.tttraj1 <- rNA
    dat$epi$incid.gcct.tttraj2 <- rNA

    dat$epi$ir100.gc.hivneg <- rNA
    dat$epi$ir100.gc.hivpos <- rNA
    dat$epi$ir100.ct.hivneg <- rNA
    dat$epi$ir100.ct.hivpos <- rNA
    dat$epi$ir100.syph.hivneg <- rNA
    dat$epi$ir100.syph.hivpos <- rNA

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

    # STI Recovery
    dat$epi$recov.rgc <- rNA
    dat$epi$recov.ugc <- rNA
    dat$epi$recov.rct <- rNA
    dat$epi$recov.uct <- rNA
    dat$epi$recov.earlysyph <- rNA
    dat$epi$recov.syphilis <- rNA

    # HIV transmissions by partner type
    dat$epi$trans.main <- rNA
    dat$epi$trans.pers <- rNA
    dat$epi$trans.inst <- rNA

    # STI treatment
    dat$epi$txGC <- rNA
    dat$epi$txGC_asympt <- rNA
    dat$epi$txCT <- rNA
    dat$epi$txCT_asympt <- rNA
    dat$epi$txsyph <- rNA
    dat$epi$txsyph_asympt <- rNA
    dat$epi$txearlysyph <- rNA
    dat$epi$txlatesyph <- rNA

    # STI testing indications
    dat$epi$stiactiveind <- rNA
    dat$epi$recentpartners <- rNA
    # dat$epi$recentSTI <- rNA
    # dat$epi$newpartner <- rNA
    # dat$epi$concurrpart <- rNA
    # dat$epi$partnersti <- rNA
    # dat$epi$uai.nmain <- rNA
    # dat$epi$uai.any <- rNA
    dat$epi$stiactiveind.prop <- rNA
    dat$epi$recentpartners.prop <- rNA

    #EPT
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

    # STI testing trajectories
    dat$epi$tt.traj.syph1.hivneg <- rep(0, nsteps)
    dat$epi$tt.traj.gc1.hivneg <- rep(0, nsteps)
    dat$epi$tt.traj.ct1.hivneg <- rep(0, nsteps)
    dat$epi$tt.traj.syph2.hivneg <- rep(0, nsteps)
    dat$epi$tt.traj.gc2.hivneg <- rep(0, nsteps)
    dat$epi$tt.traj.ct2.hivneg <- rep(0, nsteps)
    dat$epi$tt.traj.syph1.hivpos <- rep(0, nsteps)
    dat$epi$tt.traj.gc1.hivpos <- rep(0, nsteps)
    dat$epi$tt.traj.ct1.hivpos <- rep(0, nsteps)
    dat$epi$tt.traj.syph2.hivpos <- rep(0, nsteps)
    dat$epi$tt.traj.gc2.hivpos <- rep(0, nsteps)
    dat$epi$tt.traj.ct2.hivpos <- rep(0, nsteps)
    dat$epi$tt.traj.syph1 <- rep(0, nsteps)
    dat$epi$tt.traj.gc1 <- rep(0, nsteps)
    dat$epi$tt.traj.ct1 <- rep(0, nsteps)
    dat$epi$tt.traj.sti1 <- rep(0, nsteps)
    dat$epi$tt.traj.syph2 <- rep(0, nsteps)
    dat$epi$tt.traj.gc2 <- rep(0, nsteps)
    dat$epi$tt.traj.ct2 <- rep(0, nsteps)
    dat$epi$tt.traj.sti2 <- rep(0, nsteps)

    #STI Testing in last 12 months
    dat$epi$test.gc.12mo <- rNA
    dat$epi$test.ct.12mo <- rNA
    dat$epi$test.syph.12mo <- rNA

    dat$epi$test.gc.12mo.nonhivdiag <- rNA
    dat$epi$test.ct.12mo.nonhivdiag <- rNA
    dat$epi$test.syph.12mo.nonhivdiag <- rNA

    dat$epi$test.gc.12mo.hivneg <- rNA
    dat$epi$test.ct.12mo.hivneg <- rNA
    dat$epi$test.syph.12mo.hivneg <- rNA

    dat$epi$test.gc.12mo.hivdiag <- rNA
    dat$epi$test.ct.12mo.hivdiag <- rNA
    dat$epi$test.syph.12mo.hivdiag <- rNA

    dat$epi$test.gc.12mo.hivpos <- rNA
    dat$epi$test.ct.12mo.hivpos <- rNA
    dat$epi$test.syph.12mo.hivpos <- rNA

    # Incidence by risk group
    dat$epi$ir100.ct.tttraj1 <- rNA
    dat$epi$ir100.ct.tttraj2 <- rNA
    dat$epi$ir100.gc.tttraj1 <- rNA
    dat$epi$ir100.gc.tttraj2 <- rNA
    dat$epi$ir100.syph.tttraj1 <- rNA
    dat$epi$ir100.syph.tttraj2 <- rNA
    dat$epi$ir100.sti.tttraj1 <- rNA
    dat$epi$ir100.sti.tttraj2 <- rNA
    dat$epi$incid.gc.tttraj1 <- rNA
    dat$epi$incid.gc.tttraj2 <- rNA
    dat$epi$incid.ct.tttraj1 <- rNA
    dat$epi$incid.ct.tttraj2 <- rNA
    dat$epi$incid.syph.tttraj1 <- rNA
    dat$epi$incid.syph.tttraj2 <- rNA
    dat$epi$incid.sti.tttraj1 <- rNA
    dat$epi$incid.sti.tttraj2 <- rNA

    # Prevalence by risk group
    dat$epi$prev.gcct.tttraj1 <- rNA
    dat$epi$prev.gcct.tttraj2 <- rNA
    dat$epi$prev.gc.tttraj1 <- rNA
    dat$epi$prev.gc.tttraj2 <- rNA
    dat$epi$prev.ct.tttraj1 <- rNA
    dat$epi$prev.ct.tttraj2 <- rNA
    dat$epi$prev.syph.tttraj1 <- rNA
    dat$epi$prev.syph.tttraj2 <- rNA
    dat$epi$prev.primsecosyph.tttraj1 <- rNA
    dat$epi$prev.primsecosyph.tttraj2 <- rNA
    dat$epi$prev.sti.tttraj1 <- rNA
    dat$epi$prev.sti.tttraj2 <- rNA

    # Tests by risk group
    dat$epi$rCTasympttests.tttraj1 <- rNA
    dat$epi$rCTasympttests.tttraj2 <- rNA
    dat$epi$uCTasympttests.tttraj1 <- rNA
    dat$epi$uCTasympttests.tttraj2 <- rNA
    dat$epi$CTasympttests.tttraj1 <- rNA
    dat$epi$CTasympttests.tttraj2 <- rNA
    dat$epi$rGCasympttests.tttraj1 <- rNA
    dat$epi$rGCasympttests.tttraj2 <- rNA
    dat$epi$uGCasympttests.tttraj1 <- rNA
    dat$epi$uGCasympttests.tttraj2 <- rNA
    dat$epi$GCasympttests.tttraj1 <- rNA
    dat$epi$GCasympttests.tttraj2 <- rNA
    dat$epi$syphasympttests.tttraj1 <- rNA
    dat$epi$syphasympttests.tttraj2 <- rNA
    dat$epi$stiasympttests.tttraj1 <- rNA
    dat$epi$stiasympttests.tttraj2 <- rNA
    dat$epi$rCTsympttests.tttraj1 <- rNA
    dat$epi$rCTsympttests.tttraj2 <- rNA
    dat$epi$uCTsympttests.tttraj1 <- rNA
    dat$epi$uCTsympttests.tttraj2 <- rNA
    dat$epi$CTsympttests.tttraj1 <- rNA
    dat$epi$CTsympttests.tttraj2 <- rNA
    dat$epi$rGCsympttests.tttraj1 <- rNA
    dat$epi$rGCsympttests.tttraj2 <- rNA
    dat$epi$uGCsympttests.tttraj1 <- rNA
    dat$epi$uGCsympttests.tttraj2 <- rNA
    dat$epi$GCsympttests.tttraj1 <- rNA
    dat$epi$GCsympttests.tttraj2 <- rNA
    dat$epi$syphsympttests.tttraj1 <- rNA
    dat$epi$syphsympttests.tttraj2 <- rNA
    dat$epi$stisympttests.tttraj1 <- rNA
    dat$epi$stisympttests.tttraj2 <- rNA

    # Treatments by risk group
    dat$epi$txGC.tttraj1 <- rNA
    dat$epi$txGC_asympt.tttraj1 <- rNA
    dat$epi$txGC.tttraj2 <- rNA
    dat$epi$txGC_asympt.tttraj2 <- rNA
    dat$epi$txCT.tttraj1 <- rNA
    dat$epi$txCT_asympt.tttraj1 <- rNA
    dat$epi$txCT.tttraj2 <- rNA
    dat$epi$txCT_asympt.tttraj2 <- rNA
    dat$epi$txsyph.tttraj1 <- rNA
    dat$epi$txsyph_asympt.tttraj1 <- rNA
    dat$epi$txsyph.tttraj2 <- rNA
    dat$epi$txsyph_asympt.tttraj2 <- rNA
    dat$epi$txearlysyph.tttraj1 <- rNA
    dat$epi$txlatesyph.tttraj1 <- rNA
    dat$epi$txearlysyph.tttraj2 <- rNA
    dat$epi$txlatesyph.tttraj2 <- rNA
    dat$epi$txSTI <- rNA
    dat$epi$txSTI_asympt  <- rNA
    dat$epi$txSTI.tttraj1 <- rNA
    dat$epi$txSTI.tttraj2 <- rNA
    dat$epi$txSTI_asympt.tttraj1 <- rNA
    dat$epi$txSTI_asympt.tttraj2 <- rNA

    # Proportion of infections treated in past year
    dat$epi$tx.gc.prop <- rNA
    dat$epi$tx.ct.prop <- rNA
    dat$epi$tx.gcct.prop <- rNA
    dat$epi$tx.syph.prop <- rNA

    # Duration of infection
    dat$epi$gc.infect.dur <- rNA
    dat$epi$ct.infect.dur <- rNA
    dat$epi$gcct.infect.dur <- rNA
    dat$epi$syph.infect.dur <- rNA

    # UAI by concordancy
    dat$epi$num.acts.negneg <- rNA
    dat$epi$num.acts.negpos <- rNA
    dat$epi$num.acts.pospos <- rNA
    dat$epi$prop.uai.negneg <- rNA
    dat$epi$prop.uai.negpos <- rNA
    dat$epi$prop.uai.pospos <- rNA
    dat$epi$prop.acts.negneg <- rNA
    dat$epi$prop.acts.negpos <- rNA
    dat$epi$prop.acts.pospos <- rNA

  }
#
#   # Add new risk-group specific trackers if needed
#   if (is.null(dat$epi$ir100.ct.tttraj1)) {
#
#     # Denominators of STI testing trajectories
#     dat$epi$tt.traj.syph1.hivneg <- rep(0, nsteps)
#     dat$epi$tt.traj.gc1.hivneg <- rep(0, nsteps)
#     dat$epi$tt.traj.ct1.hivneg <- rep(0, nsteps)
#     dat$epi$tt.traj.syph2.hivneg <- rep(0, nsteps)
#     dat$epi$tt.traj.gc2.hivneg <- rep(0, nsteps)
#     dat$epi$tt.traj.ct2.hivneg <- rep(0, nsteps)
#     dat$epi$tt.traj.syph1.hivpos <- rep(0, nsteps)
#     dat$epi$tt.traj.gc1.hivpos <- rep(0, nsteps)
#     dat$epi$tt.traj.ct1.hivpos <- rep(0, nsteps)
#     dat$epi$tt.traj.syph2.hivpos <- rep(0, nsteps)
#     dat$epi$tt.traj.gc2.hivpos <- rep(0, nsteps)
#     dat$epi$tt.traj.ct2.hivpos <- rep(0, nsteps)
#     dat$epi$tt.traj.syph1 <- rep(0, nsteps)
#     dat$epi$tt.traj.gc1 <- rep(0, nsteps)
#     dat$epi$tt.traj.ct1 <- rep(0, nsteps)
#     dat$epi$tt.traj.syph2 <- rep(0, nsteps)
#     dat$epi$tt.traj.gc2 <- rep(0, nsteps)
#     dat$epi$tt.traj.ct2 <- rep(0, nsteps)
#     dat$epi$tt.traj.sti1 <- rep(0, nsteps)
#     dat$epi$tt.traj.sti2 <- rep(0, nsteps)
#
#     # Incidence by risk group
#     dat$epi$ir100.ct.tttraj1 <- rNA
#     dat$epi$ir100.ct.tttraj2 <- rNA
#     dat$epi$ir100.gc.tttraj1 <- rNA
#     dat$epi$ir100.gc.tttraj2 <- rNA
#     dat$epi$ir100.syph.tttraj1 <- rNA
#     dat$epi$ir100.syph.tttraj2 <- rNA
#     dat$epi$ir100.sti.tttraj1 <- rNA
#     dat$epi$ir100.sti.tttraj2 <- rNA
#     dat$epi$incid.gc.tttraj1 <- rNA
#     dat$epi$incid.gc.tttraj2 <- rNA
#     dat$epi$incid.ct.tttraj1 <- rNA
#     dat$epi$incid.ct.tttraj2 <- rNA
#     dat$epi$incid.syph.tttraj1 <- rNA
#     dat$epi$incid.syph.tttraj2 <- rNA
#     dat$epi$incid.sti.tttraj1 <- rNA
#     dat$epi$incid.sti.tttraj2 <- rNA
#
#     # Prevalence by risk group
#     dat$epi$prev.gc.tttraj1 <- rNA
#     dat$epi$prev.gc.tttraj2 <- rNA
#     dat$epi$prev.ct.tttraj1 <- rNA
#     dat$epi$prev.ct.tttraj2 <- rNA
#     dat$epi$prev.syph.tttraj1 <- rNA
#     dat$epi$prev.syph.tttraj2 <- rNA
#     dat$epi$prev.primsecosyph.tttraj1 <- rNA
#     dat$epi$prev.primsecosyph.tttraj2 <- rNA
#     dat$epi$prev.sti.tttraj1 <- rNA
#     dat$epi$prev.sti.tttraj2 <- rNA
#
#     # Tests by risk group
#     dat$epi$rCTasympttests.tttraj1 <- rNA
#     dat$epi$rCTasympttests.tttraj2 <- rNA
#     dat$epi$uCTasympttests.tttraj1 <- rNA
#     dat$epi$uCTasympttests.tttraj2 <- rNA
#     dat$epi$CTasympttests.tttraj1 <- rNA
#     dat$epi$CTasympttests.tttraj2 <- rNA
#     dat$epi$rGCasympttests.tttraj1 <- rNA
#     dat$epi$rGCasympttests.tttraj2 <- rNA
#     dat$epi$uGCasympttests.tttraj1 <- rNA
#     dat$epi$uGCasympttests.tttraj2 <- rNA
#     dat$epi$GCasympttests.tttraj1 <- rNA
#     dat$epi$GCasympttests.tttraj2 <- rNA
#     dat$epi$syphasympttests.tttraj1 <- rNA
#     dat$epi$syphasympttests.tttraj2 <- rNA
#     dat$epi$stiasympttests.tttraj1 <- rNA
#     dat$epi$stiasympttests.tttraj2 <- rNA
#     dat$epi$rCTsympttests.tttraj1 <- rNA
#     dat$epi$rCTsympttests.tttraj2 <- rNA
#     dat$epi$uCTsympttests.tttraj1 <- rNA
#     dat$epi$uCTsympttests.tttraj2 <- rNA
#     dat$epi$CTsympttests.tttraj1 <- rNA
#     dat$epi$CTsympttests.tttraj2 <- rNA
#     dat$epi$rGCsympttests.tttraj1 <- rNA
#     dat$epi$rGCsympttests.tttraj2 <- rNA
#     dat$epi$uGCsympttests.tttraj1 <- rNA
#     dat$epi$uGCsympttests.tttraj2 <- rNA
#     dat$epi$GCsympttests.tttraj1 <- rNA
#     dat$epi$GCsympttests.tttraj2 <- rNA
#     dat$epi$syphsympttests.tttraj1 <- rNA
#     dat$epi$syphsympttests.tttraj2 <- rNA
#     dat$epi$stisympttests.tttraj1 <- rNA
#     dat$epi$stisympttests.tttraj2 <- rNA
#
#     # Treatments by risk group
#     dat$epi$txGC.tttraj1 <- rNA
#     dat$epi$txGC_asympt.tttraj1 <- rNA
#     dat$epi$txGC.tttraj2 <- rNA
#     dat$epi$txGC_asympt.tttraj2 <- rNA
#     dat$epi$txCT.tttraj1 <- rNA
#     dat$epi$txCT_asympt.tttraj1 <- rNA
#     dat$epi$txCT.tttraj2 <- rNA
#     dat$epi$txCT_asympt.tttraj2 <- rNA
#     dat$epi$txsyph.tttraj1 <- rNA
#     dat$epi$txsyph_asympt.tttraj1 <- rNA
#     dat$epi$txsyph.tttraj2 <- rNA
#     dat$epi$txsyph_asympt.tttraj2 <- rNA
#     dat$epi$txearlysyph.tttraj1 <- rNA
#     dat$epi$txlatesyph.tttraj1 <- rNA
#     dat$epi$txearlysyph.tttraj2 <- rNA
#     dat$epi$txlatesyph.tttraj2 <- rNA
#     dat$epi$txSTI <- rNA
#     dat$epi$txSTI_asympt  <- rNA
#     dat$epi$txSTI.tttraj1 <- rNA
#     dat$epi$txSTI.tttraj2 <- rNA
#     dat$epi$txSTI_asympt.tttraj1 <- rNA
#     dat$epi$txSTI_asympt.tttraj2 <- rNA
#
#     # Proportion of infections treated in past year
#     dat$epi$tx.gc.prop <- rNA
#     dat$epi$tx.ct.prop <- rNA
#     dat$epi$tx.gcct.prop <- rNA
#     dat$epi$tx.syph.prop <- rNA
#
#     # Duration of infection
#     dat$epi$gc.infect.dur <- rNA
#     dat$epi$ct.infect.dur <- rNA
#     dat$epi$gcct.infect.dur <- rNA
#     dat$epi$syph.infect.dur <- rNA
#
#     dat$epi$incid.gc.hivneg <- rNA
#     dat$epi$incid.gc.hivpos <- rNA
#     dat$epi$incid.ct.hivneg <- rNA
#     dat$epi$incid.ct.hivpos <- rNA
#     dat$epi$incid.syph.hivneg <- rNA
#     dat$epi$incid.syph.hivpos <- rNA
#
#     dat$epi$ir100.gc.hivneg <- rNA
#     dat$epi$ir100.gc.hivpos <- rNA
#     dat$epi$ir100.ct.hivneg <- rNA
#     dat$epi$ir100.ct.hivpos <- rNA
#     dat$epi$ir100.syph.hivneg <- rNA
#     dat$epi$ir100.syph.hivpos <- rNA
#
#   }

  # Population sizes and HIV incidence/prevalence
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

  # PrEP
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

  # STI and Co-infection Prevalence
  dat$epi$prev.rgc[at] <- ifelse(dat$epi$num[at] > 0, sum(rGC == 1, na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.ugc[at] <- ifelse(dat$epi$num[at] > 0, sum(uGC == 1, na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.gc[at] <- ifelse(dat$epi$num[at] > 0, sum((rGC == 1 | uGC == 1), na.rm = TRUE) / dat$epi$num[at], 0)

  dat$epi$prev.gc.tttraj1[at] <- ifelse((dat$epi$tt.traj.gc1[at] == 0 | is.na(dat$epi$tt.traj.gc1[at])
                                         | is.nan(dat$epi$tt.traj.gc1[at]) | is.null(dat$epi$tt.traj.gc1[at])), 0,
                                        sum((rGC == 1 | uGC == 1) &
                                             (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) /
                                          dat$epi$tt.traj.gc1[at])
  dat$epi$prev.gc.tttraj2[at] <- ifelse((dat$epi$tt.traj.gc2[at] == 0 | is.na(dat$epi$tt.traj.gc2[at]) |
                                           is.nan(dat$epi$tt.traj.gc2[at]) | is.null(dat$epi$tt.traj.gc2[at])), 0 ,
                                        sum((rGC == 1 | uGC == 1) &
                                             (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) /
                                          dat$epi$tt.traj.gc2[at])

  dat$epi$prev.gcct[at] <- ifelse(dat$epi$num[at] > 0, sum((rGC == 1 | uGC == 1 | rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.gcct.tttraj1[at] <- ifelse((dat$epi$tt.traj.gc1[at] == 0 | is.na(dat$epi$tt.traj.gc1[at])
                                         | is.nan(dat$epi$tt.traj.gc1[at]) | is.null(dat$epi$tt.traj.gc1[at])), 0,
                                        sum((rGC == 1 | uGC == 1 | rCT == 1 | uCT == 1) &
                                              (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) /
                                          dat$epi$tt.traj.gc1[at])
  dat$epi$prev.gcct.tttraj2[at] <- ifelse((dat$epi$tt.traj.gc2[at] == 0 | is.na(dat$epi$tt.traj.gc2[at]) |
                                           is.nan(dat$epi$tt.traj.gc2[at]) | is.null(dat$epi$tt.traj.gc2[at])), 0 ,
                                        sum((rGC == 1 | uGC == 1 | rCT == 1 | uCT == 1) &
                                              (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) /
                                          dat$epi$tt.traj.gc2[at])
  # dat$epi$prev.gc.sympt[at] <- ifelse(dat$epi$num[at] > 0, sum((rGC.sympt == 1 | uGC.sympt == 1)) / dat$epi$num[at], 0)
  # dat$epi$prev.gc.dual[at] <- ifelse(dat$epi$num[at] > 0, sum((rGC == 1 & uGC == 1), na.rm = TRUE) / dat$epi$num[at], 0)

  dat$epi$prev.rgc.hivneg.only[at] <-  length(which(status == 0 & rGC == 1 & uGC == 0 & rCT == 0 & uCT == 0 & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$s.num[at]
  dat$epi$prev.ugc.hivneg.only[at] <- length(which(status == 0 & rGC == 0 & uGC == 1 & rCT == 0 & uCT == 0 & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$s.num[at]
  # dat$epi$prev.gc.hivneg.only[at] <-  length(which(status == 0 & (rGC == 1 | uGC == 1) & rCT == 0 & uCT == 0 & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$s.num[at]

  dat$epi$prev.rct[at] <- ifelse(dat$epi$num[at] > 0, sum(rCT == 1, na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.uct[at] <- ifelse(dat$epi$num[at] > 0, sum(uCT == 1, na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.ct[at] <- ifelse(dat$epi$num[at] > 0, sum((rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at], 0)
  dat$epi$prev.ct.tttraj1[at] <- ifelse((dat$epi$tt.traj.ct1[at] == 0 | is.na(dat$epi$tt.traj.ct1[at]) |
                                           is.nan(dat$epi$tt.traj.ct1[at]) | is.null(dat$epi$tt.traj.ct1[at])), 0,
                                        sum((rCT == 1 | uCT == 1) &
                                        (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE) /
                                          dat$epi$tt.traj.ct1[at])
  dat$epi$prev.ct.tttraj2[at] <- ifelse((dat$epi$tt.traj.ct2[at] == 0 | is.na(dat$epi$tt.traj.ct2[at]) |
                                           is.nan(dat$epi$tt.traj.ct2[at]) | is.null(dat$epi$tt.traj.ct2[at])), 0,
                                        sum((rCT == 1 | uCT == 1) &
                                                                   (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE) /
                                          dat$epi$tt.traj.ct2[at])
  # dat$epi$prev.ct.sympt[at] <- ifelse(dat$epi$num[at] > 0, sum((rCT.sympt == 1 | uCT.sympt == 1)) / dat$epi$num[at], 0)
  # dat$epi$prev.ct.dual[at] <- ifelse(dat$epi$num[at] > 0, sum((rCT == 1 & uCT == 1), na.rm = TRUE) / dat$epi$num[at], 0)

  dat$epi$prev.rct.hivneg.only[at] <-  length(which(status == 0 & rGC == 0 & uGC == 0 & rCT == 1 & uCT == 0 & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$s.num[at]
  dat$epi$prev.uct.hivneg.only[at] <- length(which(status == 0 & rGC == 0 & uGC == 0 & rCT == 0 & uCT == 1 & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$s.num[at]
  # dat$epi$prev.ct.hivneg.only[at] <-  length(which(status == 0 & rGC == 0 & uGC == 0 & (rCT == 1 | uCT == 1) & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$s.num[at]

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
  dat$epi$prev.stage.tert[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                        length(which(stage.syph == 6)) / length(which(syphilis == 1)), 0)
  dat$epi$prev.earlysyph[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                       length(which(stage.syph %in% c(1, 2, 3))) / length(which(syphilis == 1)), 0)
  dat$epi$prev.latesyph[at] <- ifelse(length(which(syphilis == 1)) > 0,
                                      length(which(stage.syph %in% c(4, 5, 6))) / length(which(syphilis == 1)), 0)
  dat$epi$num.newearlydiagsyph[at] <- length(which(last.diag.time.syph == at & stage.syph %in% c(1, 2, 3)))
  dat$epi$num.newlatediagsyph[at] <- length(which(last.diag.time.syph == at & stage.syph %in% c(4, 5, 6)))
  dat$epi$early.late.syphratio[at] <- ifelse(length(which(stage.syph %in% c(4, 5, 6))) > 0,
                                             length(which(stage.syph %in% c(1, 2, 3))) /
                                               length(which(stage.syph %in% c(4, 5, 6))), 0)
  dat$epi$early.late.diagsyphratio[at] <- ifelse(length(which(diag.status.syph == 1 & stage.syph %in% c(4, 5, 6))),
                                                 length(which(diag.status.syph == 1 & stage.syph %in% c(1, 2, 3))) /
                                                   length(which(diag.status.syph == 1 & stage.syph %in% c(4, 5, 6))), 0)

  dat$epi$prev.dxhiv.dxipssyph[at] <- ifelse(length(which(diag.status.syph == 1 & stage.syph %in% c(1, 2, 3))) == 0, 0,
                                                  length(which(diag.status == 1 &
                                                                 diag.status.syph == 1 &
                                                                 stage.syph %in% c(1, 2, 3))) /
                                            length(which(diag.status.syph == 1 & stage.syph %in% c(1, 2, 3))))

  dat$epi$prev.dxhiv.atdxipssyph[at] <- ifelse(length(which(diag.status.syph == 1 & stage.syph %in% c(1, 2, 3) &
                                                              dat$attr$last.diag.time.syph == at)) == 0, 0,
                                             length(which(diag.status == 1 &
                                                            diag.status.syph == 1 &
                                                            stage.syph %in% c(1, 2, 3) &
                                                            dat$attr$last.diag.time.syph == at)) /
                                               length(which(diag.status.syph == 1 & stage.syph %in% c(1, 2, 3)  &
                                                              dat$attr$last.diag.time.syph == at)))

  dat$epi$prev.syph[at] <- ifelse(dat$epi$num[at] > 0, length(which(syphilis == 1)) / dat$epi$num[at], 0)
  dat$epi$prev.syph.tttraj1[at] <- ifelse((dat$epi$tt.traj.syph1[at] == 0 | is.na(dat$epi$tt.traj.syph1[at]) |
                                             is.nan(dat$epi$tt.traj.syph1[at]) | is.null(dat$epi$tt.traj.syph1[at])), 0,
                                          sum((syphilis == 1) &
                                               (tt.traj.syph.hivneg == 1 | tt.traj.syph.hivpos == 1), na.rm = TRUE) /
                                            dat$epi$tt.traj.syph1[at])
  dat$epi$prev.syph.tttraj2[at] <- ifelse((dat$epi$tt.traj.syph2[at] == 0 | is.na(dat$epi$tt.traj.syph2[at]) |
                                             is.nan(dat$epi$tt.traj.syph2[at]) | is.null(dat$epi$tt.traj.syph2[at])), 0,
                                          sum((syphilis == 1) &
                                               (tt.traj.syph.hivneg == 2 | tt.traj.syph.hivpos == 2), na.rm = TRUE) / dat$epi$tt.traj.syph2[at])

  dat$epi$prev.primsecosyph[at] <- ifelse(dat$epi$num[at] > 0, length(which(stage.syph %in% c(1, 2, 3))) /
                                            dat$epi$num[at], 0)
  dat$epi$prev.primsecosyph.tttraj1[at] <- ifelse((dat$epi$tt.traj.syph1[at] == 0 | is.na(dat$epi$tt.traj.syph1[at]) |
                                                     is.nan(dat$epi$tt.traj.syph1[at]) | is.null(dat$epi$tt.traj.syph1[at])), 0,
                                                  length(which(stage.syph %in% c(1, 2, 3) &
                                                                   (tt.traj.syph.hivneg == 1 | tt.traj.syph.hivpos == 1))) /
                                                    dat$epi$tt.traj.syph1[at])
  dat$epi$prev.primsecosyph.tttraj2[at] <- ifelse((dat$epi$tt.traj.syph2[at] == 0 | is.na(dat$epi$tt.traj.syph2[at]) |
                                                     is.nan(dat$epi$tt.traj.syph2[at]) | is.null(dat$epi$tt.traj.syph2[at])), 0,
                                                  length(which(stage.syph %in% c(1, 2, 3) &
                                                          (tt.traj.syph.hivneg == 2 | tt.traj.syph.hivpos == 2))) /
                                                    dat$epi$tt.traj.syph2[at])

  dat$epi$prev.primsecosyph.hivneg.only[at] <-  length(which(status == 0 & rGC == 0 & uGC == 0 & rCT == 0 & uCT == 0 & stage.syph %in% c(1, 2, 3))) / dat$epi$s.num[at]

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

  # dat$epi$prev.rgc.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
  #                                       length(intersect(which(status == 1), which(rGC == 1))) / dat$epi$i.num[at], 0)
  # dat$epi$prev.ugc.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
  #                                       length(intersect(which(status == 1),  which(uGC == 1))) / dat$epi$i.num[at], 0)
  # dat$epi$prev.rct.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
  #                                       length(intersect(which(status == 1), which(rCT == 1))) / dat$epi$i.num[at], 0)
  # dat$epi$prev.uct.hivpos[at] <- ifelse(dat$epi$i.num[at] > 0,
  #                                       length(intersect(which(status == 1), which(uCT == 1))) / dat$epi$i.num[at], 0)
  #
  # dat$epi$prev.rgc.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
  #                                       length(intersect(which(status == 1), which(rGC == 1))) / dat$epi$s.num[at], 0)
  # dat$epi$prev.ugc.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
  #                                       length(intersect(which(status == 1), which(uGC == 1))) / dat$epi$s.num[at], 0)
  # dat$epi$prev.rct.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
  #                                       length(intersect(which(status == 1), which(rCT == 1))) / dat$epi$s.num[at], 0)
  # dat$epi$prev.uct.hivneg[at] <- ifelse(dat$epi$s.num[at] > 0,
  #                                       length(intersect(which(status == 1), which(uCT == 1))) / dat$epi$s.num[at], 0)

  # STI Co-infection prevalence
  # dat$epi$prev.rct.uct[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rCT == 1), which(uCT == 1)))) / dat$epi$num[at], 0)
  # dat$epi$prev.rgc.ugc[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rGC == 1), which(uGC == 1)))) / dat$epi$num[at], 0)
  # dat$epi$prev.rct.rgc[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rGC == 1), which(rCT == 1))))  / dat$epi$num[at], 0)
  # dat$epi$prev.uct.ugc[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(uCT == 1), which(uGC == 1)))) / dat$epi$num[at], 0)
  # dat$epi$prev.gc.syph[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rGC == 1 | uGC == 1), which(syphilis == 1)))) / dat$epi$num[at], 0)
  # dat$epi$prev.ct.syph[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rCT == 1 | uCT == 1), which(syphilis == 1)))) / dat$epi$num[at], 0)
  # dat$epi$prev.gc.primsecosyph[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rGC == 1 | uGC == 1), which(stage.syph %in% c(1, 2))))) / dat$epi$num[at], 0)
  # dat$epi$prev.ct.primsecosyph[at] <- ifelse(dat$epi$num[at] > 0, (length(intersect(which(rCT == 1 | uCT == 1), which(stage.syph %in% c(1, 2))))) / dat$epi$num[at], 0)

  dat$epi$prev.hivnegmultsti[at] <- sum(status == 0 &
                                    (rGC == 1 & (uGC == 1 | rCT == 1 | uCT == 1 | stage.syph %in% c(1, 2, 3))) |
                                    (uGC == 1 & (rGC == 1 | rCT == 1 | uCT == 1 | stage.syph %in% c(1, 2, 3))) |
                                    (rCT == 1 & (uGC == 1 | rGC == 1 | uCT == 1 | stage.syph %in% c(1, 2, 3))) |
                                    (uCT == 1 & (uGC == 1 | rGC == 1 | rCT == 1 | stage.syph %in% c(1, 2, 3))) |
                                    (stage.syph %in% c(1, 2, 3) & (uGC == 1 | rGC == 1 | rCT == 1 | uCT == 1))) / dat$epi$s.num[at]

  dat$epi$prev.rgc.hivpos.only[at] <-  length(which(status == 1 & rGC == 1 & uGC == 0 & rCT == 0 & uCT == 0 & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$i.num[at]
  dat$epi$prev.ugc.hivpos.only[at] <- length(which(status == 1 & rGC == 0 & uGC == 1 & rCT == 0 & uCT == 0 & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$i.num[at]
  # dat$epi$prev.gc.hivpos.only[at] <-  length(which(status == 1 & (rGC == 1 | uGC == 1) & rCT == 0 & uCT == 0 & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$i.num[at]
  dat$epi$prev.rct.hivpos.only[at] <-  length(which(status == 1 & rGC == 0 & uGC == 0 & rCT == 1 & uCT == 0 & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$i.num[at]
  dat$epi$prev.uct.hivpos.only[at] <- length(which(status == 1 & rGC == 0 & uGC == 0 & rCT == 0 & uCT == 1 & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$i.num[at]
  # dat$epi$prev.ct.hivpos.only[at] <-  length(which(status == 1 & rGC == 0 & uGC == 0 & (rCT == 1 | uCT == 1) & stage.syph %in% c(NA, 4, 5, 6))) / dat$epi$i.num[at]
  dat$epi$prev.primsecosyph.hivpos.only[at] <-  length(which(status == 1 & rGC == 0 & uGC == 0 & rCT == 0 & uCT == 0 & stage.syph %in% c(1, 2, 3))) / dat$epi$i.num[at]

  #HIV/Multiple STI
  dat$epi$prev.hivposmultsti[at] <- sum(status == 1 &
                                    ((rGC == 1 & (uGC == 1 | rCT == 1 | uCT == 1 | stage.syph %in% c(1, 2, 3))) |
                                    (uGC == 1 & (rGC == 1 | rCT == 1 | uCT == 1 | stage.syph %in% c(1, 2, 3))) |
                                    (rCT == 1 & (uGC == 1 | rGC == 1 | uCT == 1 | stage.syph %in% c(1, 2, 3))) |
                                    (uCT == 1 & (uGC == 1 | rGC == 1 | rCT == 1 | stage.syph %in% c(1, 2, 3))) |
                                    (stage.syph %in% c(1, 2, 3) & (uGC == 1 | rGC == 1 | rCT == 1 | uCT == 1)))) / dat$epi$i.num[at]

  # Site-specific STI incidence rates
  dat$epi$ir100.rgc[at] <- ifelse(sum(rGC == 0, na.rm = TRUE) > 0, (dat$epi$incid.rgc[at] / sum(rGC == 0, na.rm = TRUE)) * 5200, 0)
  dat$epi$ir100.ugc[at] <- ifelse(sum(uGC == 0, na.rm = TRUE) > 0, (dat$epi$incid.ugc[at] / sum(uGC == 0, na.rm = TRUE)) * 5200, 0)
  dat$epi$ir100.gc[at] <- ifelse((sum(rGC == 0, na.rm = TRUE) + sum(uGC == 0, na.rm = TRUE)) > 0, (dat$epi$incid.gc[at] / (sum(rGC == 0, na.rm = TRUE) + sum(uGC == 0, na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.gc.hivneg[at] <- ifelse((sum(rGC == 0 & status == 0, na.rm = TRUE) + sum(uGC == 0 & status == 0, na.rm = TRUE)) > 0, (dat$epi$incid.gc.hivneg[at] / (sum(rGC == 0 & status == 0, na.rm = TRUE) + sum(uGC == 0 & status == 0, na.rm = TRUE))) * 5200, 0)
  dat$epi$ir100.gc.hivpos[at] <- ifelse((sum(rGC == 0 & status == 1, na.rm = TRUE) + sum(uGC == 0 & status == 1, na.rm = TRUE)) > 0, (dat$epi$incid.gc.hivpos[at] / (sum(rGC == 0 & status == 1, na.rm = TRUE) + sum(uGC == 0 & status == 1, na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.gc.tttraj1[at] <- ifelse((sum(rGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) +
                                             sum(uGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE)) > 0,
                                          (dat$epi$incid.gc.tttraj1[at] /
                                             (sum(rGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) +
                                                sum(uGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.gc.tttraj2[at] <- ifelse((sum(rGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) +
                                            sum(uGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE)) > 0,
                                         (dat$epi$incid.gc.tttraj2[at] /
                                            (sum(rGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) +
                                               sum(uGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.rct[at] <- ifelse(sum(rCT == 0, na.rm = TRUE) > 0, (dat$epi$incid.rct[at] / sum(rCT == 0, na.rm = TRUE)) * 5200, 0)
  dat$epi$ir100.uct[at] <- ifelse(sum(uCT == 0, na.rm = TRUE) > 0, (dat$epi$incid.uct[at] / sum(uCT == 0, na.rm = TRUE)) * 5200, 0)
  dat$epi$ir100.ct[at] <- ifelse((sum(rCT == 0, na.rm = TRUE) + sum(uCT == 0, na.rm = TRUE)) > 0, (dat$epi$incid.ct[at] / (sum(rCT == 0, na.rm = TRUE) + sum(uCT == 0, na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.ct.hivneg[at] <- ifelse((sum(rCT == 0 & status == 0, na.rm = TRUE) + sum(uCT == 0 & status == 0, na.rm = TRUE)) > 0, (dat$epi$incid.ct.hivneg[at] / (sum(rCT == 0 & status == 0, na.rm = TRUE) + sum(uCT == 0 & status == 0, na.rm = TRUE))) * 5200, 0)
  dat$epi$ir100.ct.hivpos[at] <- ifelse((sum(rCT == 0 & status == 1, na.rm = TRUE) + sum(uCT == 0 & status == 1, na.rm = TRUE)) > 0, (dat$epi$incid.ct.hivpos[at] / (sum(rCT == 0 & status == 1, na.rm = TRUE) + sum(uCT == 0 & status == 1, na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.ct.tttraj1[at] <- ifelse((sum(rCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE) +
                                            sum(uCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE)) > 0,
                                         (dat$epi$incid.ct.tttraj1[at] /
                                            (sum(rCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE) +
                                               sum(uCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.ct.tttraj2[at] <- ifelse((sum(rCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE) +
                                            sum(uCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE)) > 0,
                                         (dat$epi$incid.ct.tttraj2[at] /
                                            (sum(rCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE) +
                                               sum(uCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.syph[at] <- ifelse(sum(syphilis == 0, na.rm = TRUE) > 0, (dat$epi$incid.syph[at] / sum(syphilis == 0 , na.rm = TRUE)) * 5200, 0)

  dat$epi$ir100.syph.tttraj1[at] <- ifelse((sum(syphilis == 0 & (tt.traj.syph.hivneg == 1 | tt.traj.syph.hivpos == 1), na.rm = TRUE)) > 0,
                                         (dat$epi$incid.syph.tttraj1[at] /
                                            (sum(syphilis == 0 & (tt.traj.syph.hivneg == 1 | tt.traj.syph.hivpos == 1), na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.syph.tttraj2[at] <- ifelse((sum(syphilis == 0 & (tt.traj.syph.hivneg == 2 | tt.traj.syph.hivpos == 2), na.rm = TRUE)) > 0,
                                           (dat$epi$incid.syph.tttraj2[at] /
                                              (sum(syphilis == 0 & (tt.traj.syph.hivneg == 2 | tt.traj.syph.hivpos == 2), na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.syph.hivneg[at] <- ifelse(sum(syphilis == 0 & status == 0, na.rm = TRUE) > 0, (dat$epi$incid.syph.hivneg[at] / sum(syphilis == 0 & status == 0 , na.rm = TRUE)) * 5200, 0)
  dat$epi$ir100.syph.hivpos[at] <- ifelse(sum(syphilis == 0 & status == 1, na.rm = TRUE) > 0, (dat$epi$incid.syph.hivpos[at] / sum(syphilis == 0 & status == 1 , na.rm = TRUE)) * 5200, 0)


  dat$epi$prev.sti[at] <- ifelse(sum(rGC == 1 | uGC == 1 | rCT == 1 | uCT == 1 | syphilis == 1 , na.rm = TRUE) > 0,
                                 sum(rGC == 1 | uGC == 1 | rCT == 1 | uCT == 1 | syphilis == 1 , na.rm = TRUE) / dat$epi$num[at], 0)

  dat$epi$prev.sti.tttraj1[at] <- ifelse((dat$epi$tt.traj.sti1[at] == 0 | is.na(dat$epi$tt.traj.sti1[at]) |
                                            is.nan(dat$epi$tt.traj.sti1[at]) | is.null(dat$epi$tt.traj.sti1[at])), 0,
                                         length(which((rGC == 1 | uGC == 1 | rCT == 1 | uCT == 1 | syphilis == 1) &
                                                        (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1 |
                                                           tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1 |
                                                           tt.traj.syph.hivneg == 1 | tt.traj.syph.hivpos == 1))) /
                                           dat$epi$tt.traj.sti1[at])

  dat$epi$prev.sti.tttraj2[at] <- ifelse((dat$epi$tt.traj.sti2[at] == 0 | is.na(dat$epi$tt.traj.sti2[at]) |
                                            is.nan(dat$epi$tt.traj.sti2[at]) | is.null(dat$epi$tt.traj.sti2[at])), 0,
                                 length(which((rGC == 1 | uGC == 1 | rCT == 1 | uCT == 1 | syphilis == 1) &
                                                (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2 |
                                                 tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2 |
                                                 tt.traj.syph.hivneg == 2 | tt.traj.syph.hivpos == 2))) /
                                dat$epi$tt.traj.sti2[at])

  dat$epi$ir100.sti[at] <- ifelse((sum(rGC == 0, na.rm = TRUE) + sum(uGC == 0, na.rm = TRUE) +
                                     sum(rCT == 0, na.rm = TRUE) + sum(uCT == 0, na.rm = TRUE) +
                                     sum(syphilis == 0, na.rm = TRUE)) > 0,
                                ((dat$epi$incid.sti[at]) /
                                (sum(rGC == 0, na.rm = TRUE) + sum(uGC == 0, na.rm = TRUE) +
                                 sum(rCT == 0, na.rm = TRUE) + sum(uCT == 0, na.rm = TRUE) +
                                 sum(syphilis == 0, na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.sti.tttraj1[at] <- ifelse((sum(rGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) +
                                             sum(uGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) +
                                             sum(rCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE) +
                                             sum(uCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE) +
                                             sum(syphilis == 0 & (tt.traj.syph.hivneg == 1 | tt.traj.syph.hivpos == 1), na.rm = TRUE)) > 0,
                                  ((dat$epi$incid.sti.tttraj1[at]) /
                                     (sum(rGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) +
                                        sum(uGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) +
                                        sum(rCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE) +
                                        sum(uCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE) +
                                        sum(syphilis == 0 & (tt.traj.syph.hivneg == 1 | tt.traj.syph.hivpos == 1), na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.sti.tttraj2[at] <- ifelse((sum(rGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) +
                                            sum(uGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) +
                                            sum(rCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE) +
                                            sum(uCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE) +
                                            sum(syphilis == 0 & (tt.traj.syph.hivneg == 2 | tt.traj.syph.hivpos == 2), na.rm = TRUE)) > 0,
                                         ((dat$epi$incid.sti.tttraj2[at]) /
                                            (sum(rGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) +
                                               sum(uGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) +
                                               sum(rCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE) +
                                               sum(uCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE) +
                                               sum(syphilis == 0 & (tt.traj.syph.hivneg == 2 | tt.traj.syph.hivpos == 2), na.rm = TRUE))) * 5200, 0)

  dat$epi$ir100.sti.prep[at] <- ifelse((sum(rGC == 0 & prepStat == 1, na.rm = TRUE) + sum(uGC == 0 & prepStat == 1, na.rm = TRUE) +
                                          sum(rCT == 0 & prepStat == 1, na.rm = TRUE) + sum(uCT == 0 & prepStat == 1, na.rm = TRUE) +
                                          sum(syphilis == 0 & prepStat == 1, na.rm = TRUE)) > 0,
                                  (dat$epi$incid.gcct.prep[at] + dat$epi$incid.syph.prep[at] /
                                  (sum(rGC == 0 & prepStat == 1, na.rm = TRUE) + sum(uGC == 0 & prepStat == 1, na.rm = TRUE) +
                                   sum(rCT == 0 & prepStat == 1, na.rm = TRUE) + sum(uCT == 0 & prepStat == 1, na.rm = TRUE) +
                                   sum(syphilis == 0 & prepStat == 1, na.rm = TRUE))) * 5200, 0)


  # GC/CT incidence
  dat$epi$ir100.gcct[at] <- ifelse((sum(rGC == 0, na.rm = TRUE) + sum(uGC == 0, na.rm = TRUE) + sum(rCT == 0, na.rm = TRUE) + sum(uCT == 0, na.rm = TRUE)) > 0,
                                   ((dat$epi$incid.gc[at] + dat$epi$incid.ct[at]) /
                                      (sum(rGC == 0, na.rm = TRUE) + sum(uGC == 0, na.rm = TRUE) + sum(rCT == 0, na.rm = TRUE) + sum(uCT == 0, na.rm = TRUE))) * 5200,
                                   0)

  dat$epi$ir100.gcct.tttraj1[at] <- ifelse((sum(rGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) +
                                            sum(uGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) +
                                            sum(rCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE) +
                                            sum(uCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE)) > 0,
                                         ((dat$epi$incid.gc.tttraj1[at] + dat$epi$incid.ct.tttraj1[at]) /
                                            (sum(rGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) +
                                               sum(uGC == 0 & (tt.traj.gc.hivneg == 1 | tt.traj.gc.hivpos == 1), na.rm = TRUE) +
                                               sum(rCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE) +
                                               sum(uCT == 0 & (tt.traj.ct.hivneg == 1 | tt.traj.ct.hivpos == 1), na.rm = TRUE))) * 5200,
                                         0)

  dat$epi$ir100.gcct.tttraj2[at] <- ifelse((sum(rGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) +
                                              sum(uGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) +
                                              sum(rCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE) +
                                              sum(uCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE)) > 0,
                                           ((dat$epi$incid.gc.tttraj1[at] + dat$epi$incid.ct.tttraj1[at]) /
                                              (sum(rGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) +
                                                 sum(uGC == 0 & (tt.traj.gc.hivneg == 2 | tt.traj.gc.hivpos == 2), na.rm = TRUE) +
                                                 sum(rCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE) +
                                                 sum(uCT == 0 & (tt.traj.ct.hivneg == 2 | tt.traj.ct.hivpos == 2), na.rm = TRUE))) * 5200,
                                           0)

  # Testing indications
  dat$epi$stiactiveind.prop[at] <- dat$epi$stiactiveind[at] / dat$epi$num[at]
  dat$epi$recentpartners.prop[at] <- dat$epi$recentpartners[at] / dat$epi$num[at]

  # Testing in last 12 months
  # Overall
  dat$epi$test.gc.12mo[at] <- length(which(tslt.rgc <= 52 | tslt.ugc <= 52)) / dat$epi$num[at]
  dat$epi$test.ct.12mo[at] <- length(which(tslt.rct <= 52 | tslt.uct <= 52)) / dat$epi$num[at]
  dat$epi$test.syph.12mo[at] <- length(which(tslt.syph <= 52)) / dat$epi$num[at]

  # Among those HIV-negative or undiagnosed
  dat$epi$test.gc.12mo.nonhivdiag[at] <- length(which((tslt.rgc <= 52 | tslt.ugc <= 52) &
                                           (is.na(diag.status) | diag.status == 0))) / length(which(is.na(diag.status) | diag.status == 0))
  dat$epi$test.ct.12mo.nonhivdiag[at] <- length(which((tslt.rct <= 52 | tslt.uct <= 52) &
                                           (is.na(diag.status) | diag.status == 0))) / length(which(is.na(diag.status) | diag.status == 0))
  dat$epi$test.syph.12mo.nonhivdiag[at] <- length(which((tslt.syph <= 52) &
                                             (is.na(diag.status) | diag.status == 0))) / length(which(is.na(diag.status) | diag.status == 0))

  # Among those diagnosed
  dat$epi$test.gc.12mo.hivdiag[at] <- length(which((tslt.rgc <= 52 | tslt.ugc <= 52) &
                                           diag.status == 1)) / length(which(diag.status == 1))
  dat$epi$test.ct.12mo.hivdiag[at] <- length(which((tslt.rct <= 52 | tslt.uct <= 52) &
                                           diag.status == 1)) / length(which(diag.status == 1))
  dat$epi$test.syph.12mo.hivdiag[at] <- length(which((tslt.syph <= 52) &
                                             diag.status == 1)) / length(which(diag.status == 1))

  # Among those HIV-negative
  dat$epi$test.gc.12mo.hivneg[at] <- length(which((tslt.rgc <= 52 | tslt.ugc <= 52) &
                                                   status == 0)) / length(which(status == 0))
  dat$epi$test.ct.12mo.hivneg[at] <- length(which((tslt.rct <= 52 | tslt.uct <= 52) &
                                                   status == 0)) / length(which(status == 0))
  dat$epi$test.syph.12mo.hivneg[at] <- length(which((tslt.syph <= 52) &
                                                     status == 0)) / length(which(status == 0))

  # Among those HIV-positive
  dat$epi$test.gc.12mo.hivpos[at] <- length(which((tslt.rgc <= 52 | tslt.ugc <= 52) &
                                                   status == 1)) / length(which(status == 1))
  dat$epi$test.ct.12mo.hivpos[at] <- length(which((tslt.rct <= 52 | tslt.uct <= 52) &
                                                    status == 1)) / length(which(status == 1))
  dat$epi$test.syph.12mo.hivpos[at] <- length(which((tslt.syph <= 52) &
                                                      status == 1)) / length(which(status == 1))

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
