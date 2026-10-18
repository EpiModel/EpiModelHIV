
# MSM -----------------------------------------------------------------

#' @title Epidemic Model Parameters
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}} for EpiModelHIV
#'
#' @param netstats Target statistics and related network initialization data
#'        from the standard ARTnet workflow.
#'
#' @param hiv.test.rate Mean probability of HIV testing per week for
#'        Black/Hispanic/White MSM (vector of length 3).
#' @param hiv.test.late.prob Proportion of Black/Hispanic/White MSM who test
#'        only during AIDS stage infection (vector of length 3).
#' @param test.window.int Length of the HIV test window period in weeks.
#' @param tt.partial.supp Proportion of Black/Hispanic/White MSM who enter
#'        partial viral suppression category after ART initiation (vector of
#'        length 3).
#' @param tt.full.supp Proportion of Black/Hispanic/White MSM who enter full
#'        viral suppression category after ART initiation (vector of length 3).
#' @param tt.durable.supp Proportion of Black/Hispanic/White MSM who enter
#'        durable viral suppression category after ART initiation (vector of
#'        length 3).
#'
#' @param tx.init.prob Probability per time step that a Black/Hispanic/White MSM
#'        who has tested positive will initiate treatment (vector of length 3).
#' @param tx.halt.partial.prob Probability per time step that
#'        Black/Hispanic/White MSM who have started treatment and assigned to
#'        the partial VL suppression category will stop treatment (vector of
#'        length 3).
#' @param tx.halt.full.rr Relative reduction in \code{tx.halt.partial.prob} for
#'        Black/Hispanic/White MSM in the full VL suppression category (vector
#'        of length 3).
#' @param tx.halt.durable.rr Relative reduction in \code{tx.halt.partial.prob}
#'        for Black/Hispanic/White MSM in the durable VL suppression category
#'        (vector of length 3).
#' @param tx.reinit.partial.prob Probability per time step that a
#'        Black/Hispanic/White MSM who has stopped treatment and assigned to the
#'        partial VL suppression category will restart treatment (vector of
#'        length 3).
#' @param tx.reinit.full.rr Relative reduction in \code{tx.reinit.partial.prob}
#'        for Black/Hispanic/White MSM in the full VL suppression category
#'        (vector of length 3).
#' @param tx.reinit.durable.rr Relative reduction in
#'        \code{tx.reinit.partial.prob} for Black/Hispanic/White MSM in the
#'        durable VL suppression category (vector of length 3).
#' @param max.time.off.tx.full.int Number of weeks off treatment for a full
#'        suppressor before onset of AIDS, including time before diagnosis.
#' @param max.time.on.tx.partial.int Number of weeks on treatment for a
#'        partial suppressor before onset of AIDS.
#' @param max.time.off.tx.partial.int Number of weeks off treatment for a
#'        partial suppressor before onset of AIDS, including time before
#'        diagnosis.
#' @param vl.acute.rise.int Number of weeks to peak viremia during acute
#'        infection.
#' @param vl.acute.peak Peak viral load (in log10 units) at the height of acute
#'        infection.
#' @param vl.acute.fall.int Number of weeks from peak viremia to set-point
#'        viral load during the acute infection period.
#' @param vl.set.point Set point viral load (in log10 units).
#' @param vl.aids.onset.int Number of weeks to AIDS for a treatment-naive
#'        patient.
#' @param vl.aids.int Duration of AIDS stage infection in weeks.
#' @param vl.aids.peak Maximum viral load during AIDS stage.
#' @param vl.full.supp Log10 viral load at full suppression on ART.
#' @param vl.part.supp Log10 viral load at partial suppression on ART.
#' @param vl.tx.down.slope Number of log10 units that viral load falls per time
#'        step from treatment initiation or re-initiation until the suppression
#'        level is reached (pre-AIDS stages).
#' @param vl.tx.aids.down.slope Number of log10 units that viral load falls per
#'        time step from treatment initiation or re-initiation until the
#'        suppression level is reached (AIDS stage).
#' @param vl.tx.up.slope Number of log10 units that viral load rises per time
#'        step from treatment halting until expected value.
#' @param aids.mr Mortality rate of persons in the AIDS stage who are currently
#'        off ART.
#'
#' @param a.rate Rate at which MSM enter the population.
#' @param arrival.age Age (in years) of new arrivals.
#'
#' @param URAI.prob Probability of transmission for a man having unprotected
#'        receptive anal intercourse with an infected man at set point viral
#'        load.
#' @param UIAI.prob Probability of transmission for an uncircumcised man having
#'        unprotected insertive anal intercourse with an infected man at set
#'        point viral load.
#' @param trans.scale Relative scalar on base infection probabilities for model
#'        calibration for Black/Hispanic/White men (vector of length 3).
#' @param acute.rr Relative risk of infection (compared to that predicted by
#'        elevated viral load) when positive partner is in the acute stage.
#' @param circ.rr Relative risk of infection from insertive anal sex when the
#'        negative insertive partner is circumcised.
#'
#' @param cond.eff Relative risk of HIV infection from anal sex when a condom is
#'        used properly (biological efficacy).
#' @param cond.fail Condom failure rates for HIV for Black/Hispanic/White MSM,
#'        as a reduction in the cond.eff parameter (vector of length 3).
#' @param circ.prob Probability that a Black/Hispanic/White new arrival in the
#'        population will be circumcised (vector of length 3).
#'
#' @param epistats GLMs for epidemiological parameter from the standard ARTnet
#'        workflow.
#' @param acts.aids.vl Viral load level after which sexual act rate goes to
#'        zero.
#' @param acts.scale Scalar for main/casual act rate for model calibration.
#' @param cond.scale Scalar for condom use probability for model calibration.
#'
#' @param riskh.start Time step at which behavioral risk history assessment
#'        occurs.
#' @param prep.start Time step at which the PrEP intervention should start.
#' @param prep.start.prob Probability of starting PrEP given current
#'        indications.
#' @param prep.adhr.dist Proportion of men who are low, medium, and high
#'        adherent to PrEP.
#' @param prep.adhr.hr The hazard ratio for infection per act associated with
#'        each level of adherence (from Grant).
#' @param prep.risk.reassess.method Interval for reassessment of risk
#'        indications of active PrEP users, either \code{"none"} for no
#'        reassessment, \code{"inst"} for weekly, or \code{"year"} for year.
#' @param prep.require.lnt If \code{TRUE}, only start on PrEP if current time
#'        step is equal to the last negative test.
#' @param prep.discont.rate Rate of random discontinuation from PrEP.
#' @param prep.tst.int Testing interval for those who are actively on PrEP. This
#'        overrides the mean testing interval parameters.
#' @param prep.risk.int Time window for assessment of risk eligibility for PrEP
#'        in weeks.
#'
#' @param rgc.tprob Probability of rectal gonorrhea infection per act.
#' @param ugc.tprob Probability of urethral gonorrhea infection per act.
#' @param rct.tprob Probability of rectal chlamydia infection per act.
#' @param uct.tprob Probability of urethral chlamydia infection per act.
#' @param rgc.sympt.prob Probability of symptoms given infection with rectal
#'        gonorrhea.
#' @param ugc.sympt.prob Probability of symptoms given infection with urethral
#'        gonorrhea.
#' @param rct.sympt.prob Probability of symptoms given infection with rectal
#'        chlamydia.
#' @param uct.sympt.prob Probability of symptoms given infection with urethral
#'        chlamydia.
#'
#' @param rgc.ntx.int Average duration in weeks of untreated rectal gonorrhea.
#' @param ugc.ntx.int Average duration in weeks of untreated urethral gonorrhea.
#' @param gc.tx.int Average duration in weeks of treated gonorrhea (both sites).
#' @param rct.ntx.int Average in weeks duration of untreated rectal chlamydia.
#' @param uct.ntx.int Average in weeks duration of untreated urethral chlamydia.
#' @param ct.tx.int Average in weeks duration of treated chlamydia (both sites).
#'
#' @param gc.sympt.prob.tx Probability of treatment for symptomatic gonorrhea
#'        for Black/Hispanic/White men (vector of length 3).
#' @param ct.sympt.prob.tx Probability of treatment for symptomatic chlamydia
#'        for Black/Hispanic/White men (vector of length 3).
#' @param gc.asympt.prob.tx Probability of treatment for asymptomatic gonorrhea
#'        for Black/Hispanic/White men (vector of length 3).
#' @param ct.asympt.prob.tx Probability of treatment for asymptomatic chlamydia
#'        for Black/Hispanic/White men (vector of length 3).
#'
#' @param prep.sti.screen.int Interval in weeks between STI screening at PrEP
#'        visits.
#' @param prep.sti.prob.tx Probability of treatment given positive screening
#'        during PrEP visit.
#' @param sti.cond.eff Relative risk of STI infection from anal sex when a
#'        condom is used properly (biological efficacy).
#' @param sti.cond.fail Condom failure rates for STI for
#'        Black/Hispanic/White MSM, as a reduction in \code{cond.eff} (vector of
#'        length 3).
#' @param hiv.rgc.rr Relative risk of HIV infection given current rectal
#'        gonorrhea.
#' @param hiv.ugc.rr Relative risk of HIV infection given current urethral
#'        gonorrhea.
#' @param hiv.rct.rr Relative risk of HIV infection given current rectal
#'        chlamydia.
#' @param hiv.uct.rr Relative risk of HIV infection given current urethral
#'        chlamydia.
#' @param hiv.dual.rr Additive proportional risk, from 0 to 1, for HIV infection
#'        given dual infection with both gonorrhea and chlamydia.
#'
#' @param truncate.plist Numeric value that truncates the cumulative partner
#'        list to retain partnerships within that specified number of previous
#'        timesteps. To truncate to currently active partnerships, set to
#'        \code{0}. For no truncation, set to \code{Inf}.
#' @param part.ident.start Time step to initiate the partner identification
#'        process. Default of \code{Inf} indicates no partner identification.
#' @param part.index.window Number of prior time steps for which a newly
#'        diagnosed index patient is eligible for partner notification. Default
#'        of 0 indicates index eligible only if newly diagnosed at current time
#'        step.
#' @param part.index.degree The minimum number of partnerships an index patient
#'        may have for partner elicitation. Default of \code{1} indicating all
#'        partnerships under consideration.
#' @param part.index.prob Probability that an index patient will initiate
#'        partner services.
#' @param part.ident.main.window The number of timesteps in the past that a main
#'        partnership qualifies for partner identification.
#' @param part.ident.casl.window The number of timesteps in the past that a
#'        casual partnership qualifies for partner identification.
#' @param part.ident.ooff.window The number of timesteps in the past that a
#'        one-off partnership qualifies for partner identification.
#' @param part.ident.main.prob Probability that a main partner is identified
#'        through partner identification.
#' @param part.ident.casl.prob Probability that a casual partner is identified
#'        through partner identification.
#' @param part.ident.ooff.prob Probability that a one-time partner is identified
#'        through partner identification.
#' @param part.hiv.test.rate Mean probability of HIV testing for a
#'        Black/Hispanic/White MSM who has been identified as the partner of an
#'        index HIV+ MSM.
#' @param part.prep.start.prob Probability of an individual identified through
#'        partner identification starting PrEP given current indications at the
#'        current time step.
#' @param part.tx.init.prob Probability that a Black/Hispanic/White MSM who has
#'        been identified as the partner of an incident HIV+ MSM will initiate
#'        treatment during the current time step.
#' @param part.tx.halt.prob Probability per time step that Black/Hispanic/White
#'        MSM who have have been identified through partner identification,
#'        started treatment and assigned to the partial VL suppression category
#'        will stop treatment (vector of length 3).
#' @param part.tx.reinit.prob Probability per time step that a
#'        Black/Hispanic/White MSM who has been identified through partner
#'        identification, stopped treatment will restart treatment (vector of
#'        length 3).
#' @param epi_trackers Factory functions for adding new epi output statistics
#'        to model (advanced).
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{param_msm}, which can be passed to
#' EpiModel function \code{netsim}.
#'
#' @keywords msm
#'
#' @export
#'
param_msm <- function(netstats,

                      # Clinical
                      hiv.test.rate = c(0.01325, 0.0125, 0.0124),
                      hiv.test.late.prob = c(0, 0, 0),
                      test.window.int = 21 / 7,
                      tt.partial.supp = c(0, 0, 0),
                      tt.full.supp = c(1, 1, 1),
                      tt.durable.supp = c(0, 0, 0),
                      tx.init.prob = c(0.092, 0.092, 0.127),
                      tx.halt.partial.prob = c(0.0102, 0.0102, 0.0071),
                      tx.halt.full.rr = c(0.45, 0.45, 0.45),
                      tx.halt.durable.rr = c(0.45, 0.45, 0.45),
                      tx.reinit.partial.prob = c(0.00066, 0.00066, 0.00291),
                      tx.reinit.full.rr = c(1, 1, 1),
                      tx.reinit.durable.rr = c(1, 1, 1),

                      # HIV natural history
                      max.time.off.tx.full.int = 52 * 15,
                      max.time.on.tx.partial.int = 52 * 10,
                      max.time.off.tx.partial.int = 52 * 10,
                      vl.acute.rise.int = 6.4,
                      vl.acute.peak = 6.886,
                      vl.acute.fall.int = 6.4,
                      vl.set.point = 4.5,
                      vl.aids.onset.int = 520,
                      vl.aids.int = 104,
                      vl.aids.peak = 7,
                      vl.full.supp = 1.5,
                      vl.part.supp = 3.5,
                      vl.tx.down.slope = 0.25,
                      vl.tx.aids.down.slope = 0.25,
                      vl.tx.up.slope = 0.25,
                      aids.mr = 1 / 104,

                      # Demographic
                      a.rate = 0.00052,
                      arrival.age = 15,

                      # HIV transmission prob
                      URAI.prob = 0.008938,
                      UIAI.prob = 0.003379,
                      trans.scale = c(1, 1, 1),
                      acute.rr = 6,
                      circ.rr = 0.4,
                      cond.eff = 0.95,
                      cond.fail = c(0.25, 0.25, 0.25),
                      circ.prob = c(0.874, 0.874, 0.918),

                      # Behavioral
                      epistats,
                      acts.aids.vl = 5.75,
                      acts.scale = 1,
                      cond.scale = 1,

                      # STI epi
                      rgc.tprob = 0.35,
                      ugc.tprob = 0.25,
                      rct.tprob = 0.20,
                      uct.tprob = 0.16,
                      rgc.sympt.prob = 0.16,
                      ugc.sympt.prob = 0.80,
                      rct.sympt.prob = 0.14,
                      uct.sympt.prob = 0.58,
                      rgc.ntx.int = 16.8,
                      ugc.ntx.int = 16.8,
                      gc.tx.int = 1.4,
                      rct.ntx.int = 32,
                      uct.ntx.int = 32,
                      ct.tx.int = 1.4,
                      gc.sympt.prob.tx = c(0.95, 0.95, 0.95),
                      ct.sympt.prob.tx = c(0.9, 0.9, 0.9),
                      gc.asympt.prob.tx = c(0.15, 0.15, 0.15),
                      ct.asympt.prob.tx = c(0.15, 0.15, 0.15),
                      sti.cond.eff = 0.9,
                      sti.cond.fail = c(0.20, 0.20, 0.20),
                      hiv.rgc.rr = 2.78,
                      hiv.ugc.rr = 1.73,
                      hiv.rct.rr = 2.78,
                      hiv.uct.rr = 1.73,
                      hiv.dual.rr = 0.2,

                      # PrEP
                      riskh.start = Inf,
                      prep.start = Inf,
                      prep.start.prob = c(0.2, 0.2, 0.2),
                      prep.adhr.dist = c(0.089, 0.127, 0.784),
                      prep.adhr.hr = c(0.69, 0.19, 0.01),
                      prep.discont.rate = rep( 1-(2^(-1/(224.4237/7))) , 3),
                      prep.tst.int = 90 / 7,
                      prep.risk.int = 182 / 7,
                      prep.sti.screen.int = 182 / 7,
                      prep.sti.prob.tx = 1,
                      prep.risk.reassess.method = "year",
                      prep.require.lnt = TRUE,

                      # Partner notification
                      truncate.plist = 0,
                      part.ident.start = Inf,
                      part.index.window = 0,
                      part.index.degree = 1,
                      part.index.prob = 1,
                      part.ident.main.window = 12,
                      part.ident.casl.window = 12,
                      part.ident.ooff.window = 12,
                      part.ident.main.prob = 1,
                      part.ident.casl.prob = 1,
                      part.ident.ooff.prob = 1,
                      part.hiv.test.rate = c(1, 1, 1),
                      part.prep.start.prob = c(0.5, 0.5, 0.5),
                      part.tx.init.prob = c(0.6, 0.6, 0.8),
                      part.tx.halt.prob = c(0.00102, 0.00102, 0.00071),
                      part.tx.reinit.prob = c(0.5, 0.5, 0.5),

                      # Epi Trackers
                      epi_trackers = list(),

                      ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  class(p) <- "param.net"
  return(p)
}


#' @title Epidemic Model Initial Conditions
#'
#' @description Sets the initial conditions for a stochastic epidemic models
#'              simulated with \code{\link{netsim}}.
#'
#' @param prev.ugc Initial prevalence of urethral gonorrhea.
#' @param prev.rgc Initial prevalence of rectal gonorrhea.
#' @param prev.uct Initial prevalence of urethral chlamydia.
#' @param prev.rct Initial prevalence of rectal chlamydia.
#' @param ... Additional arguments passed to function.
#'
#' @return
#' A list object of class \code{init_msm}, which can be passed to EpiModel
#' function \code{\link{netsim}}.
#'
#' @keywords msm
#'
#' @export
init_msm <- function(prev.ugc = 0,
                     prev.rgc = 0,
                     prev.uct = 0,
                     prev.rct = 0,
                     ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  class(p) <- "init.net"
  return(p)
}


#' @title Epidemic Model Control Settings
#'
#' @description Sets the controls for stochastic network models simulated with
#'              \code{\link{netsim}}.
#'
#' @param simno Unique ID for the simulation run, used for file naming purposes
#'        if used in conjunction with the \code{EpiModelHPC} package.
#' @param nsims Number of simulations.
#' @param ncores Number of cores per run, if parallel processing is used.
#' @param nsteps Number of time steps per simulation.
#' @param start Starting time step for simulation, with default to 1 to run new
#'        simulation. This may also be set to 1 greater than the final time
#'        step of a previous simulation to resume the simulation with different
#'        parameters.
#' @param initialize.FUN Module function to use for initialization of the
#'        epidemic model.
#' @param param_updater.FUN Module function to update the parameters of the
#'        model at specified time steps
#' @param aging.FUN Module function for aging.
#' @param departure.FUN Module function for general and disease-related
#'        departures.
#' @param arrival.FUN Module function for entries into the sexually active
#'        population.
#' @param partident.FUN Module function for partner identification process.
#' @param hivtest.FUN Module function for HIV diagnostic disease testing.
#' @param hivtx.FUN Module function for ART initiation and adherence.
#' @param prep.FUN Module function for PrEP initiation and utilization.
#' @param hivprogress.FUN Module function for HIV disease progression.
#' @param hivvl.FUN Module function for HIV viral load evolution.
#' @param resim_nets.FUN Module function for network resimulation at each time
#'        step.
#' @param acts.FUN Module function to simulate the number of sexual acts within
#'        partnerships.
#' @param condoms.FUN Module function to simulate condom use within acts.
#' @param position.FUN Module function to simulate sexual position within acts.
#' @param hivtrans.FUN Module function to stochastically simulate HIV
#'        transmission over acts given individual and dyadic attributes.
#' @param stitrans.FUN Module function to simulate GC/CT transmission over
#'        current edgelist.
#' @param stirecov.FUN Module function to simulate recovery from GC/CT,
#'        heterogeneous by disease, site, symptoms, and treatment status.
#' @param stitx.FUN Module function to simulate treatment of GC/CT.
#' @param prev.FUN Module function to calculate prevalence summary statistics.
#' @param verbose.FUN Module function to print model progress to the console or
#'        external text files.
#' @param save.nwstats Calculate and save network statistics as defined in the
#'        \code{simnet} modules.
#' @param tergmLite Logical indicating usage of \code{tergmLite} (non-modifiable
#'        currently for \code{EpiModelHIV}).
#' @param verbose If \code{TRUE}, print out simulation progress to the console
#'        if in interactive mode or text files if in batch mode.
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{control_msm}, which can be passed to the
#' EpiModel function \code{netsim}.
#'
#' @keywords msm
#'
#' @export
control_msm <- function(simno = 1,
                        nsims = 1,
                        ncores = 1,
                        nsteps = 100,
                        start = 1,
                        initialize.FUN = initialize_msm,
                        param_updater.FUN = param_updater,
                        aging.FUN = aging_msm,
                        departure.FUN = departure_msm,
                        arrival.FUN = arrival_msm,
                        partident.FUN = partident_msm,
                        hivtest.FUN = hivtest_msm,
                        hivtx.FUN = hivtx_msm,
                        hivprogress.FUN = hivprogress_msm,
                        hivvl.FUN = hivvl_msm,
                        resim_nets.FUN = simnet_msm,
                        acts.FUN = acts_msm,
                        condoms.FUN = condoms_msm,
                        position.FUN = position_msm,
                        prep.FUN = prep_msm,
                        hivtrans.FUN = hivtrans_msm,
                        stitrans.FUN = stitrans_msm,
                        stirecov.FUN = stirecov_msm,
                        stitx.FUN = stitx_msm,
                        prev.FUN = prevalence_msm,
                        verbose.FUN = verbose.net,
                        save.nwstats = FALSE,
                        tergmLite = TRUE,
                        verbose = TRUE,
                        ...) {

  formal.args <- formals(sys.function())
  dot.args <- list(...)
  p <- get_args(formal.args, dot.args)

  p$skip.check <- TRUE
  p$save.transmat <- FALSE

  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
                                  USE.NAMES = FALSE) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)

  p$save.other <- c("attr", "temp", "el", "p")

  p$save.network <- FALSE
  p$verbose.int <- 1

  class(p) <- "control.net"
  return(p)
}
