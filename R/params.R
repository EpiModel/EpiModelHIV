
#' @title Epidemic Model Parameters
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with [`netsim`] for EpiModelHIV
#'
#' @param netstats Target statistics and related network initialization data
#'        from the standard ARTnet workflow.
#'
#' @param hiv.test.rate Mean probability of HIV testing per time step for
#'        Black/Hispanic/White MSM (vector of length 3).
#' @param hiv.test.late.prob Proportion of Black/Hispanic/White MSM who test
#'        only during AIDS stage infection (vector of length 3).
#' @param test.window.int Length of the HIV test window period in time steps.
#'        (i.e. the number of time steps between infection and test positivity)
#' @param tt.partial.supp.prob Proportion of Black/Hispanic/White MSM who enter
#'        partial viral suppression category after ART initiation (vector of
#'        length 3).
#' @param tt.full.supp.prob Proportion of Black/Hispanic/White MSM who enter full
#'        viral suppression category after ART initiation (vector of length 3).
#' @param tt.durable.supp.prob Proportion of Black/Hispanic/White MSM who enter
#'        durable viral suppression category after ART initiation (vector of
#'        length 3).
#'
#' @param tx.init.rate Probability per time step that a Black/Hispanic/White MSM
#'        who has tested positive will initiate treatment (vector of length 3).
#' @param tx.halt.partial.rate Probability per time step that Black/Hispanic/White MSM who have
#'        started treatment and assigned to the partial VL suppression category will stop treatment
#'        (vector of length 3).
#' @param tx.halt.full.or Odds ratio comparing the odds of stopping treatment for
#'        Black/Hispanic/White MSM in the full VL suppression category vs. in the partial VL
#'        suppression category (vector of length 3).
#' @param tx.halt.durable.or Odds ratio comparing the odds of stopping treatment for
#'        Black/Hispanic/White MSM in the durable VL suppression category vs. in the partial VL
#'        suppression category (vector of length 3).
#' @param tx.reinit.partial.rate Probability per time step that a Black/Hispanic/White MSM who has
#'        stopped treatment and assigned to the partial VL suppression category will restart
#'        treatment (vector of length 3).
#' @param tx.reinit.full.or Odds ratio comparing the odds of re-starting treatment for
#'        Black/Hispanic/White MSM in the full VL suppression category vs. in the partial VL
#'        suppression category (vector of length 3).
#' @param tx.reinit.durable.or Odds ratio comparing the odds of re-starting treatment for
#'        Black/Hispanic/White MSM in the durable VL suppression category vs. in the partial VL
#'        suppression category (vector of length 3).
#' @param max.time.off.tx.full.int Number of time steps off treatment for a full
#'        suppressor before onset of AIDS, including time before diagnosis.
#' @param max.time.on.tx.partial.int Number of time steps on treatment for a partial suppressor
#'        before onset of AIDS.
#' @param max.time.off.tx.partial.int Number of time steps off treatment for a
#'        partial suppressor before onset of AIDS, including time before
#'        diagnosis.
#' @param vl.acute.rise.int Number of time steps to peak viremia during acute infection.
#' @param vl.acute.peak Peak viral load (in log10 units) at the height of acute infection.
#' @param vl.acute.fall.int Number of time steps from peak viremia to set-point viral load during
#'        the acute infection period.
#' @param vl.set.point Set point viral load (in log10 units).
#' @param vl.aids.onset.int Number of time steps to AIDS for a treatment-naive patient.
#' @param vl.aids.int Duration of AIDS stage infection in time steps.
#' @param vl.aids.peak Maximum viral load during AIDS stage.
#' @param vl.full.supp Log10 viral load at full suppression on ART.
#' @param vl.part.supp Log10 viral load at partial suppression on ART.
#' @param vl.tx.down.rate Number of log10 units that viral load falls per time step from treatment
#'        initiation or re-initiation until the suppression level is reached (pre-AIDS stages).
#' @param vl.tx.aids.down.rate Number of log10 units that viral load falls per time step from
#'        treatment initiation or re-initiation until the suppression level is reached (AIDS stage).
#' @param vl.tx.up.rate Number of log10 units that viral load rises per time step from treatment
#'        halting until expected value.
#' @param aids.on.tx.mort.rate Mortality rate per time step for persons in the AIDS stage who are
#'        currently on ART.
#' @param aids.off.tx.mort.rate Mortality rate per time step for persons in the AIDS stage who are
#'        currently off ART.
#'
#' @param a.rate Rate by time step at which MSM enter the population.
#' @param arrival.age Age (in years) of new arrivals.
#'
#' @param hiv.urai.prob Probability of infection for a man having unprotected receptive anal
#'        intercourse with an infected man at set point viral load.
#' @param hiv.uiai.prob Probability of infection for an uncircumcised man having unprotected
#'        insertive anal intercourse with an infected man at set point viral load.
#' @param hiv.trans.scale Relative scalar on base infection probabilities for model
#'        calibration for Black/Hispanic/White men (vector of length 3).
#' @param hiv.acute.or Odds ratio comparing the odds of infection when the positive partner
#'        is in the acute stage to the odds of infection predicted by elevated viral load alone.
#' @param hiv.circ.or Odds ratio comparing the odds of infection from insertive anal sex when the
#'        negative insertive partner is circumcised vs. not circumcised.
#'
#' @param hiv.cond.eff.or Odds ratio comparing the odds of HIV infection from anal sex when a condom is
#'        used properly vs. not used (biological efficacy).
#' @param hiv.cond.fail.or Absolute reduction applied to the `hiv.cond.eff.or` parameter to account for
#'        condom failures for Black/Hispanic/White MSM (vector of length 3).
#' @param hiv.circ.prob Probability that a Black/Hispanic/White new arrival in the
#'        population will be circumcised (vector of length 3).
#'
#' @param epistats GLMs for epidemiological parameter from the standard ARTnet workflow.
#' @param acts.aids.vl Viral load level after which sexual act rate goes to zero.
#' @param acts.scale Scalar for main/casual act rate for model calibration.
#' @param cond.scale Scalar for condom use probability for model calibration.
#'
#' @param riskh.start Time step at which behavioral risk history assessment occurs.
#' @param prep.start Time step at which the PrEP intervention should start.
#' @param prep.start.prob Probability of starting daily PrEP given current indications.
#' @param prep.adhr.dist Proportion of men who are low, medium, and high adherent to PrEP.
#' @param prep.adhr.or Odds ratio for infection per act associated with each level of adherence.
#' @param prep.risk.reassess.int Interval for reassessment of risk indications of active PrEP users.
#'        Default of `52` for yearly reassessment when time unit is weekly; use `Inf` for no
#'        reassessment.
#' @param prep.require.lnt If `TRUE`, only start on PrEP if current time step is equal to the last
#'        negative test.
#' @param prep.discont.int Median interval until PrEP discontinuation in time steps.
#' @param prep.tst.int Testing interval by time steps for those who are actively on PrEP. This
#'        overrides the mean testing interval parameters.
#' @param prep.risk.int Time window for assessment of risk eligibility for PrEP
#'        in time steps.
#' @param prep.edp.start Time step at which event-driven PrEP is introduced.
#' @param prep.adhr.edp.dist Distribution of event-driven PrEP adherence classe (None, Poor, Good, or Excellent).
#' @param prep.adhr.edp.rr Relative risk for infection per act associated with event-driven
#'        PrEP adherence classes.
#' @param edp.class.static If \code{TRUE}, then assigned EDP adherence classes are static.
#'        Default is \code{FALSE}, suggesting time-varying adherence classes for event-driven PrEP.
#' @param prep.switch.allow If \code{TRUE}, then PrEP users can switch between daily
#'        PrEP and event-driven PrEP regimens. If \code{FALSE}, then switching
#'        between regimens is not allowed unless PrEP was previously discontinued.
#' @param prep.switch.prob Probability of switching from daily PrEP to event-driven PrEP.
#' @param edp.switch.prob Probability of switching from event-driven PrEP to daily PrEP.
#' @param edp.start.scenario Eligibility requirements to starting event-driven PrEP.
#'        Default is `1`, where indications for EDP are same as those for daily PrEP.
#'        `2` = indications for EDP based on any sexual activity.
#'        `3` = indications for EDP based on any sexual activity and NO DO PrEP indications.
#'        `4` = indications for EDP based on prior low adherence to daily PrEP.
#' @param edp.prep.start.prob Probability of starting event-driven PrEP given current indications.
#'
#' @param rgc.prob Probability of rectal gonorrhea infection per act.
#' @param ugc.prob Probability of urethral gonorrhea infection per act.
#' @param rct.prob Probability of rectal chlamydia infection per act.
#' @param uct.prob Probability of urethral chlamydia infection per act.
#' @param rgc.sympt.prob Probability of symptoms given infection with rectal gonorrhea.
#' @param ugc.sympt.prob Probability of symptoms given infection with urethral gonorrhea.
#' @param rct.sympt.prob Probability of symptoms given infection with rectal chlamydia.
#' @param uct.sympt.prob Probability of symptoms given infection with urethral chlamydia.
#'
#' @param rgc.ntx.int Average duration in time steps of untreated rectal gonorrhea.
#' @param ugc.ntx.int Average duration in time steps of untreated urethral gonorrhea.
#' @param gc.tx.int Average duration in time steps of treated gonorrhea (both sites).
#' @param rct.ntx.int Average duration in time steps of untreated rectal chlamydia.
#' @param uct.ntx.int Average duration in time steps of untreated urethral chlamydia.
#' @param ct.tx.int Average duration in time steps of treated chlamydia (both sites).
#'
#' @param gc.sympt.tx.prob Probability of treatment for symptomatic gonorrhea
#'        for Black/Hispanic/White men (vector of length 3).
#' @param ct.sympt.tx.prob Probability of treatment for symptomatic chlamydia
#'        for Black/Hispanic/White men (vector of length 3).
#' @param gc.asympt.tx.prob Probability of treatment for asymptomatic gonorrhea
#'        for Black/Hispanic/White men (vector of length 3).
#' @param ct.asympt.tx.prob Probability of treatment for asymptomatic chlamydia
#'        for Black/Hispanic/White men (vector of length 3).
#'
#' @param prep.sti.screen.int Interval in time steps between STI screening at PrEP visits.
#' @param prep.sti.tx.prob Probability of treatment given positive screening during PrEP visit.
#' @param sti.cond.eff.or Odds ratio comparing the odds of STI infection from anal sex
#'        when a condom is used vs. not used properly (biological efficacy).
#' @param sti.cond.fail.or Absolute reduction applied to the `sti.cond.eff.or` parameter to account for
#'        condom failures for Black/Hispanic/White MSM (vector of length 3).
#' @param hiv.rgc.or Odds ratio comparing the odds of HIV infection given current rectal gonorrhea
#'        to the odds of HIV infection in the absence of current rectal gonorrhea.
#' @param hiv.ugc.or Odds ratio comparing the odds of HIV infection given current urethral gonorrhea
#'        to the odds of HIV infection in the absence of current urethral gonorrhea.
#' @param hiv.rct.or Odds ratio comparing the odds of HIV infection given current rectal chlamydia
#'        to the odds of HIV infection in the absence of current rectal chlamydia.
#' @param hiv.uct.or Odds ratio comparing the odds of HIV infection given current urethral chlamydia
#'        to the odds of HIV infection in the absence of current urethral chlamydia.
#' @param hiv.dual.rr Additive proportional risk, from 0 to 1, for HIV infection
#'        given dual infection with both gonorrhea and chlamydia.
#'
#' @param part.ident.start Time step to initiate the partner identification process. Default of
#'        `Inf` indicates no partner identification.
#' @param part.index.window.int Number of prior time steps for which a newly diagnosed index patient
#'        is eligible for partner notification. Default of 0 indicates index eligible only if newly
#'        diagnosed at current time step.
#' @param part.index.degree The minimum number of partnerships an index patient may have for partner
#'        elicitation. Default of `1` indicating all partnerships under consideration.
#' @param part.index.prob Probability that an index patient will initiate partner services.
#' @param part.ident.main.window.int The number of time steps in the past that a main
#'        partnership qualifies for partner identification.
#' @param part.ident.casl.window.int The number of time steps in the past that a
#'        casual partnership qualifies for partner identification.
#' @param part.ident.ooff.window.int The number of time steps in the past that a
#'        one-off partnership qualifies for partner identification.
#' @param part.ident.main.prob Probability that a main partner is identified
#'        through partner identification.
#' @param part.ident.casl.prob Probability that a casual partner is identified
#'        through partner identification.
#' @param part.ident.ooff.prob Probability that a one-time partner is identified
#'        through partner identification.
#' @param part.hiv.test.rate Mean probability per time step of HIV testing for a
#'        Black/Hispanic/White MSM who has been identified as the partner of an
#'        index HIV+ MSM.
#' @param part.prep.start.prob Probability of an individual identified through partner
#'        identification starting PrEP given current indications at the current time step.
#' @param part.tx.init.rate Probability that a Black/Hispanic/White MSM who has
#'        been identified as the partner of an incident HIV+ MSM will initiate
#'        treatment during the current time step.
#' @param part.tx.reinit.rate Probability per time step that a Black/Hispanic/White MSM who has been
#'        identified through partner identification, stopped treatment will restart treatment
#'        (vector of length 3).
#' @param ... Additional arguments passed to the function.
#'
#' @section Modifying Time-Specific Parameters:
#' The default time steps of these parameters are weekly time steps. Users must modify parameters
#' ending in `.rate`, `.int`, and `.start` when using non-weekly time steps. For example, when using
#' daily time steps, `.int` parameters should be multiplied by 7 and `.rate` parameters should be
#' divided by 7.
#'
#' @return
#' A list object of class `param.net`, which can be passed to EpiModel function [`netsim`].
#'
#' @export
#'
param_msm <- function(netstats,

                      # Clinical
                      hiv.test.rate = c(0.004766828, 0.004172928, 0.006405915),
                      hiv.test.late.prob = c(0, 0, 0),
                      test.window.int = 3,
                      tt.partial.supp.prob = c(0, 0, 0),
                      tt.full.supp.prob = c(1, 1, 1),
                      tt.durable.supp.prob = c(0, 0, 0),
                      tx.init.rate = c(0.298201, 0.367702, 0.3579208),
                      tx.halt.partial.rate = c(0.005033337, 0.004675912, 0.00313423),
                      tx.halt.full.or = c(0.45, 0.45, 0.45),
                      tx.halt.durable.or = c(0.45, 0.45, 0.45),
                      tx.reinit.partial.rate = c(0.00255, 0.00255, 0.00255),
                      tx.reinit.full.or = c(1, 1, 1),
                      tx.reinit.durable.or = c(1, 1, 1),

                      # HIV natural history
                      max.time.off.tx.full.int = 780,
                      max.time.on.tx.partial.int = 520,
                      max.time.off.tx.partial.int = 520,
                      vl.acute.rise.int = 3,
                      vl.acute.peak = 6.886,
                      vl.acute.fall.int = 3,
                      vl.set.point = 4.5,
                      vl.aids.onset.int = 520,
                      vl.aids.int = 104,
                      vl.aids.peak = 7,
                      vl.full.supp = 1.5,
                      vl.part.supp = 3.5,
                      vl.tx.down.rate = 0.25,
                      vl.tx.aids.down.rate = 0.25,
                      vl.tx.up.rate = 0.25,
                      aids.on.tx.mort.rate = 0,
                      aids.off.tx.mort.rate = 0.0006338922,

                      # Demographic
                      a.rate = 0.00049,
                      arrival.age = 15,

                      # HIV infection prob
                      hiv.urai.prob = 0.008938,
                      hiv.uiai.prob = 0.003379,
                      hiv.trans.scale = c(2.8598, 0.4550235, 0.3042151),
                      hiv.acute.or = 6,
                      hiv.circ.or = 0.4,
                      hiv.cond.eff.or = 0.95,
                      hiv.cond.fail.or = c(0.25, 0.25, 0.25),
                      hiv.circ.prob = c(0.874, 0.874, 0.918),

                      # Behavioral
                      epistats,
                      acts.aids.vl = 5.75,
                      acts.scale = 1,
                      cond.scale = 1,

                      # STI epi
                      rgc.prob = 0.2325359,
                      ugc.prob = 0.1951024,
                      rct.prob = 0.2108605,
                      uct.prob = 0.1761155,
                      rgc.sympt.prob = 0.1,
                      ugc.sympt.prob = 0.93,
                      rct.sympt.prob = 0.1,
                      uct.sympt.prob = 0.95,
                      rgc.ntx.int = 26,
                      ugc.ntx.int = 26,
                      gc.tx.int = 1,
                      rct.ntx.int = 32,
                      uct.ntx.int = 32,
                      ct.tx.int = 1,
                      gc.sympt.tx.prob = c(0.9, 0.9, 0.9),
                      ct.sympt.tx.prob = c(0.9, 0.9, 0.9),
                      gc.asympt.tx.prob = c(0.1, 0.1, 0.1),
                      ct.asympt.tx.prob = c(0.1, 0.1, 0.1),
                      sti.cond.eff.or = 0.95,
                      sti.cond.fail.or = c(0.39, 0.30, 0.21),
                      hiv.rgc.or = 2.78,
                      hiv.ugc.or = 1.73,
                      hiv.rct.or = 2.78,
                      hiv.uct.or = 1.73,
                      hiv.dual.rr = 0.2,

                      # PrEP
                      riskh.start = Inf,
                      prep.start = Inf,
                      prep.start.prob = c(0.005791293, 0.004454244, 0.006968663),
                      prep.adhr.dist = c(0.089, 0.127, 0.784),
                      prep.adhr.or = c(0.69, 0.19, 0.01),
                      prep.discont.int = c(33.42, 57.48, 57.39),
                      prep.tst.int = 12.8571,
                      prep.risk.int = 26,
                      prep.sti.screen.int = 26,
                      prep.sti.tx.prob = 1,
                      prep.risk.reassess.int = 52,
                      prep.require.lnt = FALSE,
                      prep.switch.allow = TRUE,
                      prep.switch.prob = 0.25,

                      # EDP-specific parameters
                      prep.edp.start = 400,
                      prep.adhr.edp.dist = c(0.11, 0.06, 0.09, 0.74),
                      prep.adhr.edp.rr = c(1, 0.24, 0.14, 0.03),
                      edp.class.static = FALSE,
                      edp.switch.prob = 0.25,
                      edp.start.scenario = 1,
                      edp.prep.start.prob = c(0.005, 0.005, 0.005),

                      # Partner notification
                      part.ident.start = Inf,
                      part.index.window.int = 0,
                      part.index.degree = 1,
                      part.index.prob = 0.666,
                      part.ident.main.window.int = 52,
                      part.ident.casl.window.int = 52,
                      part.ident.ooff.window.int = 52,
                      part.ident.main.prob = 0.157,
                      part.ident.casl.prob = 0.08518719,
                      part.ident.ooff.prob = 0.04448852,
                      part.hiv.test.rate = c(0.394, 0.394, 0.394),
                      part.prep.start.prob = c(0, 0, 0),
                      part.tx.init.rate = c(0.387, 0.387, 0.387),
                      part.tx.reinit.rate = c(0, 0, 0),
                      ...) {

  if ("prep.discont.rate" %in% names(list(...))) {
    warning(paste0("The value you specified for prep.discont.rate will be ignored; ",
                   "use prep.discont.int instead."))
  }

  if ("aids.mort.rate" %in% names(list(...))) {
    warning(paste0("aids.mort.rate is deprecated; ",
                   "use aids.on.tx.mort.rate & aids.off.tx.mort.rate instead."))
  }

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  class(p) <- "param.net"
  return(p)
}


#' @title Epidemic Model Initial Conditions
#'
#' @description Sets the initial conditions for a stochastic epidemic models
#'              simulated with [`netsim`].
#'
#' @param prev.ugc Initial prevalence of urethral gonorrhea.
#' @param prev.rgc Initial prevalence of rectal gonorrhea.
#' @param prev.uct Initial prevalence of urethral chlamydia.
#' @param prev.rct Initial prevalence of rectal chlamydia.
#' @param ... Additional arguments passed to function.
#'
#' @return
#' A list object of class `init.net`, which can be passed to EpiModel function [`netsim`].
#'
#' @export
#'
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
#' @description Sets the controls for stochastic network models simulated with [`netsim`].
#'
#' @param simno Unique ID for the simulation run, used for file naming purposes
#'        if used in conjunction with the `EpiModelHPC` package.
#' @param nsims Number of simulations.
#' @param ncores Number of cores per run, if parallel processing is used.
#' @param nsteps Number of time steps per simulation.
#' @param start Starting time step for simulation, with default to `1` to run new
#'        simulation. This may also be set to 1 greater than the final time
#'        step of a previous simulation to resume the simulation with different
#'        parameters.
#' @param cumulative.edgelist If `TRUE`, calculates a cumulative edgelist
#'        within the network simulation module. This is used when tergmLite is
#'        used and the entire networkDynamic object is not used.
#' @param truncate.el.cuml Number of time steps of the cumulative edgelist to
#'        retain. See help file for `EpiModel::update_cumulative_edgelist`
#'        for options.
#' @param initialize.FUN Module function to use for initialization of the epidemic model.
#' @param aging.FUN Module function for aging.
#' @param departure.FUN Module function for general and disease-related departures.
#' @param arrival.FUN Module function for entries into the sexually active population.
#' @param partident.FUN Module function for partner identification process.
#' @param hivtest.FUN Module function for HIV diagnostic disease testing.
#' @param hivtx.FUN Module function for ART initiation and adherence.
#' @param prep.FUN Module function for PrEP initiation and utilization.
#' @param hivprogress.FUN Module function for HIV disease progression.
#' @param hivvl.FUN Module function for HIV viral load evolution.
#' @param resim_nets.FUN Module function for network resimulation at each time step.
#' @param acts.FUN Module function to simulate the number of sexual acts within partnerships.
#' @param condoms.FUN Module function to simulate condom use within acts.
#' @param position.FUN Module function to simulate sexual position within acts.
#' @param hivtrans.FUN Module function to stochastically simulate HIV
#'        transmission over acts given individual and dyadic attributes.
#' @param stitrans.FUN Module function to simulate GC/CT transmission.
#' @param stirecov.FUN Module function to simulate recovery from GC/CT,
#'        heterogeneous by disease, site, symptoms, and treatment status.
#' @param stitx.FUN Module function to simulate treatment of GC/CT.
#' @param prev.FUN Module function to calculate prevalence summary statistics.
#' @param verbose.FUN Module function to print model progress to the console or
#'        external text files.
#' @param cleanup.FUN Module function to tidy the `dat` object at the end
#'        of each step.
#' @param save.nwstats Calculate and save network statistics as defined in the
#'        \code{simnet} modules.
#' @param nwstats.formula.1 Right-hand side formula for network statistics to
#'        monitor for network 1, with default `"formation"` equal to the
#'        model formula.
#' @param nwstats.formula.2 Right-hand side formula for network statistics to
#'        monitor for network 2, with default `"formation"` equal to the
#'        model formula.
#' @param nwstats.formula.3 Right-hand side formula for network statistics to
#'        monitor for network 3, with default `"formation"` equal to the
#'        model formula.
#' @param tergmLite Logical indicating usage of `tergmLite` (non-modifiable
#'        currently for `EpiModelHIV`).
#' @param tergmLite.track.duration If `TRUE`, track duration information
#'        for models in `tergmLite` simulations.
#' @param verbose If `TRUE`, print out simulation progress to the console
#'        if in interactive mode or text files if in batch mode.
#' @param set.control.ergm Control arguments passed to `ergm`'s `simulate_formula.network`.
#' @param set.control.tergm Control arguments passed to `tergm`'s `simulate_formula.network`.
#' @param set.control.stergm Deprecated; use `set.control.tergm` instead.
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class `control.net`, which can be passed to the EpiModel function [`netsim`].
#'
#' @export
#'
control_msm <- function(simno = 1,
                        nsims = 1,
                        ncores = 1,
                        nsteps = 100,
                        start = 1,
                        cumulative.edgelist = TRUE,
                        truncate.el.cuml = 0,
                        initialize.FUN = initialize_msm,
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
                        cleanup.FUN = cleanup_msm,
                        save.nwstats = FALSE,
                        nwstats.formula.1 = "formation",
                        nwstats.formula.2 = "formation",
                        nwstats.formula.3 = "formation",
                        tergmLite = TRUE,
                        tergmLite.track.duration = FALSE,
                        set.control.ergm = control.simulate.formula(
                          MCMC.burnin = 2e5
                        ),
                        set.control.tergm = control.simulate.formula.tergm(
                          MCMC.burnin.min = 5000
                        ),
                        set.control.stergm = NULL,
                        verbose = TRUE,
                        ...) {
  if (!missing(set.control.stergm)) {
    warning("set.control.stergm is deprecated and will be removed in a future
             version; use set.control.tergm instead.")
  }

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
  p[["f.names"]] <- c(p[["bi.mods"]], p[["user.mods"]])
  p$save.other <- c("attr", "temp", "el")

  p$save.network <- FALSE
  if (is.null(p$verbose.int)) {
    p$verbose.int <- 1
  }

  p$save.diss.stats <- FALSE

  p <- set.control.class("control.net", p)
  return(p)
}
