
# MSM -----------------------------------------------------------------

#' @title Epidemic Model Parameters
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}} for EpiModelHIV
#'
#' @param netstats Target statistics and related network initialization data from
#'        the standard ARTnet workflow.
#'
#' @param hiv.test.rate Mean probability of HIV testing per week for
#'        black/hispanic/white MSM (vector of length 3).
#' @param hiv.test.late.prob Proportion of black/hispanic/white MSM who test only
#'        during AIDS stage infection (vector of length 3).
#' @param test.window.int Length of the HIV test window period in weeks.
#' @param tt.part.supp Proportion of black/hispanic/white MSM who enter partial viral
#'        suppression category after ART initiation (vector of length 3).
#' @param tt.full.supp Proportion of black/hispanic/white MSM who enter full viral
#'        suppression category after ART initiation (vector of length 3).
#' @param tt.dur.supp Proportion of black/hispanic/white MSM who enter durable viral
#'        suppression category after ART initiation (vector of length 3).
#'
#' @param tx.init.prob Probability per time step that a black/hispanic/white MSM who has
#'        tested positive will initiate treatment (vector of length 3).
#' @param tx.halt.part.prob Probability per time step that black/hispanic/white
#'        MSM who have started treatment and assigned to the partial VL suppression
#'        category will stop treatment (vector of length 3).
#' @param tx.halt.full.rr Relative reduction in \code{tx.halt.part.prob} for
#'        black/hispanic/white MSM in the full VL suppression category (vector of length 3).
#' @param tx.halt.dur.rr Relative reduction in \code{tx.halt.part.prob} for
#'        black/hispanic/white MSM in the durable VL suppression category (vector of length 3).
#' @param tx.reinit.part.prob Probability per time step that a black/hispanic/white
#'        MSM who has stopped treatment and assigned to the partial VL suppression
#'        category will restart treatment (vector of length 3).
#' @param tx.reinit.full.rr Relative reduction in \code{tx.reinit.part.prob} for
#'        black/hispanic/white MSM in the full VL suppression category (vector of length 3).
#' @param tx.reinit.dur.rr Relative reduction in \code{tx.reinit.part.prob} for
#'        black/hispanic/white MSM in the durable VL suppression category (vector of length 3).
#' @param max.time.off.tx.full.int Number of weeks off treatment for a full
#'        suppressor before onset of AIDS, including time before diagnosis.
#' @param max.time.on.tx.part.int Number of weeks on treatment for a
#'        partial suppressor beofre onset of AIDS.
#' @param max.time.off.tx.part.int Nnumber of weeks off treatment for a
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
#' @param vl.tx.aids.down.slope Number of log10 units that viral load falls per time
#'        step from treatment initiation or re-initiation until the suppression
#'        level is reached (AIDS stage).
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
#'        calibration for black/hispanic/white men (vector of length 3).
#' @param acute.rr Relative risk of infection (compared to that predicted by
#'        elevated viral load) when positive partner is in the acute stage.
#' @param circ.rr Relative risk of infection from insertive anal sex when the
#'        negative insertive partner is circumcised.
#'
#' @param cond.eff Relative risk of HIV infection from anal sex when a condom is
#'        used properly (biological efficacy).
#' @param cond.fail Condom failure rates for HIV for black/hispanic/white MSM, as a reduction
#'        in the cond.eff parameter (vector of length 3).
#' @param circ.prob Probablity that a black/hispanic/white new arrival in the population
#'        will be circumcised (vector of length 3).
#'
#' @param epistats GLMs for epidemiological parameter from the standard ARTnet workflow.
#' @param acts.aids.vl Viral load level after which sexual act rate goes to zero.
#' @param acts.scale.main Scalar for main act rate for model calibration.
#' @param acts.scale.casl Scalar for casual act rate for model calibration.
#' @param cond.scale Scalar for condom use probability for model calibration.
#'
#' @param riskh.start Time step at which behavioral risk history assessment occurs.
#' @param prep.start Time step at which the PrEP intervention should start.
#' @param prep.start.prob Probability of starting PrEP given current indications.
#' @param prep.adhr.dist Proportion of men who are low, medium, and high
#'        adherent to PrEP.
#' @param prep.adhr.hr The hazard ratio for infection per act associated with each
#'        level of adherence (from Grant).
#' @param prep.risk.reassess.method Interval for reassessment of risk indications
#'        of active PrEP users, either \code{"none"} for no reassessment,
#'        \code{"inst"} for weekly, or \code{"year"} for year.
#' @param prep.require.lnt If \code{TRUE}, only start on PrEP if current time step is
#'        equal to the last negative test.
#'
#' @param prep.discont.rate Rate of random discontinuation from PrEP.
#'
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
#'        for black/hispanic/white men (vector of length 3).
#' @param ct.sympt.prob.tx Probability of treatment for symptomatic chlamydia
#'        for black/hispanic/white men (vector of length 3).
#' @param gc.asympt.prob.tx Probability of treatment for asymptomatic gonorrhea
#'        for black/hispanic/white men (vector of length 3).
#' @param ct.asympt.prob.tx Probability of treatment for asymptomatic chlamydia
#'        for black/hispanic/white men (vector of length 3).
#'
#' @param prep.sti.screen.int Interval in weeks between STI screening at PrEP visits.
#' @param prep.sti.prob.tx Probability of treatment given positive screening during
#'        PrEP visit.
#' @param sti.cond.eff Relative risk of STI infection from anal sex when a condom is
#'        used properly (biological efficacy).
#' @param sti.cond.fail Condom failure rates for STI for black/hispanic/white MSM, as
#'        a reduction in the cond.eff parameter (vector of length 3).
#' @param hiv.rgc.rr Relative risk of HIV infection given current rectal gonorrhea.
#' @param hiv.ugc.rr Relative risk of HIV infection given current urethral gonorrhea.
#' @param hiv.rct.rr Relative risk of HIV infection given current rectal chlamydia.
#' @param hiv.uct.rr Relative risk of HIV infection given current urethral chlamydia.
#' @param hiv.dual.rr Additive proportional risk, from 0 to 1, for HIV infection
#'        given dual infection with both gonorrhea and chlamydia.
#'
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
                      hiv.test.late.prob = c(0.25, 0.25, 0.25),
                      test.window.int = 21/7,
                      tt.part.supp = c(0.20, 0.20, 0.20),
                      tt.full.supp = c(0.40, 0.40, 0.40),
                      tt.dur.supp = c(0.40, 0.40, 0.40),
                      tx.init.prob = c(0.092, 0.092, 0.127),
                      tx.halt.part.prob = c(0.0102, 0.0102, 0.0071),
                      tx.halt.full.rr = c(0.9, 0.9, 0.9),
                      tx.halt.dur.rr = c(0.5, 0.5, 0.5),
                      tx.reinit.part.prob = c(0.00066, 0.00066, 0.00291),
                      tx.reinit.full.rr = c(1.0, 1.0, 1.0),
                      tx.reinit.dur.rr = c(1.0, 1.0, 1.0),

                      # HIV natural history
                      max.time.off.tx.full.int = 52 * 15,
                      max.time.on.tx.part.int = 52 * 10,
                      max.time.off.tx.part.int = 52 * 10,
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
                      aids.mr = 1/104,

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
                      acts.scale.main = 1,
                      acts.scale.casl = 1,
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
                      prep.start.prob = 0.2,
                      prep.adhr.dist = c(0.089, 0.127, 0.784),
                      prep.adhr.hr = c(0.69, 0.19, 0.01),
                      prep.discont.rate = 1 - (2^(-1/(224.4237/7))),
                      prep.tst.int = 90/7,
                      prep.risk.int = 182/7,
                      prep.sti.screen.int = 182/7,
                      prep.sti.prob.tx = 1,
                      prep.risk.reassess.method = "year",
                      prep.require.lnt = TRUE,
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
init_msm <- function(prev.ugc = 0.005,
                     prev.rgc = 0.005,
                     prev.uct = 0.013,
                     prev.rct = 0.013,
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
#' @param ncores Number of cores per run, if parallelization is used within the
#'        \code{EpiModelHPC} package.
#' @param nsteps Number of time steps per simulation.
#' @param start Starting time step for simulation, with default to 1 to run new
#'        simulation. This may also be set to 1 greater than the final time
#'        step of a previous simulation to resume the simulation with different
#'        parameters.
#' @param initialize.FUN Module function to use for initialization of the epidemic
#'        model.
#' @param aging.FUN Module function for aging.
#' @param departure.FUN Module function for general and disease-realted depatures.
#' @param arrival.FUN Module function for entries into the sexually active population.
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
#' @param hivtrans.FUN Module function to stochastically simulate HIV transmission
#'        over acts given individual and dyadic attributes.
#' @param stitrans.FUN Module function to simulate GC/CT transmission over current
#'        edgelist.
#' @param stirecov.FUN Module function to simulate recovery from GC/CT, heterogeneous
#'        by disease, site, symptoms, and treatment status.
#' @param stitx.FUN Module function to simulate treatment of GC/CT.
#' @param prev.FUN Module function to calculate prevalence summary statistics.
#' @param verbose.FUN Module function to print model progress to the console or
#'        external text files.
#' @param save.nwstats Calculate and save network statistics as defined in the
#'        \code{simnet} modules.
#' @param save.clin.hist Save individual-level clinical history matrices.
#' @param truncate.plist Truncate the cumulative partnership list to only include
#'        active partnerships.
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
                        save.clin.hist = FALSE,
                        truncate.plist = TRUE,
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
  p$tergmLite <- TRUE

  class(p) <- "control.net"
  return(p)
}



# HET -----------------------------------------------------------------


#' @title Parameters for Stochastic Network Model of HIV-1 Infection in
#'        Sub-Saharan Africa
#'
#' @description Sets the simulation parameters for the stochastic
#'              network model of HIV-1 Infection among Heterosexuals in
#'              Sub-Saharan Africa for the \code{EpiModelHIV} package.
#'
#' @param time.unit Unit of time relative to one day.
#'
#' @param acute.stage.mult Acute stage multiplier for increased infectiousness
#'        above impact of heightened viral load.
#' @param aids.stage.mult AIDS stage multiplier for increased infectiousness in
#'        AIDS above impact of heightened viral load.
#'
#' @param vl.acute.topeak Time in weeks to peak viremia during acute infection.
#' @param vl.acute.toset Time in weeks to viral set point following peak viremia.
#' @param vl.acute.peak Log 10 viral load at acute peak.
#' @param vl.setpoint Log 10 viral load at set point.
#' @param vl.aidsmax Maximum log 10 viral load during AIDS.
#'
#' @param cond.prob Probability of condoms per act with partners.
#' @param cond.eff Efficacy of condoms per act in HIV prevention.
#'
#' @param act.rate.early Daily per-partnership act rate in early disease.
#' @param act.rate.late Daily per-partnership act rate in late disease.
#' @param act.rate.cd4 CD4 count at which the \code{act.rate.late} applies.
#' @param acts.rand If \code{TRUE}, will draw number of total and unprotected
#'        acts from a binomial distribution parameterized by the \code{act.rate}.
#'
#' @param circ.prob.birth Proportion of men circumcised at birth.
#' @param circ.eff Efficacy of circumcision per act in HIV prevention.
#'
#' @param tx.elig.cd4 CD4 count at which a person becomes eligible for treatment.
#' @param tx.init.cd4.mean Mean CD4 count at which person presents for care.
#' @param tx.init.cd4.sd SD of CD4 count at which person presents for care.
#' @param tx.adhere.full Proportion of people who start treatment who are fully
#'        adherent.
#' @param tx.adhere.part Of the not fully adherent proportion, the percent of time
#'        they are on medication.
#' @param tx.vlsupp.time Time in weeks from treatment initiation to viral suppression.
#' @param tx.vlsupp.level Log 10 viral load level at suppression.
#' @param tx.cd4.recrat.feml Rate of CD4 recovery under treatment for males.
#' @param tx.cd4.recrat.male Rate of CD4 recovery under treatment for females.
#' @param tx.cd4.decrat.feml Rate of CD4 decline under periods of non-adherence
#'        for females.
#' @param tx.cd4.decrat.male Rate of CD4 decline under periods of non-adherence
#'        for males.
#' @param tx.coverage Proportion of treatment-eligible persons who have initiated
#'        treatment.
#' @param tx.prev.eff Proportional amount by which treatment reduces infectivity
#'        of infected partner.
#'
#' @param b.rate General entry rate per day for males and females specified.
#' @param b.rate.method Method for assigning birth rates, with options of "totpop"
#'        for births as a function of the total population size, "fpop" for births
#'        as a function of the female population size, and "stgrowth" for a constant
#'        stable growth rate.
#' @param b.propmale Proportion of entries assigned as male. If NULL, then set
#'        adaptively based on the proportion at time 1.
#'
#' @param ds.exit.age Age at which the age-specific ds.rate is set to 1, with NA
#'        value indicating no censoring.
#' @param ds.rate.mult Simple multiplier for background death rates.
#' @param di.cd4.aids CD4 count at which late-stage AIDS occurs and the risk of
#'        mortality is governed by \code{di.cd4.rate}.
#' @param di.cd4.rate Mortality in late-stage AIDS after hitting a nadir CD4 of
#'        \code{di.cd4.aids}.
#' @param ... additional arguments to be passed into model.
#'
#' @details This function sets the parameters for the models.
#'
#' @keywords het
#'
#' @export
#'
param_het <- function(time.unit = 7,

                      acute.stage.mult = 5,
                      aids.stage.mult = 1,

                      vl.acute.topeak = 14,
                      vl.acute.toset = 107,
                      vl.acute.peak = 6.7,
                      vl.setpoint = 4.5,
                      vl.aidsmax = 7,

                      cond.prob = 0.09,
                      cond.eff = 0.78,

                      act.rate.early = 0.362,
                      act.rate.late = 0.197,
                      act.rate.cd4 = 50,
                      acts.rand = TRUE,

                      circ.prob.birth = 0.9,
                      circ.eff = 0.53,

                      tx.elig.cd4 = 350,
                      tx.init.cd4.mean = 120,
                      tx.init.cd4.sd = 40,
                      tx.adhere.full = 0.76,
                      tx.adhere.part = 0.50,
                      tx.vlsupp.time = 365/3,
                      tx.vlsupp.level = 1.5,
                      tx.cd4.recrat.feml = 11.6/30,
                      tx.cd4.recrat.male = 9.75/30,
                      tx.cd4.decrat.feml = 11.6/30,
                      tx.cd4.decrat.male = 9.75/30,
                      tx.coverage = 0.3,
                      tx.prev.eff = 0.96,

                      b.rate = 0.03/365,
                      b.rate.method = "totpop",
                      b.propmale = NULL,

                      ds.exit.age = 55,
                      ds.rate.mult = 1,
                      di.cd4.aids = 50,
                      di.cd4.rate = 2/365,
                      ...) {

  ## Process parameters
  p <- list()
  formal.args <- formals(sys.function())
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg))
  }
  dot.args <- list(...)
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }


  ## trans.rate multiplier
  p$trans.rate <- p$trans.rate * p$trans.rate.mult


  ## Death rate transformations
  # ltGhana <- EpiModelHIV::ltGhana
  ltGhana <- 1
  ds.rates <- ltGhana[ltGhana$year == 2011, ]
  ds.rates$mrate <- ds.rates$mrate / 365
  if (is.numeric(ds.exit.age)) {
    ds.rates$mrate[ds.rates$agStart >= ds.exit.age] <- 1
  }
  ds.rates$reps <- ds.rates$agEnd - ds.rates$agStart + 1
  ds.rates$reps[ds.rates$agStart == 100] <- 1
  male <- rep(ds.rates$male, ds.rates$reps)
  mrate <- rep(ds.rates$mrate, ds.rates$reps)
  mrate <- pmin(1, mrate * ds.rate.mult)
  age <- rep(0:100, 2)
  ds.rates <- data.frame(male = male, age, mrate = mrate)
  ds.rates <- ds.rates[ds.rates$age != 0, ]
  p$ds.rates <- ds.rates

  ## Time unit scaling
  if (time.unit > 1) {

    ## Rates multiplied by time unit
    p$act.rate.early <- act.rate.early * time.unit
    p$act.rate.late <- act.rate.late * time.unit
    p$b.rate <- b.rate * time.unit
    p$ds.rates$mrate <- ifelse(p$ds.rates$mrate < 1,
                               p$ds.rates$mrate * time.unit,
                               p$ds.rates$mrate)

    p$dx.prob.feml <- p$dx.prob.feml * time.unit
    p$dx.prob.male <- p$dx.prob.male * time.unit
    p$tx.cd4.recrat.feml <- tx.cd4.recrat.feml * time.unit
    p$tx.cd4.recrat.male <- tx.cd4.recrat.male * time.unit
    p$tx.cd4.decrat.feml <- tx.cd4.decrat.feml * time.unit
    p$tx.cd4.decrat.male <- tx.cd4.decrat.male * time.unit
    p$di.cd4.rate <- di.cd4.rate * time.unit

    ## Intervals divided by time unit
    p$vl.acute.topeak <- vl.acute.topeak / time.unit
    p$vl.acute.toset <- vl.acute.toset / time.unit

    p$tx.vlsupp.time <- tx.vlsupp.time / time.unit

  }

  p$model <- "a2"

  class(p) <- "param.net"
  return(p)
}


#' @title Initial Conditions for Stochastic Network Model of HIV-1 Infection in
#'        Sub-Saharan Africa
#'
#' @description This function sets the initial conditions for the stochastic
#'              network models in the \code{epimethods} package.
#'
#' @param i.prev.male Prevalence of initially infected males.
#' @param i.prev.feml Prevalence of initially infected females.
#' @param ages.male initial ages of males in the population.
#' @param ages.feml initial ages of females in the population.
#' @param inf.time.dist Probability distribution for setting time of infection
#'        for nodes infected at T1, with options of \code{"geometric"} for randomly
#'        distributed on a geometric distribution with a probability of the
#'        reciprocal of the average length of infection, \code{"uniform"} for a
#'        uniformly distributed time over that same interval, or \code{"allacute"} for
#'        placing all infections in the acute stage at the start.
#' @param max.inf.time Maximum infection time in days for infection at initialization,
#'        used when \code{inf.time.dist} is \code{"geometric"} or \code{"uniform"}.
#' @param ... additional arguments to be passed into model.
#'
#' @details This function sets the initial conditions for the models.
#'
#' @keywords het
#'
#' @export
#'
init_het <- function(i.prev.male = 0.05,
                     i.prev.feml = 0.05,
                     ages.male = seq(18, 55, 7/365),
                     ages.feml = seq(18, 55, 7/365),
                     inf.time.dist = "geometric",
                     max.inf.time = 5 * 365,
                     ...) {

  ## Process parameters
  p <- list()
  formal.args <- formals(sys.function())
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg))
  }
  dot.args <- list(...)
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }


  ## Parameter checks
  if (!(inf.time.dist %in% c("uniform", "geometric", "allacute"))) {
    stop("inf.time.dist must be \"uniform\" or \"geometric\" or \"allacute\" ")
  }

  class(p) <- "init.net"
  return(p)
}


#' @title Control Settings for Stochastic Network Model of HIV-1 Infection in
#'        Sub-Saharan Africa
#'
#' @description This function sets the control settings for the stochastic
#'              network models in the \code{epimethods} package.
#'
#' @param simno Simulation ID number.
#' @param nsteps Number of time steps to simulate the model over in whatever unit
#'        implied by \code{time.unit}.
#' @param start Starting time step for simulation
#' @param nsims Number of simulations.
#' @param ncores Number of parallel cores to use for simulation jobs, if using
#'        the \code{EpiModel.hpc} package.
#' @param par.type Parallelization type, either of \code{"single"} for multi-core
#'        or \code{"mpi"} for multi-node MPI threads.
#' @param initialize.FUN Module to initialize the model at time 1.
#' @param aging.FUN Module to age active nodes.
#' @param cd4.FUN CD4 progression module.
#' @param vl.FUN HIV viral load progression module.
#' @param dx.FUN HIV diagnosis module.
#' @param tx.FUN HIV treatment module.
#' @param deaths.FUN Module to simulate death or exit.
#' @param births.FUN Module to simulate births or entries.
#' @param resim_nets.FUN Module to resimulate the network at each time step.
#' @param trans.FUN Module to simulate disease infection.
#' @param prev.FUN Module to calculate disease prevalence at each time step,
#'        with the default function of \code{\link{prevalence_het}}.
#' @param verbose.FUN Module to print simulation progress to screen, with the
#'        default function of \code{verbose.net}.
#' @param module.order A character vector of module names that lists modules the
#'        order in which they should be evaluated within each time step. If
#'        \code{NULL}, the modules will be evaluated as follows: first any
#'        new modules supplied through \code{...} in the order in which they are
#'        listed, then the built-in modules in their order of the function listing.
#'        The \code{initialize.FUN} will always be run first and the
#'        \code{verbose.FUN} always last.
#' @param save.nwstats Save out network statistics.
#' @param save.other Other list elements of dat to save out.
#' @param verbose If \code{TRUE}, print progress to console.
#' @param skip.check If \code{TRUE}, skips the error check for parameter values,
#'        initial conditions, and control settings before running the models.
#' @param ... Additional arguments passed to the function.
#'
#' @details This function sets the parameters for the models.
#'
#' @keywords het
#'
#' @export
#'
control_het <- function(simno = 1,
                        nsteps = 100,
                        start = 1,
                        nsims = 1,
                        ncores = 1,
                        par.type = "single",
                        initialize.FUN = initialize_het,
                        aging.FUN = aging_het,
                        cd4.FUN = cd4_het,
                        vl.FUN = vl_het,
                        dx.FUN = dx_het,
                        tx.FUN = tx_het,
                        deaths.FUN = deaths_het,
                        births.FUN = births_het,
                        resim_nets.FUN = simnet_het,
                        trans.FUN = trans_het,
                        prev.FUN = prevalence_het,
                        verbose.FUN = verbose.net,
                        module.order = NULL,
                        save.nwstats = FALSE,
                        save.other = c("el", "attr"),
                        verbose = TRUE,
                        skip.check = TRUE,
                        ...) {

  p <- list()
  formal.args <- formals(sys.function())
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg))
  }
  dot.args <- list(...)
  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }

  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
                                  USE.NAMES = FALSE) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names.dot.args, value = TRUE)

  p$save.transmat <- FALSE
  p$save.network <- FALSE

  class(p) <- "control.net"
  return(p)
}
