
#' @title Mardham2 Parameters
#'
#' @description Gather transmission related parameters before simulating from a
#'              stochastic epidemic model.
#'
#' @param nwstats Target statistics for the network model. An object of class
#'        "nwstats", which emerges from \code{\link{calc_nwstats.mard}}.
#' @param last.neg.test.time.range.B Time range for last negative test for
#'        black men.
#' @param mean.test.interv.B The mean intertest interval for Black MSM who test.
#' @param last.neg.test.time.range.W Time range for last negative test for
#'        white men.
#' @param mean.test.interv.W The mean intertest interval for White MSM who test
#' @param testing.pattern Whether HIV testing is memoryless (occurs with fixed
#'        probabilty, without regard to time since previous test) or occurs in
#'        fixed intervals.  Options are "memoryless" or "interval".
#' @param test.window.period The length of the HIV test's window period.
#' @param tt.traj.freq.B The proportion of Black MSM who enter into the four
#'        testing/treatment trajectories.  Values represent, in order, the
#'        propotions who are "NN" (never test or treat), "YN" (test, but then
#'        never initiate treatment), "YP" (test, and then when on treatment
#'        achieve partial suppression), "YF" (test, and then when on treatment
#'        achieve full suppression).
#' @param tt.traj.freq.W The proportion of White MSM who enter into the four
#'        testing/treatment trajectories.  Values represent, in order, the
#'        propotions who are "NN" (never test or treat), "YN" (test, but then
#'        never initiate treatment), "YP" (test, and then when on treatment
#'        achieve partial suppression), "YF" (test, and then when on treatment
#'        achieve full suppression).
#' @param prob.tx.init.B The probability per time step that a Black MSM who has
#'        tested positive will initiate treatment for the very first time.
#' @param prob.tx.init.W The probability per time step that a White MSM who has
#'        tested positive will initiate treatment for the very first time.
#' @param prob.tx.halt.B The probability per time step that a Black MSM who is
#'        currently on treatment will halt treatment.
#' @param prob.tx.halt.W The probability per time step that a White MSM who is
#'        currently on treatment will halt treatment.
#' @param prob.tx.reinit.B The probability per time step that a Black MSM who is
#'        not currently on treatment but who has been in the past will
#'        re-initiate treatment.
#' @param prob.tx.reinit.W The probability per time step that a White MSM who is
#'        not currently on treatment but who has been in the past will
#'        re-initiate treatment.
#' @param max.time.off.tx.full The number of time steps off treatment that a man
#'        who is a full suppressor while on treatment wil last before the onset
#'        of the AIDS stage.  Note that this includes the time before diagnosis.
#' @param max.time.on.tx.part The number of time steps on treatment that a man
#'        who is a partial suppressor wil last before the onset of the AIDS
#'        stage.
#' @param max.time.off.tx.part The number of time steps off treatment that a man
#'        who is a partial suppressor while on treatment wil last before the
#'        onset of the AIDS stage.  Note that this includes the time before
#'        diagnosis.
#' @param vl.acute.rise.dur The number of time steps over which viral load rises
#'        during the initial phase of acute infection.
#' @param vl.acute.peak Peak viral load (in log10 units) at the height of acute
#'        infection.
#' @param vl.acute.fall.dur The number of time steps over which viral load falls
#'        during the latter phase of acute infection.
#' @param vl.set.point Set point viral load for all men (in log10 units).
#' @param vl.aids.onset The number of time steps of infection after which the
#'        AIDS stage will initiate for a man who has never been on treatment.
#' @param vl.aids.dur The duration of the AIDS stage.
#' @param vl.fatal The viral load obtained at the end of the AIDS stage, at
#'        which death from AIDS occurs.  Note that the code can handle
#'        param$vl.fatal being lower than param$vl.acute.peak, as it tests for
#'        both a sufficiently high viral load and presence in the AIDS stage to
#'        deterine AIDS deaths.
#' @param vl.full.supp Log10 viral load when fully suppressed.
#' @param vl.part.supp Log10 viral load when partially suppressed
#' @param full.supp.down.slope For men who are full suppressors, the number of
#'        log10 units that viral load falls per time step from treatment
#'        initiation or re-initiation until the new level is hit.
#' @param full.supp.up.slope For men who are full suppressors, the number of
#'        log10 units that viral load rises per time step from treatment halting
#'        until the new level is hit.
#' @param part.supp.down.slope For men who are partial suppressors, the number
#'        of log10 units that viral load falls per time step from treatment
#'        initiation or re-initiation until the new level is hit.
#' @param part.supp.up.slope For men who are partial suppressors, the number of
#'        log10 units that viral load rises per time step from treatment halting
#'        until the new level is hit.
#' @param b.rate.B The rate at which Black MSM enter the population.
#' @param b.rate.W The rate at which White MSM enter the population.
#' @param birth.age The age (in years) of new arrivals.
#' @param betabase.URAI Transmissibility for a man having unprotected receptive
#'        anal intercourse with a positive man at set point viral load.
#' @param betabase.UIAI Transmissibility for an uncircumcised man having
#'        unprotected insertive anal intercourse with a positive man at set
#'        point viral load.
#' @param betamult.acute The factor by which acute infection is
#'        extra-infectious, above and beyond that predicted by elevated viral
#'        load.
#' @param betamult.circ The factor by which to multiply infectivity from
#'        insertive anal sex when the negative insertive partner is circumcised.
#' @param betamult.condom The factor by which to multiply infectivity for anal
#'        sex when a condom is used.
#' @param disc.main.outset.B The probability that an HIV positive Black MSM will
#'        disclose his status at the start of a main partnership.
#' @param disc.main.outset.W The probability that an HIV positive White MSM will
#'        disclose his status at the start of a main partnership.
#' @param disc.main.at.diag.B The probability that a Black MSM already in a main
#'        partnership will disclose at the time of diagnosis.
#' @param disc.main.at.diag.W The probability that a White MSM already in a main
#'        partnership will disclose at the time of diagnosis.
#' @param disc.main.post.diag.B The probability that an HIV positive Black MSM
#'        in a main partnership will disclose his status, assuming he didn't
#'        at the start of the partnership or at diagnosis.
#' @param disc.main.post.diag.W The probability that an HIV positive White MSM
#'        in a main partnership will disclose his status, assuming he didn't
#'        at the start of the partnership or at diagnosis.
#' @param disc.pers.outset.B The probability that an HIV positive Black MSM will
#'        disclose his status at the start of a casual partnership.
#' @param disc.pers.outset.W The probability that an HIV positive White MSM will
#'        disclose his status at the start of a casual partnership.
#' @param disc.pers.at.diag.B The probability that a Black MSM already in a
#'        casual partnership will disclose at the time of diagnosis.
#' @param disc.pers.at.diag.W The probability that a White MSM already in a
#'        casual partnership will disclose at the time of diagnosis.
#' @param disc.pers.post.diag.B The probability that an HIV positive Black MSM
#'        in a casual partnership will disclose his status, assuming he
#'        didn't at the start of the partnership or at diagnosis.
#' @param disc.pers.post.diag.W The probability that an HIV positive White MSM
#'        in a casual partnership will disclose his status, assuming he
#'        didn't at the start of the partnership or at diagnosis.
#' @param disc.inst.B The probability that an HIV positive Black MSM will
#'        disclose his status to a one-off partner.
#' @param disc.inst.W The probability that an HIV positive White MSM will
#'        disclose his status to a one-off partner.
#' @param circ.prev.B The probablity that a Black new arrival in the population
#'        will be circumcised.
#' @param circ.prev.W The probablity that a White new arrival in the population
#'        will be circumcised.
#' @param ccr5.freq.B A vector of length two of frequencies of the Delta 32
#'        mutation (homozygous and heterozygous, respecitively) in the CCR5 gene
#'        among Black MSM.
#' @param ccr5.freq.W A vector of length two of frequencies of the Delta 32
#'        mutation (homozygous and heterozygous, respecitively) in the CCR5 gene
#'        among White MSM.
#' @param ccr5.heteroz.rr Relative risk of infection for men who are heterzygous
#'        in the CCR5 mutation.
#' @param num.inst.ai.classes The number of quantiles into which men should be
#'        divided in determining their levels of one-off AI.
#' @param base.exp.ai.main.BB Expected coital frequency in Black-Black main
#'        partnerships (acts per week).
#' @param base.exp.ai.main.BW Expected coital frequency in Black-White main
#'        partnerships (acts per week).
#' @param base.exp.ai.main.WW Expected coital frequency in White-White main
#'        partnerships (acts per week).
#' @param base.exp.ai.pers.BB Expected coital frequency in Black-Black casual
#'        partnerships (acts per week).
#' @param base.exp.ai.pers.BW Expected coital frequency in Black-White casual
#'        partnerships (acts per week).
#' @param base.exp.ai.pers.WW Expected coital frequency in White-White casual
#'        partnerships (acts per week).
#' @param incr.exp.ai.full.supp.main Percent increase in expected coital
#'        frequency in main partnerships associated with full suppression of
#'        the disease.
#' @param incr.exp.ai.part.supp.main Percent increase in expected coital
#'        frequency in main partnerships associated with partial suppression of
#'        the disease.
#' @param redux.exp.ai.diag.main Percent reduction in expected coital frequency
#'        in main partnerships associated with diagnosis.
#' @param redux.exp.ai.discl.main Percent reduction in expected coital frequency
#'        in main partnerships associated with disclosure of diagnosis.
#' @param incr.exp.ai.full.supp.pers Percent increase in expected coital
#'        frequency in casual partnerships associated with full suppression of
#'        the disease.
#' @param incr.exp.ai.part.supp.pers Percent increase in expected coital
#'        frequency in casual partnerships associated with partial suppression of
#'        the disease.
#' @param redux.exp.ai.diag.pers Percent reduction in expected coital frequency
#'        in casual partnerships associated with diagnosis.
#' @param redux.exp.ai.discl.pers Percent reduction in expected coital frequency
#'        in casual partnerships associated with disclosure of diagnosis.
#' @param cprob.main.BB The probability of condom use in a Black-Black main
#'        partnership.
#' @param cprob.main.BW The probability of condom use in a Black-White main
#'        partnership.
#' @param cprob.main.WW The probability of condom use in a White-White main
#'        partnership.
#' @param cprob.pers.BB The probability of condom use in a Black-Black casual
#'        partnership.
#' @param cprob.pers.BW The probability of condom use in a Black-White casual
#'        partnership.
#' @param cprob.pers.WW The probability of condom use in a White-White casual
#'        partnership.
#' @param cprob.inst.BB The probability of condom use in a Black-Black one-off
#'        partnership.
#' @param cprob.inst.BW The probability of condom use in a Black-White one-off
#'        partnership.
#' @param cprob.inst.WW The probability of condom use in a White-White one-off
#'        partnership.
#' @param beta.cond.fsupp.main Beta multiplier for the log odds of using a
#'        condom in a main partnership if the HIV positive man is fully
#'        suppressed.
#' @param beta.cond.psupp.main Beta multiplier for the log odds of using a
#'        condom in a main partnership if the HIV positive man is partially
#'        suppressed.
#' @param beta.cond.diag.main Beta multiplier for the log odds of using a
#'        condom in a main partnership if the HIV positive man has been
#'        diagnosed.
#' @param beta.cond.discl.main Beta multiplier for the log odds of using a
#'        condom in a main partnership if the HIV positive man has disclosed
#'        his status.
#' @param beta.cond.fsupp.pers Beta multiplier for the log odds of using a
#'        condom in a casual partnership if the HIV positive man is fully
#'        suppressed.
#' @param beta.cond.psupp.pers Beta multiplier for the log odds of using a
#'        condom in a casual partnership if the HIV positive man is partially
#'        suppressed.
#' @param beta.cond.diag.pers Beta multiplier for the log odds of using a
#'        condom in a casual partnership if the HIV positive man has been
#'        diagnosed.
#' @param beta.cond.discl.pers Beta multiplier for the log odds of using a
#'        condom in a casual partnership if the HIV positive man has disclosed
#'        his status.
#' @param beta.cond.fsupp.inst Beta multiplier for the log odds of using a
#'        condom in a one-off partnership if the HIV positive man is fully
#'        suppressed.
#' @param beta.cond.psupp.inst Beta multiplier for the log odds of using a
#'        condom in a one-off partnership if the HIV positive man is partially
#'        suppressed.
#' @param beta.cond.diag.inst Beta multiplier for the log odds of using a
#'        condom in a one-off partnership if the HIV positive man has been
#'        diagnosed.
#' @param beta.cond.discl.inst Beta multiplier for the log odds of using a
#'        condom in a one-off partnership if the HIV positive man has disclosed
#'        his status.
#' @param vv.prob.iev.BB The probability that in a Black-Black partnership of
#'        two versatile men, they will engage in intra-event versatility
#'        ("flipping") given that they're having AI.
#' @param vv.prob.iev.BW The probability that in a Black-White partnership of
#'        two versatile men, they will engage in intra-event versatility
#'        ("flipping") given that they're having AI.
#' @param vv.prob.iev.WW The probability that in a White-White partnership of
#'        two versatile men, they will engage in intra-event versatility
#'        ("flipping") given that they're having AI.
#' @param ... Additional arguments passed to the function.
#'
#' @return A list object of class \code{param.mard}, which can be passed to
#' EpiModel function \code{netsim}.
#'
#' @export
param.mard <- function(nwstats,
                       last.neg.test.time.range.B = round(301 / 7),
                       mean.test.interv.B = round(301 / 7),
                       last.neg.test.time.range.W = round(315 / 7),
                       mean.test.interv.W = round(315 / 7),
                       testing.pattern = "interval",
                       test.window.period = 3,

                       tt.traj.freq.B = c(0.077, 0.000, 0.356, 0.567),
                       tt.traj.freq.W = c(0.052, 0.000, 0.331, 0.617),

                       prob.tx.init.B = 0.092,
                       prob.tx.init.W = 0.127,
                       prob.tx.halt.B = 0.0102,
                       prob.tx.halt.W = 0.0071,
                       prob.tx.reinit.B = 0.00066,
                       prob.tx.reinit.W = 0.00291,

                       max.time.off.tx.full = 520,
                       max.time.on.tx.part = 52 * 15,
                       max.time.off.tx.part = 520,
                       vl.acute.rise.dur = 3,
                       vl.acute.peak = 6.886,
                       vl.acute.fall.dur = 3,
                       vl.set.point = 4.5,
                       vl.aids.onset = 520,
                       vl.aids.dur = 52 * 2,
                       vl.fatal = 7,
                       vl.full.supp = 1.5,
                       vl.part.supp = 3.5,
                       full.supp.down.slope = 0.25,
                       full.supp.up.slope = 0.25,
                       part.supp.down.slope = 0.25,
                       part.supp.up.slope = 0.25,

                       b.rate.B = 1e-3,
                       b.rate.W = 1e-3,
                       birth.age = 18,

                       betabase.URAI = 0.0082 * 1.09,
                       betabase.UIAI = 0.0031 * 1.09,
                       betamult.acute = 4,
                       betamult.circ = 0.4,
                       betamult.condom = 0.25,

                       disc.main.outset.B = (0.685 + 0.889) / 2,
                       disc.main.outset.W = (0.685 + 0.889) / 2,
                       disc.main.at.diag.B = 1,
                       disc.main.at.diag.W = 1,
                       disc.main.post.diag.B = 0,
                       disc.main.post.diag.W = 0,
                       disc.pers.outset.B = (0.527 + 0.828) / 2,
                       disc.pers.outset.W = (0.527 + 0.828) / 2,
                       disc.pers.at.diag.B = 1,
                       disc.pers.at.diag.W = 1,
                       disc.pers.post.diag.B = 0,
                       disc.pers.post.diag.W = 0,
                       disc.inst.B = (0.445 + 0.691) / 2,
                       disc.inst.W = (0.445 + 0.691) / 2,

                       circ.prev.B = 0.874,
                       circ.prev.W = 0.918,

                       ccr5.freq.B = c(0.0105, 0.105),
                       ccr5.freq.W = c(0.0105, 0.105),
                       ccr5.heteroz.rr = 0.3,

                       num.inst.ai.classes = 1,
                       base.exp.ai.main.BB = 1.19,
                       base.exp.ai.main.BW = 1.79,
                       base.exp.ai.main.WW = 1.56,
                       base.exp.ai.pers.BB = 0.75,
                       base.exp.ai.pers.BW = 1.13,
                       base.exp.ai.pers.WW = 0.98,
                       incr.exp.ai.full.supp.main = 0,
                       incr.exp.ai.part.supp.main = 0,
                       redux.exp.ai.diag.main = 0,
                       redux.exp.ai.discl.main = 0,
                       incr.exp.ai.full.supp.pers = 0,
                       incr.exp.ai.part.supp.pers = 0,
                       redux.exp.ai.diag.pers = 0,
                       redux.exp.ai.discl.pers = 0,

                       cprob.main.BB = 0.38,
                       cprob.main.BW = 0.10,
                       cprob.main.WW = 0.15,
                       cprob.pers.BB = 0.39,
                       cprob.pers.BW = 0.11,
                       cprob.pers.WW = 0.16,
                       cprob.inst.BB = 0.49,
                       cprob.inst.BW = 0.15,
                       cprob.inst.WW = 0.22,

                       beta.cond.fsupp.main = 0.0,
                       beta.cond.psupp.main = 0.0,
                       beta.cond.diag.main = -0.67,
                       beta.cond.discl.main = -0.85,
                       beta.cond.fsupp.pers = 0.0,
                       beta.cond.psupp.pers = 0.0,
                       beta.cond.diag.pers = -0.67,
                       beta.cond.discl.pers = -0.85,
                       beta.cond.fsupp.inst = 0.0,
                       beta.cond.psupp.inst = 0.0,
                       beta.cond.diag.inst = -0.67,
                       beta.cond.discl.inst = -0.85,

                       vv.prob.iev.BB = 0.42,
                       vv.prob.iev.BW = 0.56,
                       vv.prob.iev.WW = 0.49,
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

  if (!(testing.pattern %in% c("memoryless", "interval"))) {
    stop("testing.pattern must be \"memoryless\" or \"interval\" ",
          call. = FALSE)
  }

  ## Pull from nwstats
  p$tUnit <- nwstats$tUnit

  p$role.freq.B <- nwstats$role.freq.B
  p$role.freq.W <- nwstats$role.freq.W

  p$inst.trans.matrix <- matrix(1, nrow = 1)
  p$role.trans.matrix <- matrix(c(1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 1),
                                nrow = 3)

  p$modes <- 1

  p$asmr.B <- nwstats$asmr.B
  p$asmr.W <- nwstats$asmr.W

  p$nwstats <- NULL

  class(p) <- "param.mard"
  return(p)
}


#' @title Mardham2 Initial Conditions
#'
#' @description Set the initial conditions for a stochastic epidemic
#'              network simulation.
#'
#' @param nwstats An list of target statistics created by
#'        \code{\link{calc_nwstats.mard}}.
#' @param prev.B Disease prevalence among Black MSM.
#' @param prev.W Disease prevalence among White MSM.
#' @param init.prev.age.slope.B Slope of initial prevalence by age for Black MSM.
#' @param init.prev.age.slope.W Slope of initial prevalence by age for White MSM.
#' @param ... Additional arguments passed to the function.
#'
#'@details \code{nwstats} must be supplied by the user. Prevalence defaults
#'          to 27.5 percent for both black and white men.
#'
#' @return A list object of class \code{init.mard}, which can be passed to
#'         EpiModel function \code{netsim}.
#'
#' @export
init.mard <- function(nwstats,
                      prev.B = 0.275,
                      prev.W = 0.275,
                      init.prev.age.slope.B = 0.05 / 12,
                      init.prev.age.slope.W = 0.05 / 12,
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

  p$num.B <- nwstats$num.B
  p$num.W <- nwstats$num.W

  p$ages <- nwstats$ages

  p$nwstats <- NULL

  class(p) <- "init.mard"
  return(p)
}


#' @title Mardham2 Control Settings
#'
#' @description Set control options for a stochastic epidemic network model.
#'
#'
#' @details Arguments that end in \code{.FUN} default to the corresponding
#'          \code{.mard} function of the same name.
#'
#' @param simno Unique identifier for the simulation run.
#' @param nsims Number of simulations to run.
#' @param ncores Number of cores to use.
#' @param nsteps Number of time steps in each simulation.
#' @param start Time step at which to start simulation.
#' @param initialize.FUN Function to use for initialization.
#' @param aging.FUN Aging function.
#' @param deaths.FUN Deaths function.
#' @param births.FUN Births function.
#' @param test.FUN Disease testing function.
#' @param tx.FUN Treatment function.
#' @param progress.FUN Disease progress function.
#' @param vl.FUN Viral load function.
#' @param aiclass.FUN AI class function.
#' @param roleclass.FUN Role class function.
#' @param edgescorr.FUN Edges coefficient correction function.
#' @param resimnets.FUN Network resimulation function.
#' @param disclose.FUN Disclosure function.
#' @param acts.FUN Acts function.
#' @param condoms.FUN Condom use function.
#' @param position.FUN Position function.
#' @param trans.FUN Transmission function.
#' @param getprev.FUN Prevalence summary function.
#' @param verbose.FUN Verbose simulation progress function.
#' @param delete.nodes Logical. Should inactive nodes be deleted?
#' @param save.dal Logical. Should the discordant act list be saved at each time
#'        step?
#' @param save.nwstats Logical. Should the network statistics be saved at each
#'        time step?
#' @param save.network Logical. Should the entire networks be saved at each time
#'        step?
#' @param save.other Other list elements of \code{dat} to save.
#' @param verbose Logical. Should simulation progress be printed from within the
#'        time loop?
#' @param verbose.int Integer giving the interval between time steps at which
#'        progress is printed.
#' @param ... Additional arguments passed to the function.
#'
#' @return A list object of class \code{control.mard}, which can be passed to
#'         EpiModel function \code{netsim}.
#'
#' @export
control.mard <- function(simno = 1,
                         nsims = 1,
                         ncores = 1,
                         nsteps = 52 * 50,
                         start = 1,
                         initialize.FUN = initialize.mard,
                         aging.FUN = aging.mard,
                         deaths.FUN = deaths.mard,
                         births.FUN = births.mard,
                         test.FUN = test.mard,
                         tx.FUN = tx.mard,
                         progress.FUN = progress.mard,
                         vl.FUN = update_vl.mard,
                         aiclass.FUN = update_aiclass.mard,
                         roleclass.FUN = update_roleclass.mard,
                         edgescorr.FUN = edges_correct.mard,
                         resimnets.FUN = simnet.mard,
                         disclose.FUN = disclose.mard,
                         acts.FUN = acts.mard,
                         condoms.FUN = condoms.mard,
                         position.FUN = position.mard,
                         trans.FUN = trans.mard,
                         getprev.FUN = prevalence.mard,
                         verbose.FUN = verbose.mard,
                         delete.nodes = TRUE,
                         save.dal = FALSE,
                         save.nwstats = TRUE,
                         save.network = FALSE,
                         save.other = "attr",
                         verbose = TRUE,
                         verbose.int = 1,
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

  p$skip.check <- TRUE
  p$save.transmat <- FALSE

  p$bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  p$user.mods <- grep(".FUN", names.dot.args, value = TRUE)

  class(p) <- "control.mard"
  return(p)
}
