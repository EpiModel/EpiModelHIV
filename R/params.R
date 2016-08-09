
# MSM -----------------------------------------------------------------

#' @title Epidemic Model Parameters
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}} for EpiModelHIV
#'
#' @param nwstats Target statistics for the network model. An object of class
#'        \code{nwstats} output from \code{\link{calc_nwstats_msm}}.
#' @param race.method Number of races in the model, with options of 1 or 2. If
#'        1, then race-specific parameters will be averaged.
#' @param last.neg.test.B.int Time range in days for last negative test for
#'        black men.
#' @param mean.test.B.int Mean intertest interval in days for black MSM who test.
#' @param last.neg.test.W.int Time range in days for last negative test for
#'        white men.
#' @param mean.test.W.int Mean intertest interval in days for white MSM who test.
#' @param testing.pattern Method for HIV testing, with options \code{"memoryless"}
#'        for constant hazard without regard to time since previous test, or
#'        \code{"interval"} deterministic fixed intervals.
#' @param test.window.int Length of the HIV test window period in days.
#' @param tt.traj.B.prob Proportion of black MSM who enter one of four
#'        testing/treatment trajectories: never test or treat, test and never
#'        initiate treatment, test and treated with partial viral suppression,
#'        and test and treated with full suppression.
#' @param tt.traj.W.prob Proportion of white MSM who enter into the four
#'        testing/treatment trajectories, as defined above.
#' @param tx.init.B.prob Probability per time step that a black MSM who has
#'        tested positive will initiate treatment.
#' @param tx.init.W.prob Probability per time step that a white MSM who has
#'        tested positive will initiate treatment.
#' @param tx.halt.B.prob Probability per time step that a black MSM who is
#'        currently on treatment will halt treatment.
#' @param tx.halt.W.prob Probability per time step that a white MSM who is
#'        currently on treatment will halt treatment.
#' @param tx.reinit.B.prob Probability per time step that a black MSM who is
#'        not currently on treatment but who has been in the past will
#'        re-initiate treatment.
#' @param tx.reinit.W.prob Probability per time step that a white MSM who is
#'        not currently on treatment but who has been in the past will
#'        re-initiate treatment.
#' @param max.time.off.tx.full.int Number of days off treatment for a full
#'        suppressor before onset of AIDS, including time before diagnosis.
#' @param max.time.on.tx.part.int Number of days on treatment for a
#'        partial suppressor beofre onset of AIDS.
#' @param max.time.off.tx.part.int Nnumber of days off treatment for a
#'        partial suppressor before onset of AIDS, including time before
#'        diagnosis.
#' @param vl.acute.rise.int Number of days to peak viremia during acute
#'        infection.
#' @param vl.acute.peak Peak viral load (in log10 units) at the height of acute
#'        infection.
#' @param vl.acute.fall.int Number of days from peak viremia to set-point
#'        viral load during the acute infection period.
#' @param vl.set.point Set point viral load (in log10 units).
#' @param vl.aids.onset.int Number of days to AIDS for a treatment-naive
#'        patient.
#' @param vl.aids.int Duration of AIDS stage infection in days.
#' @param vl.fatal Viral load in AIDS at which death occurs.
#' @param vl.full.supp Log10 viral load at full suppression on ART.
#' @param vl.part.supp Log10 viral load at partial suppression on ART.
#' @param full.supp.down.slope For full suppressors, number of log10 units that
#'        viral load falls per time step from treatment initiation or re-initiation
#'        until the level in \code{vl.full.supp}.
#' @param full.supp.up.slope For full suppressors, number of log10 units that
#'        viral load rises per time step from treatment halting until expected
#'        value.
#' @param part.supp.down.slope For partial suppressors, number of log10 units
#'        that viral load falls per time step from treatment initiation or
#'        re-initiation until the level in \code{vl.part.supp}.
#' @param part.supp.up.slope For partial suppressors, number of log10 units that
#'        viral load rises per time step from treatment halting until expected value.
#' @param b.B.rate Rate at which black MSM enter the population.
#' @param b.W.rate Rate at which white MSM enter the population.
#' @param birth.age Age (in years) of new arrivals.
#' @param b.method Method for calculating the number of expected births at each
#'        time step, with \code{"fixed"} based on the number of persons at the
#'        initial time step and \code{"varying"} based on the current time step.
#' @param URAI.prob Probability of transmission for a man having unprotected
#'        receptive anal intercourse with an infected man at set point viral
#'        load.
#' @param UIAI.prob Probability of transmission for an uncircumcised man having
#'        unprotected insertive anal intercourse with an infected man at set
#'        point viral load.
#' @param acute.rr Relative risk of infection (compared to that predicted by
#'        elevated viral load) when positive partner is in the acute stage.
#' @param circ.rr Relative risk of infection from insertive anal sex when the
#'        negative insertive partner is circumcised.
#' @param condom.rr Relative risk of infection from anal sex when a condom is
#'        used.
#' @param disc.outset.main.B.prob Probability that an HIV-infected black MSM will
#'        disclose his status at the start of a main partnership.
#' @param disc.outset.main.W.prob Probability that an HIV-infected white MSM will
#'        disclose his status at the start of a main partnership.
#' @param disc.at.diag.main.B.prob Probability that a black MSM already in a main
#'        partnership will disclose at the time of diagnosis.
#' @param disc.at.diag.main.W.prob Probability that a white MSM already in a main
#'        partnership will disclose at the time of diagnosis.
#' @param disc.post.diag.main.B.prob Probability that an HIV-infected black MSM
#'        in a main partnership will disclose his status, assuming he didn't
#'        at the start of the partnership or at diagnosis.
#' @param disc.post.diag.main.W.prob Probability that an HIV-infected white MSM
#'        in a main partnership will disclose his status, assuming he didn't
#'        at the start of the partnership or at diagnosis.
#' @param disc.outset.pers.B.prob Probability that an HIV-infected black MSM will
#'        disclose his status at the start of a casual partnership.
#' @param disc.outset.pers.W.prob Probability that an HIV-infected white MSM will
#'        disclose his status at the start of a casual partnership.
#' @param disc.at.diag.pers.B.prob Probability that a black MSM already in a
#'        casual partnership will disclose at the time of diagnosis.
#' @param disc.at.diag.pers.W.prob Probability that a white MSM already in a
#'        casual partnership will disclose at the time of diagnosis.
#' @param disc.post.diag.pers.B.prob Probability that an HIV-infected black MSM
#'        in a casual partnership will disclose his status, assuming he
#'        didn't at the start of the partnership or at diagnosis.
#' @param disc.post.diag.pers.W.prob Probability that an HIV-infected white MSM
#'        in a casual partnership will disclose his status, assuming he
#'        didn't at the start of the partnership or at diagnosis.
#' @param disc.inst.B.prob Probability that an HIV-infected black MSM will
#'        disclose his status to a one-off partner.
#' @param disc.inst.W.prob Probability that an HIV-infected white MSM will
#'        disclose his status to a one-off partner.
#' @param circ.B.prob Probablity that a black new arrival in the population
#'        will be circumcised.
#' @param circ.W.prob Probablity that a white new arrival in the population
#'        will be circumcised.
#' @param ccr5.B.prob Vector of length two of frequencies of the Delta 32
#'        mutation (homozygous and heterozygous, respectively) in the CCR5 gene
#'        among black MSM.
#' @param ccr5.W.prob Vector of length two of frequencies of the Delta 32
#'        mutation (homozygous and heterozygous, respectively) in the CCR5 gene
#'        among white MSM.
#' @param ccr5.heteroz.rr Relative risk of infection for men who are heterozygous
#'        in the CCR5 mutation.
#' @param num.inst.ai.classes Number of quantiles into which men should be
#'        divided in determining their levels of one-off anal intercourse.
#' @param base.ai.main.BB.rate Expected coital frequency in black-black main
#'        partnerships (acts per day).
#' @param base.ai.main.BW.rate Expected coital frequency in black-white main
#'        partnerships (acts per day).
#' @param base.ai.main.WW.rate Expected coital frequency in white-white main
#'        partnerships (acts per day).
#' @param base.ai.pers.BB.rate Expected coital frequency in black-black casual
#'        partnerships (acts per day).
#' @param base.ai.pers.BW.rate Expected coital frequency in black-white casual
#'        partnerships (acts per day).
#' @param base.ai.pers.WW.rate Expected coital frequency in white-white casual
#'        partnerships (acts per day).
#' @param ai.scale General relative scaler for all act rates for model
#'        calibration.
#' @param cond.main.BB.prob Probability of condom use in a black-black main
#'        partnership.
#' @param cond.main.BW.prob Probability of condom use in a black-white main
#'        partnership.
#' @param cond.main.WW.prob Probability of condom use in a white-white main
#'        partnership.
#' @param cond.pers.always.prob Fraction of men in casual partnerships who always
#'        use condoms in those partnerships.
#' @param cond.pers.BB.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a black-black casual partnerships.
#' @param cond.pers.BW.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a black-white casual partnerships.
#' @param cond.pers.WW.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a white-white casual partnerships.
#' @param cond.inst.always.prob Fraction of men in instant partnerships who always
#'        use condoms in those partnerships.
#' @param cond.inst.BB.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a black-black one-off partnerships.
#' @param cond.inst.BW.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a black-white one-off partnerships.
#' @param cond.inst.WW.prob Of men who are not consistent condom users, per-act
#'        probability of condom use in a white-white one-off partnerships.
#' @param cond.always.prob.corr Correlation coefficient for probability of always
#'        using condoms in both casual and one-off
#' @param cond.rr.BB Condom probability scaler for black-black partnerships for
#'        model calibration purposes.
#' @param cond.rr.BW Condom probability scaler for black-white partnerships for
#'        model calibration purposes.
#' @param cond.rr.WW Condom probability scaler for white-white partnerships for
#'        model calibration purposes.
#' @param cond.diag.main.beta Beta multiplier for the log odds of using a
#'        condom in a main partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.main.beta Beta multiplier for the log odds of using a
#'        condom in a main partnership if the HIV-infected man has disclosed.
#' @param cond.diag.pers.beta Beta multiplier for the log odds of using a
#'        condom in a casual partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.pers.beta Beta multiplier for the log odds of using a
#'        condom in a casual partnership if the HIV-infected man has disclosed
#'        his status.
#' @param cond.diag.inst.beta Beta multiplier for the log odds of using a
#'        condom in a one-off partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.inst.beta Beta multiplier for the log odds of using a
#'        condom in a one-off partnership if the HIV-infected man has disclosed
#'        his status.
#' @param vv.iev.BB.prob Probability that in a black-black partnership of
#'        two versatile men, they will engage in intra-event versatility
#'        ("flipping") given that they're having AI.
#' @param vv.iev.BW.prob Probability that in a black-white partnership of
#'        two versatile men, they will engage in intra-event versatility
#'        ("flipping") given that they're having AI.
#' @param vv.iev.WW.prob Probability that in a white-white partnership of
#'        two versatile men, they will engage in intra-event versatility
#'        ("flipping") given that they're having AI.
#' @param prep.start Time step at which the PrEP intervention should start.
#' @param prep.elig.model Modeling approach for determining who is eligible for
#'        PrEP. Current options are limited to: \code{"all"} for all persons who
#'        have never been on PrEP and are disease-susceptible.
#' @param prep.class.prob The probability of adherence class in non-adherent,
#'        low adherence, medium adherence, or high adherence groups (from Liu).
#' @param prep.class.hr The hazard ratio for infection per act associated with each
#'        level of adherence (from Grant).
#' @param prep.coverage The proportion of the eligible population who are start
#'        PrEP once they become eligible.
#' @param prep.cov.method The method for calculating PrEP coverage, with options
#'        of \code{"curr"} to base the numerator on the number of people currently
#'        on PrEP and \code{"ever"} to base it on the number of people ever on
#'        PrEP.
#' @param prep.cov.rate The rate at which persons initiate PrEP conditional on
#'        their eligibility, with 1 equal to instant start.
#' @param prep.tst.int Testing interval for those who are actively on PrEP. This
#'        overrides the mean testing interval parameters.
#' @param prep.risk.int Time window for assessment of risk eligibility for PrEP
#'        in days.
#' @param prep.risk.reassess If \code{TRUE}, reassess eligibility for PrEP at
#'        each testing visit.
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{param_msm}, which can be passed to
#' EpiModel function \code{netsim}.
#'
#' @keywords msm
#'
#' @export
param_msm <- function(nwstats,
                      race.method = 1,
                      last.neg.test.B.int = 301,
                      last.neg.test.W.int = 315,
                      mean.test.B.int = 301,
                      mean.test.W.int = 315,
                      testing.pattern = "memoryless",
                      test.window.int = 21,

                      tt.traj.B.prob = c(0.077, 0.000, 0.356, 0.567),
                      tt.traj.W.prob = c(0.052, 0.000, 0.331, 0.617),

                      tx.init.B.prob = 0.092,
                      tx.init.W.prob = 0.127,
                      tx.halt.B.prob = 0.0102,
                      tx.halt.W.prob = 0.0071,
                      tx.reinit.B.prob = 0.00066,
                      tx.reinit.W.prob = 0.00291,

                      max.time.off.tx.full.int = 520 * 7,
                      max.time.on.tx.part.int = 52 * 15 * 7,
                      max.time.off.tx.part.int = 520 * 7,
                      vl.acute.rise.int = 45,
                      vl.acute.peak = 6.886,
                      vl.acute.fall.int = 45,
                      vl.set.point = 4.5,
                      vl.aids.onset.int = 520 * 7,
                      vl.aids.int = 52 * 2 * 7,
                      vl.fatal = 7,
                      vl.full.supp = 1.5,
                      vl.part.supp = 3.5,
                      full.supp.down.slope = 0.25,
                      full.supp.up.slope = 0.25,
                      part.supp.down.slope = 0.25,
                      part.supp.up.slope = 0.25,

                      b.B.rate = 1e-3 / 7,
                      b.W.rate = 1e-3 / 7,
                      birth.age = 18,
                      b.method = "fixed",

                      URAI.prob = 0.0082 * 1.09,
                      UIAI.prob = 0.0031 * 1.09,
                      acute.rr = 6,
                      circ.rr = 0.4,
                      condom.rr = 0.295,

                      disc.outset.main.B.prob = 0.685,
                      disc.outset.main.W.prob = 0.889,
                      disc.at.diag.main.B.prob = 1,
                      disc.at.diag.main.W.prob = 1,
                      disc.post.diag.main.B.prob = 0,
                      disc.post.diag.main.W.prob = 0,
                      disc.outset.pers.B.prob = 0.527,
                      disc.outset.pers.W.prob = 0.828,
                      disc.at.diag.pers.B.prob = 1,
                      disc.at.diag.pers.W.prob = 1,
                      disc.post.diag.pers.B.prob = 0,
                      disc.post.diag.pers.W.prob = 0,
                      disc.inst.B.prob = 0.445,
                      disc.inst.W.prob = 0.691,

                      circ.B.prob = 0.874,
                      circ.W.prob = 0.918,

                      ccr5.B.prob = c(0, 0.034),
                      ccr5.W.prob = c(0.021, 0.176),
                      ccr5.heteroz.rr = 0.3,

                      num.inst.ai.classes = 1,
                      base.ai.main.BB.rate = 0.17,
                      base.ai.main.BW.rate = 0.26,
                      base.ai.main.WW.rate = 0.23,
                      base.ai.pers.BB.rate = 0.11,
                      base.ai.pers.BW.rate = 0.16,
                      base.ai.pers.WW.rate = 0.14,
                      ai.scale = 1,

                      cond.main.BB.prob = 0.38,
                      cond.main.BW.prob = 0.10,
                      cond.main.WW.prob = 0.15,
                      cond.pers.always.prob = 0.216,
                      cond.pers.BB.prob = 0.26,
                      cond.pers.BW.prob = 0.26,
                      cond.pers.WW.prob = 0.26,
                      cond.inst.always.prob = 0.326,
                      cond.inst.BB.prob = 0.27,
                      cond.inst.BW.prob = 0.27,
                      cond.inst.WW.prob = 0.27,
                      cond.always.prob.corr = 0.5,
                      cond.rr.BB = 1,
                      cond.rr.BW = 1,
                      cond.rr.WW = 1,
                      cond.diag.main.beta = -0.67,
                      cond.discl.main.beta = -0.85,
                      cond.diag.pers.beta = -0.67,
                      cond.discl.pers.beta = -0.85,
                      cond.diag.inst.beta = -0.67,
                      cond.discl.inst.beta = -0.85,

                      vv.iev.BB.prob = 0.42,
                      vv.iev.BW.prob = 0.56,
                      vv.iev.WW.prob = 0.49,

                      prep.start = Inf,
                      prep.elig.model = "base",
                      prep.class.prob = c(0.211, 0.07, 0.1, 0.619),
                      prep.class.hr = c(1, 0.69, 0.19, 0.05),
                      prep.coverage = 0,
                      prep.cov.method = "curr",
                      prep.cov.rate = 1,
                      prep.tst.int = 90,
                      prep.risk.int = 182,
                      prep.risk.reassess = TRUE,
                      ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  if (!(testing.pattern %in% c("memoryless", "interval"))) {
    stop("testing.pattern must be \"memoryless\" or \"interval\" ",
          call. = FALSE)
  }

  if (race.method == 1) {
    p$last.neg.test.B.int = (last.neg.test.B.int + last.neg.test.W.int)/2
    p$last.neg.test.W.int = (last.neg.test.B.int + last.neg.test.W.int)/2
    p$mean.test.B.int = (mean.test.W.int + mean.test.B.int)/2
    p$mean.test.W.int = (mean.test.W.int + mean.test.B.int)/2
    p$tt.traj.B.prob = (tt.traj.B.prob + tt.traj.W.prob)/2
    p$tt.traj.W.prob = (tt.traj.B.prob + tt.traj.W.prob)/2
    p$tx.init.B.prob = (tx.init.B.prob + tx.init.W.prob)/2
    p$tx.init.W.prob = (tx.init.B.prob + tx.init.W.prob)/2
    p$tx.halt.B.prob = (tx.halt.B.prob + tx.halt.W.prob)/2
    p$tx.halt.W.prob = (tx.halt.B.prob + tx.halt.W.prob)/2
    p$tx.reinit.B.prob = (tx.reinit.B.prob + tx.reinit.W.prob)/2
    p$tx.reinit.W.prob = (tx.reinit.B.prob + tx.reinit.W.prob)/2
    p$disc.outset.main.B.prob = (disc.outset.main.B.prob + disc.outset.main.W.prob)/2
    p$disc.outset.main.W.prob = (disc.outset.main.B.prob + disc.outset.main.W.prob)/2
    p$disc.outset.pers.B.prob = (disc.outset.pers.B.prob + disc.outset.pers.W.prob)/2
    p$disc.outset.pers.W.prob = (disc.outset.pers.B.prob + disc.outset.pers.W.prob)/2
    p$disc.inst.B.prob = (disc.inst.B.prob + disc.inst.W.prob)/2
    p$disc.inst.W.prob = (disc.inst.B.prob + disc.inst.W.prob)/2
    p$circ.B.prob = (circ.B.prob + circ.W.prob)/2
    p$circ.W.prob = (circ.B.prob + circ.W.prob)/2
    p$ccr5.B.prob = (ccr5.B.prob + ccr5.W.prob)/2
    p$ccr5.W.prob = (ccr5.B.prob + ccr5.W.prob)/2
    p$base.ai.main.BB.rate = (base.ai.main.BB.rate + base.ai.main.BW.rate +
                                base.ai.main.WW.rate)/3
    p$base.ai.main.BW.rate = (base.ai.main.BB.rate + base.ai.main.BW.rate +
                                base.ai.main.WW.rate)/3
    p$base.ai.main.WW.rate = (base.ai.main.BB.rate + base.ai.main.BW.rate +
                                base.ai.main.WW.rate)/3
    p$base.ai.pers.BB.rate = (base.ai.pers.BB.rate + base.ai.pers.BW.rate +
                                base.ai.pers.WW.rate)/3
    p$base.ai.pers.BW.rate = (base.ai.pers.BB.rate + base.ai.pers.BW.rate +
                                base.ai.pers.WW.rate)/3
    p$base.ai.pers.WW.rate = (base.ai.pers.BB.rate + base.ai.pers.BW.rate +
                                base.ai.pers.WW.rate)/3
    p$cond.main.BB.prob = (cond.main.BB.prob + cond.main.BW.prob + cond.main.WW.prob)/3
    p$cond.main.BW.prob = (cond.main.BB.prob + cond.main.BW.prob + cond.main.WW.prob)/3
    p$cond.main.WW.prob = (cond.main.BB.prob + cond.main.BW.prob + cond.main.WW.prob)/3
    p$cond.pers.BB.prob = (cond.pers.BB.prob + cond.pers.BW.prob + cond.pers.WW.prob)/3
    p$cond.pers.BW.prob = (cond.pers.BB.prob + cond.pers.BW.prob + cond.pers.WW.prob)/3
    p$cond.pers.WW.prob = (cond.pers.BB.prob + cond.pers.BW.prob + cond.pers.WW.prob)/3
    p$cond.inst.BB.prob = (cond.inst.BB.prob + cond.inst.BW.prob + cond.inst.WW.prob)/3
    p$cond.inst.BW.prob = (cond.inst.BB.prob + cond.inst.BW.prob + cond.inst.WW.prob)/3
    p$cond.inst.WW.prob = (cond.inst.BB.prob + cond.inst.BW.prob + cond.inst.WW.prob)/3
    p$vv.iev.BB.prob = (vv.iev.BB.prob + vv.iev.BW.prob + vv.iev.WW.prob)/3
    p$vv.iev.BW.prob = (vv.iev.BB.prob + vv.iev.BW.prob + vv.iev.WW.prob)/3
    p$vv.iev.WW.prob = (vv.iev.BB.prob + vv.iev.BW.prob + vv.iev.WW.prob)/3
  }

  p$time.unit <- nwstats$time.unit

  intvars <- grep(names(p), pattern = ".int", fixed = TRUE)
  p[intvars] <- lapply(p[intvars], FUN = function(x) round(x / p$time.unit))

  ratevars <- grep(names(p), pattern = ".rate", fixed = TRUE)
  p[ratevars] <- lapply(p[ratevars], FUN = function(x) x * p$time.unit)

  p$role.B.prob <- nwstats$role.B.prob
  p$role.W.prob <- nwstats$role.W.prob

  p$inst.trans.matrix <- matrix(1, nrow = 1)
  p$role.trans.matrix <- matrix(c(1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 1),
                                nrow = 3)


  p$riskh.start <- max(1, prep.start - prep.risk.int - 1)

  p$method <- nwstats$method
  p$modes <- 1

  p$asmr.B <- nwstats$asmr.B
  p$asmr.W <- nwstats$asmr.W

  p$nwstats <- NULL

  class(p) <- "param.net"
  return(p)
}


#' @title Epidemic Model Initial Conditions
#'
#' @description Sets the initial conditions for a stochastic epidemic models
#'              simulated with \code{\link{netsim}}.
#'
#' @param nwstats Target statistics for the network model. An object of class
#'        \code{nwstats} output from \code{\link{calc_nwstats_msm}}.
#' @param prev.B Initial disease prevalence among black MSM.
#' @param prev.W Initial disease prevalence among white MSM.
#' @param ... Additional arguments passed to function.
#'
#' @return
#' A list object of class \code{init_msm}, which can be passed to EpiModel
#' function \code{\link{netsim}}.
#'
#' @keywords msm
#'
#' @export
init_msm <- function(nwstats,
                     prev.B = 0.15,
                     prev.W = 0.15,
                     ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  p$num.B <- nwstats$num.B
  p$num.W <- nwstats$num.W

  p$ages <- nwstats$ages

  p$init.prev.age.slope.B <- prev.B / 12
  p$init.prev.age.slope.W <- prev.W / 12

  p$nwstats <- NULL

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
#' @param deaths.FUN Module function for general and disease-realted deaths.
#' @param births.FUN Module function for births or entries into the population.
#' @param test.FUN Module function for diagnostic disease testing.
#' @param tx.FUN Module function for ART initiation and adherence.
#' @param prep.FUN Module function for PrEP initiation and utilization.
#' @param progress.FUN Module function for HIV disease progression.
#' @param vl.FUN Module function for HIV viral load evolution.
#' @param aiclass.FUN Module function for one-off AI risk class transitions.
#' @param roleclass.FUN Module function for transitions in sexual roles.
#' @param resim_nets.FUN Module function for network resimulation at each time
#'        step.
#' @param disclose.FUN Module function for HIV status disclosure.
#' @param acts.FUN Module function to simulate the number of sexual acts within
#'        partnerships.
#' @param condoms.FUN Module function to simulate condom use within acts.
#' @param riskhist.FUN Module function to calculate risk history for uninfected
#'        persons in the population.
#' @param position.FUN Module function to simulate sexual position within acts.
#' @param trans.FUN Module function to stochastically simulate disease transmission
#'        over acts given individual and dyadic attributes.
#' @param prev.FUN Module function to calculate prevalence summary statistics.
#' @param verbose.FUN Module function to print model progress to the console or
#'        external text files.
#' @param save.nwstats Calculate and save network statistics as defined in the
#'        \code{simnet} modules.
#' @param verbose If \code{TRUE}, print out simulation progress to the console
#'        if in interactive mode or text files if in batch mode.
#' @param verbose.int Integer specifying the interval between time steps at which
#'        progress is printed.
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
                        aging.FUN = aging_msm,
                        deaths.FUN = deaths_msm,
                        births.FUN = births_msm,
                        test.FUN = test_msm,
                        tx.FUN = tx_msm,
                        prep.FUN = prep_msm,
                        progress.FUN = progress_msm,
                        vl.FUN = vl_msm,
                        aiclass.FUN = NULL,
                        roleclass.FUN = NULL,
                        resim_nets.FUN = simnet_msm,
                        disclose.FUN = disclose_msm,
                        acts.FUN = acts_msm,
                        condoms.FUN = condoms_msm,
                        riskhist.FUN = riskhist_msm,
                        position.FUN = position_msm,
                        trans.FUN = trans_msm,
                        prev.FUN = prevalence_msm,
                        verbose.FUN = verbose_msm,
                        save.nwstats = FALSE,
                        verbose = TRUE,
                        verbose.int = 1,
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

  p$save.other = c("attr", "temp", "el", "p")

  p$save.network = FALSE

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
#' @param vl.acute.topeak Time in days to peak viremia during acute infection.
#' @param vl.acute.toset Time in days to viral set point following peak viremia.
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
#' @param tx.vlsupp.time Time in days from treatment initiation to viral suppression.
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
  ltGhana <- EpiModelHIV::ltGhana
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
#'        default function of \code{\link{verbose_het}}.
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
#' @param verbose.int Interval for printing progress to console.
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
                        verbose.FUN = verbose_het,
                        module.order = NULL,
                        save.nwstats = FALSE,
                        save.other = c("el", "attr"),
                        verbose = TRUE,
                        verbose.int = 1,
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
