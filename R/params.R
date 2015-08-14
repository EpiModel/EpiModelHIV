
#' @title Epidemic Model Parameters for MARDHAM Models
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}} for MARDHAM.
#'
#' @param nwstats Target statistics for the network model. An object of class
#'        \code{nwstats} output from \code{\link{calc_nwstats.mard}}.
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
#'        The former method is compatible with \code{mardham1}.
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
#' @param prep.efficacy The per-contact efficacy of PrEP to prevent infection if
#'        used (parameter not currently used).
#' @param prep.class.prob The frequency of being in a low, medium, or high class
#'        of adherence to PrEP.
#' @param prep.class.effect The functional effectiveness of PrEP conditional on
#'        PrEP class.
#' @param prep.coverage The proportion of the eligible population who are start
#'        PrEP once they become eligible.
#' @param prep.cov.method The method for calculating PrEP coverage, with options
#'        of \code{"curr"} to base the numerator on the number of people currently
#'        on PrEP and \code{"ever"} to base it on the number of people ever on
#'        PrEP.
#' @param prep.cov.rate The rate at which persons initiate PrEP conditional on
#'        their eligibility, with 1 equal to instant start.
#' @param prep.rcomp The relative change in rate of UAI across all partnership
#'        types given current PrEP use, where 1 is no risk compensation.
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{param.mard}, which can be passed to
#' EpiModel function \code{netsim}.
#'
#' @export
param.mard <- function(nwstats,
                       race.method = 1,
                       last.neg.test.B.int = 301,
                       last.neg.test.W.int = 315,
                       mean.test.B.int = 301,
                       mean.test.W.int = 315,
                       testing.pattern = "interval",
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
                       condom.rr = 0.25,

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
                       cond.pers.always.prob = 0.38,
                       cond.pers.BB.prob = 0.26,
                       cond.pers.BW.prob = 0.26,
                       cond.pers.WW.prob = 0.26,
                       cond.inst.always.prob = 0.55,
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

                       prep.start = 1,
                       prep.elig.model = "all",
                       prep.efficacy = 0.92,
                       prep.class.prob = c(0.50, 0.25, 0.25),
                       prep.class.effect = c(0, 0.75, 0.90),
                       prep.coverage = 0,
                       prep.cov.method = "curr",
                       prep.cov.rate = 1,
                       prep.rcomp = 1,
                       ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  if (!(testing.pattern %in% c("memoryless", "interval"))) {
    stop("testing.pattern must be \"memoryless\" or \"interval\" ",
          call. = FALSE)
  }

  ## Pull from nwstats
  p$time.unit <- nwstats$time.unit

  intvars <- grep(names(p), pattern = ".int", fixed = TRUE)
  p[intvars] <- lapply(p[intvars], FUN = function(x) {round(x / p$time.unit)})

  ratevars <- grep(names(p), pattern = ".rate", fixed = TRUE)
  p[ratevars] <- lapply(p[ratevars], FUN = function(x) {x * p$time.unit})

  p$role.B.prob <- nwstats$role.B.prob
  p$role.W.prob <- nwstats$role.W.prob

  p$inst.trans.matrix <- matrix(1, nrow = 1)
  p$role.trans.matrix <- matrix(c(1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 1),
                                nrow = 3)

  p$method <- nwstats$method
  p$modes <- 1

  p$asmr.B <- nwstats$asmr.B
  p$asmr.W <- nwstats$asmr.W

  if (race.method == 1) {
    p$last.neg.test.B.int = (last.neg.test.B.int + p$last.neg.test.W.int)/2
    p$last.neg.test.W.int = (last.neg.test.B.int + p$last.neg.test.W.int)/2
    p$mean.test.W.int = (mean.test.W.int + mean.test.B.int)/2
    p$mean.test.B.int = (mean.test.W.int + mean.test.B.int)/2
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

  p$nwstats <- NULL

  class(p) <- c("param.mard", "param.net")
  return(p)
}


#' @title Epidemic Model Initial Conditions for MARDHAM Models
#'
#' @description Sets the initial conditions for a stochastic epidemic models
#'              simulated with \code{\link{netsim}}.
#'
#' @param nwstats Target statistics for the network model. An object of class
#'        \code{nwstats} output from \code{\link{calc_nwstats.mard}}.
#' @param prev.B Initial disease prevalence among black MSM.
#' @param prev.W Initial disease prevalence among white MSM.
#' @param ... Additional arguments passed to function.
#'
#' @return
#' A list object of class \code{init.mard}, which can be passed to EpiModel
#' function \code{\link{netsim}}.
#'
#' @export
init.mard <- function(nwstats, prev.B = 0.15, prev.W = 0.15, ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  p$num.B <- nwstats$num.B
  p$num.W <- nwstats$num.W

  p$ages <- nwstats$ages

  p$init.prev.age.slope.B <- prev.B / 12
  p$init.prev.age.slope.W <- prev.W / 12

  p$nwstats <- NULL

  class(p) <- c("init.mard", "init.net")
  return(p)
}


#' @title Epidemic Model Control Settings for MARDHAM Models
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
#' @param resim.int Interval unit for resimulation of network, relative to the
#'        base time unit in the model.
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
#' @param edgescorr.FUN Module function for the edges coefficient adjustment
#'        to preserve mean degree under varying population sizes.
#' @param resimnets.FUN Module function for network resimulation at each time
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
#' @param getprev.FUN Module function to calculate prevalence summary statistics.
#' @param verbose.FUN Module function to print model progress to the console or
#'        external text files.
#' @param delete.nodes If \code{TRUE}, dead nodes will be removed from the network
#'        object and only active nodes will be retained.
#' @param save.nwstats If \code{TRUE}, the network statistics will be saved.
#' @param save.network If \code{TRUE}, the \code{network} objects will be saved
#'        out at the end of simulation (necessary for restarting a simulation).
#' @param save.other Character vector containing other list elements of \code{dat}
#'        to save.
#' @param prevfull If \code{TRUE}, save extended summary statistics for prevalence
#'        and incidence (defined in the prevalence module).
#' @param verbose If \code{TRUE}, print out simulation progress to the console
#'        if in interactive mode or text files if in batch mode.
#' @param verbose.int Integer specifying the interval between time steps at which
#'        progress is printed.
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{control.mard}, which can be passed to the
#' EpiModel function \code{netsim}.
#'
#' @export
control.mard <- function(simno = 1,
                         nsims = 1,
                         ncores = 1,
                         nsteps = 100,
                         start = 1,
                         resim.int = 1,
                         initialize.FUN = initialize.mard,
                         aging.FUN = aging.mard,
                         deaths.FUN = deaths.mard,
                         births.FUN = births.mard,
                         test.FUN = test.mard,
                         tx.FUN = tx.mard,
                         prep.FUN = prep.mard,
                         progress.FUN = progress.mard,
                         vl.FUN = update_vl.mard,
                         aiclass.FUN = update_aiclass.mard,
                         roleclass.FUN = update_roleclass.mard,
                         edgescorr.FUN = edges_correct.mard,
                         resimnets.FUN = simnet.mard,
                         disclose.FUN = disclose.mard,
                         acts.FUN = acts.mard,
                         condoms.FUN = condoms.mard,
                         riskhist.FUN = riskhist.mard,
                         position.FUN = position.mard,
                         trans.FUN = trans.mard,
                         getprev.FUN = prevalence.mard,
                         verbose.FUN = verbose.mard,
                         delete.nodes = TRUE,
                         save.nwstats = TRUE,
                         save.network = FALSE,
                         save.other = "attr",
                         prevfull = FALSE,
                         verbose = TRUE,
                         verbose.int = 1,
                         ...) {

  formal.args <- formals(sys.function())
  dot.args <- list(...)
  p <- get_args(formal.args, dot.args)

  p$skip.check <- TRUE
  p$save.transmat <- FALSE

  p$bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)

  class(p) <- c("control.mard", "control.net")
  return(p)
}
