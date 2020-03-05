
# MSM -----------------------------------------------------------------

#' @title Epidemic Model Parameters
#'
#' @description Sets the epidemic parameters for a seattle based stochastic network model
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
#'
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
#'
#' @param rcomp.prob Level of risk compensation from 0 to 1, where 0 is no risk
#'        compensation, 0.5 is a 50% reduction in the probability of condom use
#'        per act, and 1 is a complete cessation of condom use following PrEP
#'        initiation.
#' @param rcomp.adh.groups PrEP adherence groups for whom risk compensation
#'        occurs, as a vector with values 0, 1, 2, 3 corresponding to non-adherent,
#'        low adherence, medium adherence, and high adherence to PrEP.
#' @param rcomp.main.only Logical, if risk compensation is limited to main
#'        partnerships only, versus all partnerships.
#' @param rcomp.discl.only Logical, if risk compensation is limited known-discordant
#'        partnerships only, versus all partnerships.
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
#' @param rgc.asympt.int Average duration in days of asymptomatic rectal gonorrhea.
#' @param ugc.asympt.int Average duration in days of asymptomatic urethral gonorrhea.
#' @param gc.tx.int Average duration in days of treated gonorrhea (both sites).
#' @param gc.ntx.int Average duration in days of untreated, symptomatic gonorrhea (both sites).
#'        If \code{NA}, uses site-specific durations for asymptomatic infections.
#' @param rct.asympt.int Average in days duration of asymptomatic rectal chlamydia.
#' @param uct.asympt.int Average in days duration of asymptomatic urethral chlamydia.
#' @param ct.tx.int Average in days duration of treated chlamydia (both sites).
#' @param ct.ntx.int Average in days duration of untreated, symptomatic chlamydia (both sites).
#'        If \code{NA}, uses site-specific durations for asymptomatic infections.
#' @param gc.prob.cease Probability of ceasing sexual activity during symptomatic
#'        infection with gonorrhea.
#' @param ct.prob.cease Probability of ceasing sexual activity during symptomatic
#'        infection with chlamydia.
#' @param gc.sympt.prob.tx Probability of treatment for symptomatic gonorrhea.
#' @param ct.sympt.prob.tx Probability of treatment for symptomatic chlamydia.
#' @param gc.asympt.prob.tx Probability of treatment for asymptomatic gonorrhea.
#' @param ct.asympt.prob.tx Probability of treatment for asymptomatic chlamydia.
#' @param prep.sti.screen.int Interval in days between STI screening at PrEP visits.
#' @param prep.sti.prob.tx Probability of treatment given positive screening during
#'        PrEP visit.
#' @param prep.continue.stand.tx Logical, if \code{TRUE} will continue standard
#'        STI treatment of symptomatic cases even after PrEP initiation.
#' @param sti.cond.rr Relative risk of STI infection (in either direction) given
#'        a condom used by the insertive partner.
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
param_seattle_msm <- function(nwstats,
                      race.method = 1,
                      
                      #Set for initialized testing histories
                      last.neg.test.B.int = 151,
                      last.neg.test.W.int = 372,
                      mean.test.B.int = 151,
                      mean.test.W.int = 372,
                      testing.pattern = "memoryless",
                      test.window.int = 15,
                      
                      #FOR HOME TESTING INTERVENTION
                      ##FOR HOME TESTING INTERVENTION
                      ################################
                      
                      # Population fractions for opp only, regular, and risk testers
                      tt.B.NO.Reg.Rk = c(0.141,0.666,0.193),
                      tt.W.NO.Reg.Rk = c(0.141,0.666,0.193),
                      ai.class.high.frac = .4,
                      
                      #Opportunistic testing intervals
                      mean.test.opp.NO.W.int = 183,
                      mean.test.opp.NO.B.int = 183,
                      mean.test.opp.ReT.W.int = 183,
                      mean.test.opp.ReT.B.int = 183,
                      mean.test.opp.Risk.W.int = 183,
                      mean.test.opp.Risk.B.int = 183,
                      
                      #Probability of opportunistic test uptake by tester type
                      mean.test.B.opp.prob = 0.764,
                      mean.test.W.opp.prob = 0.764,
                      mean.test.B.opp.reg.prob = 0.0953,
                      mean.test.W.opp.reg.prob = 0.0953,
                      mean.test.B.opp.risk.prob = 0.0953,
                      mean.test.W.opp.risk.prob = 0.0953,
                      
                      #Regular tester time intervals
                      mean.test.B.low.int = 372,
                      mean.test.B.high.int = 151,
                      mean.test.W.low.int = 372,
                      mean.test.W.high.int = 151,
                      
                      #Regular tester time intervals with suppliment test
                      mean.test.B.low.int.supp2 = 184,
                      mean.test.B.high.int.supp2 = 105,
                      mean.test.W.low.int.supp2 = 184,
                      mean.test.W.high.int.supp2 = 105,
                      
                      #Risk based time interval
                      mean.test.risk.CAI.nonmain.B.int = 39,
                      mean.test.risk.CAI.nonmain.W.int = 39,
                      mean.test.risk.AI.known.sd.B.int = 43,
                      mean.test.risk.AI.known.sd.W.int = 43,
                      mean.test.risk.newmain.B.int = 57,
                      mean.test.risk.newmain.W.int = 57,
                      
                     #Risk based testing probabilities
                     mean.test.risk.CAI.nonmain.B.prob = 0.359,
                     mean.test.risk.CAI.nonmain.W.prob = 0.359,
                     mean.test.risk.AI.known.sd.B.prob = 0.538,
                     mean.test.risk.AI.known.sd.W.prob = 0.538,
                     mean.test.risk.newmain.B.prob = 0.375,
                     mean.test.risk.newmain.W.prob = 0.375,
                     
                     #Maximum average number of risk-based tests per year
                     max.risk.test.time=11,

                      #Home testing intervention criteria
                      Opportunity.replace = FALSE,
                      Opportunity.supp = FALSE,
                      Regular.replace = FALSE,
                      Regular.supp = FALSE,
                      Regular.supp2 = FALSE,
                      Risk.replace = FALSE,
                      Risk.supp = FALSE,
                      Never.test.supp = FALSE,
                      
                     hometest.rep.opp.frac = .25,
                     hometest.rep.reg.frac = .25,
                     hometest.rep.risk.frac = .25,
                     

                     supp.opp.home.test.prob = 1/52,
                     supp.reg.home.test.prob = 1/52,
                     supp.risk.home.test.prob = .1,
                     supp.nevertest.home.test.prob = 1/52,
                     hometest.window.int = 90,
                      
                      ################################

                      tt.traj.B.prob = c(0.025,0.000,0.064,0.911),
                      tt.traj.W.prob = c(0.025,0.000,0.064,0.911),

                      tx.init.B.prob = 0.1467,
                      tx.init.W.prob = 0.1467,
                      tx.halt.B.prob = 0.00362,
                      tx.halt.W.prob = 0.00362,
                      tx.reinit.B.prob = 0.00796,
                      tx.reinit.W.prob = 0.00796,

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

                      disc.outset.main.B.prob = 0.9668,
                      disc.outset.main.W.prob = 0.9668,
                      disc.at.diag.main.B.prob = 1,
                      disc.at.diag.main.W.prob = 1,
                      disc.post.diag.main.B.prob = 0,
                      disc.post.diag.main.W.prob = 0,
                      disc.outset.pers.B.prob = 0.8918,
                      disc.outset.pers.W.prob = 0.8918,
                      disc.at.diag.pers.B.prob = 1,
                      disc.at.diag.pers.W.prob = 1,
                      disc.post.diag.pers.B.prob = 0,
                      disc.post.diag.pers.W.prob = 0,
                      disc.inst.B.prob = 0.8089,
                      disc.inst.W.prob = 0.8089,
                     
                      circ.B.prob = 0.6525,
                      circ.W.prob = 0.6525,

                      ccr5.B.prob = c(0.0148,0.1313),
                      ccr5.W.prob = c(0.0148,0.1313),
                      ccr5.heteroz.rr = 0.3,

                      num.inst.ai.classes = 5,
                      base.ai.main.BB.rate = 1.379 / 7,
                      base.ai.main.BW.rate = 1.379 / 7,
                      base.ai.main.WW.rate = 1.379 / 7,
                      base.ai.pers.BB.rate = 0.883 / 7,
                      base.ai.pers.BW.rate = 0.883 / 7,
                      base.ai.pers.WW.rate = 0.883 / 7,
                      ai.scale = 2.00,

                      cond.main.BB.prob = 0.062,
                      cond.main.BW.prob = 0.062,
                      cond.main.WW.prob = 0.062,
                      cond.pers.BB.prob = 0.304,
                      cond.pers.BW.prob = 0.304,
                      cond.pers.WW.prob = 0.304,
                      cond.inst.BB.prob = 0.413,
                      cond.inst.BW.prob = 0.413,
                      cond.inst.WW.prob = 0.413,
                     cond.pers.always.prob = 0,
                     cond.inst.always.prob = 0,
                     cond.always.prob.corr = 0,
                      cond.rr.BB = 1,
                      cond.rr.BW = 1,
                      cond.rr.WW = 1,
                      cond.diag.main.beta = -0.67,
                      cond.discl.main.beta = -0.85,
                      cond.diag.pers.beta = -0.67,
                      cond.discl.pers.beta = -0.85,
                      cond.diag.inst.beta = -0.67,
                      cond.discl.inst.beta = -0.85,

                      vv.iev.BB.prob = 0.479,
                      vv.iev.BW.prob = 0.479,
                      vv.iev.WW.prob = 0.479,

                      prep.start = Inf,
                      prep.elig.model = "base",
                      prep.class.prob = c(0.144, 0.041, 0.053, 0.762),
                      prep.class.hr = c(1, 0.69, 0.19, 0.05),
                      prep.coverage = 0,
                      prep.cov.method = "curr",
                      prep.cov.rate = 1,
                      prep.tst.int = 90,
                      prep.risk.int = 182,
                      prep.risk.reassess = TRUE,
                      riskh.start = 2,

                      rcomp.prob = 0,
                      rcomp.adh.groups = 0:3,
                      rcomp.main.only = FALSE,
                      rcomp.discl.only = FALSE,

                      rgc.tprob = 0.357698,
                      ugc.tprob = 0.248095,
                      rct.tprob = 0.321597,
                      uct.tprob = 0.212965,

                      rgc.sympt.prob = 0.076975,
                      ugc.sympt.prob = 0.824368,
                      rct.sympt.prob = 0.103517,
                      uct.sympt.prob = 0.885045,

                      rgc.asympt.int = 35.11851 * 7,
                      ugc.asympt.int = 35.11851 * 7,
                      gc.tx.int = 2 * 7,
                      gc.ntx.int = NA,

                      rct.asympt.int = 44.24538 * 7,
                      uct.asympt.int = 44.24538 * 7,
                      ct.tx.int = 2 * 7,
                      ct.ntx.int = NA,

                      gc.prob.cease = 0,
                      ct.prob.cease = 0,

                      gc.sympt.prob.tx = 0.90,
                      ct.sympt.prob.tx = 0.85,
                      gc.asympt.prob.tx = 0,
                      ct.asympt.prob.tx = 0,

                      prep.sti.screen.int = 182,
                      prep.sti.prob.tx = 1,
                      prep.continue.stand.tx = TRUE,

                      sti.cond.rr = 0.3,

                      hiv.rgc.rr = 2.780673,
                      hiv.ugc.rr = 1.732363,
                      hiv.rct.rr = 2.780673,
                      hiv.uct.rr = 1.732363,
                      hiv.dual.rr = 0.2,
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

