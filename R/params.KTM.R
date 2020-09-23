
# SHAMP  -----------------------------------------------------------------

#' @title Epidemic Model Parameters for SHAMP
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}} for EpiModelHIV
#'
#' @param last.neg.test.R.S.int Time range in days for last negative test for
#'        race group R where R is a charater (B, BI, H, HI or W) and sex group S
#'        where s is a charter (f, msf, msm, msmf).
#' @param mean.test.B.int Mean intertest interval in days for race group R where R 
#'        is a charater (B, BI, H, HI or W) and sex group S
#'        where s is a charter (f, msf, msm, msmf) who test.
#' @param test.window.int.ab Length of the HIV antibody test window period in weeks.
#' @param test.window.int.rna Length of the RNA HIV test window period in weeks.
#' @param tt.traj.S.prob Proportion of race group R, where R is a charater (B, BI, H, HI or W), 
#'        and sex group S, where s is a charter (f, msf, msm, msmf), who enter one of four
#'        testing/treatment trajectories: never test or treat, test and never
#'        initiate treatment, test and treated with partial viral suppression,
#'        and test and treated with full suppression.
#' @param tx.init.R.S.prob Probability per time step that a person in race group R, where R is a 
#'        charater (B, BI, H, HI or W), and sex group S, where s is a charter (f, msf, msm, msmf),
#'        who has tested positive will initiate treatment for the first 8 weeks following diagnosis.
#' @param tx.halt.R.S.prob Probability per time step that a person in race group R, where R is a 
#'        charater (B, BI, H, HI or W), and sex group S, where s is a charter (f, msf, msm, msmf),
#'        who is currently on treatment will halt treatment.
#' @param tx.reinit.R.S.prob Probability per time step that a person in race group R, where R is a 
#'        charater (B, BI, H, HI or W), and sex group S, where s is a charter (f, msf, msm, msmf),
#'        who is not currently on treatment but who has been in the past will
#'        re-initiate treatment or that someone who did not initiate treatment within 8 weeks of diagnosis will initiate treatment.
#' @param max.time.off.tx.full.int Number of days off treatment for a full
#'        suppressor before onset of AIDS, including time before diagnosis.
#' @param max.time.on.tx.part.int Number of days on treatment for a
#'        partial suppressor before onset of AIDS.
#' @param max.time.off.tx.part.int Nnumber of days off treatment for a
#'        partial suppressor before onset of AIDS, including time before
#'        diagnosis.
#' @param vl.acute.rise.int Number of days to peak viremia during acute infection.
#' @param vl.acute.peak Peak viral load (in log10 units) at the height of acute infection.
#' @param vl.acute.fall.int Number of days from peak viremia to set-point
#'        viral load during the acute infection period.
#' @param vl.set.point Set point viral load (in log10 units).
#' @param vl.aids.onset.int Number of days to AIDS for a treatment-naive patient.
#' @param vl.aids.int Duration of AIDS stage infection in days.
#' @param vl.fatal Viral load in AIDS at which death occurs.
#' @param vl.full.supp Log10 viral load at full suppression on ART.
#' @param vl.part.supp Log10 viral load at partial suppression on ART.
#' @param full.supp.down.slope For full suppressors, number of log10 units that
#'        viral load falls per time step from treatment initiation or re-initiation
#'        until the level in \code{vl.full.supp}.
#' @param full.supp.up.slope For full suppressors, number of log10 units that
#'        viral load rises per time step from treatment halting until expected value.
#' @param part.supp.down.slope For partial suppressors, number of log10 units
#'        that viral load falls per time step from treatment initiation or
#'        re-initiation until the level in \code{vl.part.supp}.
#' @param part.supp.up.slope For partial suppressors, number of log10 units that
#'        viral load rises per time step from treatment halting until expected value.
#' @param b.R.S.rate Rate at which a person in race group R, where R is a 
#'        charater (B, BI, H, HI or W), and sex group S, where s is a charter (f, msf, msm, msmf),
#'        enters the population.
#' @param birth.age Age (in years) of new arrivals.
#' @param b.method Method for calculating the number of expected births at each
#'        time step, with \code{"fixed"} based on the number of persons at the
#'        initial time step and \code{"varying"} based on the current time step.
#' @param age.adj This value is added to the age for females in age.adj
#'        so that sex age asymmetry can be captured with an absdiff term  
#' @param msm.frac The fraction of males entering the population at each time step 
#'        as a male that has intercourse exclusivly with other males "msm" - 
#'        males not designated as msm or msmf are considered "msf".
#' @param msmf.frac  The fraction of males entering the population at each time step as
#'        as a male that has intercourse with both males and females "msmf" - 
#'        males not designated as msm or msmf are considered "msf".
#' @param URAI.prob Probability of transmission for a man having unprotected
#'        receptive anal intercourse with an infected man at set point viral
#'        load.
#' @param UIAI.prob Probability of transmission for an uncircumcised man having
#'        unprotected insertive anal intercourse with an infected man at set
#'        point viral load.
#' @param URVI.prob Probability of transmission for a woman having unprotected
#'        receptive vaginal intercourse with an infected man at set point viral
#'        load.
#' @param UIVI.prob Probability of transmission for an uncircumcised man having
#'        unprotected insertive vaginal intercourse with an infected women at set
#'        point viral load.
#' @param acute.rr Relative risk of infection (compared to that predicted by
#'        elevated viral load) when positive partner is in the acute stage.
#' @param circ.rr Relative risk of infection from insertive anal sex when the
#'        negative insertive partner is circumcised.
#' @param condom.rr Relative risk of infection from anal sex when a condom is
#'        used.
#' @param disc.outset.main.R.S.prob Probability that an HIV-infected person in race group R,
#'        where R is a charater (B, BI, H, HI or W), and sex group S, where s is a charter 
#'        (f, msf, msm, msmf), will disclose their status at the start of a main partnership.
#' @param disc.at.diag.main.R.S.prob Probability that a a person in race group R, where R is a 
#'        charater (B, BI, H, HI or W), and sex group S, where s is a charter (f, msf, msm, msmf),
#'        already in a main partnership will disclose at the time of diagnosis.
#' @param disc.post.diag.main.R.S.prob Probability that an HIV-infected a person in race group R, 
#'        where R is a charater (B, BI, H, HI or W), and sex group S, where s is a charter 
#'        (f, msf, msm, msmf), in a main partnership will disclose their status, assuming they didn't
#'        at the start of the partnership or at diagnosis.
#' @param disc.outset.pers.R.S.prob Probability that an HIV-infected a person in race group R, 
#'        where R is a charater (B, BI, H, HI or W), and sex group S, where s is a charter
#'        (f, msf, msm, msmf), will disclose their status at the start of a casual partnership.
#' @param disc.at.diag.pers.R.S.prob Probability that a person in race group R, where R is a
#'        charater (B, BI, H, HI or W), and sex group S, where s is a charter (f, msf, msm, msmf),
#'        already in a casual partnership will disclose at the time of diagnosis.
#' @param disc.post.diag.pers.R.S.prob Probability that an HIV-infected person in race group R,
#'        where R is a charater (B, BI, H, HI or W), and sex group S, where s is a charter 
#'        (f, msf, msm, msmf), in a casual partnership will disclose their status, assuming they
#'        didn't at the start of the partnership or at diagnosis.
#' @param disc.inst.R.S.prob Probability that an HIV-infected person in race group R, 
#'        where R is a charater (B, BI, H, HI or W), and sex group S, where s is a charter
#'        (f, msf, msm, msmf), will disclose ther status to a one-off partner.
#' @param circ.R.prob Probablity that a male in race group R, where R is a 
#'        charater (B, BI, H, HI or W), newly arriving in the population
#'        will be circumcised.
#' @param ccr5.R.S.prob Vector of length two of frequencies of the Delta 32
#'        mutation (homozygous and heterozygous, respectively) in the CCR5 gene
#'        among people in race group R, where R is a charater (B, BI, H, HI or W), 
#'        and sex group S, where s is a charter (f or m).
#' @param ccr5.heteroz.rr Relative risk of infection for people who are heterozygous
#'        in the CCR5 mutation.
#' @param num.inst.ai.classes Number of quantiles into which men should be
#'        divided in determining their levels of one-off anal intercourse.
#' @param base.ai.main.R.rate Expected coital frequency (for anal intercourse) for males 
#'        (msm and msmf) in race group R, where R is a charater (B, BI, H, HI or W), 
#'        in main partnerships (acts per day).
#' @param base.vi.main.R.rate Expected coital frequency (for vaginal intercourse) for 
#'        males and females in race group R, where R is a charater (B, BI, H, HI or W), 
#'        in main partnerships (acts per day).
#' @param base.ai.pers.R.rate Expected coital frequency (for anal intercourse) for males 
#'        (msm and msmf) in race group R, where R is a charater (B, BI, H, HI or W), 
#'        in casual partnerships (acts per day).
#' @param base.vi.pers.R.rate Expected coital frequency (for vaginal intercourse) for 
#'        males and females in race group R, where R is a charater (B, BI, H, HI or W), 
#'        in casual partnerships (acts per day).
#' @param ai.scale General relative scaler for all ai act rates for model
#'        calibration.
#' @param vi.scale General relative scaler for all vi act rates for model
#'        calibration.
#' @param cond.main.R.prob.msm Probability of condom use by an male of race group R, 
#'        where R is a charater (B, BI, H, HI or W), in a male-male main partnership.
#' @param cond.pers.always.prob.msm Fraction of men in male-male casual partnerships who always
#'        use condoms in those partnerships.
#' @param cond.pers.R.prob.msm Of men who are not consistent condom users, per-act
#'        probability of condom use in a casual male-male partnership by 
#'        males of race group R, where R is a charater (B, BI, H, HI or W).
#' @param cond.inst.always.prob.msm Fraction of men in male-male one-time partnerships who always
#'        use condoms in those partnerships.
#' @param cond.inst.R.prob.msm Of men who are not consistent condom users, per-act
#'        probability of condom use in a one-time male-male partnership by 
#'        males of race group R, where R is a charater (B, BI, H, HI or W).
#' @param cond.always.prob.corr.msm Correlation coefficient for probability of always
#'        using condoms in both casual and one-off male-male partnerships
#' @param cond.main.R.prob.het Probability of condom use by an person of race group R, 
#'        where R is a charater (B, BI, H, HI or W), in a heterosexual main partnership.
#' @param cond.pers.always.prob.het Fraction of people in heterosexual casual partnerships who always
#'        use condoms in those partnerships.
#' @param cond.pers.R.prob.het Of people who are not consistent condom users, per-act
#'        probability of condom use in a casual heterosexual partnership by 
#'        people of race group R, where R is a charater (B, BI, H, HI or W).
#' @param cond.inst.always.prob.het Fraction of people in heterosexual one-time partnerships who always
#'        use condoms in those partnerships.
#' @param cond.inst.R.prob.het Of people who are not consistent condom users, per-act
#'        probability of condom use in a one-time heterosexual partnership by 
#'        males of race group R, where R is a charater (B, BI, H, HI or W).
#' @param cond.always.prob.corr.het Correlation coefficient for probability of always
#'        using condoms in both casual and one-off heterosexual partnerships
#' @param cond.diag.main.beta.msm Beta multiplier for the log odds of using a
#'        condom in a male-male main partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.main.beta.msm Beta multiplier for the log odds of using a
#'        condom in a male-male main partnership if the HIV-infected man has disclosed.
#' @param cond.diag.pers.beta.msm Beta multiplier for the log odds of using a
#'        condom in a male-male casual partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.pers.beta.msm Beta multiplier for the log odds of using a
#'        condom in a male-male casual partnership if the HIV-infected man has disclosed
#'        his status.
#' @param cond.diag.inst.beta.msm Beta multiplier for the log odds of using a
#'        condom in a male-male one-off partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.inst.beta.msm Beta multiplier for the log odds of using a
#'        condom in a male-male one-off partnership if the HIV-infected man has disclosed
#'        his status.
#' @param cond.diag.main.beta.het Beta multiplier for the log odds of using a
#'        condom in a heterosexual main partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.main.beta.het Beta multiplier for the log odds of using a
#'        condom in a heterosexual main partnership if the HIV-infected man has disclosed.
#' @param cond.diag.pers.beta.het Beta multiplier for the log odds of using a
#'        condom in a heterosexual casual partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.pers.beta.het Beta multiplier for the log odds of using a
#'        condom in a heterosexual casual partnership if the HIV-infected man has disclosed
#'        his status.
#' @param cond.diag.inst.beta.het Beta multiplier for the log odds of using a
#'        condom in a heterosexual one-off partnership if the HIV-infected man has been
#'        diagnosed.
#' @param cond.discl.inst.beta.het Beta multiplier for the log odds of using a
#'        condom in a heterosexual one-off partnership if the HIV-infected man has disclosed
#'        his status.
#' @param vv.iev.R.prob Probability that a person of group R, where R is a charater (B, BI, H, HI or W),
#'        in a male male relationship will engage in intra-event versatility
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
#'        
#' @param immig.depart.R.S The probability per time step that a person of group R,
#' where R is a charater (BI or HI), and sex S, were S is a charater (f or m), leaves the population
#' being simulated by returning to their home country.
#' @param immig.return.R.S The probability per time step that a person of group R,
#' where R is a charater (BI or HI), and sex S, were S is a charater (f or m), who has left
#' the population being simulated by returning to their home country returns to the population
#' being simulated.
#' @param immig.aq.prob.R.S The probability per time step that a person of group R,
#' where R is a charater (BI or HI), and sex S, were S is a charater (f or m), who has left
#' the population being simulated aquires HIV while in their home country.
#' 
#' @param pos.entry.BI The probability that a BI entering into the population that is etering at an
#' age older than the birth age (replacing deaths from mortality proportional to the demographic
#' distribution and approximating in-migration) is HIV positive
#' 
#' @param pos.entry.HI The probability that a HI entering into the population that is etering at an
#' age older than the birth age (replacing deaths from mortality proportional to the demographic
#' distribution and approximating in-migration) is HIV positive
#' 
#' @param conc_dur_dx If TRUE the simulation will store cel.temp (a current edgelist with partnership start time) 
#' and cel.complete (a complete list of all partnerships including start and end times.  These are used for the analysis of concurrency
#' and multiplexity.
#' @param death_stats  IF TRUE the data.frame death.stats stores all care continuum information for individuals at the time of death.  
#' HIV related death data are store in death.stats, care continuum information for those that age out positive is store in death.stats$age.out.
#'
#'
#' @param Ecohab.window The number of time steps individuals in a cohab relationship will have the nodal attribute Ecohab=1.
#' @param p.growth If TRUE the number of additional births per time step can be set with 'p.growth.size' and the number of time steps with 'p.growth.nsteps'.
#' @param p.growth.nsteps The number of time steps for which there are 'p.growth.size' additional entries into the population if p.growth = TRUE.
#' @param p.growth.size The number of additioan entries made at each time step 2 to p.growth.nsteps if p.growth =TRUE.
#' @param add.demog.groups if TRUE will include the calcualtion of additional attributes that are conmbinations of attributes as well as
#' attributes that change based on concurrency status. These were used as attributes during the estimation process but no included in our final model. additioanl attribute combinations as an options 

#'        
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{param_KTM}, which can be passed to
#' EpiModel function \code{netsim}.
#'
#' @keywords shamp
#'
#' @export
param_KTM <- function(race.method = 1,
                      time.unit = 7,
                      method = 1,
                      age.unit = 52,
                      
                      ##For initializing status
                      last.neg.test.B.f.int = 1042,
                      last.neg.test.B.msf.int = 1446,
                      mean.test.B.f.int = 1042,
                      mean.test.B.msf.int = 1446,
                      
                      ##Kenya testing
                      ##Sex and age specific background testing rates to match KAIS 2012 ever tested and testing in the last year.
                      ##Rate are different for never testers and those that have tested before
                      test.prob.nevtest.m = c(0.004656910, 0.012052342, 0.006761833),
                      test.prob.nevtest.f = c(0.006597747, 0.035640496, 0.022500000),
                      test.prob.tested.m = c(0.014126, 0.019322, 0.01983, 0.01352),
                      test.prob.tested.f = c(0.015517, 0.021565, 0.020283, 0.015),
                      
                      #Two test types avaialble: antibody is standard, rna for the Kenya TM plus intervention.
                      test.window.int.ab = 3*7,
                      test.window.int.rna = 1*7,
                      
                      #Assuming 50:50 partial and full viral suppression
                      #From KARPR report 2018 77% of PLHIV achieve viral suppression
                      tt.traj.B.f.prob = c(0,0,.23,.77), 
                      tt.traj.B.msf.prob = c(0,0,.23,.77) ,
                     
                      #KTM INTERVENTION NONE TMP
                      intervention_TM ="NONE",
                      #From Sander 2015 65.6% within 3 weeks
                      seek.hc.AHI.prob = 1-(1-.3)^3,
                      
                      
                      ##sym.seek.prob can take on 4 values based on estimates of illness and treatment seeking.
                      ## (0.002606, 0.002994, 0.008179, 0.009395)  
                      sym.seek.prob = 0.002606,
                      
                      
                      #Probability of being tested when presenting with symptoms
                      sym.test.prob.bl = .28,
                      sym.test.prob.tm = .9,
                      
                      #NEED TO GET THE DURATION OF FOLLOWUP AND THE PROBABILITY OF PARTNER TESTING
                      #PArtner probability if based on being found and testing within 6 weeks.
                      partner.test.prob.bl = .023,
                      partner.test.prob.tm = .169,
                      p.prev.follow.win = 52,
                      p.acute.follow.win = 12,
                      PS.time = 6,
                      tx.init.PS.prob = .159,
                    
                      

                      #Art params worksheet 
                      tx.init.B.f.prob = 0.15912,
                      tx.init.BI.f.prob = 0.15912,
                      tx.init.H.f.prob = 0.15912,
                      tx.init.HI.f.prob = 0.15912,
                      tx.init.W.f.prob = 0.15912,
                      tx.init.B.msf.prob = 0.15912,
                      tx.init.BI.msf.prob = 0.15912,
                      tx.init.H.msf.prob = 0.15912,
                      tx.init.HI.msf.prob = 0.15912,
                      tx.init.W.msf.prob = 0.15912,
                      
                      #NO MSM or MSMF in Kenya TM
                      tx.init.B.msm.prob = 0.5000,
                      tx.init.BI.msm.prob = 0.5000,
                      tx.init.H.msm.prob = 0.5000,
                      tx.init.HI.msm.prob = 0.5000,
                      tx.init.W.msm.prob = 0.5000,
                      tx.init.B.msmf.prob = 0.5000,
                      tx.init.BI.msmf.prob = 0.5000,
                      tx.init.H.msmf.prob = 0.5000,
                      tx.init.HI.msmf.prob = 0.5000,
                      tx.init.W.msmf.prob = 0.5000,

                      #Art params worksheet 
                      tx.halt.B.f.prob = 0.0020245,
                      tx.halt.BI.f.prob = 0.0020245,
                      tx.halt.H.f.prob = 0.0020245,
                      tx.halt.HI.f.prob = 0.0020245,
                      tx.halt.W.f.prob = 0.0020245,
                      tx.halt.B.msf.prob = 0.0020245,
                      tx.halt.BI.msf.prob = 0.0020245,
                      tx.halt.H.msf.prob = 0.0020245,
                      tx.halt.HI.msf.prob = 0.0020245,
                      tx.halt.W.msf.prob = 0.0020245,
                      
                      #NO MSM or MSMF in Kenya TM
                      tx.halt.B.msm.prob = 0.5000,
                      tx.halt.BI.msm.prob = 0.5000,
                      tx.halt.H.msm.prob = 0.5000,
                      tx.halt.HI.msm.prob = 0.5000,
                      tx.halt.W.msm.prob = 0.5000,
                      tx.halt.B.msmf.prob = 0.5000,
                      tx.halt.BI.msmf.prob = 0.5000,
                      tx.halt.H.msmf.prob = 0.5000,
                      tx.halt.HI.msmf.prob = 0.5000,
                      tx.halt.W.msmf.prob = 0.5000,

                      ##re-initiation is the same as uptake
                      tx.reinit.B.f.prob = 0.00404902,
                      tx.reinit.BI.f.prob = 0.00404902,
                      tx.reinit.H.f.prob = 0.00404902,
                      tx.reinit.HI.f.prob = 0.00404902,
                      tx.reinit.W.f.prob = 0.00404902,
                      tx.reinit.B.msf.prob = 0.00404902,
                      tx.reinit.BI.msf.prob = 0.00404902,
                      tx.reinit.H.msf.prob = 0.00404902,
                      tx.reinit.HI.msf.prob = 0.00404902,
                      tx.reinit.W.msf.prob = 0.00404902,
                      
                      ##No MSM or MSMF in Kenya
                      tx.reinit.B.msm.prob = 1,
                      tx.reinit.BI.msm.prob = 1,
                      tx.reinit.H.msm.prob = 1,
                      tx.reinit.HI.msm.prob = 1,
                      tx.reinit.W.msm.prob = 1,
                      tx.reinit.B.msmf.prob = 1,
                      tx.reinit.BI.msmf.prob = 1,
                      tx.reinit.H.msmf.prob = 1,
                      tx.reinit.HI.msmf.prob = 1,
                      tx.reinit.W.msmf.prob = 1,

                      ##Lit updated
                      max.time.off.tx.full.int = 3297,
                      #NEED TO CONFIRM LIFE EXPECTANCY FOR PART
                      max.time.on.tx.part.int = 52 * 15 * 7,
                      max.time.off.tx.part.int = 3297,
                      vl.acute.rise.int = 12,
                      vl.acute.peak = 6.76,
                      vl.acute.fall.int = 20,
                      vl.set.point = 4.0,
                      vl.aids.onset.int = 3297,
                      vl.aids.int = 280,
                      vl.fatal = 7,
                      #Full is 550 copies
                      vl.full.supp = 2.75,
                      #Partial is 3000 copies
                      vl.part.supp = 3.5,
                      full.supp.down.slope = 0.25,
                      full.supp.up.slope = 0.25,
                      part.supp.down.slope = 0.25,
                      part.supp.up.slope = 0.25,
                      


                      #Births are set for consistent population demographics.
                      #b.B.f.rate = 1e-3 / 7,
                      #b.BI.f.rate = 1e-3 / 7,
                      #b.H.f.rate = 1e-3 / 7,
                      #b.HI.f.rate = 1e-3 / 7,
                      #b.W.f.rate = 1e-3 / 7,
                      #b.B.m.rate = 1e-3 / 7,
                      #b.BI.m.rate = 1e-3 / 7,
                      #b.H.m.rate = 1e-3 / 7,
                      #b.HI.m.rate = 1e-3 / 7,
                      #b.W.m.rate = 1e-3 / 7,

                      birth.age = 18,
                      exit.age = 76,
                      age.adj = 7.925357,
                      
                      b.method = "fixed",
                      #NA no msmf
                      msm.frac=0.0,
                      msmf.frac.B=0,
                      msmf.frac.BI=0,
                      msmf.frac.H=0,
                      msmf.frac.HI=0,
                      msmf.frac.W=0,
                      percent.male=0.4868,

                      URAI.prob = 138/10000,
                      UIAI.prob = 11/10000,
                      
                      #URVI.prob = 8/10000,
                      #UIVI.prob = 4/10000,
                      
                      URVI.prob = 0.0182,
                      UIVI.prob = 0.0182,
                      
                      
                      #VI range 1-12.
                      VI.foi.scale = 1,
                      AI.foi.scale = 1,
                      
                      acute.rr = 4,
                      circ.rr = .4,
                      condom.rr = 0.25,

                      ##For now aasume no disclosure (may use for calibration)
                      disc.outset.main.B.f.prob = .03,
                      disc.outset.main.BI.f.prob = 0,
                      disc.outset.main.H.f.prob = 0,
                      disc.outset.main.HI.f.prob = 0,
                      disc.outset.main.W.f.prob = 0,
                      disc.outset.main.B.msf.prob = 0,
                      disc.outset.main.BI.msf.prob = 0,
                      disc.outset.main.H.msf.prob = 0,
                      disc.outset.main.HI.msf.prob = 0,
                      disc.outset.main.W.msf.prob = 0,
                      disc.outset.main.B.msm.prob = 0,
                      disc.outset.main.BI.msm.prob = 0,
                      disc.outset.main.H.msm.prob = 0,
                      disc.outset.main.HI.msm.prob = 0,
                      disc.outset.main.W.msm.prob = 0,
                      disc.outset.main.B.msmf.prob = 0,
                      disc.outset.main.BI.msmf.prob = 0,
                      disc.outset.main.H.msmf.prob = 0,
                      disc.outset.main.HI.msmf.prob = 0,
                      disc.outset.main.W.msmf.prob = 0,
                      
                      disc.at.diag.main.B.f.prob = .03,
                      disc.at.diag.main.BI.f.prob = 0,
                      disc.at.diag.main.H.f.prob = 0,
                      disc.at.diag.main.HI.f.prob = 0,
                      disc.at.diag.main.W.f.prob = 0,
                      disc.at.diag.main.B.msf.prob = 0,
                      disc.at.diag.main.BI.msf.prob = 0,
                      disc.at.diag.main.H.msf.prob = 0,
                      disc.at.diag.main.HI.msf.prob = 0,
                      disc.at.diag.main.W.msf.prob = 0,
                      disc.at.diag.main.B.msm.prob = 0,
                      disc.at.diag.main.BI.msm.prob = 0,
                      disc.at.diag.main.H.msm.prob = 0,
                      disc.at.diag.main.HI.msm.prob = 0,
                      disc.at.diag.main.W.msm.prob = 0,
                      disc.at.diag.main.B.msmf.prob = 0,
                      disc.at.diag.main.BI.msmf.prob = 0,
                      disc.at.diag.main.H.msmf.prob = 0,
                      disc.at.diag.main.HI.msmf.prob = 0,
                      disc.at.diag.main.W.msmf.prob = 0,
                      
                      disc.post.diag.main.B.f.prob = 0,
                      disc.post.diag.main.BI.f.prob = 0,
                      disc.post.diag.main.H.f.prob = 0,
                      disc.post.diag.main.HI.f.prob = 0,
                      disc.post.diag.main.W.f.prob = 0,
                      disc.post.diag.main.B.msf.prob = 0,
                      disc.post.diag.main.BI.msf.prob = 0,
                      disc.post.diag.main.H.msf.prob = 0,
                      disc.post.diag.main.HI.msf.prob = 0,
                      disc.post.diag.main.W.msf.prob = 0,
                      disc.post.diag.main.B.msm.prob = 0,
                      disc.post.diag.main.BI.msm.prob = 0,
                      disc.post.diag.main.H.msm.prob = 0,
                      disc.post.diag.main.HI.msm.prob = 0,
                      disc.post.diag.main.W.msm.prob = 0,
                      disc.post.diag.main.B.msmf.prob = 0,
                      disc.post.diag.main.BI.msmf.prob = 0,
                      disc.post.diag.main.H.msmf.prob = 0,
                      disc.post.diag.main.HI.msmf.prob = 0,
                      disc.post.diag.main.W.msmf.prob = 0,
                      
                      
                      disc.outset.pers.B.f.prob = .02,
                      disc.outset.pers.BI.f.prob = 0,
                      disc.outset.pers.H.f.prob = 0,
                      disc.outset.pers.HI.f.prob = 0,
                      disc.outset.pers.W.f.prob = 0,
                      disc.outset.pers.B.msf.prob = 0,
                      disc.outset.pers.BI.msf.prob = 0,
                      disc.outset.pers.H.msf.prob = 0,
                      disc.outset.pers.HI.msf.prob = 0,
                      disc.outset.pers.W.msf.prob = 0,
                      disc.outset.pers.B.msm.prob = 0,
                      disc.outset.pers.BI.msm.prob = 0,
                      disc.outset.pers.H.msm.prob = 0,
                      disc.outset.pers.HI.msm.prob = 0,
                      disc.outset.pers.W.msm.prob = 0,
                      disc.outset.pers.B.msmf.prob = 0,
                      disc.outset.pers.BI.msmf.prob = 0,
                      disc.outset.pers.H.msmf.prob = 0,
                      disc.outset.pers.HI.msmf.prob = 0,
                      disc.outset.pers.W.msmf.prob = 0,
                      
                      disc.at.diag.pers.B.f.prob = .02,
                      disc.at.diag.pers.BI.f.prob = 0,
                      disc.at.diag.pers.H.f.prob = 0,
                      disc.at.diag.pers.HI.f.prob = 0,
                      disc.at.diag.pers.W.f.prob = 0,
                      disc.at.diag.pers.B.msf.prob = 0,
                      disc.at.diag.pers.BI.msf.prob = 0,
                      disc.at.diag.pers.H.msf.prob = 0,
                      disc.at.diag.pers.HI.msf.prob = 0,
                      disc.at.diag.pers.W.msf.prob = 0,
                      disc.at.diag.pers.B.msm.prob = 0,
                      disc.at.diag.pers.BI.msm.prob = 0,
                      disc.at.diag.pers.H.msm.prob = 0,
                      disc.at.diag.pers.HI.msm.prob = 0,
                      disc.at.diag.pers.W.msm.prob = 0,
                      disc.at.diag.pers.B.msmf.prob = 0,
                      disc.at.diag.pers.BI.msmf.prob = 0,
                      disc.at.diag.pers.H.msmf.prob = 0,
                      disc.at.diag.pers.HI.msmf.prob = 0,
                      disc.at.diag.pers.W.msmf.prob = 0,
                      
                      disc.post.diag.pers.B.f.prob = 0,
                      disc.post.diag.pers.BI.f.prob = 0,
                      disc.post.diag.pers.H.f.prob = 0,
                      disc.post.diag.pers.HI.f.prob = 0,
                      disc.post.diag.pers.W.f.prob = 0,
                      disc.post.diag.pers.B.msf.prob = 0,
                      disc.post.diag.pers.BI.msf.prob = 0,
                      disc.post.diag.pers.H.msf.prob = 0,
                      disc.post.diag.pers.HI.msf.prob = 0,
                      disc.post.diag.pers.W.msf.prob = 0, 
                      disc.post.diag.pers.B.msm.prob = 0,
                      disc.post.diag.pers.BI.msm.prob = 0,
                      disc.post.diag.pers.H.msm.prob = 0,
                      disc.post.diag.pers.HI.msm.prob = 0,
                      disc.post.diag.pers.W.msm.prob = 0, 
                      disc.post.diag.pers.B.msmf.prob = 0,
                      disc.post.diag.pers.BI.msmf.prob = 0,
                      disc.post.diag.pers.H.msmf.prob = 0,
                      disc.post.diag.pers.HI.msmf.prob = 0,
                      disc.post.diag.pers.W.msmf.prob = 0, 
                      
                      disc.inst.B.f.prob = .01,
                      disc.inst.BI.f.prob = 0,
                      disc.inst.H.f.prob = 0,
                      disc.inst.HI.f.prob = 0,
                      disc.inst.W.f.prob = 0,
                      disc.inst.B.msf.prob = 0,
                      disc.inst.BI.msf.prob = 0,
                      disc.inst.H.msf.prob = 0,
                      disc.inst.HI.msf.prob = 0,
                      disc.inst.W.msf.prob = 0,
                      disc.inst.B.msm.prob = 0,
                      disc.inst.BI.msm.prob = 0,
                      disc.inst.H.msm.prob = 0,
                      disc.inst.HI.msm.prob = 0,
                      disc.inst.W.msm.prob = 0,
                      disc.inst.B.msmf.prob = 0,
                      disc.inst.BI.msmf.prob = 0,
                      disc.inst.H.msmf.prob = 0,
                      disc.inst.HI.msmf.prob = 0,
                      disc.inst.W.msmf.prob = 0,

                      #From Kenya 2014 DHS
                      circ.B.prob = 0.925,
                      circ.BI.prob = 0.500,
                      circ.H.prob = 0.500,
                      circ.HI.prob = 0.500,
                      circ.W.prob = 0.500,
                      
                      #No CCR5 in Kenya 
                      ccr5.B.f.prob = c(0, 0.0),
                      ccr5.BI.f.prob = c(0, 0.0),
                      ccr5.H.f.prob = c(0, 0.0),
                      ccr5.HI.f.prob = c(0, 0.0),
                      ccr5.W.f.prob = c(0.0, 0.0),
                      ccr5.B.m.prob = c(0, 0.0),
                      ccr5.BI.m.prob = c(0, 0.0),
                      ccr5.H.m.prob = c(0, 0.0),
                      ccr5.HI.m.prob = c(0, 0.0),
                      ccr5.W.m.prob = c(0.0, 0.0),
                      ccr5.heteroz.rr = 0,

                      num.inst.ai.classes = 1,
                      #not modeling msm but will include home testing SEattle estimates as placeholders
                      base.ai.main.B.rate = 1.379/7,
                      base.ai.main.BI.rate = 1.379/7,
                      base.ai.main.H.rate = 1.379/7,
                      base.ai.main.HI.rate = 1.379/7,
                      base.ai.main.W.rate = 1.379/7,
                      base.ai.pers.B.rate = 0.883/7,
                      base.ai.pers.BI.rate = 0.883/7,
                      base.ai.pers.H.rate = 0.883/7,
                      base.ai.pers.HI.rate = 0.883/7,
                      base.ai.pers.W.rate = 0.883/7,
                      
                      #Waiting on KENYA TM
                      base.vi.main.B.rate = 0.4818 ,
                      base.vi.main.BI.rate = 0 ,
                      base.vi.main.H.rate = 0,
                      base.vi.main.HI.rate = 0,
                      base.vi.main.W.rate = 0,
                      base.vi.pers.B.rate = 0.2311,
                      base.vi.pers.BI.rate = 0,
                      base.vi.pers.H.rate = 0,
                      base.vi.pers.HI.rate = 0,
                      base.vi.pers.W.rate = 0,
                      ai.scale = 1,
                      vi.scale = 1,

                      
                      # NA for MSM
                      cond.main.B.prob.msm = 0.500,
                      cond.main.BI.prob.msm = 0.500,
                      cond.main.H.prob.msm = 0.500,
                      cond.main.HI.prob.msm = 0.500,
                      cond.main.W.prob.msm = 0.500,
                      cond.pers.always.prob.msm = 0.0001,
                      cond.pers.B.prob.msm = 0.500,
                      cond.pers.BI.prob.msm = 0.500,
                      cond.pers.H.prob.msm = 0.500,
                      cond.pers.HI.prob.msm = 0.500,
                      cond.pers.W.prob.msm = 0.500,
                      cond.inst.always.prob.msm = 0.0001,
                      cond.inst.B.prob.msm = 0.500,
                      cond.inst.BI.prob.msm = 0.500,
                      cond.inst.H.prob.msm = 0.500,
                      cond.inst.HI.prob.msm = 0.500,
                      cond.inst.W.prob.msm = 0.500,
                      cond.always.prob.corr.msm = 0.001,
                      
                      #From KAIS II 2014
                      cond.main.always.prob.het = .01,
                      cond.main.B.prob.het=.04,
                      
                      cond.pers.always.prob.het = .2,
                      cond.pers.B.prob.het=.09,
                      
                      cond.inst.always.prob.het = .46,
                      cond.inst.B.prob.het =.11,
                      
                      cond.always.prob.corr.het = .5,
                      
                      #NO MSM
                      cond.diag.main.beta.msm = -0.67,
                      cond.discl.main.beta.msm = -0.85,
                      cond.diag.pers.beta.msm = -0.67,
                      cond.discl.pers.beta.msm = -0.85,  
                      cond.diag.inst.beta.msm = -0.67,
                      cond.discl.inst.beta.msm = -0.85,
                      
                      
                      cond.diag.main.beta.het = -0.67,
                      cond.discl.main.beta.het = -0.85,
                      cond.diag.pers.beta.het = -0.67,
                      cond.discl.pers.beta.het = -0.85,  
                      cond.diag.inst.beta.het = -0.67,
                      cond.discl.inst.beta.het = -0.85,
                      
                      ##  ROLE NA
                      role.B.msm.prob = c(0.242, 0.321, 0.437),
                      role.BI.msm.prob = c(0.242, 0.321, 0.437),
                      role.H.msm.prob = c(0.242, 0.321, 0.437),
                      role.HI.msm.prob = c(0.242, 0.321, 0.437),
                      role.W.msm.prob = c(0.242, 0.321, 0.437),
                      role.B.msmf.prob = c(0.242, 0.321, 0.437),
                      role.BI.msmf.prob = c(0.242, 0.321, 0.437),
                      role.H.msmf.prob = c(0.242, 0.321, 0.437),
                      role.HI.msmf.prob = c(0.242, 0.321, 0.437),
                      role.W.msmf.prob = c(0.242, 0.321, 0.437),

                      ##NA position not used for HET
                      vv.iev.B.prob = 0.47874,
                      vv.iev.BI.prob = 0.47874,
                      vv.iev.H.prob = 0.47874,
                      vv.iev.HI.prob = 0.47874,
                      vv.iev.W.prob = 0.47874,

                      ##PrEP.
                      prep.start = Inf,
                      prep.elig.model = "NONE",
                #need a real values for these
                      prep.class.prob = c(0, 0, 0, 1),
                      ##46% reduction in aquisition risk from nam aidsmap May 2020 Roger Pebody / Roger Chou et al 2019 JAMA
                      prep.class.hr = c(1, 0, 0, 0.54),
                      prep.coverage = 0,
                      prep.cov.method = "curr",
                      prep.cov.rate = 1,
                      prep.tst.int = 90,
                      prep.risk.int = 182,
                      prep.risk.reassess = TRUE,
                #need a real values for these
                      prep.start.prob = .236,
                      prep.stop.prob = .0414,
                      prep.window = 6,
                      
                      

                      conc_dur_dx = TRUE,
                      cel.complete.trim = TRUE,
                      cel.complete.lag = 52,
                      kill.poi.rels = TRUE,
                      death_stats = FALSE,
                      
                      p.growth = FALSE,
                      p.growth.size = 0,
                      
                      ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))



  p$time.unit <- time.unit

  intvars <- grep(names(p), pattern = ".int", fixed = TRUE)
  p[intvars] <- lapply(p[intvars], FUN = function(x) round(x / p$time.unit))

  ratevars <- grep(names(p), pattern = ".rate", fixed = TRUE)
  p[ratevars] <- lapply(p[ratevars], FUN = function(x) x * p$time.unit)

  
  p$role.B.msm.prob <- role.B.msm.prob
  p$role.BI.msm.prob <- role.BI.msm.prob
  p$role.H.msm.prob <- role.H.msm.prob
  p$role.HI.msm.prob <- role.HI.msm.prob
  p$role.W.msm.prob <- role.W.msm.prob
  
  p$role.B.msmf.prob <- role.B.msmf.prob
  p$role.BI.msmf.prob <- role.BI.msmf.prob
  p$role.H.msmf.prob <- role.H.msmf.prob
  p$role.HI.msmf.prob <- role.HI.msmf.prob
  p$role.W.msmf.prob <- role.W.msmf.prob

  p$inst.trans.matrix <- matrix(1, nrow = 1)
  p$role.trans.matrix <- matrix(c(1, 0, 0,
                                  0, 1, 0,
                                  0, 0, 1),
                                nrow = 3)


  p$riskh.start <- max(1, prep.start - prep.risk.int - 1)

  p$method <- method
  p$modes <- 1

  
#  p$asmr.B.f <- c(rep(0, 17),
#                1-(1-c(rep(0.001941328, 7),
#                       rep(0.003499821, 50)))^(1 / age.unit),1)
#  
#  p$asmr.B.m <- c(rep(0, 17),
#                1-(1-c(rep(0.001941328, 7),
#                       rep(0.003499821, 50)))^(1 / age.unit),1)
#  
  p$asmr.B.f <- c(rep(0, 17),
                1-(1-c(rep(0.001941328, 7),
                       rep(0.003499821, 15),
                       rep(0.003499821, 36)))^(1 / age.unit),1)
  
  p$asmr.B.m <- c(rep(0, 17),
                1-(1-c(rep(0.001941328, 7),
                       rep(0.003499821, 15),
                       rep(0.003499821, 36)))^(1 / age.unit),1)
  

  class(p) <- "param.net"
  return(p)
}


#' @title Epidemic Model Initial Conditions for SHAMP
#'
#' @description Sets the initial conditions for a stochastic epidemic models
#'              simulated with \code{\link{netsim}}.
#'
#' @param prev.R.S Initial disease prevalence among persons of group R,
#'  where R is a charater (B, BI, H, HI or W), and sex group S
#'        where s is a charter (f, msf, msm, msmf)..
#' @param ... Additional arguments passed to function.
#'
#' @return
#' A list object of class \code{init_shamp}, which can be passed to EpiModel
#' function \code{\link{netsim}}.
#'
#' @keywords SHAMP
#'
#' @export
init_KTM <- function(prev.B.f = 0.02,
                     prev.B.msf = 0.02,

                     ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  
  p$init.prev.age.slope.B.f <- prev.B.f / 21
  p$init.prev.age.slope.B.msf <- prev.B.msf / 21

 

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
#' @param ConcDurDx.FUN Module that analyzes concurrency durations by type (Dianistics only, do not use for general simualtions)
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
#' A list object of class \code{control_KTM}, which can be passed to the
#' EpiModel function \code{netsim}.
#'
#' @keywords SHAMP 
#'
#' @export
control_KTM <- function(simno = 1,
                        nsims = 1,
                        ncores = 1,
                        nsteps = 100,
                        start = 1,
                        #tergmLite = TRUE,
                        initialize.FUN = initialize_KTM,
                        aging.FUN = aging_KTM,
                        deaths.FUN = deaths_KTM,
                        births.FUN = births_KTM,
                        test.FUN = test_KTM,
                        pservices.FUN = pservices_KTM,
                        tx.FUN = tx_KTM,
                        prep.FUN = prep_shamp,
                        progress.FUN = progress_msm,
                        vl.FUN = vl_shamp,
                        #aiclass.FUN = NULL,
                        #roleclass.FUN = NULL,
                        resim_nets.FUN = simnet_shamp,
                        disclose.FUN = disclose_shamp,
                        ConcDurDx.FUN = ConcDurDx_shamp,
                        acts.FUN = acts_KTM,
                        condoms.FUN = condoms_KTM,
                        riskhist.FUN = riskhist_shamp,
                        position.FUN = position_shamp,
                        trans.FUN = trans_KTM,
                        prev.FUN = prevalence_KTM,
                        verbose.FUN = verbose_KTM,
                        save.nwstats = FALSE,
                        save.other = "attr",
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

  p$save.other = c("attr", "temp", "el", "p" ,"trans.el", "death.stats", "cel.temp", "cel.complete")

  p$save.network = FALSE

  class(p) <- "control.net"
  return(p)
}

