
############# Read in the excel parameters sheet ==============

library(tidyverse)

# A. Read in the Excel sheet.
parameter_sheet<- readxl::read_excel( 
  "data/Model Parameters.Sara.work.xlsx",
  sheet = 1) %>% 
  select(`Parameter name`, `SHAMP Estimates`) %>% # only grab the parameter name and SHAMP estimates
  rename(param = 1,
         est = 2)   # rename the parameters for easier viewing
  
  # B. Clean the parameters sheet to only have the estimate and non-missing parameters
  parameter_sheet<-parameter_sheet %>% mutate(est = as.numeric(est)) %>% 
  filter(is.na(est)==F,# only have non-missing estimates or parameters
         is.na(param)==F)

param_list<-list()
# C. Recursively save all of the parameters as single objects in the global env.
param_list$est<-parameter_sheet$est
param_list$params<-parameter_sheet$param

# create a temporary list to store all the parameters
temp_list<-list()

# take out each estimate, and apply it to the empty list
for (i in 1:length(param_list$params)){
  temp_list[i]<-param_list$est[i]
}

# apply the name of the parameter to each element in the list
names(temp_list)<-c(param_list$params)

# finally, bring all the parameters to the environment
list2env(temp_list,globalenv())

# remove the temporary list
rm(temp_list)

#### END - reading in the parameter sheet

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
#' @param testing.pattern Method for HIV testing, with options \code{"memoryless"}
#'        for constant hazard without regard to time since previous test, or
#'        \code{"interval"} deterministic fixed intervals.
#' @param test.window.int Length of the HIV test window period in days.
#' @param aids.sypm.test.prob probability per time steap of testing given stage 4 symptoms
#' @param tt.traj.R.S.prob Proportion of race group R, where R is a charater (B, BI, H, HI or W), 
#'        and sex group S, where s is a charter (f, msf, msm, msmf), who enter one of four
#'        testing/treatment trajectories: never test or treat, test and never
#'        initiate treatment, test and treated with partial viral suppression,
#'        and test and treated with full suppression.
#' @param tx.init.R.S.prob Probability per time step that a person in race group R, where R is a 
#'        charater (B, BI, H, HI or W), and sex group S, where s is a charter (f, msf, msm, msmf),
#'        who has tested positive will initiate treatment.
#' @param tx.halt.R.S.prob Probability per time step that a person in race group R, where R is a 
#'        charater (B, BI, H, HI or W), and sex group S, where s is a charter (f, msf, msm, msmf),
#'        who is currently on treatment will halt treatment.
#' @param tx.reinit.R.S.prob Probability per time step that a person in race group R, where R is a 
#'        charater (B, BI, H, HI or W), and sex group S, where s is a charter (f, msf, msm, msmf),
#'        who is not currently on treatment but who has been in the past will
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
#' @param b.R.S.rate Rate at which a person in race group R, where R is a 
#'        charater (B, BI, H, HI or W), and sex group S, where s is a charter (f, msf, msm, msmf),
#'        enters the population.
#' @param birth.age Age (in years) of new arrivals.
#' @param b.method Method for calculating the number of expected births at each
#'        time step, with \code{"fixed"} based on the number of persons at the
#'        initial time step and \code{"varying"} based on the current time step.
#' @param age.adj This value is added to the squareroot of age for females in sqrt.age.adj
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
#' A list object of class \code{param_shamp}, which can be passed to
#' EpiModel function \code{netsim}.
#'
#' @keywords shamp
#'
#' @export
param_shamp <- function(race.method = 1,
                      time.unit = 7,
                      method = 1,
                      age.unit = 52,
                      
                      ##Martina's NSFG estimates
                      last.neg.test.B.f.int = last.neg.test.B.f.int,
                      last.neg.test.BI.f.int = last.neg.test.BI.f.int,
                      last.neg.test.H.f.int = last.neg.test.H.f.int,
                      last.neg.test.HI.f.int = last.neg.test.HI.f.int,
                      last.neg.test.W.f.int = last.neg.test.W.f.int,
                      last.neg.test.B.msf.int = last.neg.test.B.msf.int,
                      last.neg.test.BI.msf.int = last.neg.test.BI.msf.int,
                      last.neg.test.H.msf.int = last.neg.test.H.msf.int,
                      last.neg.test.HI.msf.int = last.neg.test.HI.msf.int,
                      last.neg.test.W.msf.int = last.neg.test.W.msf.int,
                      last.neg.test.B.msm.int = last.neg.test.B.msm.int,
                      last.neg.test.BI.msm.int = last.neg.test.BI.msm.int,
                      last.neg.test.H.msm.int = last.neg.test.H.msm.int,
                      last.neg.test.HI.msm.int = last.neg.test.HI.msm.int,
                      last.neg.test.W.msm.int = last.neg.test.W.msm.int,
                      last.neg.test.B.msmf.int = last.neg.test.B.msmf.int,
                      last.neg.test.BI.msmf.int = last.neg.test.BI.msmf.int,
                      last.neg.test.H.msmf.int = last.neg.test.H.msmf.int,
                      last.neg.test.HI.msmf.int = last.neg.test.HI.msmf.int,
                      last.neg.test.W.msmf.int = last.neg.test.W.msmf.int,
                      mean.test.B.f.int = last.neg.test.B.f.int,
                      mean.test.BI.f.int = last.neg.test.BI.f.int,
                      mean.test.H.f.int = last.neg.test.H.f.int,
                      mean.test.HI.f.int = last.neg.test.HI.f.int,
                      mean.test.W.f.int = last.neg.test.W.f.int,
                      mean.test.B.msf.int = last.neg.test.B.msf.int,
                      mean.test.BI.msf.int = last.neg.test.BI.msf.int,
                      mean.test.H.msf.int = last.neg.test.H.msf.int,
                      mean.test.HI.msf.int = last.neg.test.HI.msf.int,
                      mean.test.W.msf.int = last.neg.test.W.msf.int,
                      mean.test.B.msm.int = last.neg.test.B.msm.int,
                      mean.test.BI.msm.int = last.neg.test.BI.msm.int,
                      mean.test.H.msm.int = last.neg.test.H.msm.int,
                      mean.test.HI.msm.int = last.neg.test.HI.msm.int,
                      mean.test.W.msm.int = last.neg.test.W.msm.int,
                      mean.test.B.msmf.int = last.neg.test.B.msmf.int,
                      mean.test.BI.msmf.int = last.neg.test.BI.msmf.int,
                      mean.test.H.msmf.int = last.neg.test.H.msmf.int,
                      mean.test.HI.msmf.int = last.neg.test.HI.msmf.int,
                      mean.test.W.msmf.int = last.neg.test.W.msmf.int,
                      testing.pattern = "memoryless",
                      test.window.int = test.window.int,
                      aids.sypm.test.prob = 1/26,
                      
                      #change trajectories based on NSFG NEver test fractions among 35-40
                      #set MF ratio to observed population ratio but within race
                      
                      #tt.traj.B.f.prob = c(0.741, 0.000, 0.053, 0.206),
                      #tt.traj.BI.f.prob = c(0.741, 0.000, 0.053, 0.206),
                      #tt.traj.H.f.prob = c(0.741, 0.000, 0.041, 0.218),
                      #tt.traj.HI.f.prob = c(0.741, 0.000, 0.041, 0.218),
                      #tt.traj.W.f.prob = c(0.741, 0.000, 0.057, 0.202),
                      #tt.traj.B.msf.prob = c(0.741, 0.000, 0.066, 0.193),
                      #tt.traj.BI.msf.prob = c(0.741, 0.000, 0.066, 0.193),
                      #tt.traj.H.msf.prob = c(0.741, 0.000, 0.067, 0.192),
                      #tt.traj.HI.msf.prob = c(0.741, 0.000, 0.067, 0.192),
                      #tt.traj.W.msf.prob = c(0.741, 0.000, 0.052, 0.207),
                      #tt.traj.B.msm.prob = c(.025, 0, .065, .911),
                      #tt.traj.BI.msm.prob = c(.025, 0, .065, .911),
                      #tt.traj.H.msm.prob = c(.025, 0, .065, .911),
                      #tt.traj.HI.msm.prob = c(.025, 0, .065, .911),
                      #tt.traj.W.msm.prob = c(.025, 0, .065, .911),
                      #tt.traj.B.msmf.prob = c(0.741, 0.000, 0.066, 0.193),
                      #tt.traj.BI.msmf.prob = c(0.741, 0.000, 0.066, 0.193),
                      #tt.traj.H.msmf.prob = c(0.741, 0.000, 0.041, 0.218),
                      #tt.traj.HI.msmf.prob = c(0.741, 0.000, 0.041, 0.218),
                      #tt.traj.W.msmf.prob = c(0.741, 0.000, 0.052, 0.207),
                      
                      tt.traj.B.f.prob = c(1-b_f/100, 0, 1-(1-b_f/100)-( b_f/100*het_f_b/100) , b_f/100*het_f_b/100),
                      tt.traj.BI.f.prob = c(1-bi_f/100, 0, 1-(1-bi_f/100)-( b_f/100*het_f_b/100) , bi_f/100*het_f_b/100),
                      tt.traj.H.f.prob =  c(1-h_f/100, 0, 1-(1-h_f/100)-( h_f/100*het_f_h/100) , h_f/100*het_f_h/100),
                      tt.traj.HI.f.prob = c(1-hi_f/100, 0, 1-(1-hi_f/100)-( hi_f/100*het_f_h/100) , hi_f/100*het_f_h/100),
                      tt.traj.W.f.prob = c(1-w_f/100, 0, 1-(1-w_f/100)-( w_f/100*het_f_w/100) , w_f/100*het_f_w/100),
                      
                      
                      tt.traj.B.msf.prob = c(1-b_m/100, 0, 1-(1-b_m/100)-( b_m/100*het_m_b/100) , b_m/100*het_m_b/100),
                      tt.traj.BI.msf.prob = c(1-bi_m/100, 0, 1-(1-bi_m/100)-( b_m/100*het_m_b/100) , bi_m/100*het_m_b/100),
                      tt.traj.H.msf.prob =  c(1-h_m/100, 0, 1-(1-h_m/100)-( h_m/100*het_m_h/100) , h_m/100*het_m_h/100),
                      tt.traj.HI.msf.prob = c(1-hi_m/100, 0, 1-(1-hi_m/100)-( hi_m/100*het_m_h/100) , hi_m/100*het_m_h/100),
                      tt.traj.W.msf.prob = c(1-w_m/100, 0, 1-(1-w_m/100)-( w_m/100*het_m_w/100) , w_m/100*het_m_w/100),
                      
                      tt.traj.B.msmf.prob = c(1-b_m/100, 0, 1-(1-b_m/100)-( b_m/100*het_m_b/100) , b_m/100*het_m_b/100),
                      tt.traj.BI.msmf.prob = c(1-bi_m/100, 0, 1-(1-bi_m/100)-( b_m/100*het_m_b/100) , bi_m/100*het_m_b/100),
                      tt.traj.H.msmf.prob =  c(1-h_m/100, 0, 1-(1-h_m/100)-( h_m/100*het_m_h/100) , h_m/100*het_m_h/100),
                      tt.traj.HI.msmf.prob = c(1-hi_m/100, 0, 1-(1-hi_m/100)-( hi_m/100*het_m_h/100) , hi_m/100*het_m_h/100),
                      tt.traj.W.msmf.prob = c(1-w_m/100, 0, 1-(1-w_m/100)-( w_m/100*het_m_w/100) , w_m/100*het_m_w/100),
                      
                      tt.traj.B.msm.prob = c(1-b_m/100, 0, 1-(1-b_m/100)-( b_m/100*msm_b/100) , b_m/100*msm_b/100),
                      tt.traj.BI.msm.prob = c(1-bi_m/100, 0, 1-(1-bi_m/100)-( b_m/100*msm_b/100) , bi_m/100*msm_b/100),
                      tt.traj.H.msm.prob =  c(1-h_m/100, 0, 1-(1-h_m/100)-( h_m/100*msm_h/100) , h_m/100*msm_h/100),
                      tt.traj.HI.msm.prob = c(1-hi_m/100, 0, 1-(1-hi_m/100)-( hi_m/100*msm_h/100) , hi_m/100*msm_h/100),
                      tt.traj.W.msm.prob = c(1-w_m/100, 0, 1-(1-w_m/100)-( w_m/100*msm_w/100) , w_m/100*msm_w/100),
                      # 
                      # tt.traj.B.msm.prob = c(.025, 0, .065, .911),
                      # tt.traj.BI.msm.prob = c(.025, 0, .065, .911),
                      # tt.traj.H.msm.prob = c(.025, 0, .065, .911),
                      # tt.traj.HI.msm.prob = c(.025, 0, .065, .911),
                      # tt.traj.W.msm.prob = c(.025, 0, .065, .911),
                      # 
                      
                      tx.init.B.f.prob = tx.init.B.f.prob,
                      tx.init.BI.f.prob = tx.init.BI.f.prob,
                      tx.init.H.f.prob = tx.init.H.f.prob,
                      tx.init.HI.f.prob = tx.init.HI.f.prob,
                      tx.init.W.f.prob = tx.init.W.f.prob,
                      tx.init.B.msf.prob = tx.init.B.msf.prob,
                      tx.init.BI.msf.prob = tx.init.BI.msf.prob,
                      tx.init.H.msf.prob = tx.init.H.msf.prob,
                      tx.init.HI.msf.prob = tx.init.HI.msf.prob,
                      tx.init.W.msf.prob = tx.init.W.msf.prob,
                      tx.init.B.msm.prob = tx.init.B.msm.prob,
                      tx.init.BI.msm.prob = tx.init.BI.msm.prob,
                      tx.init.H.msm.prob = tx.init.H.msm.prob,
                      tx.init.HI.msm.prob = tx.init.HI.msm.prob,
                      tx.init.W.msm.prob = tx.init.W.msm.prob,
                      tx.init.B.msmf.prob = tx.init.B.msmf.prob,
                      tx.init.BI.msmf.prob = tx.init.BI.msmf.prob,
                      tx.init.H.msmf.prob = tx.init.H.msmf.prob,
                      tx.init.HI.msmf.prob = tx.init.HI.msmf.prob,
                      tx.init.W.msmf.prob = tx.init.W.msmf.prob,

                      #No halting in the basic SHAMP model
                      tx.halt.B.f.prob = halt,
                      tx.halt.BI.f.prob = halt,
                      tx.halt.H.f.prob = halt,
                      tx.halt.HI.f.prob = halt,
                      tx.halt.W.f.prob = halt,
                      tx.halt.B.msf.prob = halt,
                      tx.halt.BI.msf.prob = halt,
                      tx.halt.H.msf.prob = halt,
                      tx.halt.HI.msf.prob = halt,
                      tx.halt.W.msf.prob = halt,
                      tx.halt.B.msm.prob = halt,
                      tx.halt.BI.msm.prob = halt,
                      tx.halt.H.msm.prob = halt,
                      tx.halt.HI.msm.prob = halt,
                      tx.halt.W.msm.prob = halt,
                      tx.halt.B.msmf.prob = halt,
                      tx.halt.BI.msmf.prob = halt,
                      tx.halt.H.msmf.prob = halt,
                      tx.halt.HI.msmf.prob = halt,
                      tx.halt.W.msmf.prob = halt,

                      ##With no halting there is also no re-initiation
                      tx.reinit.B.f.prob = reinit,
                      tx.reinit.BI.f.prob = reinit,
                      tx.reinit.H.f.prob = reinit,
                      tx.reinit.HI.f.prob = reinit,
                      tx.reinit.W.f.prob = reinit,
                      tx.reinit.B.msf.prob = reinit,
                      tx.reinit.BI.msf.prob = reinit,
                      tx.reinit.H.msf.prob = reinit,
                      tx.reinit.HI.msf.prob = reinit,
                      tx.reinit.W.msf.prob = reinit,
                      tx.reinit.B.msm.prob = reinit,
                      tx.reinit.BI.msm.prob = reinit,
                      tx.reinit.H.msm.prob = reinit,
                      tx.reinit.HI.msm.prob = reinit,
                      tx.reinit.W.msm.prob = reinit,
                      tx.reinit.B.msmf.prob = reinit,
                      tx.reinit.BI.msmf.prob = reinit,
                      tx.reinit.H.msmf.prob = reinit,
                      tx.reinit.HI.msmf.prob = reinit,
                      tx.reinit.W.msmf.prob = reinit,

                      max.time.off.tx.full.int = max.time.off.tx.full.int,
                      max.time.on.tx.part.int = max.time.on.tx.part.int,
                      max.time.off.tx.part.int = max.time.off.tx.part.int,
                      vl.acute.rise.int = vl.acute.rise.int,
                      vl.acute.peak = vl.acute.peak,
                      vl.acute.fall.int = vl.acute.fall.int,
                      vl.set.point = vl.set.point,
                      vl.aids.onset.int = vl.aids.onset.int,
                      vl.aids.int = vl.aids.int,
                      vl.fatal = vl.fatal,
                      vl.full.supp = vl.full.supp,
                      vl.part.supp = vl.part.supp,
                      full.supp.down.slope = full.supp.down.slope,
                      full.supp.up.slope = full.supp.up.slope,
                      part.supp.down.slope = part.supp.down.slope,
                      part.supp.up.slope = part.supp.up.slope,

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

                      birth.age = birth.age,
                      exit.age = exit.age,
                      b.method = "fixed",
                      msm.frac=msm.frac,
                      msmf.frac.B=msmf.frac.B,
                      msmf.frac.BI=msmf.frac.BI,
                      msmf.frac.H=msmf.frac.H,
                      msmf.frac.HI=msmf.frac.HI,
                      msmf.frac.W=msmf.frac.W,

                      URAI.prob = URAI.prob/10000,
                      UIAI.prob = UIAI.prob/10000,
                      
                      URVI.prob = URVI.prob/10000,
                      UIVI.prob = UIVI.prob/10000,
                      
                      #VI range 1-12.
                      VI.foi.scale = vi.scale,
                      AI.foi.scale = ai.scale,
                      
                      acute.rr = acute.rr,
                      circ.rr = circ.rr,
                      condom.rr = condom.rr,

                      ##For now assume no disclosure (may use for calibration)
                      # for sara: set these all to the same thing as diag thing
                      disc.outset.main.B.f.prob=disc.outset.main.B.f.prob,
                      disc.outset.main.BI.f.prob=disc.outset.main.BI.f.prob,
                      disc.outset.main.H.f.prob=disc.outset.main.H.f.prob,
                      disc.outset.main.HI.f.prob=disc.outset.main.HI.f.prob,
                      disc.outset.main.W.f.prob=disc.outset.main.W.f.prob,
                      disc.outset.main.B.msf.prob=disc.outset.main.B.msf.prob,
                      disc.outset.main.BI.msf.prob=disc.outset.main.BI.msf.prob,
                      disc.outset.main.H.msf.prob=disc.outset.main.H.msf.prob,
                      disc.outset.main.HI.msf.prob=disc.outset.main.HI.msf.prob,
                      disc.outset.main.W.msf.prob=disc.outset.main.W.msf.prob,
                      disc.outset.main.B.msm.prob=disc.outset.main.B.msm.prob,
                      disc.outset.main.BI.msm.prob=disc.outset.main.BI.msm.prob,
                      disc.outset.main.H.msm.prob=disc.outset.main.H.msm.prob,
                      disc.outset.main.HI.msm.prob=disc.outset.main.HI.msm.prob,
                      disc.outset.main.W.msm.prob=disc.outset.main.W.msm.prob,
                      disc.outset.main.B.msmf.prob=disc.outset.main.B.msmf.prob,
                      disc.outset.main.BI.msmf.prob=disc.outset.main.BI.msmf.prob,
                      disc.outset.main.H.msmf.prob=disc.outset.main.H.msmf.prob,
                      disc.outset.main.HI.msmf.prob=disc.outset.main.HI.msmf.prob,
                      disc.outset.main.W.msmf.prob=disc.outset.main.W.msmf.prob,
                      
                      disc.at.diag.main.B.f.prob=disc.outset.main.B.f.prob,
                      disc.at.diag.main.BI.f.prob=disc.outset.main.BI.f.prob,
                      disc.at.diag.main.H.f.prob=disc.outset.main.H.f.prob,
                      disc.at.diag.main.HI.f.prob=disc.outset.main.HI.f.prob,
                      disc.at.diag.main.W.f.prob=disc.outset.main.W.f.prob,
                      disc.at.diag.main.B.msf.prob=disc.outset.main.B.msf.prob,
                      disc.at.diag.main.BI.msf.prob=disc.outset.main.BI.msf.prob,
                      disc.at.diag.main.H.msf.prob=disc.outset.main.H.msf.prob,
                      disc.at.diag.main.HI.msf.prob=disc.outset.main.HI.msf.prob,
                      disc.at.diag.main.W.msf.prob=disc.outset.main.W.msf.prob,
                      disc.at.diag.main.B.msm.prob=disc.outset.main.B.msm.prob,
                      disc.at.diag.main.BI.msm.prob=disc.outset.main.BI.msm.prob,
                      disc.at.diag.main.H.msm.prob=disc.outset.main.H.msm.prob,
                      disc.at.diag.main.HI.msm.prob=disc.outset.main.HI.msm.prob,
                      disc.at.diag.main.W.msm.prob=disc.outset.main.W.msm.prob,
                      disc.at.diag.main.B.msmf.prob=disc.outset.main.B.msmf.prob,
                      disc.at.diag.main.BI.msmf.prob=disc.outset.main.BI.msmf.prob,
                      disc.at.diag.main.H.msmf.prob=disc.outset.main.H.msmf.prob,
                      disc.at.diag.main.HI.msmf.prob=disc.outset.main.HI.msmf.prob,
                      disc.at.diag.main.W.msmf.prob=disc.outset.main.W.msmf.prob,
                      
                      disc.post.diag.main.B.f.prob=disc.outset.main.B.f.prob,
                      disc.post.diag.main.BI.f.prob=disc.outset.main.BI.f.prob,
                      disc.post.diag.main.H.f.prob=disc.outset.main.H.f.prob,
                      disc.post.diag.main.HI.f.prob=disc.outset.main.HI.f.prob,
                      disc.post.diag.main.W.f.prob=disc.outset.main.W.f.prob,
                      disc.post.diag.main.B.msf.prob=disc.outset.main.B.msf.prob,
                      disc.post.diag.main.BI.msf.prob=disc.outset.main.BI.msf.prob,
                      disc.post.diag.main.H.msf.prob=disc.outset.main.H.msf.prob,
                      disc.post.diag.main.HI.msf.prob=disc.outset.main.HI.msf.prob,
                      disc.post.diag.main.W.msf.prob=disc.outset.main.W.msf.prob,
                      disc.post.diag.main.B.msm.prob=disc.outset.main.B.msm.prob,
                      disc.post.diag.main.BI.msm.prob=disc.outset.main.BI.msm.prob,
                      disc.post.diag.main.H.msm.prob=disc.outset.main.H.msm.prob,
                      disc.post.diag.main.HI.msm.prob=disc.outset.main.HI.msm.prob,
                      disc.post.diag.main.W.msm.prob=disc.outset.main.W.msm.prob,
                      disc.post.diag.main.B.msmf.prob=disc.outset.main.B.msmf.prob,
                      disc.post.diag.main.BI.msmf.prob=disc.outset.main.BI.msmf.prob,
                      disc.post.diag.main.H.msmf.prob=disc.outset.main.H.msmf.prob,
                      disc.post.diag.main.HI.msmf.prob=disc.outset.main.HI.msmf.prob,
                      disc.post.diag.main.W.msmf.prob=disc.outset.main.W.msmf.prob,
                      
                      
                      disc.outset.pers.B.f.prob=disc.outset.main.B.f.prob,
                      disc.outset.pers.BI.f.prob=disc.outset.main.BI.f.prob,
                      disc.outset.pers.H.f.prob=disc.outset.main.H.f.prob,
                      disc.outset.pers.HI.f.prob=disc.outset.main.HI.f.prob,
                      disc.outset.pers.W.f.prob=disc.outset.main.W.f.prob,
                      disc.outset.pers.B.msf.prob=disc.outset.main.B.msf.prob,
                      disc.outset.pers.BI.msf.prob=disc.outset.main.BI.msf.prob,
                      disc.outset.pers.H.msf.prob=disc.outset.main.H.msf.prob,
                      disc.outset.pers.HI.msf.prob=disc.outset.main.HI.msf.prob,
                      disc.outset.pers.W.msf.prob=disc.outset.main.W.msf.prob,
                      disc.outset.pers.B.msm.prob=disc.outset.main.B.msm.prob,
                      disc.outset.pers.BI.msm.prob=disc.outset.main.BI.msm.prob,
                      disc.outset.pers.H.msm.prob=disc.outset.main.H.msm.prob,
                      disc.outset.pers.HI.msm.prob=disc.outset.main.HI.msm.prob,
                      disc.outset.pers.W.msm.prob=disc.outset.main.W.msm.prob,
                      disc.outset.pers.B.msmf.prob=disc.outset.main.B.msmf.prob,
                      disc.outset.pers.BI.msmf.prob=disc.outset.main.BI.msmf.prob,
                      disc.outset.pers.H.msmf.prob=disc.outset.main.H.msmf.prob,
                      disc.outset.pers.HI.msmf.prob=disc.outset.main.HI.msmf.prob,
                      disc.outset.pers.W.msmf.prob=disc.outset.main.W.msmf.prob,
                      
                      disc.at.diag.pers.B.f.prob=disc.outset.main.B.f.prob,
                      disc.at.diag.pers.BI.f.prob=disc.outset.main.BI.f.prob,
                      disc.at.diag.pers.H.f.prob=disc.outset.main.H.f.prob,
                      disc.at.diag.pers.HI.f.prob=disc.outset.main.HI.f.prob,
                      disc.at.diag.pers.W.f.prob=disc.outset.main.W.f.prob,
                      disc.at.diag.pers.B.msf.prob=disc.outset.main.B.msf.prob,
                      disc.at.diag.pers.BI.msf.prob=disc.outset.main.BI.msf.prob,
                      disc.at.diag.pers.H.msf.prob=disc.outset.main.H.msf.prob,
                      disc.at.diag.pers.HI.msf.prob=disc.outset.main.HI.msf.prob,
                      disc.at.diag.pers.W.msf.prob=disc.outset.main.W.msf.prob,
                      disc.at.diag.pers.B.msm.prob=disc.outset.main.B.msm.prob,
                      disc.at.diag.pers.BI.msm.prob=disc.outset.main.BI.msm.prob,
                      disc.at.diag.pers.H.msm.prob=disc.outset.main.H.msm.prob,
                      disc.at.diag.pers.HI.msm.prob=disc.outset.main.HI.msm.prob,
                      disc.at.diag.pers.W.msm.prob=disc.outset.main.W.msm.prob,
                      disc.at.diag.pers.B.msmf.prob=disc.outset.main.B.msmf.prob,
                      disc.at.diag.pers.BI.msmf.prob=disc.outset.main.BI.msmf.prob,
                      disc.at.diag.pers.H.msmf.prob=disc.outset.main.H.msmf.prob,
                      disc.at.diag.pers.HI.msmf.prob=disc.outset.main.HI.msmf.prob,
                      disc.at.diag.pers.W.msmf.prob=disc.outset.main.W.msmf.prob,
                      
                      disc.post.diag.pers.B.f.prob=disc.outset.main.B.f.prob,
                      disc.post.diag.pers.BI.f.prob=disc.outset.main.BI.f.prob,
                      disc.post.diag.pers.H.f.prob=disc.outset.main.H.f.prob,
                      disc.post.diag.pers.HI.f.prob=disc.outset.main.HI.f.prob,
                      disc.post.diag.pers.W.f.prob=disc.outset.main.W.f.prob,
                      disc.post.diag.pers.B.msf.prob=disc.outset.main.B.msf.prob,
                      disc.post.diag.pers.BI.msf.prob=disc.outset.main.BI.msf.prob,
                      disc.post.diag.pers.H.msf.prob=disc.outset.main.H.msf.prob,
                      disc.post.diag.pers.HI.msf.prob=disc.outset.main.HI.msf.prob,
                      disc.post.diag.pers.W.msf.prob=disc.outset.main.W.msf.prob,
                      disc.post.diag.pers.B.msm.prob=disc.outset.main.B.msm.prob,
                      disc.post.diag.pers.BI.msm.prob=disc.outset.main.BI.msm.prob,
                      disc.post.diag.pers.H.msm.prob=disc.outset.main.H.msm.prob,
                      disc.post.diag.pers.HI.msm.prob=disc.outset.main.HI.msm.prob,
                      disc.post.diag.pers.W.msm.prob=disc.outset.main.W.msm.prob,
                      disc.post.diag.pers.B.msmf.prob=disc.outset.main.B.msmf.prob,
                      disc.post.diag.pers.BI.msmf.prob=disc.outset.main.BI.msmf.prob,
                      disc.post.diag.pers.H.msmf.prob=disc.outset.main.H.msmf.prob,
                      disc.post.diag.pers.HI.msmf.prob=disc.outset.main.HI.msmf.prob,
                      disc.post.diag.pers.W.msmf.prob=disc.outset.main.W.msmf.prob,
                      
                      disc.inst.B.f.prob=disc.outset.main.B.f.prob,
                      disc.inst.BI.f.prob=disc.outset.main.BI.f.prob,
                      disc.inst.H.f.prob=disc.outset.main.H.f.prob,
                      disc.inst.HI.f.prob=disc.outset.main.HI.f.prob,
                      disc.inst.W.f.prob=disc.outset.main.W.f.prob,
                      disc.inst.B.msf.prob=disc.outset.main.B.msf.prob,
                      disc.inst.BI.msf.prob=disc.outset.main.BI.msf.prob,
                      disc.inst.H.msf.prob=disc.outset.main.H.msf.prob,
                      disc.inst.HI.msf.prob=disc.outset.main.HI.msf.prob,
                      disc.inst.W.msf.prob=disc.outset.main.W.msf.prob,
                      disc.inst.B.msm.prob=disc.outset.main.B.msm.prob,
                      disc.inst.BI.msm.prob=disc.outset.main.BI.msm.prob,
                      disc.inst.H.msm.prob=disc.outset.main.H.msm.prob,
                      disc.inst.HI.msm.prob=disc.outset.main.HI.msm.prob,
                      disc.inst.W.msm.prob=disc.outset.main.W.msm.prob,
                      disc.inst.B.msmf.prob=disc.outset.main.B.msmf.prob,
                      disc.inst.BI.msmf.prob=disc.outset.main.BI.msmf.prob,
                      disc.inst.H.msmf.prob=disc.outset.main.H.msmf.prob,
                      disc.inst.HI.msmf.prob=disc.outset.main.HI.msmf.prob,
                      disc.inst.W.msmf.prob=disc.outset.main.W.msmf.prob,
                      
                      circ.B.prob = circ.B.prob,
                      circ.BI.prob = circ.BI.prob,
                      circ.H.prob = circ.H.prob,
                      circ.HI.prob = circ.HI.prob,
                      circ.W.prob = circ.W.prob,
                      
                      ccr5.B.f.prob = c(ccr5.nonW, 0.0),
                      ccr5.BI.f.prob = c(ccr5.nonW, 0.0),
                      ccr5.H.f.prob = c(ccr5.nonW, 0.0),
                      ccr5.HI.f.prob = c(ccr5.nonW),
                      ccr5.W.f.prob = c(ccr5.W.f, 0.0),
                      ccr5.B.m.prob = c(ccr5.nonW, 0.0),
                      ccr5.BI.m.prob = c(ccr5.nonW),
                      ccr5.H.m.prob = c(ccr5.nonW, 0.0),
                      ccr5.HI.m.prob = c(ccr5.nonW, 0.0),
                      ccr5.W.m.prob = c(ccr5.W.m, 0.0),
                      ccr5.heteroz.rr = ccr5.heteroz.rr,

                      num.inst.ai.classes = 1,
                      
                      #not modeling msm but will include home testing SEattle estimates as placeholders
                      base.ai.main.B.rate = base.ai.main.B.rate,
                      base.ai.main.BI.rate = base.ai.main.BI.rate,
                      base.ai.main.H.rate = base.ai.main.H.rate,
                      base.ai.main.HI.rate = base.ai.main.HI.rate,
                      base.ai.main.W.rate = base.ai.main.W.rate,
                      base.ai.pers.B.rate = base.ai.pers.B.rate,
                      base.ai.pers.BI.rate = base.ai.pers.BI.rate,
                      base.ai.pers.H.rate = base.ai.pers.H.rate,
                      base.ai.pers.HI.rate = base.ai.pers.HI.rate,
                      base.ai.pers.W.rate = base.ai.pers.W.rate,
                      
                      
                      base.vi.main.B.rate = base.ai.pers.W.rate ,
                      base.vi.main.BI.rate = base.vi.main.BI.rate,
                      base.vi.main.H.rate = base.vi.main.H.rate,
                      base.vi.main.HI.rate = base.vi.main.HI.rate,
                      base.vi.main.W.rate = base.vi.main.W.rate,
                      base.vi.pers.B.rate =base.vi.pers.B.rate,
                      base.vi.pers.BI.rate = base.vi.pers.BI.rate,
                      base.vi.pers.H.rate = base.vi.pers.H.rate,
                      base.vi.pers.HI.rate = base.vi.pers.HI.rate,
                      base.vi.pers.W.rate = base.vi.pers.W.rate,
                      ai.scale = ai.scale,
                      vi.scale = vi.scale,

                      
                      # For MSM pulled data from the mobile study but we are not modeling MSM for now.
                      cond.main.B.prob.msm = cond.main.B.prob.msm,
                      cond.main.BI.prob.msm = cond.main.BI.prob.msm,
                      cond.main.H.prob.msm = cond.main.H.prob.msm,
                      cond.main.HI.prob.msm = cond.main.HI.prob.msm,
                      cond.main.W.prob.msm = cond.main.W.prob.msm,
                      cond.pers.always.prob.msm = cond.pers.always.prob.msm,
                      cond.pers.B.prob.msm = cond.pers.B.prob.msm,
                      cond.pers.BI.prob.msm = cond.pers.BI.prob.msm,
                      cond.pers.H.prob.msm = cond.pers.H.prob.msm,
                      cond.pers.HI.prob.msm = cond.pers.HI.prob.msm,
                      cond.pers.W.prob.msm = cond.pers.W.prob.msm,
                      cond.inst.always.prob.msm = cond.inst.always.prob.msm,
                      cond.inst.B.prob.msm = cond.inst.B.prob.msm,
                      cond.inst.BI.prob.msm = cond.inst.BI.prob.msm ,
                      cond.inst.H.prob.msm = cond.inst.H.prob.msm ,
                      cond.inst.HI.prob.msm = cond.inst.HI.prob.msm ,
                      cond.inst.W.prob.msm = cond.inst.W.prob.msm ,
                      cond.always.prob.corr.msm = cond.always.prob.corr.msm,
                      
                      cond.main.B.prob.het = cond.main.B.prob.het,
                      cond.main.BI.prob.het = cond.main.BI.prob.het,
                      cond.main.H.prob.het = cond.main.H.prob.het,
                      cond.main.HI.prob.het = cond.main.HI.prob.het,
                      cond.main.W.prob.het = cond.main.W.prob.het,
                      cond.pers.always.prob.het = cond.pers.always.prob.het,
                      cond.pers.B.prob.het = cond.pers.B.prob.het,
                      cond.pers.BI.prob.het = cond.pers.BI.prob.het,
                      cond.pers.H.prob.het = cond.pers.H.prob.het,
                      cond.pers.HI.prob.het = cond.pers.HI.prob.het,
                      cond.pers.W.prob.het = cond.pers.W.prob.het,
                      cond.inst.always.prob.het = cond.inst.always.prob.het,
                      cond.inst.B.prob.het = cond.inst.B.prob.het,
                      cond.inst.BI.prob.het = cond.inst.BI.prob.het,
                      cond.inst.H.prob.het = cond.inst.H.prob.het,
                      cond.inst.HI.prob.het = cond.inst.HI.prob.het,
                      cond.inst.W.prob.het = cond.inst.W.prob.het,
                      cond.always.prob.corr.het = cond.always.prob.corr.het,
                      
                      #WAITING ON APRIL AS POSIBLE DATA SOURCE.
                      #Currect placeholder from home testing project
                      cond.diag.main.beta.msm = -0.67,
                      cond.discl.main.beta.msm = -0.85,
                      cond.diag.pers.beta.msm = -0.67,
                      cond.discl.pers.beta.msm = -0.85,  
                      cond.diag.inst.beta.msm = -0.67,
                      cond.discl.inst.beta.msm = -0.85,
                      
                      cond.diag.main.beta.het = 0,
                      cond.discl.main.beta.het = 0,
                      cond.diag.pers.beta.het = 0,
                      cond.discl.pers.beta.het = 0,  
                      cond.diag.inst.beta.het = 0,
                      cond.discl.inst.beta.het = 0,
                      
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

                      vv.iev.B.prob = 0.47874,
                      vv.iev.BI.prob = 0.47874,
                      vv.iev.H.prob = 0.47874,
                      vv.iev.HI.prob = 0.47874,
                      vv.iev.W.prob = 0.47874,

                      ##No PrEP in the SHAMP model.
                      prep.start = Inf,
                      prep.elig.model = "base",
                      prep.class.prob = c(0, 0, 0, 0),
                      prep.class.hr = c(1, 0.69, 0.19, 0.05),
                      prep.coverage = 0,
                      prep.cov.method = "curr",
                      prep.cov.rate = 1,
                      prep.tst.int = 90,
                      prep.risk.int = 182,
                      prep.risk.reassess = TRUE,
                      
                      #max 59253
                      msm.foi.scale = 1,
                      
                      #max 992063
                      fa.BI.m.foi.scale = 1,
                      
                      #max 496031
                      fa.BI.f.foi.scale = 1,
                      
                      #max 10101010
                      fa.HI.m.foi.scale = 1,
                      
                      #max 5050505
                      fa.HI.f.foi.scale = 1,
                      
                      depart.scale = 1,
                      return.scale = 1,
                      immig.simple=TRUE,
                      
                      
                      immig.depart.BI.f = (prob.return.home/(agewindow*weeks.in.year)),
                      immig.depart.HI.f = (prob.return.home/(agewindow*weeks.in.year)),
                      immig.depart.BI.m = (prob.return.home/(agewindow*weeks.in.year)),
                      immig.depart.BI.m = (prob.return.home/(agewindow*weeks.in.year)),
                      immig.return.BI.f = immig.return.BI.f,
                      immig.return.HI.f = immig.return.HI.f,
                      immig.return.BI.m = immig.return.BI.m,
                      immig.return.HI.m = immig.return.HI.m,
                      immig.aq.prob.BI.f = prop.sex.home*prop.sex.home.BI*((base.vi.main.BI.rate+base.vi.pers.BI.rate)/2)*URVI.prob,
                      immig.aq.prob.HI.f = prop.sex.home*prop.sex.home.HI*((base.vi.main.HI.rate+base.vi.pers.HI.rate)/2)*URVI.prob,
                      immig.aq.prob.BI.m = prop.sex.home*prop.sex.home.BI*((base.vi.main.BI.rate+base.vi.pers.BI.rate)/2)*UIVI.prob,
                      immig.aq.prob.HI.m = prop.sex.home*prop.sex.home.HI*((base.vi.main.HI.rate+base.vi.pers.HI.rate)/2)*UIVI.prob,  
                      #pos.entry.BI = 0,
                      #pos.entry.HI = 0,
                      
                      msm.aq = (((1-condom.use)*prev.hiv.msm*(((transmission.receptive*1.09)+(transmission.insertive*1.09)) /2)*(contacts.per.week/weeks.in.year)) + condom.use*prev.hiv.msm*(((transmission.receptive*1.09)+(transmission.insertive*1.09))/2)*(contacts.per.week/weeks.in.year))/2,
                      
                      msm.aq.prob.B=msm.aq,
                      msm.aq.prob.BI=msm.aq,
                      msm.aq.prob.H=msm.aq,
                      msm.aq.prob.HI=msm.aq,
                      msm.aq.prob.W= msm.aq,
                      
                      ##
                      
                      demog.list=data.params[[1]]$demog.list,
                      demog.dist=data.params[[1]]$demog.dist,
                      sex.groups=data.params[[1]]$sex.groups,
                      race.groups=data.params[[1]]$race.groups,
                      age.groups=data.params[[1]]$age.groups,
                      age.adj=data.params[[1]]$age.adj,
                      conc_dur_dx = TRUE,
                      death_stats = TRUE,
                      add.demog.groups = FALSE,
                      
                      Ecohab.window = round(5*52),
                      p.growth = FALSE,
                      p.growth.nsteps = 0,
                      p.growth.size = 0,
                      
                      ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  #Method for HIV testing with options \code{"memoryless"}  for constant hazard without regard to time since previous test or 
  #\code{"interval"} deterministic fixed intervals.
  if (!(testing.pattern %in% c("memoryless", "interval"))) {
    stop("testing.pattern must be \"memoryless\" or \"interval\" ",
          call. = FALSE)
  }

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

  for(r in c("B", "BI", "H", "HI", "W")){
    
    Ps <-  c(rep(0, 17),
             1-(1-c(rep(asmr.1830.f, 12),
                    rep(asmr.3040.f, 10),
                    rep(asmr.4050.f, 6)))^(1 / age.unit),
             1)
    
  
    p[[paste0("asmr.", r ,".f")]] <- Ps
  }
  
  for(r in c("B", "BI", "H", "HI", "W")){
    Ps <- c(rep(0, 17),
            1-(1-c(rep(asmr.1830.m, 12),
                   rep(asmr.3040.m, 10),
                   rep(asmr.4050.m, 6)))^(1 / age.unit),
            1)
    p[[paste0("asmr.", r ,".m")]] <- Ps
  }

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
init_shamp <- function(prev.B.f = prev.B.f,
                     prev.BI.f =prev.BI.f,
                     prev.H.f =prev.H.f,
                     prev.HI.f =prev.HI.f,
                     prev.W.f = prev.W.f,
                     prev.B.msf = prev.B.msf,
                     prev.BI.msf =prev.BI.msf,
                     prev.H.msf =prev.H.msf,
                     prev.HI.msf =prev.HI.msf,
                     prev.W.msf = prev.W.msf,
                     prev.B.msm = prev.B.msm,
                     prev.BI.msm =prev.BI.msm,
                     prev.H.msm =prev.H.msm,
                     prev.HI.msm =prev.HI.msm,
                     prev.W.msm = prev.W.msm,
                     prev.B.msmf = prev.B.msmf,
                     prev.BI.msmf =prev.BI.msmf,
                     prev.H.msmf =prev.H.msmf,
                     prev.HI.msmf =prev.HI.msmf,
                     prev.W.msmf = prev.W.msmf,
                     ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  
  p$init.prev.age.slope.B.f <- prev.B.f / 28
  p$init.prev.age.slope.BI.f <- prev.BI.f / 28
  p$init.prev.age.slope.H.f <- prev.H.f / 28
  p$init.prev.age.slope.HI.f <- prev.HI.f / 28
  p$init.prev.age.slope.W.f <- prev.W.f / 28
  p$init.prev.age.slope.B.msf <- prev.B.msf / 28
  p$init.prev.age.slope.BI.msf <- prev.BI.msf / 28
  p$init.prev.age.slope.H.msf <- prev.H.msf / 28
  p$init.prev.age.slope.HI.msf <- prev.HI.msf / 28
  p$init.prev.age.slope.W.msf <- prev.W.msf / 28
  p$init.prev.age.slope.B.msm <- prev.B.msm / 28
  p$init.prev.age.slope.BI.msm <- prev.BI.msm / 28
  p$init.prev.age.slope.H.msm <- prev.H.msm / 28
  p$init.prev.age.slope.HI.msm <- prev.HI.msm / 28
  p$init.prev.age.slope.W.msm <- prev.W.msm / 28
  p$init.prev.age.slope.B.msmf <- prev.B.msmf / 28
  p$init.prev.age.slope.BI.msmf <- prev.BI.msmf / 28
  p$init.prev.age.slope.H.msmf <- prev.H.msmf / 28
  p$init.prev.age.slope.HI.msmf <- prev.HI.msmf / 28
  p$init.prev.age.slope.W.msmf <- prev.W.msmf / 28
 

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
#' @param heatbath.FUN Module function to stochastically simulate disease aquisition by msmf via 
#'        contact with an unsimulated msm population. 
#'        over acts given individual and dyadic attributes.
#' @param immigration.FUN Module function to stochastically simulate disease aquisition by migrants
#'        via conatct with unsimulated populations in home counties during return trips.        
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
#' A list object of class \code{control_shamp}, which can be passed to the
#' EpiModel function \code{netsim}.
#'
#' @keywords SHAMP 
#'
#' @export
control_shamp <- function(simno = 1,
                        nsims = 1,
                        ncores = 1,
                        nsteps = 100,
                        start = 1,
                        initialize.FUN = initialize_shamp,
                        aging.FUN = aging_shamp,
                        deaths.FUN = deaths_shamp,
                        births.FUN = births_shamp,
                        demogupdate.FUN = demogupdate_shamp,
                        test.FUN = test_shamp,
                        tx.FUN = tx_shamp,
                        #prep.FUN = prep_shamp,
                        progress.FUN = progress_msm,
                        vl.FUN = vl_shamp,
                        #aiclass.FUN = NULL,
                        #roleclass.FUN = NULL,
                        resim_nets.FUN = simnet_shamp,
                        disclose.FUN = disclose_shamp,
                        ConcDurDx.FUN = ConcDurDx_shamp,
                        acts.FUN = acts_shamp,
                        condoms.FUN = condoms_shamp,
                        riskhist.FUN = riskhist_shamp,
                        position.FUN = position_shamp,
                        trans.FUN = trans_shamp,
                        heatbath.FUN = heatbath_msmf_shamp,
                        immigration.FUN = immigration_shamp,
                        prev.FUN = prevalence_shamp,
                        verbose.FUN = verbose_shamp,
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

