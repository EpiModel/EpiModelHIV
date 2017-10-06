##Running network simulations with EPIModelHIV and ergm.ego.

#Load required packages.
library(EpiModelHIV)
library(EpiModelHPC)
#devtools::install_github("statnet/ergm.ego-private",ref="3.7-compat", auth_token ="")
#devtools::install_github("statnet/tergmLite")

library(ergm.ego)
library(tergmLite)
library(parallel)
np = detectCores()

#Load the ego and alter data frames and look at them
load("~/SHAMP/egonet/data/all.egodata.rda")
str(all.egodata)


data.params<-list ()


new_data<-input_shamp(all.egodata, data.params, immigration=TRUE, msm.msmf=FALSE)
data.params<-as.list(new_data[1])


##Make the three ergm.ego objects.
ego.obj_c<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersCohab,egoIDcol="ego", egoWt=new_data[[2]]$egos$weight)
ego.obj_p<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersPers,egoIDcol="ego", egoWt=new_data[[2]]$egos$weight)
ego.obj_i<-as.egodata(new_data[[2]]$egos,alters=new_data[[2]]$altersOT,egoIDcol="ego", egoWt=new_data[[2]]$egos$weight)


save(ego.obj_c, file = "~/EpiModelHIV_SHAMP2/data/ego.obj_c.rda")
save(ego.obj_p, file = "~/EpiModelHIV_SHAMP2/data/ego.obj_p.rda")
save(ego.obj_i, file = "~/EpiModelHIV_SHAMP2/data/ego.obj_i.rda")



###cohab partnership network.

summary(ego.obj_c ~edges + 
                  nodefactor("race.sex",base=0) +
                  nodematch("race",diff=TRUE) +
                  nodefactor("agecat", base=1) + 
                  absdiff("sqrt.age.adj") + 
                  nodefactor("deg.pers.c",base=1))

##~bd(maxout=1) still produced 4% with 2 ties and a handful with more
##Try offset concurrent() ??              
fit.c<-ergm.ego(ego.obj_c ~edges + 
                    nodefactor("race",base=5) +
                    nodematch("race",diff=TRUE) +
                    nodefactor("agecat", base=1) + 
                    absdiff("sqrt.age.adj") + 
                    nodefactor("deg.pers.c",base=1) + 
                    offset(nodematch("sex", diff=FALSE)) + 
                    offset(concurrent()),
                    offset.coef = c(-Inf, -Inf),
                    constraints=~bd(maxout=1),
                    control=control.ergm.ego(ppopsize=50000, stats.est="asymptotic",
                                           ergm.control = control.ergm(MCMC.interval=7500,
                                                                       MCMC.samplesize=7500,
                                                                       MCMC.burnin = 7500,
                                                                       MPLE.max.dyad.types = 1e7,
                                                                       init.method = "zeros",
                                                                       MCMLE.maxit = 400,
                                                                       parallel = np, 
                                                                       parallel.type="PSOCK")))


summary(fit.c)
test <- simulate(fit.c)
degreedist(test)

####Casual partnership network.

summary(ego.obj_p ~edges + 
          nodefactor("race",base=5) + 
          nodematch("race",diff=TRUE) +
          nodefactor("agecat", base=1) + 
          absdiff("sqrt.age.adj") + 
          nodefactor("deg.cohab.c",base=1) +
          concurrent(by="sex") + 
          concurrent(by="race"))

fit.p<-ergm.ego(ego.obj_p ~edges + 
                    nodefactor("race",base=5) + 
                    nodematch("race",diff=TRUE) +
                    nodefactor("agecat", base=1) + 
                    absdiff("sqrt.age.adj") + 
                    nodefactor("deg.cohab.c",base=1) +
                    concurrent(by="sex") + 
                    concurrent(by="race") +
                    offset(nodematch("sex", diff=FALSE)),
                  offset.coef = -Inf,
                  constraints=~bd(maxout=3),
                  control=control.ergm.ego(ppopsize=50000, stats.est="asymptotic",
                                           ergm.control = control.ergm(MCMC.interval=7000,
                                                                       MCMC.samplesize=7000,
                                                                       MCMC.burnin = 7000,
                                                                       MPLE.max.dyad.types = 1e7,
                                                                       init.method = "zeros",
                                                                       MCMLE.maxit = 350,
                                                                       parallel = np, 
                                                                       parallel.type="PSOCK")))

summary(fit.p)
test <- simulate(fit.p)
degreedist(test)


######One time partnerships.

summary(ego.obj_i ~edges + 
          nodefactor("race", base=0) +
          nodematch("race",diff=TRUE) +
          nodefactor("agecat", base=0) + 
          nodefactor("deg.cohab.c",base=0) +
          nodefactor("deg.pers.c",base=0))

fit.i<-ergm.ego(ego.obj_i ~edges + 
                    nodefactor("race", base=5) +
                    nodematch("race",diff=TRUE) +
                    nodefactor("agecat", base=c(3,4)) + 
                    nodefactor("deg.cohab.c",base=1) +
                    nodefactor("deg.pers.c",base=1) +
                    offset(nodematch("sex", diff=FALSE)),
                  offset.coef = -Inf,
                  verbose=TRUE,
                  control=control.ergm.ego(ppopsize=50000, stats.est="asymptotic",
                                           ergm.control = control.ergm(SAN.maxit=50,
                                                                       SAN.burnin.times=100,
                                                                       MCMC.interval=7000,
                                                                       MCMC.samplesize=7000,
                                                                       MCMC.burnin = 7000,
                                                                       MPLE.max.dyad.types = 1e4,
                                                                       init.method = "MPLE",
                                                                       MCMLE.maxit = 350,
                                                                       parallel = np, 
                                                                       parallel.type="PSOCK")))

summary(fit.i)



modelfits <- list(fit.c, fit.p, fit.i)
save(modelfits, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/modelfits.rda")

load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/modelfits.rda")


#Create additional required elements for est
# Time unit for simulation, relative to 1 day
time.unit <- 7
method<-1

#Simulation size.
sim.size<-47733


# Mean durations
diss_c = ~offset(edges) 

diss_p = ~offset(edges)

#Mortality
ages <- 18:59
age.unit <- 52

asmr.B.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 10),
                   rep(0.003065837, 10)))^(1 / age.unit),
            1)
asmr.BI.f <- c(rep(0, 17),
             1-(1-c(rep(0.000405376, 12),
                    rep(0.000661066, 10),
                    rep(0.001378053, 10),
                    rep(0.003065837, 10)))^(1 / age.unit),
             1)

asmr.H.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 10),
                   rep(0.003065837, 10)))^(1/age.unit),
            1)
asmr.HI.f <- c(rep(0, 17),
             1-(1-c(rep(0.000405376, 12),
                    rep(0.000661066, 10),
                    rep(0.001378053, 10),
                    rep(0.003065837, 10)))^(1/age.unit),
             1)
asmr.W.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 10),
                   rep(0.003065837, 10)))^(1/age.unit),
            1)


asmr.B.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 10),
                     rep(0.005400669, 10)))^(1 / age.unit),
              1)
asmr.BI.m <- c(rep(0, 17),
               1-(1-c(rep(0.000853417, 12),
                      rep(0.001084014, 10),
                      rep(0.001982864, 10),
                      rep(0.005400669, 10)))^(1 / age.unit),
               1)

asmr.H.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 10),
                     rep(0.005400669, 10)))^(1/age.unit),
              1)
asmr.HI.m <- c(rep(0, 17),
               1-(1-c(rep(0.000853417, 12),
                      rep(0.001084014, 10),
                      rep(0.001982864, 10),
                      rep(0.005400669, 10)))^(1/age.unit),
               1)
asmr.W.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 10),
                     rep(0.005400669, 10)))^(1/age.unit),
              1)




exp.mort <- (mean(asmr.B.f[ages]) + mean(asmr.BI.f[ages]) + mean(asmr.H.f[ages])
             + mean(asmr.HI.f[ages]) + mean(asmr.W.f[ages]) + mean(asmr.B.m[ages]) 
             + mean(asmr.BI.m[ages]) + mean(asmr.H.m[ages])
             + mean(asmr.HI.m[ages]) + mean(asmr.W.m[ages]) ) / 10

coef.diss_c <- dissolution_coefs(dissolution = diss_c,
                               duration = data.params[[1]]$durs_c / time.unit,
                               d.rate = exp.mort)

coef.diss_p <- dissolution_coefs(dissolution = diss_p,
                               duration = data.params[[1]]$durs_p / time.unit,
                               d.rate = exp.mort)


coef.diss_i <-dissolution_coefs(~offset(edges), 1)

target.stats_c<-as.numeric(fit.c$target.stats[2:length(fit.c$target.stats)])
target.stats_p<-as.numeric(fit.p$target.stats[2:length(fit.p$target.stats)])
target.stats_i<-as.numeric(fit.i$target.stats[2:length(fit.i$target.stats)])

target.stats.names_c<- names(fit.c$target.stats[2:length(fit.c$target.stats)])
target.stats.names_p<- names(fit.p$target.stats[2:length(fit.p$target.stats)])
target.stats.names_i<- names(fit.i$target.stats[2:length(fit.i$target.stats)])


coef.form.crude_c<-fit.c$coef[2:length(fit.c$coef)]
coef.form_c<-coef.form.crude_c
coef.form_c[1]<-coef.form_c[1]+fit.c$coef[1]
coef.form_c[1]<- coef.form_c[1]- coef.diss_c$coef.adj
constraints_c <- ~.

coef.form.crude_p<-fit.p$coef[2:length(fit.p$coef)]
coef.form_p<-coef.form.crude_p
coef.form_p[1]<-coef.form_p[1]+fit.p$coef[1]
coef.form_p[1]<- coef.form_p[1]- coef.diss_p$coef.adj
constraints_p <- ~.

coef.form.crude_i<-fit.i$coef[2:length(fit.i$coef)]
coef.form_i<-coef.form.crude_i
coef.form_i[1]<-coef.form_i[1]+fit.i$coef[1]
coef.form_i[1]<- coef.form_i[1]
coef.form_i[1]<-coef.form_i[1]+ -log(52)
constraints_i <- ~.



nw<-network.initialize(sim.size, directed = FALSE, hyper = FALSE, loops = FALSE,
                       multiple = FALSE, bipartite = FALSE)

fit.1 <- list(fit= fit.c, formation=fit.c$formula, target.stats= target.stats_c,
            target.stats.names= target.stats.names_c, coef.form = coef.form_c,
            coef.form.crude= coef.form.crude_c, coef.diss=coef.diss_c, constraints= constraints_c,
            edapprox=TRUE)

fit.2 <- list(fit= fit.p, formation=fit.p$formula, target.stats= target.stats_p,
              target.stats.names= target.stats.names_p, coef.form = coef.form_p,
              coef.form.crude= coef.form.crude_p, coef.diss=coef.diss_p, constraints= constraints_p,
              edapprox=TRUE)

fit.3 <- list(fit= fit.i, formation=fit.i$formula, target.stats= target.stats_i,
              target.stats.names= target.stats.names_i, coef.form = coef.form_i,
              coef.form.crude= coef.form.crude_i, coef.diss=coef.diss_i, constraints= constraints_i,
              edapprox=TRUE)


est <- list(fit.1, fit.2, fit.3)
#save(est, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fit.rda")
save(est, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fitsmall.rda")
save(data.params, file = "~/EpiModelHIV_shamp_modeling/scenarios/est/data.params.rda")

param <- param_shamp(data.params,temp.adjust=100)
init <- init_shamp()
control <- control_shamp(nsteps = 520)




#load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fit.rda")
load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fitsmall.rda")
load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/data.params.rda")


sim<-netsim(est, param, init, control)




netsim(est, param, init, control)

###################################################################################
########testing#########

test <- simulate(est[[2]]$fit)

degreedist(test)

summary(test ~edges + degree
          nodefactor("race",base=5) + 
          nodematch("race",diff=TRUE) +
          nodefactor("agecat", base=1) + 
          absdiff("sqrt.age.adj") + 
          nodefactor("deg.cohab.c",base=1))


param <- param_shamp(data.params, msm.temp.adjust = 100, fa.temp.adjust = 500, depart.adjust = 20, return.adjust = .5)
init <- init_shamp()
control <- control_shamp(nsteps = 1040)
sim1<-netsim(est, param, init, control)

save(sim1, file = "~/EpiModelHIV_shamp_modeling/scenarios/sim100_500_20_5.rda")

time<-1:1040
plot(time,sim1$epi$i.prev[,1] , main = "Prev: temp.adjust - 10, depart.adjust - 5")

