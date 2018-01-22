##Running network simulations with EPIModelHIV and ergm.ego.

#Load required packages.
library(EpiModelHIV)
library(EpiModelHPC)
#devtools::install_github("statnet/ergm.ego-private",ref="3.7-compat", auth_token ="")
#devtools::install_github("statnet/tergmLite")

library(latticeExtra)
library(ergm.ego)
library(tergmLite)
library(parallel)
library(rlist)
np = detectCores()

#Load the ego and alter data frames and look at them
load("~/SHAMP/egonet/data/all.impute.egodata.rda")
str(all.impute.egodata)

data.params<-list ()

new_data<-input_shamp(all.impute.egodata, data.params, immigration=TRUE, msm.msmf=FALSE)

##Make smaller data set
data.params<-as.list(new_data[1])
object.size(new_data)
new_data[2] <- make_small(new_data[2])
object.size(new_data)

##Make the three ergm.ego objects.
ego.obj_c<-as.egodata(new_data[[2]][[1]]$egos,alters=new_data[[2]][[1]]$altersCohab,egoIDcol="ego", egoWt=new_data[[2]][[1]]$egos$weight)
ego.obj_p<-as.egodata(new_data[[2]][[1]]$egos,alters=new_data[[2]][[1]]$altersPers,egoIDcol="ego", egoWt=new_data[[2]][[1]]$egos$weight)
ego.obj_i<-as.egodata(new_data[[2]][[1]]$egos,alters=new_data[[2]][[1]]$altersOT,egoIDcol="ego", egoWt=new_data[[2]][[1]]$egos$weight)

# Location of EpiModelHIV folders
  epifold <- '~/EpiModelHIV_SHAMP'
# epifold <- '~/'
# epifold <- '~/Dropbox/CFAR/'

save(ego.obj_c, file = file.path(epifold, "EpiModelHIV_shamp_modeling/scenarios/base/data/ego.obj_c.rda"))
save(ego.obj_p, file = file.path(epifold, "EpiModelHIV_shamp_modeling/scenarios/base/data/ego.obj_p.rda"))
save(ego.obj_i, file = file.path(epifold, "EpiModelHIV_shamp_modeling/scenarios/base/data/ego.obj_i.rda"))

save(ego.obj_c, file = file.path(epifold, "EpiModelHIV/data/ego.obj_c.rda"))
save(ego.obj_p, file = file.path(epifold, "EpiModelHIV/data/ego.obj_p.rda"))
save(ego.obj_i, file = file.path(epifold, "EpiModelHIV/data/ego.obj_i.rda"))


###cohab partnership network.

summary(ego.obj_c ~edges + 
          nodemix("sex",base=0))

summary(ego.obj_c ~edges + 
                  nodefactor("race.sex",base=0) +
                  nodematch("race",diff=TRUE) +
                  nodefactor("agecat", base=1) + 
                  absdiff("sqrt.age.adj") + 
                  nodefactor("deg.pers.c",base=1))

##~bd(maxout=1) still produced 4% with 2 ties and a handful with more
##Try offset concurrent() ??    
#nodefactor("race",base=5) +
  
fit.c.temp <- ergm.ego(ego.obj_c ~ edges + 
                nodefactor("race",base=5) + 
                nodefactor("agecat", base=1) + 
                offset(nodematch("sex", diff=FALSE)), 
                offset.coef = -Inf,
                control=control.ergm.ego(ppopsize=20000, 
                                         stats.est="asymptotic",
                                         ergm.control = 
                                           control.ergm(MCMC.interval=7500,
                                                        MCMC.samplesize=7500,
                                                        MCMC.burnin = 7500,
                                                        MPLE.max.dyad.types = 1e7,
                                                        init.method = "zeros",
                                                        MCMLE.maxit = 400,
                                                        parallel = np, 
                                                        parallel.type="PSOCK")))



fit.c <- fit.c.temp
object.size(fit.c)

fit.c$egodata <- fit.c$network <- fit.c$newnetworks <- fit.c$sample <- fit.c$constrained <-NULL 
object.size(fit.c)
summary(fit.c)

####Casual partnership network.
##For concurrent and m.deg use 7 catagory base (BI males for deg cohab).
##Black male and Black female for concurrent.
##will need to create male and female comparison groups that are all other race groups.

##Create a new dat object that is a time ordered list of 
##transmission events with all seeds coming from UID 0.
##All MSM come from 
summary(ego.obj_p ~edges + 
          nodefactor("race") + 
          nodematch("race",diff=TRUE) +
          nodefactor("agecat", base=0) + 
          absdiff("sqrt.age.adj") + 
          nodefactor("cross.net.group", base=0) +
          concurrent(by="pers.conc.group"))

fit.p.temp<-ergm.ego(ego.obj_p ~edges + 
                    nodefactor("race3",base=3) + 
                    nodefactor("agecat", base=1) + 
                    offset(nodematch("sex", diff=FALSE)),
                  offset.coef = -Inf,
                  constraints=~bd(maxout=3),
                  control=control.ergm.ego(ppopsize=20000, stats.est="asymptotic",
                                           ergm.control = control.ergm(MCMC.interval=7500,
                                                                       MCMC.samplesize=7500,
                                                                       MCMC.burnin = 7500,
                                                                       MPLE.max.dyad.types = 1e7,
                                                                       init.method = "zeros",
                                                                       MCMLE.maxit = 400,
                                                                       parallel = np, 
                                                                       parallel.type="PSOCK")))
##Creates a sigularity?
# nodefactor("cross.net.group",base=4) +

fit.p <- fit.p.temp
object.size(fit.p)
fit.p$egodata <- fit.p$newnetworks  <- fit.p$network <- fit.p$sample <- fit.p$constrained <-NULL 
object.size(fit.p)

summary(fit.p)



######One time partnerships.

summary(ego.obj_i ~edges + 
          nodefactor("race", base=0) +
          nodematch("race",diff=TRUE) +
          nodefactor("agecat", base=0) + 
          nodefactor("deg.cohab.c",base=0) +
          nodefactor("deg.pers.c",base=0))

fit.i.temp<-ergm.ego(ego.obj_i ~edges + 
                    nodefactor("agecat", base=1) + 
                    offset(nodematch("sex", diff=FALSE)),
                  offset.coef = -Inf,
                  control=control.ergm.ego(ppopsize=20000, stats.est="asymptotic",
                                           ergm.control = control.ergm(MCMC.interval=7000,
                                                                       MCMC.samplesize=7000,
                                                                       MCMC.burnin = 7000,
                                                                       MPLE.max.dyad.types = 1e4,
                                                                       init.method = "MPLE",
                                                                       MCMLE.maxit = 350,
                                                                       parallel = np, 
                                                                       parallel.type="PSOCK")))

fit.i <- fit.i.temp
object.size(fit.i)

fit.i$egodata <- fit.i$newnetworks <- fit.i$network <- fit.i$sample <- fit.i$constrained <-NULL 
object.size(fit.i)

summary(fit.i)


fullmodelfits <- list(fit.c.temp, fit.p.temp, fit.i.temp)
modelfits <- list(fit.c, fit.p, fit.i)
save(fullmodelfits, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/fullmodelfits.rda")
save(modelfits, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/modelfits.rda")

load(file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/fullmodelfits.rda")
load(file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/modelfits.rda")
fit.c <- list.extract(modelfits,1)
fit.p <- list.extract(modelfits,2)
fit.i <- list.extract(modelfits,3)

  
#Create additional required elements for est
# Time unit for simulation, relative to 1 day
time.unit <- 7
method<-1

#Simulation size.
sim.size<-19149


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


formation_c <- fit.c$formula
coef.form.crude_c <- fit.c$coef[2:length(fit.c$coef)]
coef.form_c <- coef.form.crude_c
coef.form_c[1]<-coef.form_c[1]+fit.c$coef[1]
coef.form.crude_c <- coef.form_c
coef.form_c[1]<- coef.form_c[1] - coef.diss_c$coef.adj
constraints_c <- ~.

formation_p <- fit.p$formula
coef.form.crude_p <- fit.p$coef[2:length(fit.p$coef)]
coef.form_p <- coef.form.crude_p
coef.form_p[1] <- coef.form_p[1]+fit.p$coef[1]
coef.form.crude_p <- coef.form_p
coef.form_p[1]<- coef.form_p[1]- coef.diss_p$coef.adj
constraints_p <- ~.

formation_i <- fit.i$formula
coef.form.crude_i <- fit.i$coef[2:length(fit.i$coef)]
coef.form_i <- coef.form.crude_i
coef.form_i[1] <- coef.form.crude_i[1]+fit.i$coef[1]
coef.form.crude_i <- coef.form_i
coef.form_i[1] <- coef.form_i[1]+ -log(52)
constraints_i <- ~.



nw<-network.initialize(sim.size, directed = FALSE, hyper = FALSE, loops = FALSE,
                       multiple = FALSE, bipartite = FALSE)

fit.1<-fit.2<-fit.3<-list()
fit.1$fit <- fit.c
fit.1$formation <- formation_c
#fit.1$fit$coef <- coef.form.crude_c
fit.1$target.stats <- target.stats_c
fit.1$target.stats.names <- target.stats.names_c 
fit.1$coef.form <- coef.form_c
fit.1$coef.form.crude <- coef.form.crude_c
#fit.1$fit$ergm.formula <- fit.c$formula
#fit.1$fit$ergm.offset.coef <- fit.1$fit$offset.coef
#fit.1$fit$target.stats<-fit.1$fit$target.stats[2:length(fit.1$fit$target.stats)]
#fit.1$fit$MCMCtheta[2]<-fit.1$fit$MCMCtheta[1] + fit.1$fit$MCMCtheta[2]
#fit.1$fit$MCMCtheta <- fit.1$fit$MCMCtheta[2:length(fit.1$fit$MCMCtheta)]
fit.1$coef.diss <- coef.diss_c 
fit.1$constraints <- constraints_c
fit.1$edapprox <- TRUE



fit.2$fit <- fit.p
#fit.2$fit$coef <- coef.form.crude_p
fit.2$formation <- formation_p
fit.2$target.stats <- target.stats_p
fit.2$target.stats.names <- target.stats.names_p 
fit.2$coef.form <- coef.form_p
fit.2$coef.form.crude <- coef.form.crude_p 
#fit.2$fit$ergm.formula <- fit.p$formula
#fit.2$fit$ergm.offset.coef <- fit.2$fit$offset.coef
#fit.2$fit$target.stats<-fit.2$fit$target.stats[2:length(fit.2$fit$target.stats)]
#fit.2$fit$MCMCtheta[2]<-fit.2$fit$MCMCtheta[1] + fit.2$fit$MCMCtheta[2]
#fit.2$fit$MCMCtheta <- fit.2$fit$MCMCtheta[2:length(fit.2$fit$MCMCtheta)]
fit.2$coef.diss <- coef.diss_p 
fit.2$constraints <- constraints_p
fit.2$edapprox <- TRUE

fit.3$fit <- fit.c
fit.3$formation <- formation_i
#fit.3$fit$coef <- coef.form.crude_i
fit.3$target.stats <- target.stats_i
fit.3$target.stats.names <- target.stats.names_i 
fit.3$coef.form <- coef.form_i
fit.3$coef.form.crude <- coef.form.crude_i
#fit.3$fit$ergm.formula <- fit.i$formula
#fit.3$fit$ergm.offset.coef <- fit.3$fit$offset.coef
#fit.3$fit$target.stats<-fit.3$fit$target.stats[2:length(fit.3$fit$target.stats)]
#fit.3$fit$MCMCtheta[2]<-fit.3$fit$MCMCtheta[1] + fit.3$fit$MCMCtheta[2]
#fit.3$fit$MCMCtheta <- fit.3$fit$MCMCtheta[2:length(fit.3$fit$MCMCtheta)]
fit.3$coef.diss <- coef.diss_i 
fit.3$constraints <- constraints_i
fit.3$edapprox <- TRUE



est <- list(fit.1, fit.2, fit.3)





save(est, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/fit.rda")
save(data.params, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/data.params.rda")





load(file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/fit.rda")
load(file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/data.params.rda")


param <- param_shamp(data.params)

init <- init_shamp()
control <- control_shamp(nsteps = 1040, save.other = c("attr", "trans.el"), verbose = TRUE)

sim1 <- netsim(est, param, init, control)
save(sim1, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/sim1.rda")

sim2 <- netsim(est, param, init, control)
save(sim2, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/sim2.rda")

sim3 <- netsim(est, param, init, control)
save(sim3, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/sim3.rda")




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

save(sim1, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/sim100_500_20_5.rda")

time<-1:1040
plot(time,sim1$epi$i.prev[,1] , main = "Prev: temp.adjust - 10, depart.adjust - 5")

