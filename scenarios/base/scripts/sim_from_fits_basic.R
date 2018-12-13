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


load(file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/data.params.rda")
fit.c<- readRDS("/net/proj/SHAMPnetdat/model-fits/current/cohab-main-A.rds")
fit.p<- readRDS("/net/proj/SHAMPnetdat/model-fits/current/pers-main-A.rds")
fit.i<- readRDS("/net/proj/SHAMPnetdat/model-fits/current/ot-main-A.rds")


fit.c$egodata <- fit.c$newnetworks <- fit.c$network <- fit.c$constrained <-NULL 
fit.p$egodata <- fit.p$newnetworks <- fit.p$network <- fit.p$constrained <-NULL 
fit.i$egodata <- fit.i$newnetworks <- fit.i$network <- fit.i$constrained <-NULL 


#Create additional required elements for est
# Time unit for simulation, relative to 1 day
time.unit <- 7
method<-1

#Simulation size.
sim.size<-fit.c$ppopsize


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



est_50K_base <- list(fit.1, fit.2, fit.3)

save(est_50K_base, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/fit_50K_base.rda")


load(file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/fit_50K_base.rda")
load(file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/data.params.rda")





param <- param_shamp(data.params)

init <- init_shamp()
control <- control_shamp(nsteps = 104, save.other = c("attr", "trans.el"), verbose = TRUE)

sim_50K_base <- netsim(est_50K_base, param, init, control)
save(sim_50K_base, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/sim_50K_base.rda")



