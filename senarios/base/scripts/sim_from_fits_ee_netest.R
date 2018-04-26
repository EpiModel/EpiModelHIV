##Running network simulations with EPIModelHIV and ergm.ego.

#Load required packages.
library(EpiModelHIV)
library(EpiModelHPC)
#devtools::install_github("statnet/ergm.ego-private",ref="3.7-compat", auth_token ="")
#devtools::install_github("statnet/tergmLite")

library(latticeExtra)
library(ergm.ego)
library(tergmLite)
library(rlist)
library(egonet)



load("~/SHAMP/egonet/data/nsfg.impute.egodata.rda")
str(nsfg.egodata)
data.params<-list ()
new_data<-input_shamp(nsfg.impute.egodata, data.params, immigration=TRUE, msm.msmf=FALSE)
data.params<-as.list(new_data[1])
save(data.params, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/data.params.rda")


fit.c<- readRDS("/net/proj/SHAMPnetdat/model-fits/current/cohab.rds")
fit.p<- readRDS("/net/proj/SHAMPnetdat/model-fits/current/pers.rds")
fit.i<- readRDS("/net/proj/SHAMPnetdat/model-fits/current/ot.rds")




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
ages <- 18:45
age.unit <- 52

asmr.B.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 6)))^(1 / age.unit),
            1)
asmr.BI.f <- c(rep(0, 17),
             1-(1-c(rep(0.000405376, 12),
                    rep(0.000661066, 10),
                    rep(0.001378053, 6)))^(1 / age.unit),
             1)

asmr.H.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 6)))^(1/age.unit),
            1)
asmr.HI.f <- c(rep(0, 17),
             1-(1-c(rep(0.000405376, 12),
                    rep(0.000661066, 10),
                    rep(0.001378053, 6)))^(1/age.unit),
             1)
asmr.W.f <- c(rep(0, 17),
            1-(1-c(rep(0.000405376, 12),
                   rep(0.000661066, 10),
                   rep(0.001378053, 6)))^(1/age.unit),
            1)


asmr.B.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 6)))^(1 / age.unit),
              1)
asmr.BI.m <- c(rep(0, 17),
               1-(1-c(rep(0.000853417, 12),
                      rep(0.001084014, 10),
                      rep(0.001982864, 6)))^(1 / age.unit),
               1)

asmr.H.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 6)))^(1/age.unit),
              1)
asmr.HI.m <- c(rep(0, 17),
               1-(1-c(rep(0.000853417, 12),
                      rep(0.001084014, 10),
                      rep(0.001982864, 6)))^(1/age.unit),
               1)
asmr.W.m <- c(rep(0, 17),
              1-(1-c(rep(0.000853417, 12),
                     rep(0.001084014, 10),
                     rep(0.001982864, 6)))^(1/age.unit),
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



##Use ee.netest to create netest objects
##Delete newnetwroks for space


fit.1<-ee.netest(fit.c, coef.diss_c)
fit.1$fit$newnetworks <- fit.1$fit$network <- fit.1$fit$sample <- NULL
fit.2<-ee.netest(fit.p, coef.diss_p)
fit.2$fit$newnetworks <- fit.2$fit$network <- fit.2$fit$sample <- NULL
fit.3<-ee.netest(fit.i, coef.diss_i)
fit.3$fit$newnetworks <- fit.3$fit$network <- fit.3$fit$sample <- NULL
est_f <- list(fit.1, fit.2, fit.3)
fit.c<-fit.p<-fit.i<-NULL

save(est_f, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/fit_f.rda")
save(data.params, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/data.params.rda")



load(file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/fit_f.rda")
load(file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/data.params.rda")



#nw<-simulate(sim.size, directed = FALSE, hyper = FALSE, loops = FALSE,
#                       multiple = FALSE, bipartite = FALSE)
nw<-fit.1$fit$newnetwork
count <- network.edgecount(nw)
delete.edges(nw,1:count)

param <- param_shamp(data.params, 
                     VI.foi.scale = 3,
                     msm.foi.scale = 1000,
                     fa.foi.scale = 1000)

init <- init_shamp()
control <- control_shamp(nsteps = 500, save.other = c("attr", "trans.el", "death.stats","cel.temp"), verbose = TRUE)

sim <- netsim(est_f, param, init, control)

 
save(sim, file = "~/EpiModelHIV_SHAMP/EpiModelHIV_shamp_modeling/scenarios/base/est/sim.rda")




x<-as.data.frame<-sim$death.stats[[1]]
x<-as.data.frame(x)




time<-1:520
plot(time,sim_fits$epi$i.prev[,1] , main = "Baseline")



pdf(file = "Baseline.pdf", height = 12, width = 9, pointsize = 16)

par(mfrow = c(2,1), mar = c(3,3,2.5,5), mgp = c(2,1,0))

plot(time,sim_fits$epi$i.prev[,1], main = "HIV Prevalence", xlab = "Weeks", ylab = "Prevalence among adult heterosexuals 18-59")
plot(time,sim_fits$epi$prop.Lhet.inf[,1], main = "Proportion of infections by class", ylim = c(0, 1), xlab = "Weeks", ylab = "Proportion of infections", col="black")
lines(time,sim_fits$epi$prop.MSM.inf[,1], col="red")
lines(time,sim_fits$epi$prop.MSMds.inf[,1], col="pink")
lines(time,sim_fits$epi$prop.FA.inf[,1], col="blue")
lines(time,sim_fits$epi$prop.FAds.inf[,1], col="purple")



legend("right", c("Local Heterosexual - Black","Directly aquired from the MSM associated force of infection - Red",
                  "Downstream MSM associated infections - Pink",
                  "Directly Aquired from the foreign associated force of infection - Blue",
                  "Downstream foreign associated infection - Purple" ), cex=.8)


dev.off()