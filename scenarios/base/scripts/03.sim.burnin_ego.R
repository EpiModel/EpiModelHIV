
library("EpiModelHPC")
library("EpiModelHIV")

rm(list = ls())

load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/fitsmall.rda")
load(file = "~/EpiModelHIV_shamp_modeling/scenarios/est/data.params.rda")

param <- param_shamp()
init <- init_shamp()
control <- control_shamp(nsteps = 5, nsims = 1, prevfull = TRUE, save.transmat = TRUE,
                        save.other = "attr", save.network = TRUE,
                        verbose = TRUE, verbose.int = 1)

netsim(est, param, init, control)

sim <- netsim(est, param, init, control)

names(sim.counts$sim$epi)

save(sim, file = "~/EpimodelHIV_shamp_modeling/scenario/est/sim.rda")




