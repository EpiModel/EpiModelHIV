
source("scripts/burnin/abcfx.R")

# rect.prev, ureth.prev, gc.incid, ct.incid, hiv.prev
targets <- c(0.17, 0.07, 43, 48, 0.26)


# Batch 1 -------------------------------------------------------------

## PRE ##

# done with abcr

## POST ##

system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataAll.b1.rda scripts/burnin/")
load("scripts/burnin/simDataAll.b1.rda")
dim(sim)

cbind(sapply(sim, function(x) length(unique(x))))

plotStats(sim, targets)
plotParam(sim)

simChosen <- rejection(sim, targets = targets, threshold = 0.01)

plotStats(simChosen, targets)
plotParam(simChosen)

# calculate standard deviations
sds <- apply(simChosen, 2, sd)

# set weight
simChosen$weight <- 1/nrow(simChosen)

save(simChosen, sds, file = "scripts/burnin/simChosen.b1.rda")


# Batch 2 -------------------------------------------------------------

## PRE ##
init.samp.size = 5e5
lastParticleSamp <- simChosen[sample(x = 1:nrow(simChosen),
                                      size = init.samp.size,
                                      prob = simChosen$weight,
                                      replace = TRUE), ]

# restrict new sample to initial priors
batch.size <- 10000
simParms <- perturbParticle(parms = lastParticleSamp,
                            pnames = names(simChosen[1:14]),
                            from = NULL,
                            sds = sds[1:14])
simParms <- simParms[sample(1:nrow(simParms), size = batch.size), ]
row.names(simParms) <- NULL
dim(simParms)
sapply(simParms, function(x) length(unique(x)))
head(simParms)

save(simParms, file = "scripts/burnin/simParms.b2.rda")
system("scp scripts/burnin/simParms.b2.rda hyak:/gscratch/csde/sjenness/sti/")
system("scp scripts/burnin/*.abcsmc.[Rs]* hyak:/gscratch/csde/sjenness/sti/")


## POST ##

system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataAll.b2.rda scripts/burnin/")
load("scripts/burnin/simDataAll.b2.rda")
dim(sim)

cbind(sapply(sim, function(x) length(unique(x))))

plotStats(sim, targets)
plotParam(sim)

simChosen <- rejection(sim, targets = targets, threshold = 0.05)

plotStats(simChosen, targets)
plotParam(simChosen)

sdsNew <- apply(simChosen, 2, sd)

# this is the batch 2 parameter set
simChosenNew <- simChosen

# loading the batch 1 parameter set
load("scripts/burnin/simChosen.b1.rda")
simChosenOld <- simChosen
sdsOld <- sds

# Weighting the new parameter set
simChosenNew$weight <- weightParticles(pnames = names(simChosen[1:14]),
                                       currentBatch = simChosenNew,
                                       lastBatch = simChosenOld,
                                       sdsUse = sdsOld[1:14])
simChosenNew$weight[simChosenNew$weight == Inf] <- 0
simChosenNew$weight <- simChosenNew$weight/sum(simChosenNew$weight)

# save the new simChosen with weights and sds
simChosen <- simChosenNew
sds <- sdsNew

save(simChosen, sds, file = "scripts/burnin/simChosen.b2.rda")


# Batch 3 -------------------------------------------------------------

rm(list=ls())
source("scripts/burnin/abcfx.R")

## PRE ##

load("scripts/burnin/simChosen.b2.rda")

init.samp.size = 5e5
lastParticleSamp <- simChosen[sample(x = 1:nrow(simChosen),
                                     size = init.samp.size,
                                     prob = simChosen$weight,
                                     replace = TRUE), ]

# restrict new sample to initial priors
batch.size <- 10000
simParms <- perturbParticle(parms = lastParticleSamp,
                            pnames = names(simChosen[1:14]),
                            from = NULL,
                            sds = sds[1:14])
simParms <- simParms[sample(1:nrow(simParms), size = batch.size), ]
row.names(simParms) <- NULL
dim(simParms)
sapply(simParms, function(x) length(unique(x)))
head(simParms)

save(simParms, file = "scripts/burnin/simParms.b3.rda")
system("scp scripts/burnin/simParms.b3.rda hyak:/gscratch/csde/sjenness/sti/")
system("scp scripts/burnin/*.abcsmc.[Rs]* hyak:/gscratch/csde/sjenness/sti/")


## POST ##

system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataAll.b3.rda scripts/burnin/")
load("scripts/burnin/simDataAll.b3.rda")
dim(sim)

cbind(sapply(sim, function(x) length(unique(x))))

plotStats(sim, targets)
plotParam(sim)

simChosen <- rejection(sim, targets = targets, threshold = 0.05)

plotStats(simChosen, targets)
plotParam(simChosen)

sdsNew <- apply(simChosen, 2, sd)

# this is the batch 2 parameter set
simChosenNew <- simChosen

# loading the batch 2 parameter set
load("scripts/burnin/simChosen.b2.rda")
simChosenOld <- simChosen
sdsOld <- sds

# Weighting the new parameter set
simChosenNew$weight <- weightParticles(pnames = names(simChosen[1:14]),
                                       currentBatch = simChosenNew,
                                       lastBatch = simChosenOld,
                                       sdsUse = sdsOld[1:14])
simChosenNew$weight[simChosenNew$weight == Inf] <- 0
simChosenNew$weight <- simChosenNew$weight/sum(simChosenNew$weight)

# save the new simChosen with weights and sds
simChosen <- simChosenNew
sds <- sdsNew

save(simChosen, sds, file = "scripts/burnin/simChosen.b3.rda")


# Batch 4 -------------------------------------------------------------

rm(list=ls())
source("scripts/burnin/abcfx.R")

## PRE ##

load("scripts/burnin/simChosen.b3.rda")

init.samp.size = 5e5
lastParticleSamp <- simChosen[sample(x = 1:nrow(simChosen),
                                     size = init.samp.size,
                                     prob = simChosen$weight,
                                     replace = TRUE), ]

# restrict new sample to initial priors
batch.size <- 10000
simParms <- perturbParticle(parms = lastParticleSamp,
                            pnames = names(simChosen[1:14]),
                            from = NULL,
                            sds = sds[1:14])
simParms <- simParms[sample(1:nrow(simParms), size = batch.size), ]
row.names(simParms) <- NULL
dim(simParms)
sapply(simParms, function(x) length(unique(x)))
head(simParms)

save(simParms, file = "scripts/burnin/simParms.b4.rda")
system("scp scripts/burnin/simParms.b4.rda hyak:/gscratch/csde/sjenness/sti/")
system("scp scripts/burnin/*.abcsmc.[Rs]* hyak:/gscratch/csde/sjenness/sti/")


## POST ##

system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataAll.b4.rda scripts/burnin/")
load("scripts/burnin/simDataAll.b4.rda")
dim(sim)

cbind(sapply(sim, function(x) length(unique(x))))
sim <- sim[sample(10000, 10000), ]

plotStats(sim, targets)
plotParam(sim)

simChosen <- rejection(sim, targets = targets, threshold = 0.01)

plotStats(simChosen, targets)
plotParam(simChosen)

sdsNew <- apply(simChosen, 2, sd)

# this is the batch 2 parameter set
simChosenNew <- simChosen

# loading the batch 2 parameter set
load("scripts/burnin/simChosen.b3.rda")
simChosenOld <- simChosen
sdsOld <- sds

# Weighting the new parameter set
simChosenNew$weight <- weightParticles(pnames = names(simChosen[1:14]),
                                       currentBatch = simChosenNew,
                                       lastBatch = simChosenOld,
                                       sdsUse = sdsOld[1:14])
simChosenNew$weight[simChosenNew$weight == Inf] <- 0
simChosenNew$weight <- simChosenNew$weight/sum(simChosenNew$weight)

# save the new simChosen with weights and sds
simChosen <- simChosenNew
sds <- sdsNew

save(simChosen, sds, file = "scripts/burnin/simChosen.b4.rda")


# Batch 5 -------------------------------------------------------------

rm(list=ls())
source("scripts/burnin/abcfx.R")

## PRE ##

simParms <- abc_pre(batch = 5, batch.size = 10000)

save(simParms, file = "scripts/burnin/simParms.b5.rda")
system("scp scripts/burnin/simParms.b5.rda hyak:/gscratch/csde/sjenness/sti/")
system("scp scripts/burnin/*.abcsmc.[Rs]* hyak:/gscratch/csde/sjenness/sti/")


## POST ##

system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataAll.b5.rda scripts/burnin/")
load("scripts/burnin/simDataAll.b5.rda")
dim(sim)

cbind(sapply(sim, function(x) length(unique(x))))

plotStats(sim, targets)
plotParam(sim)

simChosen <- rejection(sim, targets = targets, threshold = 0.01)

plotStats(simChosen, targets)
plotParam(simChosen)

sdsNew <- apply(simChosen, 2, sd)

# this is the batch 2 parameter set
simChosenNew <- simChosen

# loading the batch 2 parameter set
load("scripts/burnin/simChosen.b4.rda")
simChosenOld <- simChosen
sdsOld <- sds

# Weighting the new parameter set
simChosenNew$weight <- weightParticles(pnames = names(simChosen[1:14]),
                                       currentBatch = simChosenNew,
                                       lastBatch = simChosenOld,
                                       sdsUse = sdsOld[1:14])
simChosenNew$weight[simChosenNew$weight == Inf] <- 0
simChosenNew$weight <- simChosenNew$weight/sum(simChosenNew$weight)

# save the new simChosen with weights and sds
simChosen <- simChosenNew
sds <- sdsNew

save(simChosen, sds, file = "scripts/burnin/simChosen.b5.rda")


# Batch 6 -------------------------------------------------------------

rm(list=ls())
source("scripts/burnin/abcfx.R")

## PRE ##

simParms <- abc_pre(batch = 6, batch.size = 25000)

save(simParms, file = "scripts/burnin/simParms.b6.rda")
system("scp scripts/burnin/simParms.b6.rda hyak:/gscratch/csde/sjenness/sti/")
system("scp scripts/burnin/*.abcsmc.[Rs]* hyak:/gscratch/csde/sjenness/sti/")


## POST ##

system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataAll.b6.rda scripts/burnin/")
load("scripts/burnin/simDataAll.b6.rda")
dim(sim)

cbind(sapply(sim, function(x) length(unique(x))))

plotStats(sim, targets)
plotParam(sim)

up.good <- which(abs(sim$ureth.prev - 0.07) < 0.01)
sim.up.good <- sim[up.good, ]

plotStats(sim.up.good, targets)

simChosen <- rejection(sim, targets = targets, threshold = 0.1)

plotStats(simChosen, targets)
plotParam(simChosen)

simChosen <- abc_post(simChosen, batch = 6)

save(simChosen, sds, file = "scripts/burnin/simChosen.b6.rda")



