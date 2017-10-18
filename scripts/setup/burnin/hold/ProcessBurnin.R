
## Process burn-in
library("EpiModelHPC")

# Examine output
# scp hyak:/gscratch/csde/camp/data/*.rda data/
sim <- merge_simfiles(1000, indir = "data/", ftype = "min")
plot(sim, y = "i.prev", ylim = c(0.2, 0.3), qnts = 0.5)
abline(h = 0.26)

df <- as.data.frame(sim)
round(mean(tail(df$i.prev, 100)), 3)

# Save burn-in file for FU sims
sim <- merge_simfiles(1000, indir = "data/", ftype = "max")
sim <- get_sims(sim, sims = "mean", var = "i.prev")
tail(as.data.frame(sim)$i.prev)

save(sim, file = "est/p2.burnin.rda")
# scp est/p2.burnin.rda hyak:/gscratch/csde/camp/est/
