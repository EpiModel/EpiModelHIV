
# stiPrEP burnin analysis

rm(list=ls())

list.files("data/")
unlink("data/sim.*")
system("scp hyak:/gscratch/csde/sjenness/sti/data/sim.n300.rda scripts/burnin/data/")
load("scripts/burnin/data/sim.n300.rda")

ls()


i.prev <- as.numeric(sim$epi$i.prev[2600, ])
round(quantile(i.prev, probs = c(0.025, 0.5, 0.975)), 3)

prev.rgcct <- as.numeric(sim$epi$prev.rgcct[2600, ])
round(quantile(prev.rgcct, probs = c(0.025, 0.5, 0.975)), 3)

prev.ugcct <- as.numeric(sim$epi$prev.ugcct[2600, ])
round(quantile(prev.ugcct, probs = c(0.025, 0.5, 0.975)), 3)

ir100.gc <- as.numeric(sim$epi$ir100.gc[2600, ])
round(quantile(ir100.gc, probs = c(0.025, 0.5, 0.975)), 3)

ir100.ct <- as.numeric(sim$epi$ir100.ct[2600, ])
round(quantile(ir100.ct, probs = c(0.025, 0.5, 0.975)), 3)

ir100 <- as.numeric(sim$epi$ir100[2600, ])
round(quantile(ir100, probs = c(0.025, 0.5, 0.975)), 3)

# Summary of 500 sims
par(mfrow = c(3,2), oma=c(0,0,2,0))
plot(sim, y = "i.prev", ylim = c(0.23, 0.29))
title("HIV Prevalence")
plot(sim, y = "prev.rgcct", ylim = c(0.12, 0.17))
title("Rectal STI Prevalence")
plot(sim, y = "prev.ugcct", ylim = c(0.06, 0.11))
title("Urethral STI Prevalence")
plot(sim, y = "ir100.gc", ylim = c(18, 28))
title("Gonococcal Incidence Rate")
plot(sim, y = "ir100.ct", ylim = c(20, 30))
title("Chlamydial Incidence Rate")
plot(sim, y = "ir100", ylim = c(3, 5))
title("HIV Incidence Rate")
title(main="Statistics for 500 Simulations", outer=TRUE)


# Summary of mean simulation
# Run after burnin.process
par(mfrow = c(3,2), oma=c(0,0,2,0))
plot(sim, y = "i.prev", ylim = c(0.23, 0.29))
title("HIV Prevalence")
plot(sim, y = "prev.rgcct", ylim = c(0.12, 0.17))
title("Rectal STI Prevalence")
plot(sim, y = "prev.ugcct", ylim = c(0.06, 0.11))
title("Urethral STI Prevalence")
plot(sim, y = "ir100.gc", ylim = c(18, 28))
title("Gonococcal Incidence Rate")
plot(sim, y = "ir100.ct", ylim = c(20, 30))
title("Chlamydial Incidence Rate")
plot(sim, y = "ir100", ylim = c(3, 5))
title("HIV Incidence Rate")
title(main="Statistics for Best-fitting Simulation", outer=TRUE)

df <- as.data.frame(sim)
head(df$prev.rgcct, 100)

plot(sim, y = "prev.rgcct", mean.smooth = FALSE)

df <- tail(as.data.frame(sim), 500)
sum(df$trans.main) / sum(df$incid)
sum(df$trans.casl) / sum(df$incid)
sum(df$trans.inst) / sum(df$incid)

sum(df$trans.main.gc) / sum(df$incid.gc)
sum(df$trans.casl.gc) / sum(df$incid.gc)
sum(df$trans.inst.gc) / sum(df$incid.gc)

sum(df$trans.main.ct) / sum(df$incid.ct)
sum(df$trans.casl.ct) / sum(df$incid.ct)
sum(df$trans.inst.ct) / sum(df$incid.ct)
