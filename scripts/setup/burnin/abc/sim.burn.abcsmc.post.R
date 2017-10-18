
## abc model analysis

library("EpiModelHIV")
library("EasyABC")

system("scp scripts/burnin/abc/*.abcsmc.[Rs]* hyak:/gscratch/csde/sjenness/sti2")
system("scp scripts/burnin/abc/*.abcsmc2.[Rs]* hyak:/gscratch/csde/sjenness/sti2")
system("scp scripts/burnin/abc/*.abcsmc3.[Rs]* hyak:/gscratch/csde/sjenness/sti2")
system("scp scripts/burnin/abc/*.abcsmc4.[Rs]* hyak:/gscratch/csde/sjenness/stia")

system("scp hyak:/gscratch/csde/sjenness/stia/data/*.rda scripts/burnin/abc/")


## averaged ATL/demo fits
load("scripts/burnin/abc/smc.5pct.100sim.rda")
load("scripts/burnin/abc/smc2.5pct.100sim.rda")
load("scripts/burnin/abc/smc3.10pct.100sim.rda")
load("scripts/burnin/abc/smc3.5pct.100sim.rda")
load("scripts/burnin/abc/smc3.1pct.100sim.rda")


p <- as.data.frame(a$param)
s <- as.data.frame(a$stats)
w <- a$weights

# names(p) <- c("rgc.tprob", "ugc.tprob", "rct.tprob", "uct.tprob",
#               "rgc.sympt.prob", "ugc.sympt.prob", "rct.sympt.prob", "uct.sympt.prob",
#               "rgc.dur.asympt", "ugc.dur.asympt", "rct.dur.asympt", "uct.dur.asympt",
#               "hiv.rect.rr", "hiv.ureth.rr")

names(p) <- c("gc.dur.asympt", "ct.dur.asympt", "hiv.rect.rr", "hiv.ureth.rr")

# for averaged fits
# names(s) <- c("rect.prev", "ureth.prev", "gc.incid", "ct.incid", "hiv.incid", "hiv.prev")

# for AIDS meta
# names(s) <- c("gc.incid", "ct.incid", "hiv.incid", "hiv.prev")

names(s) <- c("gc.incid", "ct.incid", "hiv.prev")

( mean.s <- apply(s, 2, function(x) sum(x * w)) )
( mean.p <- apply(p, 2, function(x) sum(x * w)) )


tar.avg <- c(0.135, 0.046, 23.2, 26.8, 3.8, 0.26)

data.frame(mean.s, tar.atl)

mean.p

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(3,2))
for (i in 1:ncol(s)) {
  hist(s[, i], col = "bisque2", border = "white", main = names(s)[i])
  abline(v = tar[i], lwd = 2, col = "red")
}

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(4,4))
for (i in 1:ncol(p)) {
  hist(p[, i], col = "bisque2", border = "white", main = names(p)[i])
}

save(mean.p, file = "scripts/burnin/abc/abc.avg.parms.1pct.rda")
save(mean.p, file = "est/meta.parms.rda")

for (i in seq_along(mean.p)) {
  assign(names(mean.p)[i], unname(mean.p[i]))
}

