
rm(list = ls())
suppressMessages(library("EpiModelHIV"))
devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV-p")

scr.dir <- "~/Dropbox/Dev/ARTnet/"
netstats <- readRDS(file.path(scr.dir, "data/artnet.NetStats.Atlanta.rda"))
epistats <- readRDS(file.path(scr.dir, "data/artnet.EpiStats.Atlanta.rda"))
est <- readRDS(file.path(scr.dir, "data/artnet.NetEst.Atlanta.rda"))

param <- param_msm(netstats = netstats,
                   epistats = epistats,
                   hiv.test.int = c(43, 43, 45),
                   a.rate = 0.00055,
                   riskh.start = 2,
                   prep.start = 30,
                   prep.start.prob = 0.10,
                   tt.part.supp = c(0.20, 0.20, 0.20),
                   tt.full.supp = c(0.40, 0.40, 0.40),
                   tt.dur.supp = c(0.40, 0.40, 0.40),
                   tx.halt.full.rr = 0.8,
                   tx.halt.dur.rr = 0.1,
                   tx.reinit.full.rr = 2.0,
                   tx.reinit.dur.rr = 5.0,
                   hiv.rgc.rr = 2.5,
                   hiv.ugc.rr = 1.5,
                   hiv.rct.rr = 2.5,
                   hiv.uct.rr = 1.5,
                   hiv.dual.rr = 0.0,
                   rgc.tprob = 0.35,
                   ugc.tprob = 0.25,
                   rct.tprob = 0.20,
                   uct.tprob = 0.16,
                   rgc.ntx.int = 16.8,
                   ugc.ntx.int = 16.8,
                   rct.ntx.int = 32,
                   uct.ntx.int = 32,
                   acts.aids.vl = 5.75)
init <- init_msm()
control <- control_msm(simno = 1,
                       nsteps = 52 * 2,
                       nsims = 1,
                       ncores = 1,
                       save.nwstats = FALSE,
                       save.clin.hist = FALSE)

sim <- netsim(est, param, init, control)

# Explore clinical history
par(mar = c(3,3,1,1), mgp = c(2,1,0))
m1 <- sim$temp[[1]]$clin.hist[[1]]
m2 <- sim$temp[[1]]$clin.hist[[2]]
m3 <- sim$temp[[1]]$clin.hist[[3]]
a <- sim$attr[[1]]
h <- which(a$status == 1)

m1[h[1:10], 95:104]
aids <- which(a$stage == 4)
id <- h[58]
plot(m1[id, ], type = "o", ylim = c(1, 7))
data.frame(vl = m1[id, ], stage = m2[id, ], tx = m3[id, ])
a$tt.traj[id]
matplot(t(m1[h[1:500], ]), type = "l", lty = 1, ylim = c(1, 7))


df <- as.data.frame(sim)
names(df)

par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = "i.prev", mean.smooth = FALSE, ylim = c(0, 1))
plot(sim, y = "num")
plot(sim, y = "dep.gen", mean.smooth = TRUE)
plot(sim, y = "dep.AIDS", mean.smooth = FALSE)
plot(sim, y = "prepCurr")
plot(sim, y = "cc.dx", mean.smooth = FALSE)
plot(sim, y = "cc.linked", mean.smooth = FALSE, ylim = c(0.8, 1))
plot(sim, y = "cc.linked1m", mean.smooth = FALSE)
plot(sim, y = "cc.tx", mean.smooth = FALSE)
plot(sim, y = "cc.tx.any1y", mean.smooth = FALSE)
plot(sim, y = "cc.vsupp", mean.smooth = FALSE)
plot(sim, y = "cc.vsupp.tt1", mean.smooth = FALSE)
plot(sim, y = "cc.vsupp.tt2", mean.smooth = FALSE)
plot(sim, y = "cc.vsupp.tt3", mean.smooth = FALSE)
plot(sim, y = "cc.vsupp.dur1y", mean.smooth = FALSE)

plot(sim, y = "hstage.acute", mean.smooth = TRUE)
plot(sim, y = "hstage.chronic", mean.smooth = FALSE)
plot(sim, y = "hstage.aids", mean.smooth = FALSE)

plot(sim, y = "ir100.gc", mean.smooth = FALSE)
plot(sim, y = "ir100.ct", mean.smooth = FALSE)
plot(sim, y = "ir100.sti", mean.smooth = FALSE)
plot(sim, y = "prev.gc", mean.smooth = FALSE)
plot(sim, y = "prev.ct", mean.smooth = FALSE)

plot(sim, type = "formation", network = 1, plots.joined = FALSE)
plot(sim, type = "formation", network = 2, plots.joined = FALSE)
plot(sim, type = "formation", network = 3, plots.joined = FALSE)


# Testing/Timing ------------------------------------------------------

m <- microbenchmark::microbenchmark(hivvl_msm(dat, at))
print(m, unit = "ms")

dat <- initialize_msm(est, param, init, control, s = 1)

for (at in 2:200) {
  dat <- aging_msm(dat, at)
  dat <- departure_msm(dat, at)
  dat <- arrival_msm(dat, at)
  dat <- hivtest_msm(dat, at)
  dat <- hivtx_msm(dat, at)
  dat <- hivprogress_msm(dat, at)
  dat <- hivvl_msm(dat, at)
  dat <- simnet_msm(dat, at)
  dat <- acts_msm(dat, at)
  dat <- condoms_msm(dat, at)
  dat <- position_msm(dat, at)
  dat <- prep_msm(dat, at)
  dat <- hivtrans_msm(dat, at)
  dat <- stitrans_msm(dat, at)
  dat <- stirecov_msm(dat, at)
  dat <- stitx_msm(dat, at)
  dat <- prevalence_msm(dat, at)
  verbose.net(dat, "progress", at = at)
}

nrow(dat$temp$plist)
table(dat$temp$plist[, "start"])
table(dat$temp$plist[, "stop"])
head(dat$temp$plist)

plist <- as.data.frame(dat$temp$plist)
pmain <- filter(plist, ptype == 2)
table(pmain$start)
hist(pmain$start)
