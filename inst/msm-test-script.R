
rm(list = ls())
suppressMessages(library("EpiModelHIV"))
devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV-p")

scr.dir <- "~/Dropbox/Projects/SexualDistancing/"
netstats <- readRDS(file.path(scr.dir, "out/est/netstats.rds"))
epistats <- readRDS(file.path(scr.dir, "out/est/epistats.rds"))
est <- readRDS(file.path(scr.dir, "out/est/netest.rds"))

full_tx_eff <- rep(1, 3)

param <- param_msm(
  netstats = netstats,
  epistats = epistats,
  hiv.test.rate = c(0.00385, 0.00385, 0.0069),
  hiv.test.late.prob = rep(0, 3),
  tx.init.prob = c(0.1775, 0.19, 0.2521),
  tt.part.supp = 1 - full_tx_eff,
  tt.full.supp = full_tx_eff,
  tt.dur.supp = rep(0, 3),
  tx.halt.part.prob = c(0.0065, 0.0053, 0.003),
  tx.halt.full.rr = rep(0.45, 3),
  tx.halt.dur.rr = rep(0.45, 3),
  tx.reinit.part.prob = rep(0.00255, 3),
  tx.reinit.full.rr = rep(1, 3),
  tx.reinit.dur.rr = rep(1, 3),
  max.time.off.tx.full.int = 52 * 15,
  max.time.on.tx.part.int = 52 * 10,
  max.time.off.tx.part.int = 52 * 10,
  aids.mr = 1 / 250,
  trans.scale =  c(2.7, 0.35, 0.243), #c(2.21, 0.405, 0.255),
  acts.scale.main = 1.00,
  acts.scale.casl = 0.10,
  acts.aids.vl = 5.75,
  circ.prob = c(0.874, 0.874, 0.918),
  a.rate = 0.00052,
  prep.start = (52 * 60) + 1,
  riskh.start = 52 * 59,
  prep.adhr.dist = c(0.089, 0.127, 0.784),
  prep.adhr.hr = c(0.69, 0.19, 0.01),
  prep.start.prob =  0.71, # 0.00896,
  prep.discont.rate = 0.02138792, # 1 - (2^(-1/(224.4237/7)))
  ## prep.tst.int = 90/7,         # do I need that?
  ## prep.risk.int = 182/7,       # do I need that?
  ## prep.sti.screen.int = 182/7,
  ## prep.sti.prob.tx = 1,
  prep.risk.reassess.method = "year",
  prep.require.lnt = TRUE, # FALSE -> start with random PrEP initiation

  ## STI PARAMS (default: from combprev2, make it gaps)
  ## Using values in prep-race: scripts/burnin/sim.burn.R
  ## If not mentionned -> default from prep disparities
  ## for H : mean(c(B, W))
  #ok
  rgc.tprob = 0.2267303, #logistic(logit(0.19) + log(1.25)) #0.357,  # gaps appendix 9.4
  ugc.tprob = 0.19,# 0.248,  # gaps appendix 9.4
  rct.tprob = 0.2038369, #logistic(logit(0.17) + log(1.25)) #0.3216, # gaps appendix 9.3
  uct.tprob = 0.17,#0.213,  # gaps appendix 9.3
  rgc.sympt.prob = 0.1,#0.077, # gaps appendix 10.3
  ugc.sympt.prob = 0.9333333,#0.824, # gaps appendix 10.3
  rct.sympt.prob = 0.1,#0.1035,# gaps appendix 10.2
  uct.sympt.prob = 0.95,#0.885, # gaps appendix 10.2
  rgc.ntx.int = 26,#35.11851, # gaps appendix 11.2
  ugc.ntx.int = 26,#35.11851, # gaps appendix 11.2
  gc.tx.int   = 2, # gaps appendix 11.2 - mentionned, not explicit
  rct.ntx.int = 32,#44.24538, # gaps appendix 11.1
  uct.ntx.int = 32,#44.24538, # gaps appendix 11.1
  ct.tx.int   = 2, # gaps appendix 11.1 - mentionned, not explicit

  gc.sympt.prob.tx =  rep(0.9, 3),  #c(0.86, 0.91, 0.96),
  ct.sympt.prob.tx =  rep(0.9, 3),  #c(0.72, 0.785, 0.85),
  gc.asympt.prob.tx = rep(0.1, 3), #c(0.10, 0.145, 0.19),
  ct.asympt.prob.tx = rep(0.1, 3), #c(0.05, 0.525, 0.10),
  # gaps appendix 9.3 - 9.4 (not explained this way but similar result)
  sti.cond.eff = 0.95,
  sti.cond.fail = c(0.39, 0.3, 0.21),
  # gaps appendix 9.2
  hiv.rgc.rr = 2.78,
  hiv.ugc.rr = 1.73,
  hiv.rct.rr = 2.78,
  hiv.uct.rr = 1.73,
  # if both ct + gc -> log(RRgc) + 0.2 * log(RRct) | swap ct and gc if RRct > RRgc
  hiv.dual.rr = 0.2, # not mentionned in appendix

  netresim.form.rr = rep(1, 3),
  netresim.disl.rr = rep(1, 2),

)

init <- init_msm(
  prev.ugc = 0.05,
  prev.rct = 0.05,
  prev.rgc = 0.05,
  prev.uct = 0.05
)
control <- control_msm(simno = 1,
                       nsteps = 52 * 2,
                       nsims = 1,
                       ncores = 1,
                       save.nwstats = FALSE,
                       save.clin.hist = FALSE)

debug(acts_msm)
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
