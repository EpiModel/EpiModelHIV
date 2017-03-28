suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

i.prev <- as.numeric(mean(sim$epi$i.prev[2500:2600, ]))
ir100 <- as.numeric(mean(sim$epi$ir100[2500:2600, ]))
prev.hiv.syphpos <- as.numeric(mean(sim$epi$prev.hiv.syphpos[2500:2600, ]))
prev.hiv.primsecosyphpos <- as.numeric(mean(sim$epi$prev.hiv.primsecosyphpos[2500:2600, ]))

ir100.syph <- as.numeric(mean(sim$epi$ir100.syph[2500:2600, ]))
prev.syph <- as.numeric(mean(sim$epi$prev.syph[2500:2600, ]))
prev.syph.hivpos <- as.numeric(mean(sim$epi$prev.syph.hivpos[2500:2600, ]))
prev.syph.hivneg <- as.numeric(mean(sim$epi$prev.syph.hivneg[2500:2600, ]))
prev.primsecosyph.hivpos <- as.numeric(mean(sim$epi$prev.syph.hivpos[2500:2600, ]))
prev.primsecosyph.hivneg <- as.numeric(mean(sim$epi$prev.syph.hivneg[2500:2600, ]))
prev.stage.incubprim <- as.numeric(mean(sim$epi$prev.stage.incubprim[2500:2600, ]))
prev.stage.seco <- as.numeric(mean(sim$epi$prev.stage.seco[2500:2600, ]))
prev.stage.earlat <- as.numeric(mean(sim$epi$prev.stage.earlat[2500:2600, ]))
prev.stage.latelat <- as.numeric(mean(sim$epi$prev.stage.latelat[2500:2600, ]))
prev.stage.latelatelat <- as.numeric(mean(sim$epi$prev.stage.latelatelat[2500:2600, ]))
prev.stage.alllatelat <- as.numeric(mean(sim$epi$prev.stage.alllatelat[2500:2600, ]))
prev.stage.tert <- as.numeric(mean(sim$epi$prev.stage.tert[2500:2600, ]))
prev.earlysyph <- as.numeric(mean(sim$epi$prev.earlysyph[2500:2600, ]))
prev.latesyph <- as.numeric(mean(sim$epi$prev.latesyph[2500:2600, ]))

ir100.gc <- as.numeric(mean(sim$epi$ir100.gc[2500:2600, ]))
ir100.ct <- as.numeric(mean(sim$epi$ir100.ct[2600:2600, ]))

# Syphilis ratio
par(mfrow = c(1, 1))
boxplot(sim$epi$prev.syph.hivpos / sim$epi$prev.syph.hivneg)
abline(h = 3.96, lty = 2, col = "red")

# Syphilis prevalence plots
par(mfrow = c(1, 1), oma = c(0,0,2,0))
plot(sim, y = "prev.primsecosyph.hivpos", ylab = "Prevalence", ylim = c(0, 0.15), mean.col = "blue")
plot(sim, y = "prev.primsecosyph.hivneg", ylab = "Prevalence", add = TRUE, mean.col = "green")
plot(sim, y = "prev.primsecosyph", ylab = "Prevalence", add = TRUE, mean.col = "red")
abline(h = 0.103, col = "blue", lty = 2)
abline(h = 0.026, col = "green", lty = 2)
abline(h = 0.046, col = "red", lty = 2)
title("Syphilis by HIV Status")
legend("topleft", c("HIV+", "HIV-", "Overall"), col = c("blue", "green", "red"), lty = c(1, 1, 1))

plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0.00, 0.40), mean.col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, mean.col = "green")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", add = TRUE, mean.col = "purple")
abline(h = 0.20, col = "blue", lty = 2)
abline(h = 0.077, col = "green", lty = 2)
abline(h = 0.2770, col = "purple", lty = 2)
title("P, S, and Early latent Syphilis Prevalence")
legend("topright", c("Primary", "Seco", "Early latent"), col = c("blue", "green", "purple"), lty = c(1, 1, 1))

plot(sim, y = "prev.stage.alllatelat", ylab = "Prevalence", ylim = c(0.00, 0.60), mean.col = "blue")
plot(sim, y = "prev.stage.tert", ylab = "Prevalence", add = TRUE, mean.col = "green")
abline(h = 0.44, col = "blue", lty = 2)
abline(h = 0.006, col = "green", lty = 2)
title("Late Latent and Tertiary Syphilis Prevalence")
legend("right", c("Late latent", "Tertiary"), col = c("blue", "green"), lty = c(1, 1))

# All syphilis stages
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0, 1), mean.col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, mean.col = "green")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", add = TRUE, mean.col = "orange")
plot(sim, y = "prev.stage.alllatelat", ylab = "Prevalence", add = TRUE, mean.col = "purple")
plot(sim, y = "prev.stage.tert", ylab = "Prevalence", add = TRUE, mean.col = "black")
title("Stage-specific Syphilis Prevalence")
abline(h = 0.20, col = "blue", lty = 1)
abline(h = 0.077, col = "green", lty = 1)
abline(h = 0.2770, col = "orange", lty = 1)
abline(h = 0.44, col = "purple", lty = 1)
abline(h = 0.006, col = "black", lty = 1)
legend("topright", c("Primary/Incub", "Secondary", "Early Latent", "Late latent", "Tertiary"),
       col = c("blue", "green", "orange", "purple", "black"), lty = c(1, 1, 1, 1, 1))

a <- cbind(c(prev.primsecosyph.hivpos, prev.syph, prev.primsecosyph.hivneg, prev.stage.incubprim,
             prev.stage.seco, prev.stage.earlat, prev.stage.alllatelat,
             prev.stage.tert, prev.earlysyph, prev.latesyph))
rownames(a) <- c( "prev.syph.hivpos", "prev.syph", "prev.syph.hivneg", "prev.stage.incubprim",
                  "prev.stage.seco", "prev.stage.earlat", "prev.stage.alllatelat",
                  "prev.stage.tert", "prev.earlysyph", "prev.latesyph")
a

# Early vs late syphilis
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "prev.earlysyph", ylab = "Prevalence", ylim = c(0, 1), mean.col = "blue")
plot(sim, y = "prev.latesyph", ylab = "Prevalence", add = TRUE, mean.col = "green")
title("Early and Late-stage Syphilis")
abline(h = 0.554, col = "blue", lty = 1)
abline(h = 0.446, col = "green", lty = 1)
legend("topright", c("Early", "Late"), col = c("blue", "green"), lty = c(1, 1))

# STI Prevalence
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "prev.primsecosyph", ylab = "Prevalence")
abline(h = 0.046, col = "red", lty = 2)
abline(v = sim$param$stitest.start)
title("P and S Syphilis Prevalence")

plot(sim, y = "prev.ct", ylab = "Prevalence")
abline(v = sim$param$stitest.start)
abline(h = 0.040, col = "red", lty = 2)
title("CT Prevalence")

plot(sim, y = "prev.gc", ylab = "Prevalence")
abline(v = sim$param$stitest.start)
abline(h = 0.040, col = "red", lty = 2)
title("GC Prevalence")

plot(sim, y = "i.prev", ylim = c(0, 0.3), ylab = "Prevalence")
abline(v = sim$param$stitest.start)
abline(h = 0.26, col = "red",  lty = 2)
title("HIV Prevalence")
#title("Syph Tprob = XXX, Syph.HIV.RR = , HIV.Syph.RR =", outer = TRUE)

# HIV and STI Incidence
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
plot(sim, y = "ir100")
abline(h = 3.8, col = "red", lty = 2)
abline(v = sim$param$stitest.start)
title("HIV Incidence")
plot(sim, y = "ir100.gc")
abline(h = 4.2, col = "red", lty = 2)
abline(v = sim$param$stitest.start)
title("GC Incidence")
plot(sim, y = "ir100.ct")
abline(h = 6.6, col = "red", lty = 2)
abline(v = sim$param$stitest.start)
title("CT Incidence")
plot(sim, y = "ir100.syph")
abline(h = 0.9, col = "red", lty = 2)
abline(v = sim$param$stitest.start)
title("Syph Incidence")

par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "hiv_sum", mean.col = "blue")
plot(sim, y = "sti_hiv_sum", mean.col = "green", add = TRUE)
plot(sim, y = "sti_u_hiv_sum", mean.col = "red", add = TRUE)
plot(sim, y = "sti_r_hiv_sum", mean.col = "orange", add = TRUE)
plot(sim, y = "sti_syph_hiv_sum", mean.col = "purple", add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("blue", "green", "red", "orange", "purple"), 
       c("HIV Sum", "STI & HIV sum", "USTI & HIV sum", "RSTI & HIV sum", "Syph & HIV sum"))

plot(sim, y = "sti_paf", mean.col = "blue")
plot(sim, y = "sti_u_paf", mean.col = "green", add = TRUE)
plot(sim, y = "sti_r_paf", mean.col = "red", add = TRUE)
plot(sim, y = "sti_syph_paf", mean.col = "orange", add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1), col = c("blue", "green", "red", "orange"), 
       c("STI PAF", "USTI PAF", "RSTI PAF", "Syph PAF"))

plot(sim, y = "sti_u_paf", mean.col = "blue")
plot(sim, y = "sti_u_sympt_paf", mean.col = "green", add = TRUE)
plot(sim, y = "sti_u_asympt_paf", mean.col = "red", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"), 
       c("U STI PAF", "U Sympt STI PAF", "U Asympt STI PAF"))

plot(sim, y = "sti_r_paf", mean.col = "blue")
plot(sim, y = "sti_r_sympt_paf", mean.col = "green", add = TRUE)
plot(sim, y = "sti_r_asympt_paf", mean.col = "red", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"), 
       c("R STI PAF", "R Sympt STI PAF", "R Asympt STI PAF"))

plot(sim, y = "sti_syph_paf", mean.col = "blue")
plot(sim, y = "sti_syph_sympt_paf", mean.col = "green", add = TRUE)
plot(sim, y = "sti_syph_asympt_paf", mean.col = "red", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"), 
       c("Syph PAF", "Syph Sympt PAF", "Syph Asympt PAF"))

# Evaluate whether testing is leading to treatment
a <- cbind(sim$epi$txCT, sim$epi$CTsympttests, sim$epi$CTasympttests.pos,
      sim$epi$txGC, sim$epi$GCsympttests, sim$epi$GCasympttests.pos,
      sim$epi$txsyph, sim$epi$syphsympttests, sim$epi$syphasympttests.pos)
colnames(a) <- c("txCT", "CTsympttests", "CTposasympttests",
                 "txGC", "GCsympttests", "GCposasympttests",
                 "txsyph", "syphsympttests", "syphposasympttests") 
tail(a)

par(mfrow = c(1,2), oma = c(0,0,2,0))
plot(sim, y = "GCasympttests", mean.col = "blue", ylim = c(0, 550))
plot(sim, y = "CTasympttests", mean.col = "red", add = TRUE)
plot(sim, y = "syphasympttests", mean.col = "green", add = TRUE)
plot(sim, y = "GCasympttests.pos", mean.col = "blue", add = TRUE, mean.lty = 2)
plot(sim, y = "CTasympttests.pos", mean.col = "red", add = TRUE, mean.lty = 2)
plot(sim, y = "syphasympttests.pos", mean.col = "green", add = TRUE, mean.lty = 2)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(1, 1, 1, 2, 2, 2), col = c("blue", "red", "green", "blue", "red", "green"), 
       c("GC Tests", "CT Tests", "Syph Tests", "Pos GC Tests", "Pos CT Tests", "Pos Syph Tests"))
title("Asymptomatic Tests")

plot(sim, y = "txGC", mean.col = "blue", mean.lty = 2, ylim = c(0, 100))
plot(sim, y = "txCT", mean.col = "red", mean.lty = 2, add = TRUE)
plot(sim, y = "txsyph", mean.col = "green", mean.lty = 2, add = TRUE)
plot(sim, y = "GCsympttests", mean.col = "blue", add = TRUE)
plot(sim, y = "CTsympttests", mean.col = "red", add = TRUE)
plot(sim, y = "syphsympttests", mean.col = "green", add = TRUE)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(1, 1, 1, 2, 2, 2), col = c("blue", "red", "green", "blue", "red", "green"), 
       c("GC Tests", "CT Tests", "Syph Tests", "Tx GC", "Tx CT", "TX Syph"))
title("Symptomatic Tests and All Treated")

plot(sim, y = "txGC", mean.col = "blue", mean.lty = 2, ylim = c(0, 90))
plot(sim, y = "GCsympttests", mean.col = "red", add = TRUE)
plot(sim, y = "GCasympttests.pos", mean.col = "green", add = TRUE)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(2, 1, 1), col = c("blue", "red", "green"), 
       c("Tx GC", "GC Sympt Tests", "GC Pos Asympt Tests"))
title("GC")

plot(sim, y = "txCT", mean.col = "blue", mean.lty = 2, ylim = c(0, 60))
plot(sim, y = "CTsympttests", mean.col = "red", add = TRUE)
plot(sim, y = "CTasympttests.pos", mean.col = "green", add = TRUE)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(2, 1, 1), col = c("blue", "red", "green"), 
       c("Tx CT", "CT Sympt Tests", "CT Pos Asympt Tests"))
title("CT")

plot(sim, y = "prev.rct", mean.col = "purple", ylim = c(0, 0.15))
plot(sim, y = "prev.uct", mean.col = "orange", add = TRUE)
plot(sim, y = "prev.syph", mean.col = "green", add = TRUE)
plot(sim, y = "prev.rgc", mean.col = "blue", add = TRUE)
plot(sim, y = "prev.ugc", mean.col = "red", add = TRUE)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("purple", "orange", "green","blue", "red"), 
       c("rCT", "uCT", "Syph", "rGC", "uGC"))
title("Prevalence")

plot(sim, y = "ir100.rct", mean.col = "purple", ylim = c(0, 40))
plot(sim, y = "ir100.uct", mean.col = "orange", add = TRUE)
plot(sim, y = "ir100.syph", mean.col = "green", add = TRUE)
plot(sim, y = "ir100.rgc", mean.col = "blue", add = TRUE)
plot(sim, y = "ir100.ugc", mean.col = "red", add = TRUE)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("purple", "orange", "green","blue", "red"), 
       c("rCT", "uCT", "Syph", "rGC", "uGC"))
title("Incidence")

# STI Testing Indications
par(mfrow = c(1,1), oma = c(0,0,2,0))
plot(sim, y = "stiactiveind", mean.col = "purple", xlim = c(0, sim$control$nsteps - sim$param$stitest.start), ylim = c(0, 1))
plot(sim, y = "newpartner", mean.col = "orange", add = TRUE)
plot(sim, y = "recentpartners", mean.col = "green", add = TRUE)
plot(sim, y = "concurrpart", mean.col = "blue", add = TRUE)
plot(sim, y = "partnersti", mean.col = "brown", add = TRUE)
plot(sim, y = "uai.nmain", mean.col = "black", add = TRUE)
plot(sim, y = "uai.any", mean.col = "gray", add = TRUE)
plot(sim, y = "recentSTI", mean.col = "yellow", add = TRUE)
legend("topright", lty = c(rep(1, 8)), 
       col = c("purple", "orange", "green", "blue", "brown", "black", "gray", "yellow"),
       c("stiactiveind", "newpartner", "recentpartners", "concurrpart", 
         "partnersti", "uai.nmain", "uai.any", "recentSTI"))
title("STI Testing Indications")

mean.s <- rbind(ir100.gc, ir100.ct, ir100, ir100.syph, i.prev, prev.primsecosyph.hivpos, 
            prev.primsecosyph.hivneg, prev.syph, prev.hiv.primsecosyphpos, prev.earlysyph, 
            prev.latesyph)
tar.syph <- c(4.2, 6.6, 3.8, 0.9, 0.26, 0.103, 0.026, 0.046, 0.498, 0.554, 0.446) #, 0.1385, 0.1385, 0.2770, 0.2000, 0.2000, 0.0460)

b <- data.frame(mean.s, tar.syph)
rownames(b) <- c("ir100.gc", "ir100.ct", "ir100", "ir100.syph", "i.prev", "prev.syph.hivpos", "prev.syph.hivneg", "prev.syph", "prev.hiv.syphpos", "prev.earlysyph", "prev.latesyph")
b

