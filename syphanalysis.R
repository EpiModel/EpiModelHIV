suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

# Syphilis ratio
par(mfrow = c(1, 1))
boxplot(sim$epi$prev.syph.hivpos / sim$epi$prev.syph.hivneg)
abline(h = 3.96, lty = 2, col = "red")

# Syphilis prevalence plots
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "prev.syph.hivpos", ylab = "Prevalence", ylim = c(0, 0.15), mean.col = "blue")
plot(sim, y = "prev.syph.hivneg", ylab = "Prevalence", add = TRUE, mean.col = "green")
abline(h = 0.103, col = "red", lty = 2)
abline(h = 0.026, col = "red", lty = 2)
title("Syphilis by HIV Status")
plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0.10, 0.20), col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, col = "green")
abline(h = 0.1385, col = "red", lty = 2)
abline(h = 0.1385, col = "red", lty = 2)
title("P and S Syphilis Prevalence")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", ylim = c(0.15, 0.35), col = "blue")
plot(sim, y = "prev.stage.latelat", ylab = "Prevalence", add = TRUE, col = "green")
abline(h = 0.2770, col = "red", lty = 2)
abline(h = 0.20, col = "red", lty = 2)
title("Early and Late Latent Syphilis Prevalence")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", ylim = c(0, 0.35), col = "blue")
plot(sim, y = "prev.stage.latelat", ylab = "Prevalence", add = TRUE, col = "green")
abline(h = 0.2770, col = "red", lty = 2)
abline(h = 0.046, col = "red", lty = 2)
title("Late Late Latent and Tertiary Syphilis Prevalence")

# All syphilis stages
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0, 1), mean.col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, mean.col = "green")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", add = TRUE, mean.col = "orange")
plot(sim, y = "prev.stage.latelat", ylab = "Prevalence", add = TRUE, mean.col = "gray")
plot(sim, y = "prev.stage.latelatelat", ylab = "Prevalence", add = TRUE, mean.col = "purple")
plot(sim, y = "prev.stage.tert", ylab = "Prevalence", add = TRUE, mean.col = "black")
title("Stage-specific Syphilis Prevalence")
abline(h = 0.1385, col = "blue", lty = 2)
abline(h = 0.1385, col = "green", lty = 2)
abline(h = 0.2770, col = "orange", lty = 2)
abline(h = 0.20, col = "gray", lty = 2)
abline(h = 0.20, col = "yellow", lty = 2)
abline(h = 0.046, col = "black", lty = 2)
legend("topright", c("Primary/Incub", "Secondary", "Early Latent", "Late latent", "Late late latent", "Tertiary"),
       col = c("blue", "green", "orange", "gray", "purple", "black"), lty = c(1, 1, 1, 1, 1, 1))


i.prev <- as.numeric(sim$epi$i.prev[2600, ])
ir100 <- as.numeric(sim$epi$ir100[2600, ])
prev.hiv.syphpos <- as.numeric(sim$epi$prev.hiv.syphpos[2600, ])

ir100.syph <- as.numeric(sim$epi$ir100.syph[2600, ])
prev.syph <- as.numeric(sim$epi$prev.syph[2600, ])
prev.syph.hivpos <- as.numeric(sim$epi$prev.syph.hivpos[2600, ])
prev.syph.hivneg <- as.numeric(sim$epi$prev.syph.hivneg[2600, ])
prev.stage.incubprim <- as.numeric(sim$epi$prev.stage.incubprim[2600, ])
prev.stage.seco <- as.numeric(sim$epi$prev.stage.seco[2600, ])
prev.stage.earlat <- as.numeric(sim$epi$prev.stage.earlat[2600, ])
prev.stage.latelat <- as.numeric(sim$epi$prev.stage.latelat[2600, ])
prev.stage.latelatelat <- as.numeric(sim$epi$prev.stage.latelatelat[2600, ])
prev.stage.tert <- as.numeric(sim$epi$prev.stage.tert[2600, ])
prev.earlysyph <- as.numeric(sim$epi$prev.earlysyph[2600, ])
prev.latesyph <- as.numeric(sim$epi$prev.latesyph[2600, ])

ir100.gc <- as.numeric(sim$epi$ir100.gc[2600, ])
ir100.ct <- as.numeric(sim$epi$ir100.ct[2600, ])

a <- cbind(c(i.prev, ir100, prev.syph.hivpos, prev.syph, prev.syph.hivneg, prev.stage.incubprim,
             prev.stage.seco, prev.stage.earlat, prev.stage.latelat, prev.stage.latelatelat,
             prev.stage.tert, prev.earlysyph, prev.latesyph))
rownames(a) <- c("i.prev", "ir100", "prev.syph.hivpos", "prev.syph", "prev.syph.hivneg", "prev.stage.incubprim",
                 "prev.stage.seco", "prev.stage.earlat", "prev.stage.latelat", "prev.stage.latelatelat",
                 "prev.stage.tert", "prev.earlysyph", "prev.latesyph")
a

mean.s <- c(ir100.gc, ir100.ct, ir100, ir100.syph, i.prev, prev.syph.hivpos, prev.syph.hivneg, prev.syph, prev.hiv.syphpos)
tar.syph <- c(4.2, 6.6, 3.8, 0.9, 0.26, 0.103, 0.026, 0.046, 0.498) #, 0.1385, 0.1385, 0.2770, 0.2000, 0.2000, 0.0460)

b <- data.frame(mean.s, tar.syph)
rownames(b) <- c("ir100.gc", "ir100.ct", "ir100", "ir100.syph", "i.prev", "prev.syph.hivpos", "prev.syph.hivneg", "prev.syph", "prev.hiv.syphpos")
b

# STI Prevalence
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "prev.syph", ylab = "Prevalence")
title("Syphilis Prevalence")
abline(h = 0.046, col = "red", lty = 2)
plot(sim, y = "prev.ct", ylab = "Prevalence")
title("CT Prevalence")
plot(sim, y = "prev.gc", ylab = "Prevalence")
title("GC Prevalence")
plot(sim, y = "i.prev", ylim = c(0, 0.3), ylab = "Prevalence")
abline(h = 0.26, col = "red",  lty = 2)
title("HIV Prevalence")
#title("Syph Tprob = XXX, Syph.HIV.RR = , HIV.Syph.RR =", outer = TRUE)

# HIV and STI Incidence
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "ir100")
abline(h = 3.8, col = "red", lty = 2)
title("HIV Incidence")
plot(sim, y = "ir100.gc")
abline(h = 4.2, col = "red", lty = 2)
title("GC Incidence")
plot(sim, y = "ir100.ct")
abline(h = 6.6, col = "red", lty = 2)
title("CT Incidence")
plot(sim, y = "ir100.syph")
abline(h = 0.9, col = "red", lty = 2)
title("Syph Incidence")
#title("Syph Tprob = XXX, Relrisk for Syph<->HIV = XXX", outer = TRUE)
