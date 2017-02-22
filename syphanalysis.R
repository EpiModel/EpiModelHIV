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
abline(h = 0.103, col = "blue", lty = 2)
abline(h = 0.026, col = "green", lty = 2)
title("Syphilis by HIV Status")
legend("topleft", c("HIV+", "HIV-"), col = c("blue", "green"), lty = c(1, 1))

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

mean.s <- rbind(ir100.gc, ir100.ct, ir100, ir100.syph, i.prev, prev.primsecosyph.hivpos, 
            prev.primsecosyph.hivneg, prev.syph, prev.hiv.primsecosyphpos, prev.earlysyph, 
            prev.latesyph)
tar.syph <- c(4.2, 6.6, 3.8, 0.9, 0.26, 0.103, 0.026, 0.046, 0.498, 0.554, 0.446) #, 0.1385, 0.1385, 0.2770, 0.2000, 0.2000, 0.0460)

b <- data.frame(mean.s, tar.syph)
rownames(b) <- c("ir100.gc", "ir100.ct", "ir100", "ir100.syph", "i.prev", "prev.syph.hivpos", "prev.syph.hivneg", "prev.syph", "prev.hiv.syphpos", "prev.earlysyph", "prev.latesyph")
b

