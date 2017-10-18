
library("EpiModelHIV")

sim <- NULL

# All data ------------------------------------------------------------

system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataAll.b2.rda scripts/burnin/")
load("scripts/burnin/simDataAll.b2.rda")
dim(sim)

# sim
cbind(sapply(sim, function(x) length(unique(x))))

# rgc, ugc, rct, uct
# targets <- c(0.083, 0.015, 0.118, 0.027)
targets <- c(0.102, 0.111, 0.141, 0.084)

nhist <- function(...) hist(..., col = "steelblue2", border = "white")

par(mfrow = c(2,2), mar = c(3,3,3,1), mgp = c(2,1,0))
nhist(sim$rgc.prev)
  abline(v = targets[1], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(sim$rgc.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(sim$rgc.prev) - targets[1], 3)), cex = 0.8)
nhist(sim$ugc.prev)
  abline(v = targets[2], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(sim$ugc.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(sim$ugc.prev) - targets[2], 3)), cex = 0.8)
nhist(sim$rct.prev)
  abline(v = targets[3], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(sim$rct.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(sim$rct.prev) - targets[3], 3)), cex = 0.8)
nhist(sim$uct.prev)
  abline(v = targets[4], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(sim$uct.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(sim$uct.prev) - targets[4], 3)), cex = 0.8)

par(mfrow = c(3,4))
for (i in 1:12) {
  nhist(sim[[i]], main = names(sim[i]))
}
pairs(sim[, 1:12])

# Accepted data -------------------------------------------------------

rejection <- function(sim, targets, threshold) {
  diff.rgc <- abs(sim$rgc.prev - targets[1])
  diff.ugc <- abs(sim$ugc.prev - targets[2])
  diff.rct <- abs(sim$rct.prev - targets[3])
  diff.uct <- abs(sim$uct.prev - targets[4])

  choice <- which(diff.rgc <= threshold & diff.ugc <= threshold &
                    diff.rct <= threshold & diff.uct <= threshold)
  simChosen <- sim[choice, ]
  cat("\n Accepted n:", nrow(simChosen))
  return(simChosen)
}

rejection2 <- function(sim, targets, threshold = 0.01) {
  edist <- sapply(1:nrow(sim), function(x) sqrt(sum((targets - sim[x, 13:16])^2)))
  edist.quant <- quantile(edist, threshold)
  accepted <- which(edist <= edist.quant)
  post <- sim[accepted, ]
  cat("\n Accepted n:", nrow(post))
  return(post)
}

# batch 1: 0.005
threshold <- 0.005

simChosen <- rejection(sim, targets = targets, threshold = 0.015)
simChosen <- rejection2(sim, targets = targets, threshold = 0.01)

# system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataChosen.b4.rda scripts/burnin/")
# load("scripts/burnin/simDataChosen.b4.rda")
# simChosen

par(mfrow = c(2,2), mar = c(3,3,3,1), mgp = c(2,1,0))
  nhist(simChosen$rgc.prev)
  abline(v = targets[1], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(simChosen$rgc.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(simChosen$rgc.prev) - targets[1], 3)), cex = 0.8)
nhist(simChosen$ugc.prev)
  abline(v = targets[2], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(simChosen$ugc.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(simChosen$ugc.prev) - targets[2], 3)), cex = 0.8)
nhist(simChosen$rct.prev)
  abline(v = targets[3], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(simChosen$rct.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(simChosen$rct.prev) - targets[3], 3)), cex = 0.8)
nhist(simChosen$uct.prev)
  abline(v = targets[4], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(simChosen$uct.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(simChosen$uct.prev) - targets[4], 3)), cex = 0.8)


mn <- apply(simChosen, 2, mean)
sds <- apply(simChosen, 2, sd)
lo <- mn - sds
hi <- mn + sds
data.frame(mn, lo, hi)

par(mfrow = c(3,4))
for (i in 1:12) {
  nhist(simChosen[[i]], main = names(simChosen[i]))
}
for (i in 1:12) {
  plot(density(simChosen[[i]]), main = names(simChosen[i]))
}

for (i in 1:12) {
  d <- density(simChosen[[i]])
  cat(names(simChosen[i]), ":", d$x[which.max(d$y)], "\n")
}
pairs(simChosen[,1:12])

# for next batch
save(simChosen, file = "scripts/burnin/simChosen.b1.rda")

system("scp scripts/burnin/simChosen.b1.rda hyak:/gscratch/csde/sjenness/sti/")

system("scp scripts/burnin/*.abcr.* hyak:/gscratch/csde/sjenness/sti/")
