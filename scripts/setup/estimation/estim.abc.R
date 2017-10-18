
library(methods)
suppressMessages(library(EpiModel))
library(doParallel)
library(foreach)

## Environmental Arguments
args <- commandArgs(trailingOnly = TRUE)
simno <- args[1]

load("est/fit.rda")


# Casual Network ----------------------------------------------------------

# Pull casual network
est <- est[[2]]

# calibration function
f <- function(est) {

  est$coef.form[1] <- est$coef.form[1] + runif(1, -0.15, 0.15)
  est$coef.form[2] <- est$coef.form[2] + runif(1, -0.15, 0.15)
  est$coef.form[3] <- est$coef.form[3] + runif(1, -0.15, 0.15)
  est$coef.form[4] <- est$coef.form[4] + runif(1, -0.15, 0.15)

  dx <- netdx(est, nsims = 1, nsteps = 300, verbose = FALSE,
              set.control.ergm = control.simulate.ergm(MCMC.burnin = 2e6))
  stat1 <- mean(tail(dx$stats[[1]][, 1], 100))
  stat2 <- mean(tail(dx$stats[[1]][, 2], 100))
  stat3 <- mean(tail(dx$stats[[1]][, 3], 100))
  stat4 <- mean(tail(dx$stats[[1]][, 4], 100))
  return(c(p1 = est$coef.form[1],
           p2 = est$coef.form[2],
           p3 = est$coef.form[3],
           p4 = est$coef.form[4],
           out1 = stat1,
           out2 = stat2,
           out3 = stat3,
           out4 = stat4))
}

# Run parallel on each node
registerDoParallel(16)
nsims <- 500
sout <- foreach(s = 1:nsims) %dopar% {
  f(est)
}

sim <- as.data.frame(do.call("rbind", sout))
save(sim, file = paste0("data/sim.", simno, ".rda"))


# post-processing of data
fn <- list.files("est/", pattern = "sim", full.names = TRUE)

load(fn[1])
fsim <- sim
for (i in 2:length(fn)) {
  load(fn[i])
  fsim <- rbind(fsim, sim)
}
dim(fsim)

# rejection algorithm, weighted threshold
rejection <- function(sim,
                      target.stat = c(2022.5, 890, 950, 1185),
                      threshold = 0.05) {

  diff1 <- abs(sim$out1 - target.stat[1])
  diff2 <- abs(sim$out2 - target.stat[2])
  diff3 <- abs(sim$out3 - target.stat[3])
  diff4 <- abs(sim$out4 - target.stat[4])

  diff <- (3.5*diff1) + diff2 + (2*diff3) + diff4
  cutoff <- quantile(diff, threshold)

  in.threshold <- which(diff <= cutoff)

  post <- sim[in.threshold, ]
  return(post)
}

post <- rejection(fsim, threshold = 0.001)
post

# Accepted adjusted coefficients
colMeans(post)[5:8]
selection <- colMeans(post)[1:4]
selection # <- c(-13.2053497653944, -0.346004865323799, 1.04989624409589, -0.456503943369997)

# Test it
est2 <- est[[2]]
est2$coef.form[1:4] <- selection

dx <- netdx(est2, nsteps = 300, nsims = 20, ncores = 4, dynamic = TRUE,
            set.control.ergm = control.simulate.ergm(MCMC.burnin = 2e6))
dx
plot(dx)

# Write out to coefficients
load("est/fit.rda")
est[[2]]$coef.form[1:4] <- selection
save(est, file = "est/fit.rda")
