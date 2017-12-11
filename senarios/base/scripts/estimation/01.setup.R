
## Paper 1 scenario setup file

rm(list = ls())
suppressMessages(library("EpiModelHIV"))


# Time unit for simulation, relative to 1 day
time.unit <- 7

# Population size by race
num.B <- 5000
num.W <- 5000

# mean/pers degree distributions matrices.
deg.mp.B <- deg.mp.W <-
  (matrix(c(0.506, 0.151, 0.053,
            0.207, 0.061, 0.022), byrow = TRUE, nrow = 2) +
   matrix(c(0.435, 0.184, 0.095,
            0.233, 0.033, 0.020), byrow = TRUE, nrow = 2))/2

# Instant rates
mdeg.inst.B <- mdeg.inst.W <-
  (matrix(c(0.010402, 0.012954, 0.011485,
            0.007912, 0.007424, 0.007424), byrow = TRUE, nrow = 2) +
   matrix(c(0.008186, 0.012017, 0.013024,
            0.008151, 0.008341, 0.008341), byrow = TRUE, nrow = 2))/2

# Quintile distribution of overall AI rates
qnts.W <- qnts.B <- c(0.0000,
                      0.0010,
                      0.0054,
                      0.0102,
                      0.0315)

# Proportion in same-race partnerships (main, casl, inst)
prop.hom.mpi.B <- prop.hom.mpi.W <- (c(0.9484, 0.9019, 0.9085) +
                                     c(0.9154, 0.8509, 0.8944))/2

# Mean age diffs (main, casl, inst)
sqrt.adiff.BB <- c(0.417, 0.498, 0.456)
sqrt.adiff.BW <- c(0.454, 0.629, 0.585)
sqrt.adiff.WW <- c(0.520, 0.632, 0.590)

# Mean durations
rates.main <- mean(c(0.002876,
                     0.002692,
                     0.001803))
rates.pers <- mean(c(0.007617,
                     0.003501,
                     0.006931))

durs.main <- 1/rates.main
durs.pers <- 1/rates.pers

# Age-sex-specific mortality rates
ages <- 18:39
asmr.B <- c(rep(0, 17),
            1-(1-c(rep(0.00159, 7),
                   rep(0.00225, 10),
                   rep(0.00348, 5)))^(1/(365/time.unit)), 1)

asmr.W <- c(rep(0, 17),
            1-(1-c(rep(0.00103, 7),
                   rep(0.00133, 10),
                   rep(0.00214, 5)))^(1/(365/time.unit)), 1)

# I, R, V role frequencies
role.B.prob <- role.W.prob <- (c(0.242, 0.321, 0.437) +
                               c(0.228, 0.228, 0.544))/2


# Create meanstats
st <- calc_nwstats.msm(
  method = 1,
  time.unit = time.unit,
  num.B = num.B,
  num.W = num.W,
  deg.mp.B = deg.mp.B,
  deg.mp.W = deg.mp.W,
  mdeg.inst.B = mdeg.inst.B,
  mdeg.inst.W = mdeg.inst.W,
  qnts.B = qnts.B,
  qnts.W = qnts.W,
  prop.hom.mpi.B = prop.hom.mpi.B,
  prop.hom.mpi.W = prop.hom.mpi.W,
  balance = "mean",
  sqrt.adiff.BB = sqrt.adiff.BB,
  sqrt.adiff.WW = sqrt.adiff.WW,
  sqrt.adiff.BW = sqrt.adiff.BW,
  diss.main = ~offset(edges),
  diss.pers = ~offset(edges),
  durs.main = durs.main,
  durs.pers = durs.pers,
  ages = ages,
  asmr.B = asmr.B,
  asmr.W = asmr.W,
  role.B.prob = role.B.prob,
  role.W.prob = role.W.prob)

save(st, file = "est/nwstats.rda")
rm(list = ls())
