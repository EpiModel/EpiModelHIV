context("Network Initialization and Summary Statistics")

test_that("calc_nwstats.mard testing", {

  tUnit <- 7

  num.B <- 5000
  num.W <- 5000

  deg.mp.B <- matrix(c(0.506, 0.151, 0.053,
                       0.207, 0.061, 0.022),
                     byrow = TRUE, nrow = 2)
  deg.mp.W <- matrix(c(0.435, 0.184, 0.095,
                       0.233, 0.033, 0.020),
                     byrow = TRUE, nrow = 2)

  mdeg.inst.B <- matrix(c(2.39, 3.27, 2.43,
                          2.17, 1.85, 1.85) * tUnit/100,
                        byrow = TRUE, nrow = 2)
  mdeg.inst.W <- matrix(c(1.82, 1.62, 2.66,
                          1.00, 2.23, 2.23) * tUnit/100,
                        byrow = TRUE, nrow = 2)

  prop.hom.mpi.B <- c(0.9484, 0.9019, 0.9085)
  prop.hom.mpi.W <- c(0.9154, 0.8509, 0.8944)

  sqrt.adiff.BB <- c(0.417, 0.498, 0.456)
  sqrt.adiff.BW <- c(0.454, 0.629, 0.585)
  sqrt.adiff.WW <- c(0.520, 0.632, 0.590)

  durs.main <- c(421, 662, 963) / tUnit
  durs.pers <- c(326, 344, 347) / tUnit

  ages <- 18:39
  asm.B.rate <- c(rep(0, 17),
              1 - (1 - c(rep(0.00159, 7),
                         rep(0.00225, 10),
                         rep(0.00348, 5))) ^ (1/(365/tUnit)),
              1)

  asm.W.rate <- c(rep(0, 17),
              1 - (1 - c(rep(0.00103, 7),
                         rep(0.00133, 10),
                         rep(0.00214, 5))) ^ (1/(365/tUnit)),
              1)

  role.B.prob <- c(0.242, 0.321, 0.437)
  role.W.prob <- c(0.228, 0.228, 0.544)

  st <- calc_nwstats.mard(
    tUnit = tUnit,
    num.B = num.B,
    num.W = num.W,
    deg.mp.B = deg.mp.B,
    deg.mp.W = deg.mp.W,
    mdeg.inst.B = mdeg.inst.B,
    mdeg.inst.W = mdeg.inst.W,
    prop.hom.mpi.B = prop.hom.mpi.B,
    prop.hom.mpi.W = prop.hom.mpi.W,
    balance = "mean",
    sqrt.adiff.BB = sqrt.adiff.BB,
    sqrt.adiff.WW = sqrt.adiff.WW,
    sqrt.adiff.BW = sqrt.adiff.BW,
    age.method = "heterogeneous",
    diss.main = ~offset(edges) + offset(nodemix("race", base = 1)),
    diss.pers = ~offset(edges) + offset(nodemix("race", base = 1)),
    durs.main = durs.main,
    durs.pers = durs.pers,
    ages = ages,
    asm.B.rate = asm.B.rate,
    asm.W.rate = asm.W.rate,
    role.B.prob = role.B.prob,
    role.W.prob = role.W.prob)

  expect_is(st, "nwstats")
})


test_that("base.nw.mard tests", {

  tUnit <- 7

  num.B <- 5000
  num.W <- 5000

  deg.mp.B <- matrix(c(0.506, 0.151, 0.053,
                       0.207, 0.061, 0.022),
                     byrow = TRUE, nrow = 2)
  deg.mp.W <- matrix(c(0.435, 0.184, 0.095,
                       0.233, 0.033, 0.020),
                     byrow = TRUE, nrow = 2)

  mdeg.inst.B <- matrix(c(2.39, 3.27, 2.43,
                          2.17, 1.85, 1.85) * tUnit/100,
                        byrow = TRUE, nrow = 2)
  mdeg.inst.W <- matrix(c(1.82, 1.62, 2.66,
                          1.00, 2.23, 2.23) * tUnit/100,
                        byrow = TRUE, nrow = 2)

  prop.hom.mpi.B <- c(0.9484, 0.9019, 0.9085)
  prop.hom.mpi.W <- c(0.9154, 0.8509, 0.8944)

  sqrt.adiff.BB <- c(0.417, 0.498, 0.456)
  sqrt.adiff.BW <- c(0.454, 0.629, 0.585)
  sqrt.adiff.WW <- c(0.520, 0.632, 0.590)

  durs.main <- c(421, 662, 963) / tUnit
  durs.pers <- c(326, 344, 347) / tUnit

  ages <- 18:39
  asm.B.rate <- c(rep(0, 17),
              1 - (1 - c(rep(0.00159, 7),
                         rep(0.00225, 10),
                         rep(0.00348, 5))) ^ (1/(365/tUnit)),
              1)

  asm.W.rate <- c(rep(0, 17),
              1 - (1 - c(rep(0.00103, 7),
                         rep(0.00133, 10),
                         rep(0.00214, 5))) ^ (1/(365/tUnit)),
              1)

  role.B.prob <- c(0.242, 0.321, 0.437)
  role.W.prob <- c(0.228, 0.228, 0.544)

  st <- calc_nwstats.mard(
    tUnit = tUnit,
    num.B = num.B,
    num.W = num.W,
    deg.mp.B = deg.mp.B,
    deg.mp.W = deg.mp.W,
    mdeg.inst.B = mdeg.inst.B,
    mdeg.inst.W = mdeg.inst.W,
    prop.hom.mpi.B = prop.hom.mpi.B,
    prop.hom.mpi.W = prop.hom.mpi.W,
    balance = "mean",
    sqrt.adiff.BB = sqrt.adiff.BB,
    sqrt.adiff.WW = sqrt.adiff.WW,
    sqrt.adiff.BW = sqrt.adiff.BW,
    age.method = "heterogeneous",
    diss.main = ~offset(edges) + offset(nodemix("race", base = 1)),
    diss.pers = ~offset(edges) + offset(nodemix("race", base = 1)),
    durs.main = durs.main,
    durs.pers = durs.pers,
    ages = ages,
    asm.B.rate = asm.B.rate,
    asm.W.rate = asm.W.rate,
    role.B.prob = role.B.prob,
    role.W.prob = role.W.prob)

  nw <- base_nw.mard(st)

  expect_is(nw, "network")
  expect_identical(list.vertex.attributes(nw),
                   c("na", "race", "role.class", "sqrt.age", "vertex.names"))
})
