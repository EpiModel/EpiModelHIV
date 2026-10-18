# 0. Estimate Models ------------------------------------------------------

test_that("network functionality looks good", {
  param <- param_msm(netstats = netstats,
                     epistats = epistats)
  init <- init_msm()
  nsteps <- 52
  control <- control_msm(nsteps = nsteps,
                         nsims = 1,
                         ncores = 1,
                         save.nwstats = TRUE,
                         tergmLite.track.duration = TRUE,
                         save.network = TRUE,
                         resim_nets.FUN = function(dat, at) {
                           dat <- simnet_msm(dat, at)
                           dat <- set_attr(dat,
                                           "deg.casl",
                                           EpiModel::get_degree(dat$el[[2]]))
                           dat
                         },
                         verbose = FALSE)
  control$save.other <- c("nw", "attr", "el", "temp")
  sim <- netsim(est, param, init, control)

  capture_output(
    print(sim)
  )
  capture_output(
    print(sim, network = 1L)
  )
  capture_output(
    print(sim, network = 2L)
  )
  capture_output(
    print(sim, network = 3L)
  )

  plot(sim)
  plot(sim, network = 1L)
  plot(sim, network = 2L)
  plot(sim, network = 3L)

  plot(sim, type = "formation")
  plot(sim, type = "formation", network = 1L)
  plot(sim, type = "formation", network = 2L)
  plot(sim, type = "formation", network = 3L)

  # check networkLite class
  expect_s3_class(sim$nw[[1]][[1]], "networkLite")
  expect_s3_class(sim$nw[[1]][[2]], "networkLite")
  expect_s3_class(sim$nw[[1]][[3]], "networkLite")

  # check consistency of edgelists and network states
  expect_equal(sim$el[[1]][[1]], as.edgelist(sim$nw[[1]][[1]]))
  expect_equal(sim$el[[1]][[2]], as.edgelist(sim$nw[[1]][[2]]))
  expect_equal(sim$el[[1]][[3]], as.edgelist(sim$nw[[1]][[3]]))

  # check consistency of edgelists and nodal attributes
  expect_equal(sim$attr[[1]][["deg.main"]], get_degree(sim$el[[1]][[1]]))
  expect_equal(sim$attr[[1]][["deg.casl"]], get_degree(sim$el[[1]][[2]]))
  expect_equal(sim$attr[[1]][["deg.tot"]], pmin(get_degree(sim$el[[1]][[1]]) +
                                                  get_degree(sim$el[[1]][[2]]), 3))

  # check degree constraints
  expect_true(all(get_degree(sim$el[[1]][[1]]) < 3))
  expect_true(all(get_degree(sim$el[[1]][[2]]) < 4))

  # check timestamp on network
  expect_equal(sim$nw[[1]][[1]] %n% "time", nsteps)
  expect_equal(sim$nw[[1]][[2]] %n% "time", nsteps)

  # check lasttoggle information (for network 1)
  expect_true(all((sim$nw[[1]][[1]] %n% "lasttoggle")[, 3] >= 0))
  expect_true(all((sim$nw[[1]][[1]] %n% "lasttoggle")[, 3] <= nsteps))

  # check lasttoggle information (for network 2)
  expect_true(all((sim$nw[[1]][[2]] %n% "lasttoggle")[, 3] >= 0))
  expect_true(all((sim$nw[[1]][[2]] %n% "lasttoggle")[, 3] <= nsteps))

  # check network size adjustment to coefficients
  expect_equal(sim$nwparam[[1]]$coef.form[1],
               est[[1]]$coef.form[1] +
                 log(network.size(est[[1]]$newnetwork)) -
                 log(network.size(sim$nw[[1]][[1]])))

  expect_equal(sim$nwparam[[2]]$coef.form[1],
               est[[2]]$coef.form[1] +
                 log(network.size(est[[2]]$newnetwork)) -
                 log(network.size(sim$nw[[1]][[2]])))

  expect_equal(sim$nwparam[[3]]$coef.form[1],
               est[[3]]$coef.form[1] +
                 log(network.size(est[[3]]$newnetwork)) -
                 log(network.size(sim$nw[[1]][[3]])))

  # check endpoints of network size history
  expect_equal(network.size(est[[1]]$newnetwork), sim$epi$num[[1]][1])
  expect_equal(network.size(sim$nw[[1]][[1]]), sim$epi$num[[1]][nsteps])

  # check nwstats dimensions and names (for network 1)
  expect_equal(NCOL(sim$stats$nwstats[[1]][[1]]), length(est[[1]]$coef.form))
  expect_equal(NROW(sim$stats$nwstats[[1]][[1]]), nsteps)
  expect_equal(colnames(sim$stats$nwstats[[1]][[1]]), names(est[[1]]$coef.form))

  # check nwstats dimensions and names (for network 2)
  expect_equal(NCOL(sim$stats$nwstats[[1]][[2]]), length(est[[2]]$coef.form))
  expect_equal(NROW(sim$stats$nwstats[[1]][[2]]), nsteps)
  expect_equal(colnames(sim$stats$nwstats[[1]][[2]]), names(est[[2]]$coef.form))

  # check nwstats dimensions and names (for network 3)
  expect_equal(NCOL(sim$stats$nwstats[[1]][[3]]), length(est[[3]]$coef.form))
  expect_equal(NROW(sim$stats$nwstats[[1]][[3]]), nsteps)
  expect_equal(colnames(sim$stats$nwstats[[1]][[3]]), names(est[[3]]$coef.form))

  # check final nwstats values (for network 1)
  expect_equal(sim$stats$nwstats[[1]][[1]][nsteps, ],
               summary(est[[1]]$formation, basis = sim$nw[[1]][[1]]))

  # check final nwstats values (for network 2)
  expect_equal(sim$stats$nwstats[[1]][[2]][nsteps, ],
               summary(est[[2]]$formation, basis = sim$nw[[1]][[2]]))

  # check final nwstats values (for network 3)
  expect_equal(sim$stats$nwstats[[1]][[3]][nsteps, ],
               summary(est[[3]]$formation, basis = sim$nw[[1]][[3]]))
})
