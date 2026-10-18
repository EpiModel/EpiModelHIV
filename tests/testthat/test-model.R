# 0. Estimate Models ------------------------------------------------------

test_that("network functionality looks good", {
  param <- param_msm(netstats = netstats,
                     epistats = epistats)
  init <- init_msm()
  nnets <- 3L

  age_offset_stat_name <- "nodefactor.age>=age.sexual.cessation.TRUE"

  stat_names <- list(names(est[[1]]$coef.form),
                     c("edges", "degree0", "degree1", "degree2", "mean.age", age_offset_stat_name),
                     names(est[[3]]$coef.form))

  for (iteration in 1:2) {
    if (iteration == 1) {
      ## set up for start time == 1
      nsteps <- 12L
      control <- control_msm(nsteps = nsteps,
                             nsims = 2,
                             ncores = 1,
                             save.nwstats = TRUE,
                             initialize.FUN = function(x, param, init, control, s) {
                               dat <- initialize_msm(x, param, init, control, s)
                               summaries <- NULL
                               for (network in seq_len(dat$num.nw)) {
                                 summaries <- c(summaries, summary(~nodefactor(~age >= age.sexual.cessation, levels = I(TRUE)),
                                                                   basis = get_network(dat, network)))
                               }
                               ## record age.sexual.cessation nodefactor summaries immediately after initialization
                               ## (see departure.FUN for more information)
                               dat[["tempstats"]] <- summaries
                               return(dat)
                             },
                             departure.FUN = function(dat, at) {
                               dat <- departure_msm(dat, at)
                               ## record age.sexual.cessation nodefactors immediately after
                               ## departure module to confirm all edges above the sexual
                               ## cessation age were dissolved
                               age.sexual.cessation <- get_param(dat, "epistats")[["age.sexual.cessation"]]
                               summaries <- NULL
                               for (network in seq_len(dat$num.nw)) {
                                 summaries <- c(summaries, summary(~nodefactor(~age >= age.sexual.cessation, levels = I(TRUE)),
                                                                   basis = get_network(dat, network)))
                               }
                               dat[["tempstats"]] <- c(dat[["tempstats"]], summaries)
                               return(dat)
                             },
                             tergmLite.track.duration = multilayer(TRUE, TRUE, FALSE),
                             ## nwstats.formula records age.sexual.cessation nodefactors immediately after
                             ## resimulation module, confirming no new edges were formed above the sexual
                             ## cessation age
                             nwstats.formula = multilayer("formation",
                                                          ~edges
                                                           + degree(0:2)
                                                           + mean.age
                                                           + nodefactor(~age >= age.sexual.cessation, levels = I(TRUE)),
                                                          "formation"),
                             save.network = TRUE,
                             save.other = c("attr", "temp", "el", "net_attr", "tempstats"),
                             verbose = FALSE)
      basis <- est
    } else {
      ## set up for re-initialization
      control$start <- nsteps + 1L
      nsteps <- 29L
      control$nsteps <- nsteps
      control$initialize.FUN <- reinit_msm
      basis <- sim
    }

    sim <- netsim(basis, param, init, control)

    capture_output(
      print(sim)
    )
    plot(sim)
    plot(sim, type = "formation")

    tergm_nets <- c(1, 2)

    for (network in seq_len(nnets)) {
      capture_output(
        print(sim, network = network)
      )
      plot(sim, network = network)
      plot(sim, type = "formation", network = network)
    }

    for (simno in c(1, 2)) {
      ## check that there are no edges involving nodes above the sexual cessation age
      ## according to the temp statistics computed immediately after departures
      expect_equal(length(sim$tempstats[[simno]]), (nsteps - control$start + 1)*nnets)
      expect_true(is.numeric(sim$tempstats[[simno]]))
      expect_true(all(sim$tempstats[[simno]] == 0))

      for (network in seq_len(nnets)) {
        expect_s3_class(sim$network[[simno]][[network]], "networkLite")

        expect_equal(sim$el[[simno]][[network]], as.edgelist(sim$network[[simno]][[network]]))

        if (network %in% tergm_nets) {
          expect_equal(sim$network[[simno]][[network]] %n% "time", nsteps)
          expect_true(all((sim$network[[simno]][[network]] %n% "lasttoggle")[, 3] >= 0))
          expect_true(all((sim$network[[simno]][[network]] %n% "lasttoggle")[, 3] <= nsteps))
        }

        expect_equal(sim$coef.form[[simno]][[network]][1],
                     est[[network]]$coef.form[1] +
                       log(network.size(est[[network]]$newnetwork)) -
                       log(network.size(sim$network[[simno]][[network]])))

        expect_equal(network.size(est[[network]]$newnetwork), sim$epi$num[[simno]][1])
        expect_equal(network.size(sim$network[[simno]][[network]]), sim$epi$num[[simno]][nsteps])

        stats_matrix <- get_nwstats(sim, network = network, mode = "list")[[simno]]
        expect_equal(NCOL(stats_matrix), length(stat_names[[network]]))
        expect_equal(NROW(stats_matrix), nsteps)
        expect_equal(colnames(stats_matrix), stat_names[[network]])

        ## check that there are no edges involving nodes above the sexual cessation age
        ## according to the summary statistics computed immediately after resimulation
        expect_true(length(stats_matrix[, age_offset_stat_name]) == nsteps)
        expect_true(all(stats_matrix[, age_offset_stat_name] == 0))

        nwL <- networkLite(sim$el[[simno]][[network]], sim$attr[[simno]])
        if (network %in% tergm_nets) {
          nwL %n% "time" <- sim$network[[simno]][[network]] %n% "time"
          nwL %n% "lasttoggle" <- sim$network[[simno]][[network]] %n% "lasttoggle"
        }
        expect_equal(unlist(stats_matrix[nsteps, ]),
                     summary(get_network_control(sim, network, "nwstats.formula"), basis = nwL))
      }

      # check consistency of edgelists and nodal attributes
      expect_equal(sim$attr[[simno]][["deg.main"]], get_degree(sim$el[[simno]][[1]]))
      expect_equal(sim$attr[[simno]][["deg.casl"]], get_degree(sim$el[[simno]][[2]]))
      expect_equal(sim$attr[[simno]][["deg.tot"]], pmin(get_degree(sim$el[[simno]][[1]]) +
                                                      get_degree(sim$el[[simno]][[2]]), 3))

      # check degree constraints
      expect_true(all(get_degree(sim$el[[simno]][[1]]) < 3))
      expect_true(all(get_degree(sim$el[[simno]][[2]]) < 4))
    }
  }
})
