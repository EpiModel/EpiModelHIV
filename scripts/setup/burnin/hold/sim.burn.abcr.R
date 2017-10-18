
# Packages ------------------------------------------------------------

library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("doParallel"))
suppressMessages(library("foreach"))


# Functions -----------------------------------------------------------

f <- function(batch) {

  suppressMessages(library("EpiModelHIV"))
  sourceDir("source/", verbose = FALSE)

  if (batch == 1) {
    rgc.tprob <- runif(1, 0.25, 0.6)
    ugc.tprob <- runif(1, 0.25, 0.6)
    rct.tprob <- runif(1, 0.25, 0.6)
    uct.tprob <- runif(1, 0.25, 0.6)

    gc.dur.ntx <- runif(1, 26, 52)
    ct.dur.ntx <- runif(1, 26, 52)

    rgc.sympt.prob <- runif(1, 0.05, 0.20)
    ugc.sympt.prob <- runif(1, 0.60, 0.90)
    rct.sympt.prob <- runif(1, 0.05, 0.20)
    uct.sympt.prob <- runif(1, 0.60, 0.95)

    gc.asympt.prob.tx <- runif(1, 0, 0.1)
    ct.asympt.prob.tx <- runif(1, 0, 0.1)
  }
  if (batch > 1) {
    load(paste0("simChosen.b", batch-1, ".rda"))
    mn <- apply(simChosen, 2, mean)
    sds <- apply(simChosen, 2, sd)
    lo <- mn - sds
    hi <- mn + sds

    rgc.tprob <- runif(1, lo[["rgc.tprob"]], hi[["rgc.tprob"]])
    ugc.tprob <- runif(1, lo[["ugc.tprob"]], hi[["ugc.tprob"]])
    rct.tprob <- runif(1, lo[["rct.tprob"]], hi[["rct.tprob"]])
    uct.tprob <- runif(1, lo[["uct.tprob"]], hi[["uct.tprob"]])

    gc.dur.ntx <- runif(1, lo[["gc.dur.ntx"]], hi[["gc.dur.ntx"]])
    ct.dur.ntx <- runif(1, lo[["ct.dur.ntx"]], hi[["ct.dur.ntx"]])

    rgc.sympt.prob <- runif(1, lo[["rgc.sympt.prob"]], hi[["rgc.sympt.prob"]])
    ugc.sympt.prob <- runif(1, lo[["ugc.sympt.prob"]], hi[["ugc.sympt.prob"]])
    rct.sympt.prob <- runif(1, lo[["rct.sympt.prob"]], hi[["rgc.sympt.prob"]])
    uct.sympt.prob <- runif(1, lo[["uct.sympt.prob"]], hi[["uct.sympt.prob"]])

    gc.asympt.prob.tx <- runif(1, lo[["gc.asympt.prob.tx"]], hi[["gc.asympt.prob.tx"]])
    ct.asympt.prob.tx <- runif(1, lo[["ct.asympt.prob.tx"]], hi[["ct.asympt.prob.tx"]])
  }

  load("est/nwstats.rda")
  param <- param_msm(nwstats = st,
                     ai.scale = 1,

                     prep.coverage = 0,

                     rcomp.prob = 0,
                     rcomp.adh.groups = 0:3,
                     rcomp.main.only = FALSE,
                     rcomp.discl.only = FALSE,

                     rgc.tprob = rgc.tprob,
                     ugc.tprob = ugc.tprob,
                     rct.tprob = rct.tprob,
                     uct.tprob = uct.tprob,

                     rgc.sympt.prob = rgc.sympt.prob,
                     ugc.sympt.prob = ugc.sympt.prob,
                     rct.sympt.prob = rct.sympt.prob,
                     uct.sympt.prob = uct.sympt.prob,

                     rgc.dur.asympt = gc.dur.ntx,
                     ugc.dur.asympt = gc.dur.ntx,
                     gc.dur.tx = 14/7,
                     gc.dur.ntx = gc.dur.ntx,

                     rct.dur.asympt = ct.dur.ntx,
                     uct.dur.asympt = ct.dur.ntx,
                     ct.dur.tx = 14/7,
                     ct.dur.ntx = ct.dur.ntx,

                     gc.prob.cease = 0,
                     ct.prob.cease = 0,

                     gc.sympt.prob.tx = 0.90,
                     ct.sympt.prob.tx = 0.85,
                     gc.asympt.prob.tx = gc.asympt.prob.tx,
                     ct.asympt.prob.tx = ct.asympt.prob.tx,

                     prep.sti.screen.int = 182,
                     prep.sti.prob.tx = 1,

                     sti.cond.rr = 0.3,

                     hiv.rgc.rr = 2.5,
                     hiv.ugc.rr = 1.25,
                     hiv.rct.rr = 2.5,
                     hiv.uct.rr = 1.25,
                     hiv.dual.rr = 0)

  init <- init_msm(nwstats = st,
                   prev.B = 0.253,
                   prev.W = 0.253,
                   prev.ugc = 0.10,
                   prev.rgc = 0.10,
                   prev.uct = 0.10,
                   prev.rct = 0.10)

  control <- control_msm(simno = 1,
                         nsteps = 10,
                         nsims = 1,
                         ncores = 1,
                         acts.FUN = acts_sti,
                         condoms.FUN = condoms_sti,
                         initialize.FUN = initialize_sti,
                         prep.FUN = prep_sti,
                         prev.FUN = prevalence_sti,
                         riskhist.FUN = riskhist_sti,
                         position.FUN = position_sti,
                         trans.FUN = trans_sti,
                         stitrans.FUN = sti_trans,
                         stirecov.FUN = sti_recov,
                         stitx.FUN = sti_tx,
                         verbose.FUN = verbose_sti,
                         verbose = FALSE,
                         module.order = c("aging.FUN", "deaths.FUN", "births.FUN",
                                          "test.FUN", "tx.FUN", "prep.FUN",
                                          "progress.FUN", "vl.FUN",
                                          "resim_nets.FUN", "disclose.FUN",
                                          "acts.FUN", "condoms.FUN", "riskhist.FUN",
                                          "position.FUN", "trans.FUN", "stitrans.FUN",
                                          "stirecov.FUN", "stitx.FUN", "prev.FUN"))


  load("est/fit.rda")
  sim <- netsim(est, param, init, control)

  df <- tail(as.data.frame(sim), 500)
  gc.prev <- mean(df$prev.gc)
  ct.prev <- mean(df$prev.ct)
  hiv.prev <- mean(df$i.prev)

  out <- data.frame(rgc.tprob = rgc.tprob,
                    ugc.tprob = ugc.tprob,
                    rct.tprob = rct.tprob,
                    uct.tprob = uct.tprob,
                    gc.dur.ntx = gc.dur.ntx,
                    ct.dur.ntx = ct.dur.ntx,
                    rgc.sympt.prob = rgc.sympt.prob,
                    ugc.sympt.prob = ugc.sympt.prob,
                    rct.sympt.prob = rct.sympt.prob,
                    uct.sympt.prob = uct.sympt.prob,
                    gc.asympt.prob.tx = gc.asympt.prob.tx,
                    ct.asympt.prob.tx = ct.asympt.prob.tx,
                    gc.prev = gc.prev,
                    ct.prev = ct.prev,
                    hiv.prev = hiv.prev)

  return(out)
}

rejection <- function(sim, targets, threshold) {
  diff.gc <- abs(sim$gc.prev - targets[1])
  diff.ct <- abs(sim$ct.prev - targets[2])

  choice <- which(diff.gc <= threshold & diff.ct <= threshold)
  simChosen <- sim[choice, , drop = FALSE]
  return(simChosen)
}


# Parameters ----------------------------------------------------------

target.n.chosen <- 250
# targets <- c(0.102, 0.111, 0.141, 0.084)
targets <- c(0.1080, 0.1043)
threshold <- 0.03
sims.per.batch <- 25

batch <- 1
out.fn.all <- paste0("data/simDataAll.b", batch, ".rda")
out.fn.chosen <- paste0("data/simDataChosen.b", batch, ".rda")


# Algorithm -----------------------------------------------------------

# if target.n.chosen is NULL, then do the rejection manually
if (is.null(target.n.chosen)) {

  # Run batches of sims
  cl <- makeCluster(parallel::detectCores())
  registerDoParallel(cl)
  nsims <- sims.per.batch
  sout <- foreach(s = 1:nsims) %dopar% {
    f(batch = batch)
  }
  stopCluster(cl)
  sim <- as.data.frame(do.call("rbind", sout))

  ind.fn <- tempfile(pattern = paste0("simDataAll.b", batch, "."),
                     tmpdir = "data/ind/", fileext = ".rda")

  # Save all sims so far
  if (!file.exists(out.fn.all)) {
    save(sim, file = out.fn.all)
    save(sim, file = ind.fn)
  } else {
    save(sim, file = ind.fn)
    simNew <- sim
    while(inherits(try(load(out.fn.all), silent = TRUE), "try-error")) {
      Sys.sleep(1)
    }
    sim <- rbind(sim, simNew)
    save(sim, file = out.fn.all)
  }

} else {

  # Load current simChosen file to get number already chosen
  if (file.exists(out.fn.chosen)) {
    while(inherits(try(load(out.fn.chosen), silent = TRUE), "try-error")) {
      Sys.sleep(1)
    }
    n.chosen <- nrow(simChosen)
  } else {
    n.chosen <- 0
  }

  # ABC-R Loop
  while (n.chosen < target.n.chosen) {

    # Run batches of sims
    cl <- makeCluster(parallel::detectCores())
    registerDoParallel(cl)
    nsims <- sims.per.batch
    sout <- foreach(s = 1:nsims) %dopar% {
      f(batch = batch)
    }
    stopCluster(cl)
    sim <- as.data.frame(do.call("rbind", sout))

    ind.fn <- tempfile(pattern = paste0("simDataAll.b", batch, "."),
                       tmpdir = "data/ind/", fileext = ".rda")

    # Rejection algorithm
    simChosen <- rejection(sim, targets = targets, threshold = threshold)

    # Save all sims so far
    if (!file.exists(out.fn.all)) {
      save(sim, file = out.fn.all)
      save(sim, file = ind.fn)
    } else {
      save(sim, file = ind.fn)
      simNew <- sim
      while(inherits(try(load(out.fn.all), silent = TRUE), "try-error")) {
        Sys.sleep(1)
      }
      sim <- rbind(sim, simNew)
      save(sim, file = out.fn.all)
    }

    # Save accepted sims
    if (!file.exists(out.fn.chosen)) {
      save(simChosen, file = out.fn.chosen)
    }
    if (file.exists(out.fn.chosen)) {
      if (nrow(simChosen) == 0) {
        load(out.fn.chosen)
      } else {
        simChosenNew <- simChosen
        while(inherits(try(load(out.fn.chosen), silent = TRUE), "try-error")) {
          Sys.sleep(1)
        }
        simChosen <- rbind(simChosen, simChosenNew)
        save(simChosen, file = out.fn.chosen)
      }
    }

    # Update n chosen within loop
    n.chosen <- nrow(simChosen)
  }

}
