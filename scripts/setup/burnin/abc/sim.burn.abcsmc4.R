
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("doParallel"))
suppressMessages(library("foreach"))
suppressMessages(library("EasyABC"))

f <- function(x) {

  set.seed(x[1])

  suppressMessages(library("EpiModelHIV"))

  data(st)
  param <- param_msm(nwstats = st,

                     rgc.tprob = x[2],
                     ugc.tprob = x[3],
                     rct.tprob = x[4],
                     uct.tprob = x[5],
                     rct.asympt.int = x[6],
                     uct.asympt.int = x[6],
                     hiv.rgc.rr = x[7],
                     hiv.ugc.rr = x[8],
                     hiv.rct.rr = x[7],
                     hiv.uct.rr = x[8])

  init <- init_msm(nwstats = st)

  control <- control_msm(simno = 1,
                         nsteps = 2600,
                         nsims = 1, ncores = 1,
                         verbose = FALSE)

  data(est)
  sim <- netsim(est, param, init, control)

  df <- tail(as.data.frame(sim), 52)

  gc.incid <- mean(df$ir100.gc)
  ct.incid <- mean(df$ir100.ct)
  hiv.prev <- mean(df$i.prev)

  out <- c(gc.incid, ct.incid, hiv.prev)

  return(out)
}

priors <- list(c("unif", 0.35, 0.40),
               c("unif", 0.30, 0.35),
               c("unif", 0.21, 0.24),
               c("unif", 0.20, 0.22),
               c("unif", 260, 290),
               c("unif", 2.5, 3),
               c("unif", 1.5, 2))

targets <- c(4.2, 6.6, 0.26)

( nsim <- as.numeric(Sys.getenv("NSIM")) )
( pacc <- as.numeric(Sys.getenv("PACC")) )

a <- ABC_sequential(method = "Lenormand",
                    model = f,
                    prior = priors,
                    nb_simul = nsim,
                    summary_stat_target = targets,
                    p_acc_min = pacc,
                    progress_bar = TRUE,
                    n_cluster = 16,
                    use_seed = TRUE,
                    verbose = FALSE)

fn <- paste0("data/smcR.", pacc*100, "pct.", nsim, "sim.rda")
save(a, file = fn)
