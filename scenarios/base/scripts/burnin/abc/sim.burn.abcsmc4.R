
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("doParallel"))
suppressMessages(library("foreach"))
suppressMessages(library("EasyABC"))

f <- function(x) {

  set.seed(x[1])

  suppressMessages(library("EpiModelHIV"))

  load(file = "~est/fit.rda")
  load(file = "~data.params.rda")
  
  param <- param_shamp(nwstats = st,

                       vi.foi.scale = x[2],
                       msm.foi.scale = x[3],
                       fa.foi.scale = x[4])
  


  init <- init_shamp(nwstats = st)

  control <- control_msm(simno = 1,
                         nsteps = 200,
                         nsims = 1, ncores = 1,
                         verbose = FALSE)

  data(est)
  sim <- netsim(est, param, init, control)

  df <- tail(as.data.frame(sim), 52)

  prop.MSM.inf <- mean(df$prop.MSM.inf)
  prop.MSMds.inf <- mean(df$prop.MSMds.inf)
  prop.FA.inf <- mean(df$prop.FA.inf)
  prop.FAds.inf <- mean(df$prop.FAds.inf)
  prop.Lhet.inf <- mean(df$prop.Lhet.inf)
  hiv.prev<- mean(df$i.prev)

  out <- c(prop.MSM.inf, prop.MSMds.inf, prop.FA.inf, prop.FAds.inf, prop.Lhet.inf, hiv.prev)

  return(out)
}

priors <- list(c("unif", 1, 12),
               c("unif", 1, 0.35),
               c("unif", 1, 0.24))

targets <- c(.0027)

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
