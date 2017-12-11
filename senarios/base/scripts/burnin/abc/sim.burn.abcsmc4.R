
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

                       msm.temp.adjust = x[2],
                       fa.temp.adjust = x[3],
                       depart.adjust = x[4],
                       return.adjust = x[5],
                       vi.scale = x[6],
                       cond.main.scale = x[7],
                       cond.pers.scale = x[8],
                       cond.inst.scale = x[9])
  


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
