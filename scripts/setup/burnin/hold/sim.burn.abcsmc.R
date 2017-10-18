
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("doParallel"))
suppressMessages(library("foreach"))

args <- commandArgs(trailingOnly = TRUE)
simno <- as.numeric(args[1])
batch <- as.numeric(args[2])
set.size <- as.numeric(args[3])

f <- function(simno, s, batch, set.size) {

  suppressMessages(library("EpiModelHIV"))
  sourceDir("source/", verbose = FALSE)

  fn <- paste0("simParms.b", batch, ".rda")
  load(fn)

  parm.row.id <- (simno - 1) * set.size + s
  if (parm.row.id > nrow(simParms)) {
    parm.row.id <- sample(nrow(simParms), 1)
  }
  parms <- as.list(simParms[parm.row.id, ])

  data(st)
  param <- param_msm(nwstats = st,
                     ai.scale = 1,

                     prep.coverage = 0,

                     rcomp.prob = 0,
                     rcomp.adh.groups = 0:3,
                     rcomp.main.only = FALSE,
                     rcomp.discl.only = FALSE,

                     rgc.tprob = parms[["rgc.tprob"]],
                     ugc.tprob = parms[["ugc.tprob"]],
                     rct.tprob = parms[["rct.tprob"]],
                     uct.tprob = parms[["uct.tprob"]],

                     rgc.sympt.prob = parms[["rgc.sympt.prob"]],
                     ugc.sympt.prob = parms[["ugc.sympt.prob"]],
                     rct.sympt.prob = parms[["rct.sympt.prob"]],
                     uct.sympt.prob = parms[["uct.sympt.prob"]],

                     rgc.dur.asympt = parms[["rgc.dur.asympt"]],
                     ugc.dur.asympt = parms[["ugc.dur.asympt"]],
                     gc.dur.tx = 2,
                     gc.dur.ntx = NULL,

                     rct.dur.asympt = parms[["rct.dur.asympt"]],
                     uct.dur.asympt = parms[["uct.dur.asympt"]],
                     ct.dur.tx = 2,
                     ct.dur.ntx = NULL,

                     gc.prob.cease = 0,
                     ct.prob.cease = 0,

                     gc.sympt.prob.tx = 0.90,
                     ct.sympt.prob.tx = 0.85,
                     gc.asympt.prob.tx = 0,
                     ct.asympt.prob.tx = 0,

                     prep.sti.screen.int = 182,
                     prep.sti.prob.tx = 1,

                     sti.cond.rr = 0.3,

                     hiv.rgc.rr = parms[["hiv.rect.rr"]],
                     hiv.ugc.rr = parms[["hiv.ureth.rr"]],
                     hiv.rct.rr = parms[["hiv.rect.rr"]],
                     hiv.uct.rr = parms[["hiv.ureth.rr"]],
                     hiv.dual.rr = 0)

  init <- init_msm(nwstats = st,
                   prev.B = 0.253,
                   prev.W = 0.253,
                   prev.ugc = 0.10,
                   prev.rgc = 0.10,
                   prev.uct = 0.10,
                   prev.rct = 0.10)

  control <- control_msm(simno = 1,
                         nsteps = 1560,
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


  data(est)
  sim <- netsim(est, param, init, control)

  df <- tail(as.data.frame(sim), 500)
  rect.prev <- mean(df$prev.rgcct)
  ureth.prev <- mean(df$prev.ugcct)
  gc.incid <- mean(df$ir100.gc)
  ct.incid <- mean(df$ir100.ct)
  hiv.prev <- mean(df$i.prev)

  out <- data.frame(rgc.tprob = parms[["rgc.tprob"]],
                    ugc.tprob = parms[["ugc.tprob"]],
                    rct.tprob = parms[["rct.tprob"]],
                    uct.tprob = parms[["uct.tprob"]],
                    rgc.dur.asympt = parms[["rgc.dur.asympt"]],
                    ugc.dur.asympt = parms[["ugc.dur.asympt"]],
                    rct.dur.asympt = parms[["rct.dur.asympt"]],
                    uct.dur.asympt = parms[["uct.dur.asympt"]],
                    rgc.sympt.prob = parms[["rgc.sympt.prob"]],
                    ugc.sympt.prob = parms[["ugc.sympt.prob"]],
                    rct.sympt.prob = parms[["rct.sympt.prob"]],
                    uct.sympt.prob = parms[["uct.sympt.prob"]],
                    hiv.rect.rr = parms[["hiv.rect.rr"]],
                    hiv.ureth.rr = parms[["hiv.ureth.rr"]],
                    rect.prev = rect.prev,
                    ureth.prev = ureth.prev,
                    gc.incid = gc.incid,
                    ct.incid = ct.incid,
                    hiv.prev = hiv.prev)

  return(out)

}


# Run batches of sims
registerDoParallel(parallel::detectCores())
sout <- foreach(s = 1:set.size) %dopar% {
  f(simno, s, batch, set.size)
}
sim <- as.data.frame(do.call("rbind", sout))

out.fn.all <- paste0("data/simDataAll.b", batch, ".rda")
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
