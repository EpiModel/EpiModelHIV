#' @import data.table
#'
#'
#'
#'
#' @export
#'
cea_msm <- function(dat, at) {

  # Function Selection ------------------------------------------------------

  if (at < dat$param$cea.start) {
    return(dat)
  }

  dat <- track_cost_msm(dat, at)
  dat <- track_util_msm(dat, at)

  return(dat)
}

#' @import data.table
#'
#'
#'
#'
#'
#'
#' @export
#'
track_cost_msm <- function(dat, at) {

  # Attributes --------------------------------------------------------------
  alive <- dat$attr$alive
  prepStartTime <- dat$attr$prepStartTime
  prepStat <- dat$attr$prepStat
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  stage <- dat$attr$stage
  tx.status <- dat$attr$tx.status
  age <- dat$attr$age

  hiv.test.neg <- dat$epi$tot.neg.tests[at]
  hiv.test.pos <- dat$epi$newDx[at]

  OptimInitStarts <- dat$epi$OptimInitStarts[at]
  OptimAdhrStarts <- dat$epi$OptimAdhrStarts[at]
  OptimRetnStarts <- dat$epi$OptimRetnStarts[at]

  OptimRetnPrev <- dat$epi$OptimRetnPrev[at]
  OptimAdhrPrev <- dat$epi$OptimAdhrPrev[at]

  OptimInitUserStarts <- dat$epi$OptimInitUserStarts[at]

  # Parameters --------------------------------------------------------------

  prep.optim.start <- dat$param$prep.optim.start
  prep.optim.end <- dat$param$prep.optim.end

  # Discount rate
  d_r = dat$param$d_r

  # Costs
  # Initial office visit cost + 3CrCl test + HepB test
  init.visit.prep.cost <- dat$param$init.visit.prep.cost
  # Follow-up office visit cost
  ret.visit.prep.cost <- dat$param$ret.visit.prep.cost
  # HIVpos post-test counseling + step1 + step2
  # not sure how to incorporate step3 into EpiModel framework
  hiv.test.pos.cost <- dat$param$hiv.test.pos.cost
  # step1
  hiv.test.neg.cost <- dat$param$hiv.test.neg.cost
  # syphillis + chlamydia + gonorrhea + 3CrCl
  biannual.test.cost <- dat$param$biannual.test.cost

  prep.cost.ann <- dat$param$prep.cost.ann

  hiv.onart.cost.ann <- dat$param$hiv.onart.cost.ann
  hiv.offart.cost.ann <- dat$param$hiv.offart.cost.ann
  aids.onart.cost.ann <- dat$param$aids.onart.cost.ann
  aids.offart.cost.ann <- dat$param$aids.offart.cost.ann
  hiv.neg.cost.ann <- dat$param$hiv.neg.cost.ann

  prep.optim.init.init.pp <- dat$param$prep.optim.init.init.pp
  prep.optim.init.monthly.fixed <- dat$param$prep.optim.init.monthly.fixed
  prep.optim.adhr.init.pp <- dat$param$prep.optim.adhr.init.pp
  prep.optim.retn.init.pp <- dat$param$prep.optim.retn.init.pp
  prep.optim.retn.monthly.fixed <- dat$param$prep.optim.retn.monthly.fixed

  # HIV Test Costs
  hiv.test.cost <- hiv.test.neg * hiv.test.neg.cost + hiv.test.pos * hiv.test.pos.cost

  # PrEP Visit Costs
  num.init.visit <- sum(prepStat == 1 & prepStartTime == at, na.rm = TRUE)
  num.ret.visit <- sum(prepStat == 1 & (prepStartTime - at) %% 13 == 0
                     & prepStartTime != at, na.rm = TRUE)
  num.biannual.tests <- sum(prepStat == 1 & (prepStartTime - at) %% 26 == 0
                            & prepStartTime != at, na.rm = TRUE)

  init.visit.cost <- num.init.visit * init.visit.prep.cost
  ret.visit.cost <- num.ret.visit * ret.visit.prep.cost + num.biannual.tests * biannual.test.cost
  if (at < dat$param$end.horizon) {
    tot.visit.cost <- init.visit.cost + ret.visit.cost
  } else {
    tot.visit.cost <- 0
  }

  # PrEP Drug Costs
  num.prep <- sum(prepStat, na.rm = TRUE)
  if (at < dat$param$end.horizon) {
    prep.cost <- num.prep * prep.cost.ann / 52
  } else {
    prep.cost <- 0
  }


  # Annual Costs of Clinical Care
  ids.hiv.neg <- which(alive == 1 & status == 0)
  ids.hiv.onart <- which(alive == 1 & status == 1 &
                         stage != 4 & tx.status == 1)
  ids.hiv.offart <- which(alive == 1 & status == 1 &
                          stage != 4 & tx.status == 0)
  ids.aids.onart <- which(alive == 1 & status == 1 &
                          stage == 4 & tx.status == 1)
  ids.aids.offart <- which(alive == 1 & status == 1 &
                           stage == 4 & tx.status == 0)

  age.grp.neg <- as.character(cut(age,
                                  breaks = c(0, 45, 65, Inf),
                                  labels = c("18-44", "45-64", "65+")))
  age.grp.pos <- as.character(cut(age,
                                  breaks = c(0, 30, 40, 50, Inf),
                                  labels = c("18-29", "30-39", "40-49", "50+")))

  hiv.neg.cost <- sum(hiv.neg.cost.ann[.(age.grp.neg[ids.hiv.neg]), cost] / 52)
  hiv.onart.cost <- sum(hiv.onart.cost.ann[.(age.grp.pos[ids.hiv.onart]), cost] / 52)
  hiv.offart.cost <- sum(hiv.offart.cost.ann[.(age.grp.pos[ids.hiv.offart]), cost] /52)
  aids.onart.cost <- sum(aids.onart.cost.ann[.(age.grp.pos[ids.aids.onart]), cost] /52)
  aids.offart.cost <- sum(aids.offart.cost.ann[.(age.grp.pos[ids.aids.offart]), cost] /52)

  cc.cost <- sum(c(hiv.neg.cost, hiv.onart.cost, hiv.offart.cost,
    aids.onart.cost, aids.offart.cost), na.rm = TRUE)

  # Intervention Costs
  # Initialize
  int.cost <- 0
  init.int.cost <- 0
  adhr.int.cost <- 0
  retn.int.cost <- 0

  if (at >= prep.optim.start & at < prep.optim.end) {

    if (dat$param$prep.optim.init == TRUE) {
      init.pp.cost <- OptimInitUserStarts * prep.optim.init.init.pp
      init.fixed.cost <- prep.optim.init.monthly.fixed * 12 / 52

      init.int.cost <- init.pp.cost + init.fixed.cost
      int.cost <- int.cost + init.pp.cost + init.fixed.cost
    }

    if (dat$param$prep.optim.adhr == TRUE) {
      adhr.pp.cost <- OptimAdhrStarts * prep.optim.adhr.init.pp

      adhr.int.cost <- adhr.pp.cost
      int.cost <- int.cost + adhr.pp.cost
    }

    # Retention Intervention Weekly Costs per Person
    if (dat$param$prep.optim.retn == TRUE) {
      retn.pp.cost <- OptimRetnStarts * prep.optim.retn.init.pp
      retn.fixed.cost <- prep.optim.retn.monthly.fixed * 12 / 52

      retn.int.cost <- retn.pp.cost + retn.fixed.cost
      int.cost <- int.cost + retn.pp.cost + retn.fixed.cost
    }
  }

  # Discounting
  t <- (at - dat$param$cea.start)
  hiv.test.cost.disc <- hiv.test.cost * (1 - d_r) ^ (t / 52)
  tot.visit.cost.disc <- tot.visit.cost * (1 - d_r) ^ (t / 52)
  prep.cost.disc <- prep.cost * (1 - d_r) ^ (t / 52)
  cc.cost.disc <- cc.cost * (1 - d_r) ^ (t / 52)
  int.cost.disc <- int.cost * (1 - d_r) ^ (t / 52)

  init.int.cost.disc <- init.int.cost * (1 - d_r) ^ (t / 52)
  adhr.int.cost.disc <- adhr.int.cost * (1 - d_r) ^ (t / 52)
  retn.int.cost.disc <- retn.int.cost * (1 - d_r) ^ (t / 52)


  # Accounting
  total.cost <- sum(c(hiv.test.cost, tot.visit.cost,
                      prep.cost, cc.cost, int.cost), na.rm = TRUE)
  total.cost.disc <- sum(c(hiv.test.cost.disc, tot.visit.cost.disc,
                           prep.cost.disc, cc.cost.disc, int.cost.disc), na.rm = TRUE)

  if (at == dat$param$cea.start) {
    cuml.total.cost <- total.cost
    cuml.total.cost.disc <- total.cost.disc
    cuml.init.int.cost <- init.int.cost
    cuml.adhr.int.cost <- adhr.int.cost
    cuml.retn.int.cost <- retn.int.cost
    cuml.init.int.cost.disc <- init.int.cost.disc
    cuml.adhr.int.cost.disc <- adhr.int.cost.disc
    cuml.retn.int.cost.disc <- retn.int.cost.disc

  } else {
    cuml.total.cost <- dat$epi$cuml.total.cost[at - 1] + total.cost
    cuml.total.cost.disc <- dat$epi$cuml.total.cost.disc[at - 1] + total.cost.disc
    cuml.init.int.cost <- dat$epi$cuml.init.int.cost[at - 1] + init.int.cost
    cuml.adhr.int.cost <- dat$epi$cuml.adhr.int.cost[at - 1] + adhr.int.cost
    cuml.retn.int.cost <- dat$epi$cuml.retn.int.cost[at - 1] + retn.int.cost
    cuml.init.int.cost.disc <- dat$epi$cuml.init.int.cost.disc[at - 1] + init.int.cost.disc
    cuml.adhr.int.cost.disc <- dat$epi$cuml.adhr.int.cost.disc[at - 1] + adhr.int.cost.disc
    cuml.retn.int.cost.disc <- dat$epi$cuml.retn.int.cost.disc[at - 1] + retn.int.cost.disc
  }

  dat$epi$cuml.total.cost[at] <- cuml.total.cost
  dat$epi$total.cost[at] <- total.cost

  dat$epi$cuml.total.cost.disc[at] <- cuml.total.cost.disc
  dat$epi$total.cost.disc[at] <- total.cost.disc

  dat$epi$hiv.test.cost[at] <- hiv.test.cost
  dat$epi$tot.visit.cost[at] <- tot.visit.cost
  dat$epi$prep.cost[at] <- prep.cost
  dat$epi$int.cost[at] <- int.cost
  dat$epi$cc.cost[at] <- cc.cost

  dat$epi$hiv.test.cost.disc[at] <- hiv.test.cost.disc
  dat$epi$tot.visit.cost.disc[at] <- tot.visit.cost.disc
  dat$epi$prep.cost.disc[at] <- prep.cost.disc
  dat$epi$int.cost.disc[at] <- int.cost.disc
  dat$epi$cc.cost.disc[at] <- cc.cost.disc

  dat$epi$cuml.init.int.cost[at] <- cuml.init.int.cost
  dat$epi$cuml.adhr.int.cost[at] <- cuml.adhr.int.cost
  dat$epi$cuml.retn.int.cost[at] <- cuml.retn.int.cost
  dat$epi$cuml.init.int.cost.disc[at] <- cuml.init.int.cost.disc
  dat$epi$cuml.adhr.int.cost.disc[at] <- cuml.adhr.int.cost.disc
  dat$epi$cuml.retn.int.cost.disc[at] <- cuml.retn.int.cost.disc

  return(dat)
}

#' @import data.table
#'
#'
#'
#'
#'
#' @export
#'
track_util_msm <- function(dat, at) {

  # Attributes --------------------------------------------------------------

  alive <- dat$attr$alive
  prepStat <- dat$attr$prepStat
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  stage <- dat$attr$stage
  tx.status <- dat$attr$tx.status
  age <- dat$attr$age
  vl <- dat$attr$vl
  cuml.time.on.tx <- dat$attr$cuml.time.on.tx

  # Parameters --------------------------------------------------------------

  d_r <- dat$param$d_r

  nohiv.prep.qaly <- dat$param$nohiv.prep.qaly
  nohiv.noprep.qaly <- dat$param$nohiv.noprep.qaly
  aids.art.vs.qaly <- dat$param$aids.art.vs.qaly
  hiv.art.vs.qaly <- dat$param$hiv.art.vs.qaly
  aids.art.novs.qaly <- dat$param$aids.art.novs.qaly
  hiv.art.novs.qaly <- dat$param$hiv.art.novs.qaly
  aids.stopart.qaly <- dat$param$aids.stopart.qaly
  hiv.stopart.qaly <- dat$param$hiv.stopart.qaly
  aids.neverart.qaly <- dat$param$aids.neverart.qaly
  hiv.neverart.qaly <- dat$param$hiv.neverart.qaly

  age.adj <- dat$param$age.adj

  # QALY Calculations
  ids.nohiv.prep <- which(alive == 1 & prepStat == 1 & status == 0)
  ids.nohiv.noprep <- which(alive == 1 & prepStat == 0 & status == 0)
  # Other EpiModelHIV functions consider 200 vl to be suppressed cutoff--we are using 50 here intentionally
  ids.aids.art.vs <- which(alive == 1 & tx.status == 1 & vl <= log10(50) & stage == 4)
  ids.hiv.art.vs <- which(alive == 1 & tx.status == 1 & vl <= log10(50) & stage != 4)
  ids.aids.art.novs <- which(alive == 1 & tx.status == 1 & vl > log10(50) & stage == 4)
  ids.hiv.art.novs <- which(alive == 1 & tx.status == 1 & vl > log10(50) & stage != 4)
  ids.aids.stopart <- which(alive == 1 & tx.status == 0 & cuml.time.on.tx > 0 & stage == 4)
  ids.hiv.stopart <- which(alive == 1 & tx.status == 0 & cuml.time.on.tx > 0 & stage != 4)
  ids.aids.neverart <- which(alive == 1 & tx.status == 0 & cuml.time.on.tx == 0 & stage == 4)
  ids.hiv.neverart <- which(alive == 1 & tx.status == 0 & cuml.time.on.tx == 0 & stage != 4)


  nohiv.prep.qaly.sum <- sum((age.adj * age[ids.nohiv.prep] + nohiv.prep.qaly) / 52)
  nohiv.noprep.qaly.sum <- sum((age.adj * age[ids.nohiv.noprep] + nohiv.noprep.qaly) / 52)
  aids.art.vs.qaly.sum <- sum((age.adj * age[ids.aids.art.vs] + aids.art.vs.qaly) / 52)
  hiv.art.vs.qaly.sum <- sum((age.adj * age[ids.hiv.art.vs] + hiv.art.vs.qaly) / 52)
  aids.art.novs.qaly.sum <- sum((age.adj * age[ids.aids.art.novs] + aids.art.novs.qaly) / 52)
  hiv.art.novs.qaly.sum <- sum((age.adj * age[ids.hiv.art.novs] + hiv.art.novs.qaly) / 52)
  aids.stopart.qaly.sum <- sum((age.adj * age[ids.aids.stopart] + aids.stopart.qaly) / 52)
  hiv.stopart.qaly.sum <- sum((age.adj * age[ids.hiv.stopart] + hiv.stopart.qaly) / 52)
  aids.neverart.qaly.sum <- sum((age.adj * age[ids.aids.neverart] + aids.neverart.qaly) / 52)
  hiv.neverart.qaly.sum <- sum((age.adj * age[ids.hiv.neverart] + hiv.neverart.qaly) / 52)



  # Discounting
  t <- (at - dat$param$cea.start)
  nohiv.prep.qaly.sum.disc <- nohiv.prep.qaly.sum * (1 - d_r) ^ (t / 52)
  nohiv.noprep.qaly.sum.disc <- nohiv.noprep.qaly.sum * (1 - d_r) ^ (t / 52)
  aids.art.vs.qaly.sum.disc <- aids.art.vs.qaly.sum * (1 - d_r) ^ (t / 52)
  hiv.art.vs.qaly.sum.disc <- hiv.art.vs.qaly.sum * (1 - d_r) ^ (t / 52)
  aids.art.novs.qaly.sum.disc <- aids.art.novs.qaly.sum * (1 - d_r) ^ (t / 52)
  hiv.art.novs.qaly.sum.disc <- hiv.art.novs.qaly.sum * (1 - d_r) ^ (t / 52)
  aids.stopart.qaly.sum.disc <- aids.stopart.qaly.sum * (1 - d_r) ^ (t / 52)
  hiv.stopart.qaly.sum.disc <- hiv.stopart.qaly.sum * (1 - d_r) ^ (t / 52)
  aids.neverart.qaly.sum.disc <- aids.neverart.qaly.sum * (1 - d_r) ^ (t / 52)
  hiv.neverart.qaly.sum.disc <- hiv.neverart.qaly.sum * (1 - d_r) ^ (t / 52)

  # Accounting
  total.qaly <- sum(c(nohiv.prep.qaly.sum, nohiv.noprep.qaly.sum, aids.art.vs.qaly.sum,
                      hiv.art.vs.qaly.sum, aids.art.novs.qaly.sum, hiv.art.novs.qaly.sum,
                      aids.stopart.qaly.sum, hiv.stopart.qaly.sum, aids.neverart.qaly.sum,
                      hiv.neverart.qaly.sum), na.rm = TRUE)
  total.qaly.disc <- sum(c(nohiv.prep.qaly.sum.disc, nohiv.noprep.qaly.sum.disc, aids.art.vs.qaly.sum.disc,
                           hiv.art.vs.qaly.sum.disc, aids.art.novs.qaly.sum.disc, hiv.art.novs.qaly.sum.disc,
                           aids.stopart.qaly.sum.disc, hiv.stopart.qaly.sum.disc, aids.neverart.qaly.sum.disc,
                           hiv.neverart.qaly.sum.disc), na.rm = TRUE)

  dat$epi$total.qaly[at] <- total.qaly
  dat$epi$total.qaly.disc[at] <- total.qaly.disc

  if (at == dat$param$cea.start) {
    dat$epi$cuml.total.qaly[at] <- total.qaly
    dat$epi$cuml.total.qaly.disc[at] <- total.qaly.disc
  } else {
    dat$epi$cuml.total.qaly[at] <- dat$epi$cuml.total.qaly[at - 1] + total.qaly
    dat$epi$cuml.total.qaly.disc[at] <- dat$epi$cuml.total.qaly.disc[at - 1] + total.qaly.disc
  }

  dat$epi$num.nohiv.prep[at] <- length(ids.nohiv.prep)
  dat$epi$num.nohiv.noprep[at] <- length(ids.nohiv.noprep)
  dat$epi$num.aids.art.vs[at] <- length(ids.aids.art.vs)
  dat$epi$num.hiv.art.vs[at] <- length(ids.hiv.art.vs)
  dat$epi$num.aids.art.novs[at] <- length(ids.aids.art.novs)
  dat$epi$num.hiv.art.novs[at] <- length(ids.hiv.art.novs)
  dat$epi$num.aids.stopart[at] <- length(ids.aids.stopart)
  dat$epi$num.hiv.stopart[at] <- length(ids.hiv.stopart)
  dat$epi$num.aids.neverart[at] <- length(ids.aids.neverart)
  dat$epi$num.hiv.neverart[at] <- length(ids.hiv.neverart)

  return(dat)
}
