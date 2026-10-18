
#' @title Condom Use Module
#'
#' @description Module function stochastically simulates potential condom use
#'              for each act on the discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' For each act on the discordant edgelist, condom use is stochastically simulated
#' based on the partnership type and racial combination of the dyad. Other
#' modifiers for the probability of condom use in that pair are diagnosis of
#' disease, and full or partial HIV viral suppression
#' given HIV antiretroviral therapy.
#'
#' @return
#' Updates the discordant edgelist with a `uai` variable indicating whether
#' condoms were used in that act.
#'
#' @export
#'
condoms_msm <- function(dat, at) {

  ## Input
  # Attributes
  race              <- get_attr(dat, "race")
  diag.status       <- get_attr(dat, "diag.status")
  prepStat          <- get_attr(dat, "prepStat")
  age               <- get_attr(dat, "age")

  mdeCurr                  <- get_attr(dat, "mdeCurr")
  mde.symseverity          <- get_attr(dat, "mde.symseverity")
  mdd.tx.status            <- get_attr(dat, "mdd.tx.status")
  mdd.currtx.tottime       <- get_attr(dat, "mdd.currtx.tottime")
  last.neg.test            <- get_attr(dat, "last.neg.test")
  mdd.txintv.status        <- get_attr(dat, "mdd.txintv.status")
  mdd.txintvinit.time      <- get_attr(dat, "mdd.txintvinit.time")
  mdd.txintv.last.endtime  <- get_attr(dat, "mdd.txintv.last.endtime")


  # Parameters
  netstats        <- get_param(dat, "netstats")
  epistats        <- get_param(dat, "epistats")
  cond.scale      <- get_param(dat, "cond.scale")

  cond.prob.mult      <- get_param(dat, "cond.prob.mult")
  minadeq.txdur       <- get_param(dat, "minadeq.txdur")
  stigcounse.efxdur   <- get_param(dat, "stigcounse.efxdur")
  mh.stigtxefx.uai    <- get_param(dat, "mh.stigtxefx.uai")


  race.flag <- netstats[["race"]]
  geog.lvl  <- netstats[["race"]]

  # Condom Use Models
  cond.mc.mod <- epistats[["cond.mc.mod"]]
  cond.oo.mod <- epistats[["cond.oo.mod"]]


  ## Process

  # Temp edgelist
  el <- dat[["temp"]][["el"]]

  race.combo <- get_race_combo(race[el[, 1]], race[el[, 2]])

  comb.age <- age[el[, 1]] + age[el[, 2]]

  hiv.concord.pos <- rep(0, nrow(el))
  cp <- which(diag.status[el[, 1]] == 1 & diag.status[el[, 2]] == 1)
  hiv.concord.pos[cp] <- 1

  any.prep <- as.numeric((prepStat[el[, 1]] + prepStat[el[, 2]]) > 0)

  ## Main/casual partnerships ##
  mc.parts <- which(el[, "ptype"] != 3)
  el.mc <- el[mc.parts, ]

  if (!is.null(geog.lvl)) {
    if (race.flag == TRUE) {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = el.mc[, "durations"],
        race.combo = race.combo[mc.parts],
        comb.age = comb.age[mc.parts],
        hiv.concord.pos = hiv.concord.pos[mc.parts],
        prep = any.prep[mc.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
    } else {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = el.mc[, "durations"],
        comb.age = comb.age[mc.parts],
        hiv.concord.pos = hiv.concord.pos[mc.parts],
        prep = any.prep[mc.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
    }
  } else {
    if (race.flag == TRUE) {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = el.mc[, "durations"],
        race.combo = race.combo[mc.parts],
        comb.age = comb.age[mc.parts],
        hiv.concord.pos = hiv.concord.pos[mc.parts],
        prep = any.prep[mc.parts]
      )
      cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
    } else {
      x <- data.frame(
        ptype = el.mc[, "ptype"],
        duration.time = el.mc[, "durations"],
        comb.age = comb.age[mc.parts],
        hiv.concord.pos = hiv.concord.pos[mc.parts],
        prep = any.prep[mc.parts]
      )
      cond.prob <- unname(predict(cond.mc.mod, newdata = x, type = "response"))
    }
  }
  el.mc <- cbind(el.mc, cond.prob)

  ## One-off partnerships ##
  oo.parts <- which(el[, "ptype"] == 3)
  el.oo <- el[oo.parts, ]

  if (!is.null(geog.lvl)) {
    if (race.flag == TRUE) {
      x <- data.frame(
        race.combo = race.combo[oo.parts],
        comb.age = comb.age[oo.parts],
        hiv.concord.pos = hiv.concord.pos[oo.parts],
        prep = any.prep[oo.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
      el.oo <- cbind(el.oo, cond.prob)
    } else {
      x <- data.frame(
        comb.age = comb.age[oo.parts],
        hiv.concord.pos = hiv.concord.pos[oo.parts],
        prep = any.prep[oo.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
      el.oo <- cbind(el.oo, cond.prob)
    }
  } else {
    if (race.flag == TRUE) {
      x <- data.frame(
        race.combo = race.combo[oo.parts],
        comb.age = comb.age[oo.parts],
        hiv.concord.pos = hiv.concord.pos[oo.parts],
        prep = any.prep[oo.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
      el.oo <- cbind(el.oo, cond.prob)
    } else {
      x <- data.frame(
        comb.age = comb.age[oo.parts],
        hiv.concord.pos = hiv.concord.pos[oo.parts],
        prep = any.prep[oo.parts],
        geogYN = 1
      )
      cond.prob <- unname(predict(cond.oo.mod, newdata = x, type = "response"))
      el.oo <- cbind(el.oo, cond.prob)
    }
  }

  ## Bind el together
  el <- rbind(el.mc, el.oo)

  #if(at >= get_param(dat, "mhefx.start") & length(which(last.neg.test==at) > 10)) browser()

  # Acts
  ai.vec <- el[, "ai"]
  pid <- rep(seq_along(ai.vec), ai.vec)
  p1 <- rep(el[, "p1"], ai.vec)
  p2 <- rep(el[, "p2"], ai.vec)
  ptype <- rep(el[, "ptype"], ai.vec)

  cond.prob <- rep(el[, "cond.prob"], ai.vec)
  cond.prob <- cond.prob * cond.scale


  # Depression effects on condom use
  #get mdeCurr combo for both partners and determine elig for cond.prob modification
  #based on mdeCurr of the partners and the network type (only for casual ptypes)
  mdeCombo    <- mdeCurr[el[, 1]] + mdeCurr[el[, 2]]

  elig.cpmod.netmde  <- rep(0, nrow(el))
  netmde      <- which(el[, "ptype"] > 1 & mdeCombo > 0)
  elig.cpmod.netmde[netmde] <- 1


  #get mdd.tx.status and tot time on tx combo for both partners and determine elig for
  #modification based on being on depression tx for at least 4 weeks
  mddtxCombo  <- mdd.tx.status[el[, 1]] + mdd.tx.status[el[, 2]]
  mddtxTimeCombo  <- pmax(mdd.currtx.tottime[el[, 1]], mdd.currtx.tottime[el[, 2]])

  elig.cpmod.mddtx  <- rep(0, nrow(el))
  mddtxtime      <- which(mddtxCombo > 0 & mddtxTimeCombo >= minadeq.txdur)
  elig.cpmod.mddtx[mddtxtime] <- 1


  # id partnerships where one member has received stigma counseling in p6mo
  #(to model sustained interv impacts on condom usage post counseling)
  txintv.lastendtime     <- at - mdd.txintv.last.endtime
  txintv.endtimeCombo    <- pmax(txintv.lastendtime[el[, 1]], txintv.lastendtime[el[, 2]])

  elig.cpmod.txintv      <- rep(0, nrow(el))
  elig.cpmod.txintv[which(txintv.endtimeCombo <= stigcounse.efxdur)] <- 1


  #get mde.symseverity combo for both partners and determine magnitude of modification
  #to be implemented
  mdeSevCombo <- rep(0, nrow(el))
  mdeSevCombo <- pmax(mde.symseverity[el[, 1]], mde.symseverity[el[, 2]])


  cpmod  <- rep(0, nrow(el))
  cpmod[which(elig.cpmod.netmde == 1)] <- mdeSevCombo[which(elig.cpmod.netmde == 1)]
  cpmod[which(elig.cpmod.netmde == 1 & mdeSevCombo == 3)] <- 2
  cpmod[which(elig.cpmod.netmde == 1 & mdeSevCombo == 4)] <- 3
  cpmod <- cpmod + 1


  cpmod.multip <- cond.prob.mult[cpmod]
  #set multiplier for those eligible for modification based on tx status to result in
  #increased prob of using condoms (mimicking increased self-efficacy and higher likelihood of
  #adopting risk reduction strategies likely with depression tx)
  cpmod.multip[which(elig.cpmod.mddtx == 1 | elig.cpmod.txintv == 1)] <- mh.stigtxefx.uai



  el <- cbind(el, mdeCombo, mddtxCombo, mddtxTimeCombo, elig.cpmod.netmde,
              elig.cpmod.mddtx, mdeSevCombo, cpmod, cpmod.multip)


  # get per act versions
  mdeCombo2 <- rep(el[, "mdeCombo"], ai.vec)
  mddtxCombo2 <- rep(el[, "mddtxCombo"], ai.vec)
  mddtxTimeCombo2 <- rep(el[, "mddtxTimeCombo"], ai.vec)
  elig.cpmod.netmde2 <- rep(el[, "elig.cpmod.netmde"], ai.vec)
  elig.cpmod.mddtx2 <- rep(el[, "elig.cpmod.mddtx"], ai.vec)
  mdeSevCombo2 <- rep(el[, "mdeSevCombo"], ai.vec)
  cpmod2 <- rep(el[, "cpmod"], ai.vec)
  cpmod.multip2 <- rep(el[, "cpmod.multip"], ai.vec)


  # apply depression multiplier to cond.prob
  if(at >= get_param(dat, "mhefx.start")) {
    cond.prob <- plogis(qlogis(cond.prob) - log(cpmod.multip2))
  }



  # UAI draw per act
  uai <- runif(length(cond.prob)) < 1 - cond.prob

  # Act list construction
  al <- cbind(p1, p2, ptype, uai, pid)


  #epi trackers
  dat <- set_epi(dat, "mdd.ai.numall", at, length(pid))
  dat <- set_epi(dat, "mdd.ai.totmain", at, length(which(ptype == 1)))
  dat <- set_epi(dat, "mdd.ai.totcas", at, length(which(ptype == 2)))
  dat <- set_epi(dat, "mdd.ai.totoof", at, length(which(ptype == 3)))

  dat <- set_epi(dat, "mdd.ai.mde1.notx", at, length(which(elig.cpmod.netmde2 == 1 &
                                                             elig.cpmod.mddtx2 == 0)))
  dat <- set_epi(dat, "mdd.ai.mde1.tx", at, length(which(elig.cpmod.netmde2 == 1 &
                                                           elig.cpmod.mddtx2 == 1)))
  dat <- set_epi(dat, "mdd.ai.mde0", at, length(which(elig.cpmod.netmde2 == 0)))

  dat <- set_epi(dat, "mdd.ai.totcas.sev1", at,
                 length(which(ptype == 2 & mdeSevCombo2 == 1)))
  dat <- set_epi(dat, "mdd.ai.totcas.sev2", at,
                 length(which(ptype == 2 & mdeSevCombo2 == 2)))
  dat <- set_epi(dat, "mdd.ai.totcas.sev3", at,
                 length(which(ptype == 2 & mdeSevCombo2 == 3)))
  dat <- set_epi(dat, "mdd.ai.totcas.sev4", at,
                 length(which(ptype == 2 & mdeSevCombo2 == 4)))

  dat <- set_epi(dat, "mdd.uai.numall", at,  length(which(as.numeric(uai) == 1)))
  dat <- set_epi(dat, "mdd.uai.totmain", at,
                 length(which(as.numeric(uai) == 1 & ptype == 1)))
  dat <- set_epi(dat, "mdd.uai.totcas", at,
                 length(which(as.numeric(uai) == 1 & ptype == 2)))
  dat <- set_epi(dat, "mdd.uai.totoof", at,
                 length(which(as.numeric(uai) == 1 & ptype == 3)))

  dat <- set_epi(dat, "mdd.uai.mde1.notx", at,
                 length(which(as.numeric(uai) == 1 & elig.cpmod.netmde2 == 1 &
                                elig.cpmod.mddtx2 == 0)))
  dat <- set_epi(dat, "mdd.uai.mde1.tx", at,
                 length(which(as.numeric(uai) == 1 & elig.cpmod.netmde2 == 1 &
                                elig.cpmod.mddtx2 == 1)))
  dat <- set_epi(dat, "mdd.uai.mde0", at,
                 length(which(as.numeric(uai) == 1 & elig.cpmod.netmde2 == 0)))


  dat <- set_epi(dat, "mdd.uai.totcasoof.sev1", at,
                 length(which(as.numeric(uai) == 1 & ptype > 1 & mdeSevCombo2 == 1)))
  dat <- set_epi(dat, "mdd.uai.totcasoof.sev2", at,
                 length(which(as.numeric(uai) == 1 & ptype > 1 & mdeSevCombo2 == 2)))
  dat <- set_epi(dat, "mdd.uai.totcasoof.sev3", at,
                 length(which(as.numeric(uai) == 1 & ptype > 1 & mdeSevCombo2 == 3)))
  dat <- set_epi(dat, "mdd.uai.totcasoof.sev4", at,
                 length(which(as.numeric(uai) == 1 & ptype > 1 & mdeSevCombo2 == 4)))


  dat[["temp"]][["al"]] <- al

  ## Output

  return(dat)
}
