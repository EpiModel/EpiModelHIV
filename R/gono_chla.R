#' @title Initialize The Gonorrhea And Chlamydia Epidemics
#'
#' @inheritParams aging_msm
#'
#' @param disease A character string containing either "gono" or "chla" to
#' select which of the Gonorrhea or Chlamydia epidemic to initialize.
#'
#' @return
#' Returns the `dat` object after the selected epidemic has been initialized.
#'
gono_chla_init <- function(dat, disease) {
  # Attributes
  role.class <- get_attr(dat, "role.class")
  n_new <- length(role.class)

  # Parameters
  rect.sympt.prob <- get_param(dat, paste0(disease, ".rect.sympt.prob"))
  uret.sympt.prob <- get_param(dat, paste0(disease, ".uret.sympt.prob"))
  screen.rate     <- get_param(dat, paste0(disease, ".screen.hivneg.rate"))

  # Init values
  prev.uret <- get_init(dat, paste0("prev.", disease, ".uret"))
  prev.rect <- get_init(dat, paste0("prev.", disease, ".rect"))

  ## Process -------------------------------------------------------------------
  uret_ids <- which(role.class %in% c(0, 2))
  rect_ids <- which(role.class %in% c(1, 2))

  # Initialize infection (both sites)
  uret_inf_ids <- sample(uret_ids, round(prev.uret * n_new), FALSE)
  rect_inf_ids <- sample(rect_ids, round(prev.rect * n_new), FALSE)

  # Initialize symptoms
  uret_sympt_ids <- uret_inf_ids[runif(length(uret_inf_ids)) < uret.sympt.prob]
  rect_sympt_ids <- rect_inf_ids[runif(length(rect_inf_ids)) < rect.sympt.prob]

  # Initialize last test
  test_last <- 1 - rgeom(n_new, screen.rate)

  ## Output --------------------------------------------------------------------
  dat <- set_attr(dat, paste0(disease, ".test.last"), test_last)
  # Uret
  dat <- set_attr(dat, paste0(disease, ".uret"),           1, uret_inf_ids)
  dat <- set_attr(dat, paste0(disease, ".uret.inf.last"),  1, uret_inf_ids)
  dat <- set_attr(dat, paste0(disease, ".uret.inf.count"), 1, uret_inf_ids)
  dat <- set_attr(dat, paste0(disease, ".uret.sympt"),     1, uret_sympt_ids)
  # Rect
  dat <- set_attr(dat, paste0(disease, ".rect"),           1, rect_inf_ids)
  dat <- set_attr(dat, paste0(disease, ".rect.inf.last"),  1, rect_inf_ids)
  dat <- set_attr(dat, paste0(disease, ".rect.inf.count"), 1, rect_inf_ids)
  dat <- set_attr(dat, paste0(disease, ".rect.sympt"),     1, rect_sympt_ids)

  return(dat)
}

#' @title Handles the Transmission of Gonorrhea And Chlamydia
#'
#' @details
#' Gonorrhea and Chlamydia are site specific diseases. It means that to get a
#' rectal infection, an individual must be in contact with a urethreal
#' infectious partner. (and vice-versa).
#'
#' Currently, the probability of infection per infectious act is governed by the
#' parameters `<disease>.uret.prob`. The probability of rectal infection is the
#' same value with a 1.25 odds ratio applied. This changed is thus applied when
#' `i_insertive == TRUE`. That is, when the infectious partner has a urethreal
#' infection and the susceptible one is at risk of a rectal infection.
#'
#' @param i_insertive Logical, is the infectious partner the insertive one? i.e.
#'                    the infection will occur from a urethreal infected partner
#'                    to a rectal susceptible one.
#'
#' @inheritParams gono_chla_init
#'
#' @return
#' Returns the `dat` object with updated attributes for Gonorrhea and
#' Chlamydia infections.
#'
gono_chla_transmission <- function(dat, disease, i_insertive) {
  # Setup the variable differently if the infected partner is insertive
  #   site_i: site of disease for the infectious partner
  #   site_s: site of disease for the susceptible partner
  if (i_insertive) {
    site_i <- ".uret"
    site_s <- ".rect"
  } else {
    site_i <- ".rect"
    site_s <- ".uret"
  }

  ## Inputs --------------------------------------------------------------------
  at <- get_current_timestep(dat)
  # Attributes
  race <- get_attr(dat, "race")
  dis.i.status <- get_attr(dat, paste0(disease, site_i))
  dis.i.inf.last <- get_attr(dat, paste0(disease, site_i, ".inf.last"))
  dis.s.status <- get_attr(dat, paste0(disease, site_s))
  # Params
  dis.prob <- get_param(dat, paste0(disease, ".uret.prob"))
  if (i_insertive) {
    dis.prob <- apply_odds_ratio(dis.prob, 1.25)
  }

  sympt.prob <- get_param(dat, paste0(disease, site_s, ".sympt.prob"))
  sti.cond.eff.or  <- get_param(dat, "sti.cond.eff.or")
  cond.fail.or   <- get_param(dat, "cond.fail.or")

  # Pull act list
  al <- dat[["temp"]][["al"]]

  ## Process -------------------------------------------------------------------
  #
  # Create an act list where P1 always transmit to P2
  #
  # Conditions:
  #   - p1 infected on `site_i`
  #   - p1 infection occured before the current step (< at)
  #   - `ins == i_insertive`: correct position for the infection to occur
  #   - p2 not infected on `site_s`
  #
  # Repeat with p1 and p2 inverted for the situation where p2 is infected

  # acts where p1 can transmit to p2
  p1i_acts <- which(dis.i.status[al[, "p1"]] == 1)
  p1i_acts <- p1i_acts[dis.i.inf.last[al[p1i_acts, "p1"]] < at]
  p1i_acts <- p1i_acts[al[p1i_acts, "ins"] == as.numeric(i_insertive)]
  p1i_acts <- p1i_acts[dis.s.status[al[p1i_acts, "p2"]] == 0]
  # acts where p2 can transmit to p1
  p2i_acts <- which(dis.i.status[al[, "p2"]] == 1)
  p2i_acts <- p2i_acts[dis.i.inf.last[al[p2i_acts, "p2"]] < at]
  p2i_acts <- p2i_acts[al[p2i_acts, "ins"] != as.numeric(i_insertive)]
  p2i_acts <- p2i_acts[dis.s.status[al[p2i_acts, "p1"]] == 0]

  # Store the acts so that infection always is from p1 to p2
  elig_acts <- matrix(
    c(
      al[p1i_acts, "p1"], al[p2i_acts, "p2"],  # 1st column: infected partner
      al[p1i_acts, "p2"], al[p2i_acts, "p1"],  # 2nd column: susceptible partner
      al[p1i_acts, "uai"], al[p2i_acts, "uai"] # 3rd column: condomless act flag
    ),
    ncol = 3
  )
  colnames(elig_acts) <- c("infectious", "susceptible", "uai")

  # Default probability
  inf_prob <- rep(dis.prob, nrow(elig_acts))

  # Condom modifier: UAI == 0 -> a condom is used
  cond_acts <- which(elig_acts[, "uai"] == 0)
  p_condom <- if (i_insertive) "infectious" else "susceptible" # which wears the condom
  race_condom <- race[elig_acts[cond_acts, p_condom]] # race stratification
  condom_or <- 1 - (sti.cond.eff.or - cond.fail.or[race_condom])
  inf_prob[cond_acts] <- apply_odds_ratio(inf_prob[cond_acts], condom_or)

  # Stocastic transmission
  new_inf_ids <- elig_acts[runif(nrow(elig_acts)) < inf_prob, "susceptible"]
  new_inf_ids <- unique(new_inf_ids) # deduplicate (multiple acts/partnerships)
  # Symptoms attribution
  sympt_ids <- new_inf_ids[runif(length(new_inf_ids)) < sympt.prob]

  ## Output --------------------------------------------------------------------
  # Attributes
  dat <- set_attr(dat, paste0(disease, site_s),              1, new_inf_ids)
  dat <- set_attr(dat, paste0(disease, site_s, ".inf.last"), at, new_inf_ids)
  dat <- set_attr(dat, paste0(disease, site_s, ".sympt"),    1, sympt_ids)
  # Update inf counter
  count_name <- paste0(disease, site_s, ".inf.count")
  count <- get_attr(dat, count_name, new_inf_ids)
  dat <- set_attr(dat, count_name, count + 1, new_inf_ids)
  # Trackers
  dat <- set_epi(dat, paste0("incid.", paste0(disease, site_s)), at,
                 length(new_inf_ids))

  return(dat)
}

#' @title Handles the Natural History of Gonorrhea And Chlamydia
#'
#' @details
#' Gonorrhea and Chlamydia natural history is currently modeled as simply
#' infected or not. This function thus handeles the recovery of Gonorrhea and
#' Chlamydia.
#'
#' The recovery is possible either through treatment or spontaneously. At each
#' time step a stochastic process determines if the infection recovers. The
#' probability of recovery per time step is derived from the parameters
#' `<disease>.ntx.recov.int` and `<disease>.tx.recov.int` for spontaneous and
#' treatment recovery respectively.
#'
#' @inheritParams gono_chla_init
#'
#' @return
#' Returns the `dat` object with updated attributes for Gonorrhea and
#' Chlamydia infections.
#'
gono_chla_progression <- function(dat, disease) {
  ## Inputs --------------------------------------------------------------------
  at <- get_current_timestep(dat)
  # Attributes
  dis.rect <- get_attr(dat, paste0(disease, ".rect"))
  dis.uret <- get_attr(dat, paste0(disease, ".uret"))
  rect.inf.last <- get_attr(dat, paste0(disease, ".rect.inf.last"))
  uret.inf.last <- get_attr(dat, paste0(disease, ".uret.inf.last"))
  dis.tx <- get_attr(dat, paste0(disease, ".tx"))

  # Parameters
  ntx.recov.rate <- 1 / get_param(dat, paste0(disease, ".ntx.recov.int"))
  tx.recov.rate <- 1 / get_param(dat, paste0(disease, ".tx.recov.int"))

  ## Process -------------------------------------------------------------------
  # Untreated, asymptomatic and symptomatic, infection occured before `at`
  elig_ids <- which(dis.tx == 0 & (
    (dis.rect == 1 & rect.inf.last < at) |
      (dis.uret == 1 & uret.inf.last < at)
  ))
  recov_ntx_ids <- elig_ids[runif(length(elig_ids)) < ntx.recov.rate]

  # Treated, asymptomatic and symptomatic
  elig_ids <- which(dis.tx == 1)
  recov_tx_ids <- elig_ids[runif(length(elig_ids)) < tx.recov.rate]

  # Simple combine as they are exclusive
  recov_ids <- c(recov_ntx_ids, recov_tx_ids)

  ## Outputs -------------------------------------------------------------------
  # reset attrs for the recovered
  dat <- set_attr(dat, paste0(disease, ".rect"),       0, recov_ids)
  dat <- set_attr(dat, paste0(disease, ".uret"),       0, recov_ids)
  dat <- set_attr(dat, paste0(disease, ".rect.sympt"), 0, recov_ids)
  dat <- set_attr(dat, paste0(disease, ".uret.sympt"), 0, recov_ids)
  dat <- set_attr(dat, paste0(disease, ".dx"),         0, recov_ids)
  dat <- set_attr(dat, paste0(disease, ".tx"),         0, recov_ids)

  return(dat)
}

#' @title Handles the Screening for Gonorrhea And Chlamydia
#'
#' @details
#' Screening for Gonorrhea and Chlamydia can occur at every timestep for all
#' nodes not diagnosed with the disease. Different probablity are used for nodes
#' with an HIV positive diagnostic, for nodes on PrEP and for nodes not
#' diagnosed with HIV and not on PrEP.
#'
#' By default, nodes are tested for urethreal infection. They then have a
#' probability of also being tested for rectal infection.
#'
#' Additionally, the PrEP starters will test for both diseases at both sites on
#' the step they initiate PrEP. (PrEP occurs before STIs)
#' NOTE:: do we actually want rectal test for all PrEP starter? Or should we use
#' the `sti.screen.rect.prep.prob`?
#'
#' @inheritParams gono_chla_init
#'
#' @return
#' Returns the `dat` object with updated testing attributes for Gonorrhea and
#' Chlamydia
#'
gono_chla_testing <- function(dat, disease) {
  # Inputs
  at <- get_current_timestep(dat)
  # Attributes
  dis.rect <- get_attr(dat, paste0(disease, ".rect"))
  dis.uret <- get_attr(dat, paste0(disease, ".uret"))
  dx <- get_attr(dat, paste0(disease, ".dx"))

  prep.any <- get_attr(dat, "prep.any")
  prep <- get_attr(dat, "prep")
  prep.start.last <- get_attr(dat, "prep.start.last")
  hiv.dx <- get_attr(dat, "diag.status")

  prep.otc <- get_attr(dat, "prep.otc")
  prep.otc.start.last <- get_attr(dat, "prep.otc.start.last")

  # parameters
  screen.hivneg.rate <- get_param(dat, paste0(disease, ".screen.hivneg.rate"))
  screen.hivpos.rate <- get_param(dat, paste0(disease, ".screen.hivpos.rate"))
  screen.prep.rate <- get_param(dat, "sti.screen.prep.rate")

  screen.prep.otc.rate <- get_param(dat, "sti.screen.prep.otc.rate")

  screen.rect.hivpos.prob <- get_param(dat, "sti.screen.rect.hivpos.prob")
  screen.rect.hivneg.prob <- get_param(dat, "sti.screen.rect.hivneg.prob")
  screen.rect.prep.prob <- get_param(dat, "sti.screen.rect.prep.prob")

  screen.rect.prep.otc.prob <- get_param(dat, "sti.screen.rect.prep.otc.prob")
  prep.otc.always.sti.tst <- get_param(dat, "prep.otc.always.sti.tst")

  ## Processs ------------------------------------------------------------------

  # Create a list to store the 3 possible screening paths
  uret_ids <- vector(mode = "list", 4)
  names(uret_ids) <- c("nprep", "prep", "hivpos", "prep.otc")
  rect_ids <- uret_ids

  # HIV Negative - not on PrEP
  elig_ids <- which(dx == 0 & hiv.dx == 0 & prep.any == 0)
  uret_ids$nprep <- elig_ids[runif(length(elig_ids)) < screen.hivneg.rate]
  # rectal testing (some proportion of the urethreal screened)
  rect_ids$nprep <- uret_ids$nprep[
    runif(length(uret_ids$nprep)) < screen.rect.hivneg.prob
  ]

  # HIV Negative - on PrEP
  elig_ids <- which(dx == 0 & prep == 1)
  uret_ids$prep <- elig_ids[runif(length(elig_ids)) < screen.prep.rate]
  # rectal testing (some proportion of the urethreal screened)
  rect_ids$prep <- uret_ids$prep[
    runif(length(uret_ids$prep)) < screen.rect.prep.prob
  ]
  # unconditionnaly add those who started PrEP this step
  prep_start_ids <- elig_ids[prep.start.last[elig_ids] == at]
  uret_ids$prep <- unique(c(uret_ids$prep, prep_start_ids))
  rect_ids$prep <- unique(c(rect_ids$prep, prep_start_ids))

  # HIV Positive
  elig_ids <- which(dx == 0 & hiv.dx == 1)
  uret_ids$hivpos <- elig_ids[runif(length(elig_ids)) < screen.hivpos.rate]
  # rectal testing (some proportion of the urethreal screened)
  rect_ids$hivpos <- uret_ids$hivpos[
    runif(length(uret_ids$hivpos)) < screen.rect.hivpos.prob
  ]

  # HIV Negative - on PrEP OTC
  elig_ids <- which(dx == 0 & prep.otc == 1)
  uret_ids$prep.otc <- elig_ids[runif(length(elig_ids)) < screen.prep.otc.rate]
  # rectal testing (some proportion of the urethreal screened)
  rect_ids$prep.otc <- uret_ids$prep.otc[
    runif(length(uret_ids$prep.otc)) < screen.rect.prep.otc.prob
  ]
  if (prep.otc.always.sti.tst == 1) {
    prep_otc_start_ids <- elig_ids[prep.otc.start.last[elig_ids] == at]
    uret_ids$prep.otc <- unique(c(uret_ids$prep.otc, prep_otc_start_ids))
    rect_ids$prep.otc <- unique(c(rect_ids$prep.otc, prep_otc_start_ids))
  }


  ## Outputs -------------------------------------------------------------------
  # Merge all 3 (no need for `union` as they are mutually exclusive)
  final_uret_ids <- unlist(uret_ids)
  final_rect_ids <- unlist(rect_ids)

  # rect_ids is a subset of uret_ids so no need to use it here
  dat <- set_attr(dat, paste0(disease, ".test.last"), at, final_uret_ids)

  # diagnosed are: screened && infected at the correct site
  dx_ids <- unique(c(
    final_uret_ids[dis.uret[final_uret_ids] == 1],
    final_rect_ids[dis.rect[final_rect_ids] == 1]
  ))

  dat <- set_attr(dat, paste0(disease, ".dx"), 1, dx_ids)
  dat <- set_attr(dat, paste0(disease, ".dx.last"), at, dx_ids)

  if (get_control(dat, "debug")) {
    for (grp in names(uret_ids)) {
      dat <- set_epi(dat, paste0("dbg_", disease, "_uret_tst_", grp), at,
                     length(uret_ids[[grp]]))
      dat <- set_epi(dat, paste0("dbg_", disease, "_rect_tst_", grp), at,
                     length(rect_ids[[grp]]))
    }
    dat <- set_epi(dat, paste0("dbg_", disease, "_dx"), at, length(dx_ids))
  }

  return(dat)
}

#' @title Handles the Treatment Attribution for Gonorrhea And Chlamydia
#'
#' @details
#' Treatment attribution for Gonorrhea and Chlamydia can happend via two
#' non-exclusive mecanismes.
#'
#' Symptomatic treatment is a one time process. It happens one step after
#' infection. A symptomatic individual have a single opportunity to get treated
#' for this episode, with probability `<disease>.sympt.tx.prob`.
#'
#' Asymptomatic treatment is also a one time process. When a node is diagnosed
#' with either Gonorrhea or Chlamydia it has a single opportunity to get
#' treated. With probability `<disease>.asympt.tx.prob`.
#'
#' @inheritParams gono_chla_init
#'
#' @return
#' Returns the `dat` object with updated treatment attributes for Gonorrhea and
#' Chlamydia
#'
gono_chla_treatment <- function(dat, disease) {
  ## Input
  at <- get_current_timestep(dat)
  # Attributes
  rect.sympt <- get_attr(dat, paste0(disease, ".rect.sympt"))
  uret.sympt <- get_attr(dat, paste0(disease, ".uret.sympt"))
  rect.inf.last <- get_attr(dat, paste0(disease, ".rect.inf.last"))
  uret.inf.last <- get_attr(dat, paste0(disease, ".uret.inf.last"))
  dx.last <- get_attr(dat, paste0(disease, ".dx.last"))
  tx <- get_attr(dat, paste0(disease, ".tx"))
  prep <- get_attr(dat, "prep")
  prep.otc <- get_attr(dat, "prep.otc")
  prep.any <- get_attr(dat, "prep.any")
  race <- get_attr(dat, "race")

  # Parameters
  sympt.tx.prob <- get_param(dat, paste0(disease, ".sympt.tx.prob"))
  asympt.tx.prob <- get_param(dat, paste0(disease, ".asympt.tx.prob"))
  prep.tx.prob <- get_param(dat, "sti.prep.tx.prob")
  prep.otc.tx.prob <- get_param(dat, "sti.prep.otc.tx.prob")

  # Process --------------------------------------------------------------------

  ## Symptomatic treatment eligible (infected on the previous step)
  # no need to check the absence of treatment as the process occurs only once
  elig_ids <- unique(c(
    which(rect.sympt == 1 & rect.inf.last == at - 1),
    which(uret.sympt == 1 & uret.inf.last == at - 1)
  ))
  rates <- sympt.tx.prob[race[elig_ids]]
  sympt_tx_ids <- elig_ids[runif(length(elig_ids)) < rates]

  ## Asymptomatic diagnosed this step
  dx_ids <- which(dx.last == at) # precalculation to speed up the subsetting

  # not on PrEP
  elig_ids <- dx_ids[prep.any[dx_ids] == 0]
  elig_ids <- elig_ids[tx[elig_ids] == 0]
  rates <- asympt.tx.prob[race[elig_ids]]
  nprep_tx_ids <- elig_ids[runif(length(elig_ids)) < rates]

  # on PrEP
  elig_ids <- dx_ids[prep[dx_ids] == 1]
  elig_ids <- elig_ids[tx[elig_ids] == 0]
  prep_tx_ids <- elig_ids[runif(length(elig_ids)) < prep.tx.prob]

  # on PrEP OTC
  elig_ids <- dx_ids[prep.otc[dx_ids] == 1]
  elig_ids <- elig_ids[tx[elig_ids] == 0]
  prep_otc_tx_ids <- elig_ids[runif(length(elig_ids)) < prep.otc.tx.prob]

  ## Combine all of the treated (potential overlap)
  tx_ids <- unique(c(sympt_tx_ids, nprep_tx_ids, prep_tx_ids, prep_otc_tx_ids))

  # Output ---------------------------------------------------------------------
  dat <- set_attr(dat, paste0(disease, ".tx"),        1, tx_ids)
  dat <- set_attr(dat, paste0(disease, ".dx"),        1, tx_ids)
  dat <- set_attr(dat, paste0(disease, ".test.last"), at, tx_ids)
  dat <- set_attr(dat, paste0(disease, ".dx.last"),   at, tx_ids)
  # NOTE:: becoming treated also sets the tests and diagnostics attributes

  return(dat)
}
