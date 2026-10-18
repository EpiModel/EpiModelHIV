
#' @title Asymptomatic STI Screening
#'
#' @description Module function for STI screening of persons with asymptomatic infection.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This STI screening module implements site-specific (urogenital and rectal) screening of
#' gonorrhea and chlamydia. Screening rates are stratified by PrEP status (for MSM not diagnosed
#' with HIV) and diagnosed HIV status. Stochastic or deterministic screening rates may be
#' implemented with the `sti.screening.dist` parameter.
#'
#' @return
#' This function returns the `dat` object with updated `diag.status.XX`, `XX.last.screen.time`, and
#' `XX.last.dxpos.time` attributes for each STI. Summary statistics are tracked for the total number
#' of screening events and positive diagnoses for the strata above.
#'
#' @keywords module msm
#'
#' @export
#'
stiscreen_msm <- function(dat, at) {

  # Setup -----

  ## Attributes
  rGC <- get_attr(dat, "rGC")
  uGC <- get_attr(dat, "uGC")
  rCT <- get_attr(dat, "rCT")
  uCT <- get_attr(dat, "uCT")

  diag.status.gc <- get_attr(dat, "diag.status.gc")
  diag.status.gc.prep <- get_attr(dat, "diag.status.gc.prep")
  diag.status.ct <- get_attr(dat, "diag.status.ct")
  diag.status.ct.prep <- get_attr(dat, "diag.status.ct.prep")

  gc.last.screen.time <- get_attr(dat, "gc.last.screen.time")
  ct.last.screen.time <- get_attr(dat, "ct.last.screen.time")
  gc.last.dxpos.time <- get_attr(dat, "gc.last.dxpos.time")
  ct.last.dxpos.time <- get_attr(dat, "ct.last.dxpos.time")

  diag.status <- get_attr(dat, "diag.status")

  prepStat <- get_attr(dat, "prepStat")
  prepStartTime <- get_attr(dat, "prepStartTime")


  ## Parameters
  sti.screening.dist <- get_param(dat, "sti.screening.dist")
  if (!(sti.screening.dist %in% c("deterministic", "stochastic")))
    stop("sti.screening.dist must be `deterministic` or `stochastic`")

  # Screening rates for urogenital infection, by HIV and PrEP status
  gc.screen.hivneg.rate <- get_param(dat, "gc.screen.hivneg.rate")
  gc.screen.hivpos.rate <- get_param(dat, "gc.screen.hivpos.rate")
  ct.screen.hivneg.rate <- get_param(dat, "ct.screen.hivneg.rate")
  ct.screen.hivpos.rate <- get_param(dat, "ct.screen.hivpos.rate")
  sti.screen.prep.start <- get_param(dat, "sti.screen.prep.start")
  sti.screen.prep.rate <- get_param(dat, "sti.screen.prep.rate")

  # Probability of rectal screening given urogenital screening
  sti.screen.rect.hivneg.prob <- get_param(dat, "sti.screen.rect.hivneg.prob")
  sti.screen.rect.hivpos.prob <- get_param(dat, "sti.screen.rect.hivpos.prob")
  sti.screen.rect.prep.prob <- get_param(dat, "sti.screen.rect.prep.prob")


  # Screening for Non-HIV+ Dx, Off-PrEP -----

  gc.hivneg.nprep.elig <- which(
    (diag.status.gc == 0 | is.na(diag.status.gc)) &
    (diag.status == 0 | is.na(diag.status)) &
    prepStat == 0)
  gc.hivneg.nprep <- sti_screen_subfx(
    elig = gc.hivneg.nprep.elig,
    test.dist = sti.screening.dist,
    screen.rate = gc.screen.hivneg.rate,
    last.screen.time = gc.last.screen.time,
    r.screen.prob = sti.screen.rect.hivneg.prob,
    u.status = uGC,
    r.status = rGC,
    at = at
  )
  diag.status.gc[gc.hivneg.nprep$dx.pos.ids] <- 1
  gc.last.screen.time[gc.hivneg.nprep$screen.ids] <- at
  gc.last.dxpos.time[gc.hivneg.nprep$dx.pos.ids] <- at

  ct.hivneg.nprep.elig <- which(
    (diag.status.ct == 0 | is.na(diag.status.ct)) &
    (diag.status == 0 | is.na(diag.status)) &
    prepStat == 0)
  ct.hivneg.nprep <- sti_screen_subfx(
    elig = ct.hivneg.nprep.elig,
    test.dist = sti.screening.dist,
    screen.rate = ct.screen.hivneg.rate,
    last.screen.time = ct.last.screen.time,
    r.screen.prob = sti.screen.rect.hivneg.prob,
    u.status = uCT,
    r.status = rCT,
    at = at
  )
  diag.status.ct[ct.hivneg.nprep$dx.pos.ids] <- 1
  ct.last.screen.time[ct.hivneg.nprep$screen.ids] <- at
  ct.last.dxpos.time[ct.hivneg.nprep$dx.pos.ids] <- at


  # Screening for Non-HIV+ Dx, On-PrEP -----
  if (at >= sti.screen.prep.start) {
    gc.hivneg.prep.elig <- which(
      (diag.status.gc == 0 | is.na(diag.status.gc)) &
      (diag.status == 0 | is.na(diag.status)) &
      prepStat == 1)
    gc.hivneg.prep <- sti_screen_subfx(
      elig = gc.hivneg.prep.elig,
      test.dist = sti.screening.dist,
      screen.rate = sti.screen.prep.rate,
      last.screen.time = gc.last.screen.time,
      r.screen.prob = sti.screen.rect.prep.prob,
      u.status = uGC,
      r.status = rGC,
      prep = TRUE,
      prep.start.time = prepStartTime,
      at = at
    )
    diag.status.gc[gc.hivneg.prep$dx.pos.ids] <- 1
    diag.status.gc.prep[gc.hivneg.prep$dx.pos.ids] <- 1
    gc.last.screen.time[gc.hivneg.prep$screen.ids] <- at
    gc.last.dxpos.time[gc.hivneg.prep$dx.pos.ids] <- at

    ct.hivneg.prep.elig <- which(
      (diag.status.ct == 0 | is.na(diag.status.ct)) &
      (diag.status == 0 | is.na(diag.status)) &
      prepStat == 1)
    ct.hivneg.prep <- sti_screen_subfx(
      elig = ct.hivneg.prep.elig,
      test.dist = sti.screening.dist,
      screen.rate = sti.screen.prep.rate,
      last.screen.time = ct.last.screen.time,
      r.screen.prob = sti.screen.rect.prep.prob,
      u.status = uCT,
      r.status = rCT,
      prep = TRUE,
      prep.start.time = prepStartTime,
      at = at
    )
    diag.status.ct[ct.hivneg.prep$dx.pos.ids] <- 1
    diag.status.ct.prep[ct.hivneg.prep$dx.pos.ids] <- 1
    ct.last.screen.time[ct.hivneg.prep$screen.ids] <- at
    ct.last.dxpos.time[ct.hivneg.prep$dx.pos.ids] <- at
  }


  # Screening for HIV+ Dx -----

  gc.hivpos.elig <- which(
    (diag.status.gc == 0 | is.na(diag.status.gc)) &
    diag.status == 1 &
    prepStat == 0)
  gc.hivpos <- sti_screen_subfx(
    elig = gc.hivpos.elig,
    test.dist = sti.screening.dist,
    screen.rate = gc.screen.hivpos.rate,
    last.screen.time = gc.last.screen.time,
    r.screen.prob = sti.screen.rect.hivpos.prob,
    u.status = uGC,
    r.status = rGC,
    at = at
  )
  diag.status.gc[gc.hivpos$dx.pos.ids] <- 1
  gc.last.screen.time[gc.hivpos$screen.ids] <- at
  gc.last.dxpos.time[gc.hivpos$dx.pos.ids] <- at

  ct.hivpos.elig <- which(
    (diag.status.ct == 0 | is.na(diag.status.ct)) &
    diag.status == 1 &
    prepStat == 0)
  ct.hivpos <- sti_screen_subfx(
    elig = ct.hivpos.elig,
    test.dist = sti.screening.dist,
    screen.rate = ct.screen.hivpos.rate,
    last.screen.time = ct.last.screen.time,
    r.screen.prob = sti.screen.rect.hivpos.prob,
    u.status = uCT,
    r.status = rCT,
    at = at
  )
  diag.status.ct[ct.hivpos$dx.pos.ids] <- 1
  ct.last.screen.time[ct.hivpos$screen.ids] <- at
  ct.last.dxpos.time[ct.hivpos$dx.pos.ids] <- at

  # Calculate sti.dx.n.window --------------------------------------------------
  # GC
  dat <- update_windowed_counter(
    dat,
    posit_ids = unique(c(
      gc.hivneg.nprep$dx.pos.ids,
      gc.hivneg.prep$dx.pos.ids,
      gc.hivpos$dx.pos.ids
    )),
    window_size = get_param(dat, "gc.dx.window"),
    counter_name = "gc.dx.count.window"
  )

  # CT
  dat <- update_windowed_counter(
    dat,
    posit_ids = unique(c(
      ct.hivneg.nprep$dx.pos.ids,
      ct.hivneg.prep$dx.pos.ids,
      ct.hivpos$dx.pos.ids
    )),
    window_size = get_param(dat, "ct.dx.window"),
    counter_name = "ct.dx.count.window"
  )

  # Both STIs
  dat <- set_attr(
    dat, "sti.dx.count.window",
    get_attr(dat, "gc.dx.count.window") + get_attr(dat, "ct.dx.count.window")
  )

  # Output -----

  # Screening episodes
  sti.screens.hivneg.nprep <- length(unique(c(
    gc.hivneg.nprep$screen.ids, ct.hivneg.nprep$screen.ids
  )))
  sti.dx.hivneg.nprep <- length(unique(c(
    gc.hivneg.nprep$dx.pos.ids, ct.hivneg.nprep$dx.pos.ids
  )))
  sti.screens.hivneg.prep <- length(unique(c(
    gc.hivneg.prep$screen.ids, ct.hivneg.prep$screen.ids
  )))
  sti.dx.hivneg.prep <- length(unique(c(
    gc.hivneg.prep$dx.pos.ids, ct.hivneg.prep$dx.pos.ids
  )))
  sti.screens.hivpos <- length(unique(c(
    gc.hivpos$screen.ids, ct.hivpos$screen.ids
  )))
  sti.dx.hivpos <- length(unique(c(
    gc.hivpos$dx.pos.ids, ct.hivpos$dx.pos.ids
  )))

  dat <- set_epi(dat, "sti.screens.hivneg.nprep", at, sti.screens.hivneg.nprep)
  dat <- set_epi(dat, "sti.dx.hivneg.nprep", at, sti.dx.hivneg.nprep)
  dat <- set_epi(dat, "sti.screens.hivneg.prep", at, sti.screens.hivneg.prep)
  dat <- set_epi(dat, "sti.dx.hivneg.prep", at, sti.dx.hivneg.prep)
  dat <- set_epi(dat, "sti.screens.hivpos", at, sti.screens.hivpos)
  dat <- set_epi(dat, "sti.dx.hivpos", at, sti.dx.hivpos)

  # Any past year screening (proportions)
  hivneg.denom <- which(diag.status == 0 & prepStat == 0)
  hivpos.denom <- which(diag.status == 1 & prepStat == 0)

  ids.gc.screened.py <- which(at - gc.last.screen.time <= 364 / get_param(dat, "time.unit"))
  ids.ct.screened.py <- which(at - ct.last.screen.time <= 364 / get_param(dat, "time.unit"))

  gc.screened.py.hivneg <- length(intersect(ids.gc.screened.py, hivneg.denom)) / length(hivneg.denom)
  gc.screened.py.hivpos <- length(intersect(ids.gc.screened.py, hivpos.denom)) / length(hivpos.denom)
  ct.screened.py.hivneg <- length(intersect(ids.ct.screened.py, hivneg.denom)) / length(hivneg.denom)
  ct.screened.py.hivpos <- length(intersect(ids.ct.screened.py, hivpos.denom)) / length(hivpos.denom)

  dat <- set_epi(dat, "gc.screened.py.hivneg", at, gc.screened.py.hivneg)
  dat <- set_epi(dat, "gc.screened.py.hivpos", at, gc.screened.py.hivpos)
  dat <- set_epi(dat, "ct.screened.py.hivneg", at, ct.screened.py.hivneg)
  dat <- set_epi(dat, "ct.screened.py.hivpos", at, ct.screened.py.hivpos)

  # Update Attributes
  dat <- set_attr(dat, "diag.status.gc", diag.status.gc)
  dat <- set_attr(dat, "diag.status.gc.prep", diag.status.gc.prep)
  dat <- set_attr(dat, "gc.last.screen.time", gc.last.screen.time)
  dat <- set_attr(dat, "gc.last.dxpos.time", gc.last.dxpos.time)

  dat <- set_attr(dat, "diag.status.ct", diag.status.ct)
  dat <- set_attr(dat, "diag.status.ct.prep", diag.status.ct.prep)
  dat <- set_attr(dat, "ct.last.screen.time", ct.last.screen.time)
  dat <- set_attr(dat, "ct.last.dxpos.time", ct.last.dxpos.time)

  return(dat)
}

sti_screen_subfx <- function(elig, test.dist, screen.rate, last.screen.time,
                             r.screen.prob, u.status, r.status,
                             prep = FALSE, prep.start.time, at) {

  if (length(elig) > 1 || prep == TRUE) {
    if (test.dist == "stochastic") {
      u.screen <- rbinom(length(elig), 1, screen.rate)
      u.screen.ids <- elig[u.screen == 1]
    }
    if (test.dist == "deterministic") {
      u.screen.ids <- which((at - last.screen.time) >= (1 / screen.rate) | is.na(last.screen.time))
    }
    r.screen.ids <- ssample(u.screen.ids, length(u.screen.ids) * r.screen.prob)

    if (prep == TRUE && sum(prep.start.time == at, na.rm = TRUE) > 0) {
      u.screen.ids <- union(u.screen.ids, which(prep.start.time == at))
      r.screen.ids <- union(r.screen.ids, which(prep.start.time == at))
    }

    screen.ids <- union(u.screen.ids, r.screen.ids)

    u.dx.pos.ids <- intersect(which(u.status == 1), u.screen.ids)
    r.dx.pos.ids <- intersect(which(r.status == 1), r.screen.ids)

    dx.pos.ids <- union(u.dx.pos.ids, r.dx.pos.ids)

    out <- list(
      screen.ids = screen.ids,
      dx.pos.ids = dx.pos.ids
    )
  } else {
    out <- list(
      screen.ids = NULL,
      dx.pos.ids = NULL
    )
  }

  return(out)
}
