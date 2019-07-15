
#' @title Verbose Module
#'
#' @description Prints simulation model progress within the time loop.
#'
#' @param x If the \code{type} is "startup", then an object of class
#'        \code{control.net}, otherwise the all master data object in \code{netsim}
#'        simulations.
#' @param type Progress type, either of "startup" for starting messages before
#'        all simulations, or "progress" for time step specific messages.
#' @param s Current simulation number, if type is "progress"
#' @param at Current time step, if type is "progress"
#'
#' @details
#' In interactive mode, this module function prints out a standard set of
#' demographic and epidemiologic summary statistics to the console window. In
#' non-interactive, batch mode these are saved onto \code{.txt} files in a
#' \code{verb/} subdirectory. This subdirectory will be created if it does not
#' exist.
#'
#' @keywords module msm
#'
#' @export
#'
verbose_msm <- function(x, type, s, at) {

  if (type == "startup") {

  }

  if (type == "progress") {
    if (x$control$verbose == TRUE) {
      if (!interactive()) {
        if (!is.null(x$control$ncores) && x$control$ncores > 1) {
          if (x$control$verbose.int > 0 && (at %% x$control$verbose.int == 0)) {
            simno <- x$control$simno
            currsim <- x$control$currsim
            if (file.exists("verb/") == FALSE) {
              dir.create("verb/")
            }
            fn <- paste0("verb/sim", simno, ".s", currsim, ".txt")
            cat("SIMNO ", paste(simno, currsim, sep = "."),
                "\n====================",
                "\nStep: ", at, " (", round(at/x$control$nsteps, 2), ")",
                "\nPop Size: ", x$epi$num[at],
                "\nTot Prev: ", round(x$epi$i.num[at] / x$epi$num[at], 3),
                "\n\n", sep = "", file = fn)
          }
        }
      } else {
        verbose.int <- x$control$verbose.int
        if (verbose.int > 0 && (at %% verbose.int == 0)) {

          nsteps <- x$control$nsteps
          time.unit <- x$param$time.unit
          prev <- round(x$epi$i.prev[at], 3)
          prev.rgc <- round(x$epi$prev.rgc[at], 3)
          prev.ugc <- round(x$epi$prev.ugc[at], 3)
          prev.rct <- round(x$epi$prev.rct[at], 3)
          prev.uct <- round(x$epi$prev.uct[at], 3)
          incid.hiv <- round(x$epi$incid[at], 3)
          incid.gc <- round(x$epi$incid.gc[at], 3)
          incid.ct <- round(x$epi$incid.ct[at], 3)
          # incid.syph <- round(x$epi$incid.syph[at], 3)
          # prev.syph <- round(x$epi$prev.syph[at], 3)
          # prev.pssyph <- round(x$epi$prev.primsecosyph[at], 3)
          # prev.hiv.pssyph <- round(x$epi$prev.hiv.primsecosyphpos[at], 3)

          cat("\014")
          cat("\nEpidemic Simulation")
          cat("\n==============================")
          cat("\nTimestep: ", at, "/", nsteps, sep = "")
          cat("\nDay: ", at * time.unit, "/", nsteps * time.unit, sep = "")
          cat("\nYear: ", round((at * time.unit) / 365, 1), "/",
              round((nsteps * time.unit) / 365, 1), sep = "")
          cat("\n------------------------------")
          cat("\nTotal Pop Size:", x$epi$num[at])
          cat("\n------------------------------")
          #cat("\nHIV Curr Incidence:", x$epi$incid[at])
          #cat("\nHIV Cuml Incidence:", sum(x$epi$incid, na.rm = TRUE))
          cat("\nHIV Prevalence: ", x$epi$i.num[at], " (", prev, ")", sep = "")
          cat("\nrGC Prevalence: ", prev.rgc, sep = "")
          cat("\nuGC Prevalence: ", prev.ugc, sep = "")
          cat("\nrCT Prevalence: ", prev.rct, sep = "")
          cat("\nuCT Prevalence: ", prev.uct, sep = "")
          # cat("\nSyphilis Prevalence: ", prev.syph, sep = "")
          # cat("\nP and S Syphilis Prevalence: ", prev.pssyph, sep = "")
          # cat("\nHIV Prevalence in P and S Syphilis: ", prev.hiv.pssyph, sep = "")
          cat("\n------------------------------")
          cat("\nHIV Incidence: ", incid.hiv, sep = "")
          cat("\nGC Incidence: ", incid.gc, sep = "")
          cat("\nCT Incidence: ", incid.ct, sep = "")
          # cat("\nSyph Incidence: ", incid.syph, sep = "")


          cat("\n==============================")

        }
      }
    }


  }
}


#' @title Verbose Module
#'
#' @description Prints simulation model progress within the time loop.
#'
#' @param x If the \code{type} is "startup", then an object of class
#'        \code{control.net}, otherwise the all master data object in \code{netsim}
#'        simulations.
#' @param type Progress type, either of "startup" for starting messages before
#'        all simulations, or "progress" for time step specific messages.
#' @param s Current simulation number, if type is "progress"
#' @param at Current time step, if type is "progress"
#'
#' @details
#' In interactive mode, this module function prints out a standard set of
#' demographic and epidemiologic summary statistics to the console window. In
#' non-interactive, batch mode these are saved onto \code{.txt} files in a
#' \code{verb/} subdirectory. This subdirectory will be created if it does not
#' exist.
#'
#' @keywords module het
#'
#' @export
#'
verbose_het <- function(x, type, s, at) {

  if (type == "startup") {

  }

  if (type == "progress") {
    if (x$control$verbose == TRUE) {

      if (!interactive()) {
        if (!is.null(x$control$ncores) && x$control$ncores > 1) {
          if (x$control$verbose.int > 0 && (at %% x$control$verbose.int == 0)) {
            simno <- x$control$simno
            currsim <- x$control$currsim
            if (file.exists("verb/") == FALSE) {
              dir.create("verb/")
            }
            fn <- paste0("verb/sim", simno, ".s", currsim, ".txt")
            cat("SIMNO ", paste(simno, currsim, sep = "."),
                "\n====================",
                "\nStep: ", at,
                "\nPop Size: ", x$epi$num[at],
                "\nPrev: ", x$epi$i.num[at],
                "\n\n", sep = "", file = fn)
          }
        }
      } else {
        verbose.int <- x$control$verbose.int
        if (verbose.int > 0 && (at %% verbose.int == 0)) {

          nsteps <- x$control$nsteps
          tunit <- x$param$time.unit
          prev <- round(x$epi$i.num[at] / x$epi$num[at], 3)

          cat("\014")
          cat("\nDisease Simulation")
          cat("\n==============================")
          cat("\nTimestep: ", at, "/", nsteps, sep = "")
          cat("\nDay: ", at*tunit, "/", nsteps * tunit, sep = "")
          cat("\nYear: ", round((at*tunit)/365, 1), "/",
              round((nsteps * tunit)/365, 1), sep = "")
          cat("\n------------------------------")
          cat("\nCurr Incidence:", x$epi$si.flow[at])
          cat("\nCuml Incidence:", sum(x$epi$si.flow, na.rm = TRUE))
          cat("\nTotal Prevalence: ", x$epi$i.num[at], " (", prev, ")", sep = "")
          cat("\n------------------------------")
          cat("\nCurr Pop Size:", x$epi$num[at])
          cat("\n  Births:", x$epi$b.flow[at])
          cat("\n  Deaths (Sus):", x$epi$ds.flow[at])
          cat("\n  Deaths (Inf):", x$epi$di.flow[at])
          cat("\n==============================")

        }
      }
    }

  }
}
