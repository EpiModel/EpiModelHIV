
#' @title Verbose Module for SHAMP
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
verbose_KTM <- function(x, type, s, at) {

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
          prev <- round(x$epi$prev.poi[at], 3)
          prev.f <- round(x$epi$prev.f.poi[at], 3)
          prev.m <- round(x$epi$prev.m.poi[at], 3)


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
          cat("\nHIV Curr Incidence:", x$epi$incid[at])
          cat("\nHIV Cuml Incidence:", sum(x$epi$incid, na.rm = TRUE))
          cat("\nHIV Prevalence: ", x$epi$i.num[at], " (", prev, ")", sep = "")
          cat("\n------------------------------")
          cat("\nprev.f Prevalence: ", prev.f, sep = "")
          cat("\nprev.m Prevalence: ", prev.m, sep = "")
          cat("\n==============================")

        }
      }
    }


  }
}

