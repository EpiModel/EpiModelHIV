
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
                "\nWht Prev: ", round(x$epi$i.num.W[at] / x$epi$num.W[at], 3),
                "\nBlk Prev: ", round(x$epi$i.num.B[at] / x$epi$num.B[at], 3),
                "\n\n", sep = "", file = fn)
          }
        }
      } else {
        verbose.int <- x$control$verbose.int
        if (verbose.int > 0 && (at %% verbose.int == 0)) {

          nsteps <- x$control$nsteps
          time.unit <- x$param$time.unit
          prev <- round(x$epi$i.prev[at], 3)
          prev.B <- round(x$epi$i.prev.B[at], 3)
          prev.W <- round(x$epi$i.prev.W[at], 3)

          cat("\014")
          cat("\nEpidemic Simulation")
          cat("\n==============================")
          cat("\nTimestep: ", at, "/", nsteps, sep = "")
          cat("\nDay: ", at * time.unit, "/", nsteps * time.unit, sep = "")
          cat("\nYear: ", round((at * time.unit) / 365, 1), "/",
              round((nsteps * time.unit) / 365, 1), sep = "")
          cat("\n------------------------------")
          cat("\nTotal Pop Size:", x$epi$num[at])
          cat("\nBlack Pop Size:", x$epi$num.B[at])
          cat("\nWhite Pop Size:", x$epi$num.W[at])
          cat("\n------------------------------")
          cat("\nCurr Incidence:", x$epi$incid[at])
          cat("\nCuml Incidence:", sum(x$epi$incid, na.rm = TRUE))
          cat("\nTotal Prevalence: ", x$epi$i.num[at], " (", prev, ")", sep = "")
          cat("\nBlack Prevalence: ", x$epi$i.num.B[at], " (", prev.B, ")", sep = "")
          cat("\nWhite Prevalence: ", x$epi$i.num.W[at], " (", prev.W, ")", sep = "")
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
