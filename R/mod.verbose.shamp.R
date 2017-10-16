
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
verbose_shamp <- function(x, type, s, at) {

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
          prev.B.f <- round(x$epi$i.prev.B.f[at], 3)
          prev.BI.f <- round(x$epi$i.prev.BI.f[at], 3)
          prev.H.f <- round(x$epi$i.prev.H.f[at], 3)
          prev.HI.f <- round(x$epi$i.prev.HI.f[at], 3)
          prev.W.f <- round(x$epi$i.prev.W.f[at], 3)
          
          prev.B.m <- round(x$epi$i.prev.B.m[at], 3)
          prev.BI.m <- round(x$epi$i.prev.BI.m[at], 3)
          prev.H.m <- round(x$epi$i.prev.H.m[at], 3)
          prev.HI.m <- round(x$epi$i.prev.HI.m[at], 3)
          prev.W.m <- round(x$epi$i.prev.W.m[at], 3)
          
          prev.B.msf <- round(x$epi$i.prev.B.msf[at], 3)
          prev.BI.msf <- round(x$epi$i.prev.BI.msf[at], 3)
          prev.H.msf <- round(x$epi$i.prev.H.msf[at], 3)
          prev.HI.msf <- round(x$epi$i.prev.HI.msf[at], 3)
          prev.W.msf <- round(x$epi$i.prev.W.msf[at], 3)
          
          prev.B.msm <- round(x$epi$i.prev.B.msm[at], 3)
          prev.BI.msm <- round(x$epi$i.prev.BI.msm[at], 3)
          prev.H.msm <- round(x$epi$i.prev.H.msm[at], 3)
          prev.HI.msm <- round(x$epi$i.prev.HI.msm[at], 3)
          prev.W.msm <- round(x$epi$i.prev.W.msm[at], 3)
          
          prev.B.msmf <- round(x$epi$i.prev.B.msmf[at], 3)
          prev.BI.msmf <- round(x$epi$i.prev.BI.msmf[at], 3)
          prev.H.msmf <- round(x$epi$i.prev.H.msmf[at], 3)
          prev.HI.msmf <- round(x$epi$i.prev.HI.msmf[at], 3)
          prev.W.msmf <- round(x$epi$i.prev.W.msmf[at], 3)

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
          cat("\nprev.B.f Prevalence: ", prev.B.f, sep = "")
          cat("\nprev.BI.f Prevalence: ", prev.BI.f, sep = "")
          cat("\nprev.H.f Prevalence: ", prev.H.f, sep = "")
          cat("\nprev.HI.f Prevalence: ", prev.HI.f, sep = "")
          cat("\nprev.W.f Prevalence: ", prev.W.f, sep = "")
          cat("\nprev.B.m Prevalence: ", prev.B.m, sep = "")
          cat("\nprev.BI.m Prevalence: ", prev.BI.m, sep = "")
          cat("\nprev.H.m Prevalence: ", prev.H.m, sep = "")
          cat("\nprev.HI.m Prevalence: ", prev.HI.m, sep = "")
          cat("\nprev.W.m Prevalence: ", prev.W.m, sep = "")
          cat("\nprev.B.msf Prevalence: ", prev.B.msf, sep = "")
          cat("\nprev.BI.msf Prevalence: ", prev.BI.msf, sep = "")
          cat("\nprev.H.msf Prevalence: ", prev.H.msf, sep = "")
          cat("\nprev.HI.msf Prevalence: ", prev.HI.msf, sep = "")
          cat("\nprev.W.msf Prevalence: ", prev.W.msf, sep = "")
          cat("\nprev.B.msm Prevalence: ", prev.B.msm, sep = "")
          cat("\nprev.BI.msm Prevalence: ", prev.BI.msm, sep = "")
          cat("\nprev.H.msm Prevalence: ", prev.H.msm, sep = "")
          cat("\nprev.HI.msm Prevalence: ", prev.HI.msm, sep = "")
          cat("\nprev.W.msm Prevalence: ", prev.W.msm, sep = "")
          cat("\nprev.B.msmf Prevalence: ", prev.B.msmf, sep = "")
          cat("\nprev.BI.msmf Prevalence: ", prev.BI.msmf, sep = "")
          cat("\nprev.H.msmf Prevalence: ", prev.H.msmf, sep = "")
          cat("\nprev.HI.msmf Prevalence: ", prev.HI.msmf, sep = "")
          cat("\nprev.W.msmf Prevalence: ", prev.W.msmf, sep = "")
          cat("\n==============================")

        }
      }
    }


  }
}

