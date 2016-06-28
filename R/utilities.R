
#' @title Apportion Least-Remainder Method
#'
#' @description Apportions a vector of values given a specified frequency
#'              distribution of those values such that the length of the output
#'              is robust to rounding and other instabilities.
#'
#' @param vector.length Length for the output vector.
#' @param values Values for the output vector.
#' @param proportions Proportion distribution with one number for each value. This
#'        must sum to 1.
#' @param shuffled If \code{TRUE}, randomly shuffle the order of the vector.
#'
#' @export
#'
apportion_lr <- function(vector.length, values,
                         proportions, shuffled = FALSE) {

  if (vector.length != round(vector.length)) {
    stop("argument vector.length must be a positive integer")
  }
  if (vector.length <= 0) {
    stop("argument vector.length must be a positive integer")
  }
  if (is.vector(values) == FALSE) {
    stop("argument values must be a vector")
  }
  if (!(length(proportions) == length(values) && round(sum(proportions), 10) == 1) &&
     (!(length(proportions) == length(values) - 1 && round(sum(proportions), 10) <= 1 &&
        round(sum(proportions), 10) >= 0))) {
    stop("error in proportions length or proportions sum")
  }

  if (length(proportions) == length(values) - 1) {
    proportions <- c(proportions, 1 - round(sum(proportions), 10))
  }
  result <- rep(NA, vector.length)
  exp.nums <- proportions * vector.length
  counts <- floor(exp.nums)
  remainders <- exp.nums - counts
  leftovers <- vector.length - sum(counts)
  if (leftovers > 0) {
    additions <- order(remainders, decreasing = TRUE)[1:leftovers]
    counts[additions]   <- counts[additions] + 1
  }
  result <- rep(values, counts)
  if (shuffled == TRUE) {
    result <- sample(result, length(result))
  }

  return(result)
}


#' @title Get Arguments from EpiModel Parameterization Functions
#'
#' @description Returns a list of argument names and values for use for parameter
#'              processing functions.
#'
#' @param formal.args The output of \code{formals(sys.function())}.
#' @param dot.args The output of \code{list(...)}.
#'
#' @export
#'
get_args <- function(formal.args, dot.args){
  p <- list()
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg, pos = parent.frame()))
  }

  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }
  return(p)
}


#' @title Proportionally Reallocate PrEP Adherence Class Probability
#'
#' @description Shifts probabilities from the high-adherence category to the lower
#'              three adherence categories while maintaining the proportional
#'              distribution of those lower categories.
#'
#' @param in.pcp Input vector of length four for the \code{prep.class.prob}
#'        parameter.
#' @param reall The pure percentage points to shift from the high adherence
#'        group to the lower three groups.
#'
#' @export
#'
reallocate_pcp <- function(in.pcp = c(0.211, 0.07, 0.1, 0.619), reall = 0) {

  dist <- in.pcp[1]/sum(in.pcp[1:3])
  dist[2] <- in.pcp[2]/sum(in.pcp[1:3])
  dist[3] <- in.pcp[3]/sum(in.pcp[1:3])

  out.pcp <- rep(NA, 4)
  out.pcp[1:3] <- in.pcp[1:3] - (dist * reall)
  out.pcp[4] <- 1 - sum(out.pcp[1:3])

  return(out.pcp)
}


#' @title Truncate Simulation Time Series
#'
#' @description Left-truncates a simulation epidemiological summary statistics and
#'              network statistics at a specified time step.
#'
#' @param x Object of class \code{netsim}.
#' @param at Time step at which to left-truncate the time series.
#'
#' @details
#' This function would be used when running a follow-up simulation from time steps
#' \code{b} to \code{c} after a burnin period from time \code{a} to \code{b},
#' where the final time window of interest for data analysis is \code{b} to \code{c}
#' only.
#'
#' @export
#'
truncate_sim <- function(x, at) {

  rows <- at:(x$control$nsteps)

  # epi
  x$epi <- lapply(x$epi, function(r) r[rows, ])

  # control settings
  x$control$start <- 1
  x$control$nsteps <- max(seq_along(rows))

  return(x)
}



#' @title Source All Files in a Directory
#'
#' @description Loops over all files in a directory to source them to the
#'              Global Environment.
#'
#' @param path Directory of files to source.
#' @param verbose Print names of sourced files to console.
#' @param ... Additional arguments passed to \code{source}.
#'
#' @export
#'
sourceDir <- function(path, verbose = TRUE, ...) {

  fn <- list.files(path, pattern = "\\.[Rr]$")
  if (length(fn) == 0) {
    stop("No R files in that path", call. = FALSE)
  }

  if (verbose == TRUE) {
    cat("\n Sourced Files in", path,
        "\n -----------------")
  }

  for (nm in fn) {
    if (verbose == TRUE) cat("\n", nm)
    source(file.path(path, nm), ...)
  }
}


nbsdtosize <- function(mu, sd) {
  mu ^ 2 / (sd ^ 2 - mu)
}

get_attr <- function(x, sim = 1) {
  if (is.null(x$attr)) {
    stop("No attr on x")
  } else {
    x$attr[[1]]
  }
}

cut_age <- function(age, breaks = c(0, 29, 39, Inf)) {
  cut(age, breaks = breaks, labels = FALSE)
}

keep.attr <- function(attrList, keep) {
  lapply(attrList, function(x) x[keep])
}
