
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
#' @keywords internal
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
