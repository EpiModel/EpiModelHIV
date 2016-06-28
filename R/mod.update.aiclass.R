
#' @title Update Role Class in One-Off Partnerships
#'
#' @description Module function for updating act class in one-off partnerships
#'              based on probabilities of transition.
#'
#' @inheritParams aging_msm
#'
#' @return
#' This function updates the individual-level attribute \code{inst.ai.class} on
#' \code{dat$attr}.
#'
#' @keywords module msm
#' 
#' @export
#'
update_aiclass_msm <- function(dat, at) {

  inst.trans.matrix <- dat$param$inst.trans.matrix
  if (sum(colSums(dat$param$inst.trans.matrix) != 1) > 0) {
    stop("Column sums in argument inst.trans.matrix must all equal 1.")
  }
  old.inst.ai.class <- dat$attr$inst.ai.class

  new.inst.ai.class <- sapply(1:length(old.inst.ai.class),
                              function(x) {
                                which(rmultinom(1, 1,
                                                inst.trans.matrix[, old.inst.ai.class[x]]) == 1)
                              })

  dat$attr$inst.ai.class <- new.inst.ai.class

  return(dat)
}
