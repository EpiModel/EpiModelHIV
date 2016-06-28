
#' @title Update Role Class in Main and Casual Partnerships
#'
#' @description Module function for updating act class in main and casual
#'              partnerships based on probabilities of transition.
#'
#' @inheritParams aging_msm
#'
#' @return
#' This function updates the individual-level attribute \code{role.class} on
#' \code{dat$attr}.
#'
#' @keywords module msm
#' 
#' @export
#'
update_roleclass_msm <- function(dat, at) {

  role.trans.matrix <- dat$param$role.trans.matrix
  if (sum(colSums(role.trans.matrix) != 1) > 0) {
    stop("Column sums in argument role.trans.matrix must all equal 1.")
  }
  old.role.class <- dat$attr$role.class

  new.role.class <- sapply(1:length(old.role.class),
                           function(x) {
                             sample(c("I", "V", "R"), size = 1,
                                     prob = role.trans.matrix[, 1] *
                                               (old.role.class[x] == "I") +
                                            role.trans.matrix[, 2] *
                                               (old.role.class[x] == "V") +
                                            role.trans.matrix[, 3] *
                                               (old.role.class[x] == "R"))
                           })

  dat$attr$role.class <- new.role.class

  return(dat)
}
