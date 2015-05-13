
#' @method print nwstats
#' @export
print.nwstats <- function(x, ...) {

  cat("Network Statistics Summary")
  cat("\n==========================")
  cat("\nMean degree frequencies by race")
  cat("\nBlack (0/1)", rowSums(x$deg.mp.B))
  cat("\nWhite (0/1):", rowSums(x$deg.mp.W))
  cat("\n\nCasual degree frequencies by race")
  cat("\nBlack (0/1/2)", colSums(x$deg.mp.B))
  cat("\nWhite (0/1/2):", colSums(x$deg.mp.W))

  cat("\n\nMain network model target statistics:\n")
  cat(x$stats.m)

  cat("\n\nCasual network model target statistics:\n")
  cat(x$stats.p)

  cat("\n\nInstant network model target statistics:\n")
  cat(x$stats.i)

  cat("\n\nMain Model ")
  print(x$coef.diss.m)
  cat("\n\nCasual Model ")
  print(x$coef.diss.p)

}
