
#' @title HIV Testing Module
#'
#' @description Module function for HIV diagnostic testing of infected persons.
#'
#' @inheritParams aging.msm
#'
#' @details
#' This testing module supports two testing parameterizations, input via the
#' \code{testing.pattern} parameter: memoryless for stochastc and
#' geometrically-distributed waiting times to test (constant hazard); and interval
#' for deterministic tested after defined waiting time intervals.
#'
#' @return
#' This function returns the \code{dat} object with updated \code{last.neg.test},
#' \code{diag.status} and \code{diag.time} attributes.
#'
#' @keywords module
#' @export
#'
test.msm <- function(dat, at) {

  ## Variables

  # Attributes
  active <- dat$attr$active
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  tt.traj <- dat$attr$tt.traj
  status <- dat$attr$status
  inf.time <- dat$attr$inf.time

  prepStat <- dat$attr$prepStat
  prep.tst.int <- dat$param$prep.tst.int

  # Parameters
  testing.pattern <- dat$param$testing.pattern
  mean.test.B.int <- dat$param$mean.test.B.int
  mean.test.W.int <- dat$param$mean.test.W.int
  twind.int <- dat$param$test.window.int

  ## Process

  if (testing.pattern == "memoryless") {
    elig.B <- which(active == 1 & race == "B" & tt.traj != "NN" &
                    (diag.status == 0 | is.na(diag.status)))
    rates.B <- rep(1/mean.test.B.int, length(elig.B))
    rates.B[which(prepStat[elig.B] == 1)] <- 1/prep.tst.int
    tst.B <- elig.B[rbinom(length(elig.B), 1, rates.B) == 1]

    elig.W <- which(active == 1 & race == "W" & tt.traj != "NN" &
                    (diag.status == 0 | is.na(diag.status)))
    rates.W <- rep(1/mean.test.W.int, length(elig.W))
    rates.W[which(prepStat[elig.W] == 1)] <- 1/prep.tst.int
    tst.W <- elig.W[rbinom(length(elig.W), 1, rates.W) == 1]
  }

  if (testing.pattern == "interval") {
    # Time since last neg test
    tsincelntst <- at - dat$attr$last.neg.test
    tsincelntst[is.na(tsincelntst)] <- at - dat$attr$arrival.time[is.na(tsincelntst)]

    tst.B.nprep <- which(active == 1 & race == "B" & tt.traj != "NN" &
                         (diag.status == 0 | is.na(diag.status)) &
                         tsincelntst >= mean.test.B.int)
    tst.B.prep <- which(active == 1 & race == "B" & tt.traj != "NN" &
                        (diag.status == 0 | is.na(diag.status)) &
                        prepStat == 1 & tsincelntst >= prep.tst.int)
    tst.B <- c(tst.B.nprep, tst.B.prep)

    tst.W.nprep <- which(active == 1 & race == "W" & tt.traj != "NN" &
                         (diag.status == 0 | is.na(diag.status)) &
                         tsincelntst >= mean.test.W.int)
    tst.W.prep <- which(active == 1 & race == "W" & tt.traj != "NN" &
                        (diag.status == 0 | is.na(diag.status)) &
                        prepStat == 1 & tsincelntst >= prep.tst.int)
    tst.W <- c(tst.W.nprep, tst.W.prep)
  }

  tst.pos.B <- tst.B[status[tst.B] == 1 & inf.time[tst.B] <= at - twind.int]
  tst.neg.B <- setdiff(tst.B, tst.pos.B)

  tst.pos.W <- tst.W[status[tst.W] == 1 & inf.time[tst.W] <= at - twind.int]
  tst.neg.W <- setdiff(tst.W, tst.pos.W)

  tst.pos <- c(tst.pos.B, tst.pos.W)
  tst.neg <- c(tst.neg.B, tst.neg.W)


  ## Output

  # Attributes
  dat$attr$last.neg.test[tst.neg] <- at
  dat$attr$diag.status[tst.pos] <- 1
  dat$attr$diag.time[tst.pos] <- at

  ## Summary statistics
  dat$epi$tst.W.inc[at] <- length(tst.W)
  dat$epi$tst.B.inc[at] <- length(tst.B)


  return(dat)
}
