
#' @export
prep.mard <- function(dat, at) {

  ## Variables
  active <- dat$attr$active
  status <- dat$attr$status
  prepElig <- dat$attr$prepElig
  prepEligTime <- dat$attr$prepEligTime
  prepStat <- dat$attr$prepStat
  prepStartTime <- dat$attr$prepStartTime
  prepEver <- dat$attr$prepEver
  prepClass <- dat$attr$prepClass

  prep.start <- dat$param$prep.start
  prep.elig.model <- dat$param$prep.elig.model
  prep.coverage <- dat$param$prep.coverage
  prep.cov.method <- dat$param$prep.cov.method
  prep.cov.rate <- dat$param$prep.cov.rate
  prep.class.prob <- dat$param$prep.class.prob


  ## Eligibility
  if (prep.elig.model == "all") {
    idsElig <- which(active == 1 & status == 0 & is.na(prepElig))
  }
  prepElig[idsElig] <- 1
  prepEligTime[idsElig] <- at


  ## Initiation
  if (prep.cov.method == "curr") {
    prepCov <- sum(prepStat == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  }
  if (prep.cov.method == "ever") {
    prepCov <- sum(prepEver == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  }
  prepCov <- ifelse(is.nan(prepCov), 0, prepCov)

  idsEligSt <- which(prepElig == 1)
  nEligSt <- length(idsEligSt)

  nStart <- max(0, min(nEligSt, round((prep.coverage - prepCov) *
                                             sum(prepElig == 1, na.rm = TRUE))))
  idsStart <- NULL
  if (nStart > 0 & at >= prep.start) {
    if (prep.cov.rate >= 1) {
      idsStart <- ssample(idsEligSt, nStart)
    } else {
      idsStart <- idsEligSt[rbinom(nStart, 1, prep.cov.rate) == 1]
    }
  }

  ## Attributes
  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1
    prepStartTime[idsStart] <- at
    prepEver[idsStart] <- 1
    prepClass[idsStart] <- sample(c("l", "m", "h"), length(idsStart),
                                  TRUE, prep.class.prob)
  }

  dat$attr$prepElig <- prepElig
  dat$attr$prepEligTime <- prepEligTime
  dat$attr$prepStat <- prepStat
  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepEver <- prepEver
  dat$attr$prepClass <- prepClass

  # Summary Statistics
  dat$epi$prepCov[at] <- prepCov
  dat$epi$prepStart[at] <- length(idsStart)

  return(dat)
}
