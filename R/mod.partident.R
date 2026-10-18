
#' @title Partner Identification Module
#'
#' @description Module function for identifying partners of incident HIV+ MSM.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This partner identification module handles the elicitation of partners of
#' incident HIV+ MSM for the purposes of HIV testing (\link[EpiModelHIV]{hivtest_msm}).
#'
#' @return
#' This function returns the \code{dat} object with at updated \code{part.ident}
#' attribute detailing time of identification.
#'
#' @keywords module msm
#'
#' @export
#'
partident_msm <- function(dat, at) {

  if (at < get_param(dat, "part.ident.start")) {
    return(dat)
  }

  # Input
  ## Attributes
  active      <- get_attr(dat, "active")
  diag.status <- get_attr(dat, "diag.status")
  diag.time   <- get_attr(dat, "diag.time")
  time.on.tx  <- get_attr(dat, "cuml.time.on.tx")
  race        <- get_attr(dat, "race")
  uid         <- get_unique_ids(dat)

  # Parameters
  part.index.window <- get_param(dat, "part.index.window")
  part.index.degree <- get_param(dat, "part.index.degree")
  part.index.prob   <- get_param(dat, "part.index.prob")

  # Concatenate Partnership Parameters
  part.ident.window <- c(
    get_param(dat, "part.ident.main.window"),
    get_param(dat, "part.ident.casl.window"),
    get_param(dat, "part.ident.ooff.window")
  )
  part.ident.prob <- c(
    get_param(dat, "part.ident.main.prob"),
    get_param(dat, "part.ident.casl.prob"),
    get_param(dat, "part.ident.ooff.prob")
  )

  ## Process

  # Who was diagnosed last step. Identification occurs before testing
  at.prev <- at - 1

  # Identify individuals with recent HIV+ diagnoses
  recent.diag <- at.prev - diag.time
  hivpos.at <- which(diag.status == 1 & recent.diag <= part.index.window)
  if (length(hivpos.at) == 0) {
    dat <- set_epi(dat, "elig_indexes", at, 0)
    dat <- set_epi(dat, "found_indexes", at, 0)
    return(dat)
  }


  # Subset those partnerships where one partner has a recent HIV+ diagnosis.
  hivpos.uid <- uid[hivpos.at]
  plist.temp <- dat$temp$plist
  p1_pos <- plist.temp[, "p1"] %in% hivpos.uid
  p2_pos <- plist.temp[, "p2"] %in% hivpos.uid
  plist.temp <- plist.temp[p1_pos | p2_pos, , drop = FALSE]

  # Remove partnerships from temp partner list that involve inactive individuals
  active1 <- plist.temp[, 1] %in% uid
  active2 <- plist.temp[, 2] %in% uid
  plist.temp <- plist.temp[active1 & active2, , drop = FALSE]

  # Remove partnerships of index patients with mean degree less than part.index.degree
  all.ids <- c(plist.temp[, 1], plist.temp[, 2])
  keep.ids <- which(table(all.ids) >= part.index.degree)
  keep.ids <- as.numeric(names(keep.ids))
  hivpos.uid <- intersect(hivpos.uid, keep.ids)

  # Each indexe is found with probability `part.index.prob`
  found.uid <- hivpos.uid[rbinom(length(hivpos.uid), 1, part.index.prob) == 1]
  p1_found <- plist.temp[, "p1"] %in% found.uid
  p2_found <- plist.temp[, "p2"] %in% found.uid
  plist.temp <- plist.temp[p1_found | p2_found, , drop = FALSE]

  # Construct new identification partnership list based on active partnerships
  # and recent HIV+ diagnoses.
  plist_colnames <- c("partnerUID", "ptype", "start", "stop")

  p <- plist.temp
  # Replace 'found.uid's in `p1` and `p2` by NA
  p[, 1:2][ p[, 1:2] %in% found.uid ] <- NA
  # If NA in `p1`, put `p2`. `p1` now only contains partners
  p[, "p1"] <- ifelse(is.na(p[, "p1"]), p[, "p2"], p[, "p1"])
  # In case of a partnership where both are identified, `p1` is still NA.
  # We also remove `p2` which is now useless
  p <- p[!is.na(p[, "p1"]), -2, drop = FALSE]

  colnames(p) <- plist_colnames
  plist.ident <- p

  # Set trackers for the indexes
  elig_indexes <- length(unique(hivpos.uid))
  dat <- set_epi(dat, "elig_indexes", at, elig_indexes)

  found_indexes <- length(unique(found.uid))
  dat <- set_epi(dat, "found_indexes", at, found_indexes)

  if (nrow(plist.ident) == 0) {
    return(dat)
  }

  # Subset plist.ident by partners who are in a qualified relationship,
  # defined by relationship time ("stop" is NA for ongoing relationships)
  rel.time <- at - plist.ident[, "stop"]
  rel.elig <- which(
    is.na(plist.ident[, "stop"])
    | rel.time <= part.ident.window[plist.ident[, "ptype"]]
  )
  plist.ident <- plist.ident[rel.elig, , drop = FALSE]

  if (nrow(plist.ident) == 0) {
    return(dat)
  }

  ## Append PartnerID to the amended plist.ident for effecient indexing
  plist.ident <- cbind(
    plist.ident,
    "partnerPID" = get_posit_ids(dat, plist.ident[, "partnerUID"])
  )

  ## Identification of partners of recent HIV+ diagnosed individuals:
  partner.pos <- unique(plist.ident[, "partnerPID"])

  if (length(partner.pos) == 0) {
    return(dat)
  }

  ## Attempt to identify eligible partners

  # Pull probability of identification based on partnership type
  part.ident.prob.type <- part.ident.prob[plist.ident[, "ptype"]]

  # Random binomial to determine if identification is successful
  part.identified <- rbinom(nrow(plist.ident), 1, part.ident.prob.type)

  # Subset temporary partnership list by those who were identified
  plist.ident <- plist.ident[part.identified == 1, , drop = FALSE]

  # Note identified partners
  identified <- unique(plist.ident[, "partnerPID"])

  # If any identified individuals have any time on ART's, update part.art
  # attribute which tracks identified partners who have been on ART's.
  # dat$attr$part.art[time.on.tx[identified] > 0] <- 1

  # Update attributes
  dat$attr$part.ident[identified] <- at
  dat$attr$part.ident.counter[identified] <- ifelse(
    is.na(dat$attr$part.ident.counter[identified]),
    1, dat$attr$part.ident.counter[identified] + 1
  )

  return(dat)
}
