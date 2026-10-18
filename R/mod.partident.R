
#' @title Partner Identification Module
#'
#' @description Module function for identifying partners of incident HIV+ MSM.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This partner identification module handles the elicitation of partners of
#' incident HIV+ MSM for the purposes of HIV testing ([`hivtest_msm`]).
#'
#' @return
#' This function returns the `dat` object with at updated `part.ident` attribute detailing time of
#' identification.
#'
#' @export
#'
partident_msm <- function(dat, at) {
  if (at < get_param(dat, "part.ident.start")) {
    return(dat)
  }

  # Input
  ## Attributes
  diag.status <- get_attr(dat, "diag.status")
  diag.time   <- get_attr(dat, "diag.time")

  # Parameters
  part.index.window.int    <- get_param(dat, "part.index.window.int")
  part.index.degree        <- get_param(dat, "part.index.degree")
  part.index.prob          <- get_param(dat, "part.index.prob")

  # Concatenate Partnership Parameters
  part.ident.int <- c(
    get_param(dat, "part.ident.main.window.int"),
    get_param(dat, "part.ident.casl.window.int"),
    get_param(dat, "part.ident.ooff.window.int")
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

  # set the epi to 0, they will be overwritten if needed. (prevent NAs in early
  # return cases)
  dat <- set_epi(dat, "elig_indexes", at, 0)
  dat <- set_epi(dat, "found_indexes", at, 0)
  dat <- set_epi(dat, "elig_partners", at, 0)
  dat <- set_epi(dat, "found_partners", at, 0)

  hivpos.at <- which(diag.status == 1 & recent.diag <= part.index.window.int)
  if (length(hivpos.at) == 0) {
    return(dat)
  }

  # Partnerships involving an HIV positive as index
  hivpos_pel <- get_partners(dat, hivpos.at, only.active.nodes = TRUE)

  # Remove partnerships of index patients with mean degree < part.index.degree
  keep.ids <- which(table(hivpos_pel[["index"]]) >= part.index.degree)
  keep.ids <- names(keep.ids)
  hivpos.uid <- intersect(keep.ids, hivpos_pel[["index"]])

  # Each index is found with probability `part.index.prob`
  found.uid <- hivpos.uid[runif(length(hivpos.uid)) < part.index.prob]
  hivpos_pel <- hivpos_pel[hivpos_pel[["index"]] %in% found.uid, ]

  # Set trackers for the indexes
  dat <- set_epi(dat, "elig_indexes", at, length(hivpos.uid))
  dat <- set_epi(dat, "found_indexes", at, length(found.uid))

  if (nrow(hivpos_pel) == 0) {
    return(dat)
  }

  # Subset plist.ident by partners who are in a qualified relationship,
  # defined by relationship time ("stop" is NA for ongoing relationships)
  rel.age <- at - hivpos_pel[["stop"]]
  rel.age <- ifelse(is.na(rel.age), 0, rel.age)
  age.thresh <- part.ident.int[hivpos_pel[["network"]]]
  rel.elig <- rel.age <= age.thresh
  hivpos_pel <- hivpos_pel[rel.elig, ]

  if (nrow(hivpos_pel) == 0) {
    return(dat)
  }

  # Set trackers for the eligible partners
  elig.partners <- unique(hivpos_pel[["partner"]])
  dat <- set_epi(dat, "elig_partners", at, length(elig.partners))

  if (length(elig.partners) == 0) {
    return(dat)
  }

  ## Attempt to identify eligible partners
  part.ident.prob.type <- part.ident.prob[hivpos_pel[["network"]]]
  found.rel <- runif(nrow(hivpos_pel)) < part.ident.prob.type
  hivpos_pel <- hivpos_pel[found.rel, ]

  identified <- unique(hivpos_pel[["partner"]])
  identified <- get_posit_ids(dat, identified)

  dat <- set_epi(dat, "found_partners", at, length(identified))

  # Update attributes
  dat <- set_attr(dat, "part.ident", at, posit_ids = identified)

  counter <- get_attr(dat, "part.ident.counter", posit_ids = identified)
  counter <- ifelse(is.na(counter), 1, counter + 1)
  dat <- set_attr(dat, "part.ident.counter", counter, posit_ids = identified)

  return(dat)
}
