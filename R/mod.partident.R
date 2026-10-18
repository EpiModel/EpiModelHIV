
#' @title Partner Identification Module
#'
#' @description Module function for identifying partners of incident HIV+ MSM.
#' #description Module function for identifying partners of incident and prevalent HIV+ MSM.
#'
#' @inheritParams aging_msm
#'
#' @details
#' This partner identification module handles the elicitation of partners of
#' incident HIV+ MSM for the purposes of HIV testing ([`hivtest_msm`]).
#'
#' This partner identification module handles the elicitation of partners of
#' incident and prevalent HIV+ MSM for the purposes of HIV testing (\link[EpiModelHIV]{hivtest_msm}).
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
  # Attributes
  diag.status   <- get_attr(dat, "diag.status")
  diag.time     <- get_attr(dat, "diag.time")
  last.pos.test <- get_attr(dat, "last.pos.test")

  # Parameters
  part.index.window.int    <- get_param(dat, "part.index.window.int")
  part.index.degree        <- get_param(dat, "part.index.degree")
  part.index.prob          <- get_param(dat, "part.index.prob")
  part.ppindex.prob        <- get_param(dat, "part.ppindex.prob")

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


## PS Process
  # Who was diagnosed last step. Identification occurs before testing
  at.prev <- at - 1

  #A. 1st generation partners to index cases
  #-----------------------------------------------------------------------------
  # set the epi to 0 (prevent NAs in early return cases),
  # they will be overwritten if needed.

  dat <- set_epi(dat, "recent.newdiagn", at, 0)
  dat <- set_epi(dat, "elig.indexes.nd", at, 0)
  dat <- set_epi(dat, "found.indexes.nd", at, 0)

  dat <- set_epi(dat, "recent.ppretested", at, 0)
  dat <- set_epi(dat, "elig.indexes.pp", at, 0)
  dat <- set_epi(dat, "found.indexes.pp", at, 0)
  dat <- set_epi(dat, "found.indexes.pp.un", at, 0)

  dat <- set_epi(dat, "elig.partners", at, 0)
  dat <- set_epi(dat, "found.partners", at, 0)
  dat <- set_epi(dat, "negunkPart.indexes", at, 0)
  dat <- set_epi(dat, "posPart.indexes", at, 0)

  dat <- set_epi(dat, "elig.indexes.posPart", at, 0)
  dat <- set_epi(dat, "found.indexes.posPart", at, 0)

  dat <- set_epi(dat, "elig.partners.gen2", at, 0)
  dat <- set_epi(dat, "found.partners.gen2", at, 0)
  dat <-  set_epi(dat, "posPart.gen2", at, 0)
  dat <-  set_epi(dat, "negunkPart.gen2", at, 0)

  dat <- set_epi(dat, "found.partners.all", at, 0)


  # identify indexes (new diagnoses and retested prev pos)
  recent.diag <- at.prev - diag.time
  recent.retests <- at.prev - last.pos.test

  hivpos.at <- which(diag.status == 1 & recent.diag <= part.index.window.int)
  retests.at <- which(diag.status == 1 & recent.retests <=  part.index.window.int)

  dat <-  set_epi(dat, "recent.newdiagn", at, length(hivpos.at))
  dat <-  set_epi(dat, "recent.ppretested", at, length(retests.at))

      if (length(hivpos.at) + length(retests.at) == 0) {
        return(dat)
      }

  # Get partnerships involving indexes
  hivpos_pel <- get_partners(dat, hivpos.at, only.active.nodes = TRUE)
  retests_pel <- get_partners(dat, retests.at, only.active.nodes = TRUE)

  # Remove partnerships with mean degree < part.index.degree
  keep.ids.nd <- which(table(hivpos_pel[["index"]]) >= part.index.degree)
  keep.ids.nd <- as.integer(names(keep.ids.nd))
  hivpos.uid <- intersect(keep.ids.nd, hivpos_pel[["index"]])

  keep.ids.pp <- which(table(retests_pel[["index"]]) >= part.index.degree)
  keep.ids.pp <- as.integer(names(keep.ids.pp))
  retests.uid <- intersect(keep.ids.pp, retests_pel[["index"]])

  dat <- set_epi(dat, "elig.indexes.nd", at, length(hivpos.uid))
  dat <- set_epi(dat, "elig.indexes.pp", at, length(retests.uid))

  # Prob. of index participating in PS interview
  found.uid.nd <- hivpos.uid[runif(length(hivpos.uid)) < part.index.prob]
  hivpos_pel <- hivpos_pel[hivpos_pel[["index"]] %in% found.uid.nd, ]

  found.uid.pp <- retests.uid[runif(length(retests.uid)) < part.ppindex.prob]
  retests_pel <- retests_pel[retests_pel[["index"]] %in% found.uid.pp, ]

  dat <- set_epi(dat, "found.indexes.nd", at, length(found.uid.nd))
  dat <- set_epi(dat, "found.indexes.pp", at, length(found.uid.pp))

      #update found pp indexes (for tx check and art restart)
      if (length(found.uid.pp)>0) {
        found.uid.pp <- get_posit_ids(dat, unique(found.uid.pp))
        dat <- set_attr(dat, "foundPrevPos.ident", at, posit_ids = found.uid.pp)
        dat <- set_epi(dat, "found.indexes.pp.un", at, length(found.uid.pp))
      }

  #merge indexes
  hivpos_pel <- rbind(hivpos_pel, retests_pel)

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

  # eligible and found partners
  elig.partners <- unique(hivpos_pel[["partner"]])

  dat <- set_epi(dat, "elig.partners", at, length(elig.partners))

      if (length(elig.partners) == 0) {
        return(dat)
      }

  ## Attempt to identify wave 1 partners
  part.ident.prob.type <- part.ident.prob[hivpos_pel[["network"]]]
  found.rel <- runif(nrow(hivpos_pel)) < part.ident.prob.type
  hivpos_pel <- hivpos_pel[found.rel, ]

  identified <- unique(hivpos_pel[["partner"]])
  identified <- get_posit_ids(dat, identified)

  dat <- set_epi(dat, "found.partners", at, length(identified))

      if (length(identified) == 0) {
        return(dat)
      }


  # get HIV status of found wave 1 partners
  diagstatus.ident <- diag.status[identified]
  hivposPart.at <- identified[diagstatus.ident == 1 & !is.na(diagstatus.ident)]
  negunkPart.at <- identified[diagstatus.ident == 0 | is.na(diagstatus.ident)]

  dat <-  set_epi(dat, "posPart.indexes", at, length(hivposPart.at))
  dat <-  set_epi(dat, "negunkPart.indexes", at, length(negunkPart.at))



  #B. 2nd generation partners (partners to found HIV+ partners)
  #-----------------------------------------------------------------------------
  if (at < get_param(dat, "second.genps.start")) {
     identified.gen2 <- NULL
  }

  if (at >= get_param(dat, "second.genps.start")) {

    if (length(hivposPart.at) > 0) {

      # get partnerships and remove partnerships with mean degree < part.index.degree
      posPart_pel <- get_partners(dat, hivposPart.at, only.active.nodes = TRUE)

      keep.ids.posPart <- which(table(posPart_pel[["index"]]) >= part.index.degree)
      keep.ids.posPart <- as.integer(names(keep.ids.posPart))
      posPart.uid <- intersect(keep.ids.posPart, posPart_pel[["index"]])

      # Probability of posPart-index participating in PS interview
      # Since this partner is already found, update partnership list
      posPart_pel <- posPart_pel[posPart_pel[["index"]] %in% posPart.uid, ]

      dat <- set_epi(dat, "elig.indexes.posPart", at, length(posPart.uid))

        if (nrow(posPart_pel) > 0) {

          # Subset plist.ident by partners who are in a qualified relationship,
          # defined by relationship time ("stop" is NA for ongoing relationships)
          rel.age <- at - posPart_pel[["stop"]]
          rel.age <- ifelse(is.na(rel.age), 0, rel.age)
          age.thresh <- part.ident.int[posPart_pel[["network"]]]
          rel.elig <- rel.age <= age.thresh
          posPart_pel <- posPart_pel[rel.elig, ]

            if (nrow(posPart_pel) > 0) {

              # get eligible Wave 2 partners
              elig.partners.gen2 <- unique(posPart_pel[["partner"]])
              dat <- set_epi(dat, "elig.partners.gen2", at, length(elig.partners.gen2))

                if (length(elig.partners.gen2) > 0) {

                  # get found wave 2 partners
                  part.ident.prob.type <- part.ident.prob[posPart_pel[["network"]]]
                  found.rel.gen2 <- runif(nrow(posPart_pel)) < part.ident.prob.type
                  posPart_pel <- posPart_pel[found.rel.gen2, ]

                  identified.gen2 <- unique(posPart_pel[["partner"]])
                  identified.gen2 <- get_posit_ids(dat, identified.gen2)
                  dat <- set_epi(dat, "found.partners.gen2", at, length(identified.gen2))

                  # get HIV status of found wave 2 partners
                  diagstatus.ident.gen2 <- diag.status[identified.gen2]
                  hivposPartgen2.at <- identified.gen2[diagstatus.ident.gen2 == 1 &
                                                     !is.na(diagstatus.ident.gen2)]
                  negunkPartgen2.at <- identified.gen2[diagstatus.ident.gen2 == 0 |
                                                     is.na(diagstatus.ident.gen2)]

                  dat <-  set_epi(dat, "posPart.gen2", at, length(hivposPartgen2.at))
                  dat <-  set_epi(dat, "negunkPart.gen2", at, length(negunkPartgen2.at))


                  #merge wave 1 and 2 partners
                  identified <- unique(c(identified, identified.gen2))
                }
            }
        }
    }
  }


  dat <- set_epi(dat, "found.partners.all", at, length(identified))


  # Update attributes
  dat <- set_attr(dat, "part.ident", at, posit_ids = identified)
  counter <- get_attr(dat, "part.ident.counter", posit_ids = identified)
  counter <- ifelse(is.na(counter), 1, counter + 1)
  dat <- set_attr(dat, "part.ident.counter", counter, posit_ids = identified)


  return(dat)

}
