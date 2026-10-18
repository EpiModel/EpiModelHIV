#' @title Treatment Resistance Function Module
#'
#' @description Module function managing the resistance to cART.
#'
#' @inheritParams aging_msm
#'
#' @export
#'
tx_resist_msm <- function(dat, at) {
  drugs <- c("tdf", "ftc")

  # Attributes
  status <- get_attr(dat, "status")
  prep.start.last <- get_attr(dat, "prep.start.last")
  prep.otc.start.last <- get_attr(dat, "prep.otc.start.last")
  prep.any <- get_attr(dat, "prep.any")
  inf.time <- get_attr(dat, "inf.time")

  all_resists_ids <- numeric(0)
  all_prep_resists_ids <- numeric(0)
  all_hiv_resists_ids <- numeric(0)

  for (drug_name in drugs) {
    # on PrEP -> get HIV
    prep_elig_ids <- which(prep.any == 1 & inf.time == at)
    resist_prob <- get_param(dat, paste0(drug_name, ".resist.prep.prob"))
    prep_resist_ids <- prep_elig_ids[runif(length(prep_elig_ids)) < resist_prob]

    # HIV+ start PrEP (cannot start PrEP if HIV diagnosed)
    hiv_elig_ids <- which(status == 1 &
                      (prep.start.last == at | prep.otc.start.last == at))
    resist_prob <- get_param(dat, paste0(drug_name, ".resist.hiv.prob"))
    hiv_resist_ids <- hiv_elig_ids[runif(length(hiv_elig_ids)) < resist_prob]

    # Outputs
    resist_ids <- c(prep_resist_ids, hiv_resist_ids)
    all_resists_ids <- union(all_resists_ids, resist_ids)
    all_prep_resists_ids <- union(all_prep_resists_ids, prep_resist_ids)
    all_hiv_resists_ids <- union(all_hiv_resists_ids, hiv_resist_ids)

    dat <- set_attr(dat, paste0(drug_name, ".resist"), 1, posit_ids = resist_ids)
    dat <- set_epi(dat, paste0(drug_name, ".resist.incid"),
                   at, length(resist_ids))

    if (get_control(dat, "debug")) {
      dat <- set_epi(dat, paste0("dbg_", drug_name, ".resist.prep_elig"),
                     at, length(prep_elig_ids))
      dat <- set_epi(dat, paste0("dbg_", drug_name, ".resist.hiv_elig"),
                     at, length(hiv_elig_ids))
    }
  }

  dat <- set_epi(dat, "any.resist.incid", at, length(all_resists_ids))
  dat <- set_epi(dat, "any.resist.prep.incid", at, length(all_prep_resists_ids))
  dat <- set_epi(dat, "any.resist.hiv.incid", at, length(all_hiv_resists_ids))
  tdf.resist <- get_attr(dat, "tdf.resist")
  ftc.resist <- get_attr(dat, "ftc.resist")
  dat <- set_epi(dat, "any.resist.prev", at, sum(tdf.resist | ftc.resist))
  dat <- set_epi(dat, "tdf.resist.prev", at, sum(tdf.resist))
  dat <- set_epi(dat, "ftc.resist.prev", at, sum(ftc.resist))

  return(dat)
}
