
#' @title Position Module for msm and heterosexual contacts.
#'
#' @description Module function for establishing sexual role or position in each
#'              act on the discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' The sexual role within each act is determined by each nodes "role identity"
#' as exclusively receptive, exclusively insertive, or versatile. For heterosexual dyads
#' the females node is always recpetive and the male is always incertive.  This function
#' determines whether the infected or the susceptible partner is the insertive
#' partner for that act. For male-male dyads the first two role identity types, that is
#' deterministic based on identity. For male-male versatile-versatile pairs, this is
#' determined stochastically for each act.
#'
#' @return
#' This function returns the updated discordant edgelist with a \code{ins}
#' attribute for values of whether the infected node is insertive or the
#' susceptible node is insertive for that act.
#'
#' @keywords module SHAMP het msm
#' 
#' @export
#'
position_shamp <- function(dat, at) {

  ## Variables
  al <- dat$temp$al
  if (nrow(al) == 0) {
    return(dat)
  }

  status <- dat$attr$status
  sex <- dat$attr$sex
  dal <- al[which(status[al[, 1]] == 1 & status[al[, 2]] == 0), ]
  dat$temp$al <- NULL

  role.class <- dat$attr$role.class
  ins.quot <- dat$attr$ins.quot
  race <- dat$attr$race


  vv.iev.B.prob <- dat$param$vv.iev.B.prob
  vv.iev.BI.prob <- dat$param$vv.iev.BI.prob
  vv.iev.H.prob <- dat$param$vv.iev.H.prob
  vv.iev.HI.prob <- dat$param$vv.iev.HI.prob
  vv.iev.W.prob <- dat$param$vv.iev.W.prob


  pos.role.class <- role.class[dal[, 1]]
  neg.role.class <- role.class[dal[, 2]]


  ins <- rep(NA, length(pos.role.class))
  ins[which(pos.role.class == "I")] <- 1  # "P"
  ins[which(pos.role.class == "R")] <- 0  # "N"
  ins[which(neg.role.class == "I")] <- 0  # "N"
  ins[which(neg.role.class == "R")] <- 1  # "P"

  vv <- which(pos.role.class == "V" & neg.role.class == "V")
  vv.race.combo <- paste0(race[dal[, 1]][vv], race[dal[, 2]][vv])
  vv.race.combo[vv.race.combo == "BIB"] <- "BBI"
  vv.race.combo[vv.race.combo == "HB"] <- "BH"
  vv.race.combo[vv.race.combo == "HIB"] <- "BHI"
  vv.race.combo[vv.race.combo == "WB"] <- "BW"
  
  vv.race.combo[vv.race.combo == "HBI"] <- "BIH"
  vv.race.combo[vv.race.combo == "HIBI"] <- "BIHI"
  vv.race.combo[vv.race.combo == "WBI"] <- "BIW"

  vv.race.combo[vv.race.combo == "HIH"] <- "HHI"
  vv.race.combo[vv.race.combo == "WH"] <- "HW"
  
  vv.race.combo[vv.race.combo == "WHI"] <- "HIW"
  
  vv.iev.prob <- (vv.race.combo == "BB") * vv.iev.B.prob +
                 (vv.race.combo == "BH") * mean(vv.iev.B.prob,vv.iev.H.prob) +
                 (vv.race.combo == "BHI") * mean(vv.iev.B.prob,vv.iev.HI.prob) +
                 (vv.race.combo == "BW") * mean(vv.iev.B.prob,vv.iev.W.prob) +
                 (vv.race.combo == "BIBI") * vv.iev.BI.prob +
                 (vv.race.combo == "BIH") * mean(vv.iev.BI.prob,vv.iev.H.prob) +
                 (vv.race.combo == "BIHI") * mean(vv.iev.BI.prob,vv.iev.HI.prob) +
                 (vv.race.combo == "BIW") * mean(vv.iev.BI.prob,vv.iev.W.prob) +
                 (vv.race.combo == "HH") * vv.iev.H.prob +
                 (vv.race.combo == "HHI") * mean(vv.iev.H.prob,vv.iev.HI.prob) +
                 (vv.race.combo == "HW") * mean(vv.iev.H.prob,vv.iev.W.prob) +
                 (vv.race.combo == "HIHI") * vv.iev.HI.prob +
                 (vv.race.combo == "HIW") * mean(vv.iev.HI.prob,vv.iev.W.prob) +
                 (vv.race.combo == "WW") * vv.iev.W.prob 
    


  iev <- rbinom(length(vv), 1, vv.iev.prob)
  ins[vv[iev == 1]] <- 2 # "B"
  vv.remaining <- vv[iev == 0]

  inspos.prob <- ins.quot[dal[, 1][vv.remaining]] /
                 (ins.quot[dal[, 1][vv.remaining]] + ins.quot[dal[, 2][vv.remaining]])
  inspos <- rbinom(length(vv.remaining), 1, inspos.prob)
  ins[vv.remaining[inspos == 1]] <- 1  # "P"
  ins[vv.remaining[inspos == 0]] <- 0  # "N"


  ## Output
  dat$temp$dal <- cbind(dal, ins)

  return(dat)
}
