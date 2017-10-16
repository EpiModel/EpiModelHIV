# SHMAP -----------------------------------------------------------------

#' @title Network Resimulation Module
#'
#' @description Module function for resimulating the main, casual, and one-off
#'              networks for one time step.
#'
#' @inheritParams aging_msm
#'
#' @keywords module shamp
#'
#' @export
#'
simnet_shamp <- function(dat, at) {
  
  ## Edges correction
  dat <- edges_correct_msm(dat, at)
  
  ## Main network
  nwparam.c <- EpiModel::get_nwparam(dat, network = 1)
  
  if (dat$param$method == 1) {
    dat$attr$deg.pers <- get_degree(dat$el[[2]])
    dat$attr$deg.pers.c <- dat$attr$deg.pers
    dat$attr$deg.pers.c<-ifelse(dat$attr$deg.pers.c > 0,1,dat$attr$deg.pers.c)
    
  } else {
    dat$attr$deg.pers <- paste0(dat$attr$race, get_degree(dat$el[[2]]))
    dat$attr$deg.pers.c <- dat$attr$deg.pers
    dat$attr$deg.pers.c<-ifelse(dat$attr$deg.pers.c > 0,1,dat$attr$deg.pers.c)
  }
  dat <- tergmLite::updateModelTermInputs(dat, network = 1)
  

  dat$el[[1]] <- tergmLite::simulate_network(p = dat$p[[1]],
                                             el = dat$el[[1]],
                                             coef.form = nwparam.c$coef.form,
                                             coef.diss = nwparam.c$coef.diss$coef.adj,
                                             save.changes = TRUE)
  
  dat$temp$new.edges <- NULL
  if (at == 2) {
    new.edges.c <- matrix(dat$el[[1]], ncol = 2)
  } else {
    new.edges.c <- attributes(dat$el[[1]])$changes
    new.edges.c <- new.edges.c[new.edges.c[, "to"] == 1, 1:2, drop = FALSE]
  }
  dat$temp$new.edges <- matrix(dat$attr$uid[new.edges.c], ncol = 2)
  
  
  ## Casual network
  nwparam.p <- EpiModel::get_nwparam(dat, network = 2)
  
  if (dat$param$method == 1) {
    dat$attr$deg.cohab <- get_degree(dat$el[[1]])
    dat$attr$deg.cohab.c <- dat$attr$deg.cohab
    dat$attr$deg.cohab.c<-ifelse(dat$attr$deg.cohab.c > 0 , 1,dat$attr$deg.cohab.c)
    
  } else {
    dat$attr$deg.cohab <- paste0(dat$attr$race, get_degree(dat$el[[1]]))
    dat$attr$deg.cohab.c <- dat$attr$deg.cohab
    dat$attr$deg.cohab.c<-ifelse(dat$attr$deg.cohab.c > 0, 1,dat$attr$deg.cohab.c)
  }
  dat <- tergmLite::updateModelTermInputs(dat, network = 2)
  
  dat$el[[2]] <- tergmLite::simulate_network(p = dat$p[[2]],
                                             el = dat$el[[2]],
                                             coef.form = nwparam.p$coef.form,
                                             coef.diss = nwparam.p$coef.diss$coef.adj,
                                             save.changes = TRUE)
  
  if (at == 2) {
    new.edges.p <- matrix(dat$el[[2]], ncol = 2)
  } else {
    new.edges.p <- attributes(dat$el[[2]])$changes
    new.edges.p <- new.edges.p[new.edges.p[, "to"] == 1, 1:2, drop = FALSE]
  }
  dat$temp$new.edges <- rbind(dat$temp$new.edges,
                              matrix(dat$attr$uid[new.edges.p], ncol = 2))
  
  
  ## One-off network
  nwparam.i <- EpiModel::get_nwparam(dat, network = 3)
  
  if (dat$param$method == 1) {
    dat$attr$deg.pers <- get_degree(dat$el[[2]])
    dat$attr$deg.pers.c <- dat$attr$deg.pers
    dat$attr$deg.pers.c<-ifelse(dat$attr$deg.pers.c > 0, 1,dat$attr$deg.pers.c)
    
  } else {
    dat$attr$deg.pers <- paste0(dat$attr$race, get_degree(dat$el[[2]]))
    dat$attr$deg.pers.c <- dat$attr$deg.pers
    dat$attr$deg.pers.c<-ifelse(dat$attr$deg.pers.c > 0, 1,dat$attr$deg.pers.c)
  }
  dat <- tergmLite::updateModelTermInputs(dat, network = 3)
  
  dat$el[[3]] <- tergmLite::simulate_ergm(p = dat$p[[3]],
                                          el = dat$el[[3]],
                                          coef = nwparam.i$coef.form)
  
  if (dat$control$save.nwstats == TRUE) {
    dat <- calc_resim_nwstats(dat, at)
  }
  
  return(dat)
}



calc_resim_nwstats <- function(dat, at) {
  
  for (nw in 1:3) {
    n <- attr(dat$el[[nw]], "n")
    edges <- nrow(dat$el[[nw]])
    meandeg <- round(edges / n, 3)
    concurrent <- round(mean(get_degree(dat$el[[nw]]) > 1), 3)
    mat <- matrix(c(edges, meandeg, concurrent), ncol = 3, nrow = 1)
    if (at == 2) {
      dat$stats$nwstats[[nw]] <- mat
      colnames(dat$stats$nwstats[[nw]]) <- c("edges", "meand", "conc")
    }
    if (at > 2) {
      dat$stats$nwstats[[nw]] <- rbind(dat$stats$nwstats[[nw]], mat)
    }
  }
  
  return(dat)
}



#' @title Adjustment for the Edges Coefficient with Changing Network Size
#'
#' @description Adjusts the edges coefficients in a dynamic network model
#'              to preserve the mean degree.
#'
#' @inheritParams aging_msm
#'
#' @details
#' In HIV/STI modeling, there is typically an assumption that changes in
#' population size do not affect one's number of partners, specified as the
#' mean degree for network models. A person would not have 10 times the number
#' of partners should he move from a city 10 times as large. This module uses
#' the adjustment of Krivitsky et al. to adjust the edges coefficients on the
#' three network models to account for varying population size in order to
#' preserve that mean degree.
#'
#' @return
#' The network model parameters stored in \code{dat$nwparam} are updated for
#' each of the three network models.
#'
#' @references
#' Krivitsky PN, Handcock MS, and Morris M. "Adjusting for network size and
#' composition effects in exponential-family random graph models." Statistical
#' Methodology. 2011; 8.4: 319-339.
#'
#' @keywords module msm shamp
#'
#' @export
#'
edges_correct_msm <- function(dat, at) {
  
  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)
  adjust <- log(old.num) - log(new.num)
  
  coef.form.c <- get_nwparam(dat, network = 1)$coef.form
  coef.form.c[1] <- coef.form.c[1] + adjust
  dat$nwparam[[1]]$coef.form <- coef.form.c
  
  coef.form.p <- get_nwparam(dat, network = 2)$coef.form
  coef.form.p[1] <- coef.form.p[1] + adjust
  dat$nwparam[[2]]$coef.form <- coef.form.p
  
  coef.form.i <- get_nwparam(dat, network = 3)$coef.form
  coef.form.i[1] <- coef.form.i[1] + adjust
  dat$nwparam[[3]]$coef.form <- coef.form.i
  
  return(dat)
}

