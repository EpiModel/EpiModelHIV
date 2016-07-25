
# MSM -----------------------------------------------------------------

#' @title Network Resimulation Module
#'
#' @description Module function for resimulating the main, casual, and one-off
#'              networks for one time step.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
simnet_msm <- function(dat, at) {

  ## Edges correction
  dat <- edges_correct_msm(dat, at)

  ## Main network
  nwparam.m <- EpiModel::get_nwparam(dat, network = 1)
  dat <- updatenwp_msm(dat, network = 1)


  dat$el[[1]] <- tergmLite::simulate_network(p = dat$p[[1]],
                                             el = dat$el[[1]],
                                             coef.form = nwparam.m$coef.form,
                                             coef.diss = nwparam.m$coef.diss$coef.adj,
                                             save.changes = TRUE)


  dat$temp$new.edges <- NULL
  if (at == 2) {
    new.edges.m <- matrix(dat$el[[1]], ncol = 2)
  } else {
    new.edges.m <- attributes(dat$el[[1]])$changes
    new.edges.m <- new.edges.m[new.edges.m[, "to"] == 1, 1:2, drop = FALSE]
  }
  dat$temp$new.edges <- matrix(dat$attr$uid[new.edges.m], ncol = 2)


  ## Casual network
  nwparam.p <- EpiModel::get_nwparam(dat, network = 2)
  dat <- updatenwp_msm(dat, network = 2)

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
  dat <- updatenwp_msm(dat, network = 3)

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
#' @keywords module msm
#'
#' @export
#'
edges_correct_msm <- function(dat, at) {

  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)
  adjust <- log(old.num) - log(new.num)

  coef.form.m <- get_nwparam(dat, network = 1)$coef.form
  coef.form.m[1] <- coef.form.m[1] + adjust
  dat$nwparam[[1]]$coef.form <- coef.form.m

  coef.form.p <- get_nwparam(dat, network = 2)$coef.form
  coef.form.p[1] <- coef.form.p[1] + adjust
  dat$nwparam[[2]]$coef.form <- coef.form.p

  coef.form.i <- get_nwparam(dat, network = 3)$coef.form
  coef.form.i[1] <- coef.form.i[1] + adjust
  dat$nwparam[[3]]$coef.form <- coef.form.i

  return(dat)
}


#' @title Update Network Data Structure and Parameters
#'
#' @description Updates the internal data structure containing the main data
#'              passed into the TERGM resimulation algorithm. This step is
#'              necessary with the new tergmLite approach.
#'
#' @param dat Data object created in initialization module.
#' @param network Integer value for network number
#'
#' @keywords module msm
#'
#' @export
#'
#'
updatenwp_msm <- function(dat, network) {

  n <- attributes(dat$el[[1]])$n
  maxdyads <- choose(n, 2)

  p <- dat$p[[network]]

  if (network == 1) {

    mf <- p$model.form
    md <- p$model.diss
    mhf <- p$MHproposal.form
    mhd <- p$MHproposal.diss

    if (!identical(mf$coef.names, c("edges",
                                    "nodefactor.deg.pers.1",
                                    "nodefactor.deg.pers.2",
                                    "absdiff.sqrt.age",
                                    "nodematch.role.class.I",
                                    "nodematch.role.class.R"))) {
      stop("updatenwp_msm will not currently work with this formula, contact SJ")
    }

    ## Update model.form ##

    # edges
    mf$terms[[1]]$maxval <- maxdyads

    # nodefactor("deg.pers")
    dat$attr$deg.pers <- get_degree(dat$el[[2]])

    nodecov <- dat$attr$deg.pers
    u <- sort(unique(nodecov))
    u <- u[-1] # remove base values here
    nodecov <- match(nodecov, u, nomatch = length(u) + 1)
    ui <- seq(along = u)
    inputs <- c(ui, nodecov)
    mf$terms[[2]]$inputs <- c(length(ui), length(mf$terms[[2]]$coef.names),
                              length(inputs), inputs)

    # absdiff("sqrt.age")
    nodecov <- dat$attr$sqrt.age
    power <- 1
    inputs <- c(power, nodecov)
    mf$terms[[3]]$inputs <- c(0, length(mf$terms[[3]]$coef.names),
                              length(inputs), inputs)

    # nodematch("role.class)
    nodecov <- dat$attr$role.class
    u <- sort(unique(nodecov))
    u <- u[1:2] # keep = 1:2
    nodecov <- match(nodecov, u, nomatch = length(u) + 1)
    dontmatch <- nodecov == (length(u) + 1)
    nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
    ui <- seq(along = u)
    inputs <- c(ui, nodecov)
    mf$terms[[4]]$inputs <- c(0, length(mf$terms[[4]]$coef.names),
                              length(inputs), inputs)

    ## combined maxval ##
    mf$maxval[1] <- maxdyads


    ## Update model.diss ##
    md$terms[[1]]$maxval <- maxdyads
    md$maxval <- maxdyads

    ## Update MHproposal.form ##
    mhf$arguments$constraints$bd$attribs <-
           matrix(rep(mhf$arguments$constraints$bd$attribs[1], n), ncol = 1)
    mhf$arguments$constraints$bd$maxout <-
            matrix(rep(mhf$arguments$constraints$bd$maxout[1], n), ncol = 1)
    mhf$arguments$constraints$bd$maxin <- matrix(rep(n - 1, n), ncol = 1)
    mhf$arguments$constraints$bd$minout <-
           mhf$arguments$constraints$bd$minin <- matrix(rep(0, n), ncol = 1)

    ## Update MHproposal.diss ##
    mhd$arguments$constraints$bd <- mhf$arguments$constraints$bd

    dat$p[[network]] <- list(model.form = mf, model.diss = md,
                             MHproposal.form = mhf, MHproposal.diss = mhd)

  }

  if (network == 2) {

    mf <- p$model.form
    md <- p$model.diss
    mhf <- p$MHproposal.form
    mhd <- p$MHproposal.diss

    if (!identical(mf$coef.names, c("edges",
                                    "nodefactor.deg.main.1",
                                    "concurrent",
                                    "absdiff.sqrt.age",
                                    "nodematch.role.class.I",
                                    "nodematch.role.class.R"))) {
      stop("updatenwp_msm will not currently work with this formula, contact SJ")
    }


    ## Update model.form ##

    # edges
    mf$terms[[1]]$maxval <- maxdyads

    # nodefactor("deg.main")
    dat$attr$deg.main <- get_degree(dat$el[[1]])

    nodecov <- dat$attr$deg.main
    u <- sort(unique(nodecov))
    u <- u[-1] # remove base values here
    nodecov <- match(nodecov, u, nomatch = length(u) + 1)
    ui <- seq(along = u)
    inputs <- c(ui, nodecov)
    mf$terms[[2]]$inputs <- c(length(ui), length(mf$terms[[2]]$coef.names),
                              length(inputs), inputs)

    # concurrent
    mf$terms[[3]]$maxval <- n

    # absdiff("sqrt.age")
    nodecov <- dat$attr$sqrt.age
    power <- 1
    inputs <- c(power, nodecov)
    mf$terms[[4]]$inputs <- c(0, length(mf$terms[[4]]$coef.names),
                              length(inputs), inputs)

    # nodematch("role.class)
    nodecov <- dat$attr$role.class
    u <- sort(unique(nodecov))
    u <- u[1:2] # keep = 1:2
    nodecov <- match(nodecov, u, nomatch = length(u) + 1)
    dontmatch <- nodecov == (length(u) + 1)
    nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
    ui <- seq(along = u)
    inputs <- c(ui, nodecov)
    mf$terms[[5]]$inputs <- c(0, length(mf$terms[[5]]$coef.names),
                              length(inputs), inputs)

    ## combined maxval ##
    mf$maxval[1] <- maxdyads
    mf$maxval[3] <- n

    ## Update model.diss ##
    md$terms[[1]]$maxval <- maxdyads
    md$maxval <- maxdyads

    ## Update MHproposal.form ##
    mhf$arguments$constraints$bd$attribs <-
           matrix(rep(mhf$arguments$constraints$bd$attribs[1], n), ncol = 1)
    mhf$arguments$constraints$bd$maxout <-
            matrix(rep(mhf$arguments$constraints$bd$maxout[1], n), ncol = 1)
    mhf$arguments$constraints$bd$maxin <- matrix(rep(n - 1, n), ncol = 1)
    mhf$arguments$constraints$bd$minout <-
           mhf$arguments$constraints$bd$minin <- matrix(rep(0, n), ncol = 1)

    ## Update MHproposal.diss ##
    mhd$arguments$constraints$bd <- mhf$arguments$constraints$bd

    dat$p[[network]] <- list(model.form = mf, model.diss = md,
                             MHproposal.form = mhf, MHproposal.diss = mhd)

  }

  if (network == 3) {

    mf <- p$model.form
    mhf <- p$MHproposal

    if (!identical(mf$coef.names, c("edges",
                                    "nodefactor.deg.main.deg.pers.0.1",
                                    "nodefactor.deg.main.deg.pers.0.2",
                                    "nodefactor.deg.main.deg.pers.1.0",
                                    "nodefactor.deg.main.deg.pers.1.1",
                                    "nodefactor.deg.main.deg.pers.1.2",
                                    "nodefactor.riskg.1",
                                    "nodefactor.riskg.2",
                                    "nodefactor.riskg.4",
                                    "nodefactor.riskg.5",
                                    "absdiff.sqrt.age",
                                    "nodematch.role.class.I",
                                    "nodematch.role.class.R"))) {
      stop("updatenwp_msm will not currently work with this formula, contact SJ")
    }


    ## Update model.form ##

    # edges
    mf$terms[[1]]$maxval <- maxdyads

    # nodefactor(c("deg.main", "deg.pers"))
    # current main degree already written in last conditional block
    dat$attr$deg.pers <- get_degree(dat$el[[2]])

    nodecov <- do.call(paste, c(sapply(c("deg.main", "deg.pers"),
                                       function(oneattr) dat$attr[[oneattr]],
                                       simplify = FALSE), sep = "."))
    u <- sort(unique(nodecov))
    u <- u[-1] # remove base values here
    nodecov <- match(nodecov, u, nomatch = length(u) + 1)
    ui <- seq(along = u)
    inputs <- c(ui, nodecov)
    mf$terms[[2]]$inputs <- c(length(ui), length(mf$terms[[2]]$coef.names),
                              length(inputs), inputs)


    # nodefactor("riskg", base = 3)
    nodecov <- dat$attr$riskg
    u <- sort(unique(nodecov))
    u <- u[-3] # remove base values here
    nodecov <- match(nodecov, u, nomatch = length(u) + 1)
    ui <- seq(along = u)
    inputs <- c(ui, nodecov)
    mf$terms[[3]]$inputs <- c(length(ui), length(mf$terms[[3]]$coef.names),
                              length(inputs), inputs)

    # absdiff("sqrt.age")
    nodecov <- dat$attr$sqrt.age
    power <- 1
    inputs <- c(power, nodecov)
    mf$terms[[4]]$inputs <- c(0, length(mf$terms[[4]]$coef.names),
                              length(inputs), inputs)

    # nodematch("role.class)
    nodecov <- dat$attr$role.class
    u <- sort(unique(nodecov))
    u <- u[1:2] # keep = 1:2
    nodecov <- match(nodecov, u, nomatch = length(u) + 1)
    dontmatch <- nodecov == (length(u) + 1)
    nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
    ui <- seq(along = u)
    inputs <- c(ui, nodecov)
    mf$terms[[5]]$inputs <- c(0, length(mf$terms[[5]]$coef.names),
                              length(inputs), inputs)

    ## combined maxval ##
    mf$maxval[1] <- maxdyads

    ## Update MHproposal ##
    # no changes

    dat$p[[network]] <- list(model.form = mf, MHproposal = mhf)
  }

  return(dat)
}



# HET -----------------------------------------------------------------


#' @title Network Resimulation Module
#'
#' @description Module function to resimulate the dynamic network forward one
#'              time step conditional on current network structure and vertex
#'              attributes.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
simnet_het <- function(dat, at) {

  # Update edges coefficients
  dat <- edges_correct_het(dat, at)

  # Update internal ergm data
  dat <- update_nwp_het(dat)

  # Pull network parameters
  nwparam <- get_nwparam(dat)

  # Simulate edgelist
  dat$el <- tergmLite::simulate_network(p = dat$p,
                                        el = dat$el,
                                        coef.form = nwparam$coef.form,
                                        coef.diss = nwparam$coef.diss$coef.adj)

  return(dat)
}

update_nwp_het <- function(dat) {

  mf <- dat$p$model.form
  md <- dat$p$model.diss
  mhf <- dat$p$MHproposal.form
  mhd <- dat$p$MHproposal.diss

  n <- attributes(dat$el)$n
  maxdyads <- choose(n, 2)

  ## 1. Update model.form ##

  # edges
  # inputs <- c(0, 1, 0) # not changed
  mf$terms[[1]]$maxval <- maxdyads

  # nodematch
  nodecov <- dat$attr$male
  u <- sort(unique(nodecov))
  nodecov <- match(nodecov, u, nomatch = length(u) + 1)
  inputs <- nodecov
  mf$terms[[2]]$inputs <- c(0, 1, length(inputs), inputs)

  ## Update combined maxval here
  mf$maxval <- c(maxdyads, Inf)


  ## 2. Update model.diss ##
  md$terms[[1]]$maxval <- maxdyads
  md$maxval <- maxdyads


  ## 3. Update MHproposal.form ##
  mhf$arguments$constraints$bd$attribs <-
    matrix(rep(mhf$arguments$constraints$bd$attribs[1], n), ncol = 1)
  mhf$arguments$constraints$bd$maxout <-
    matrix(rep(mhf$arguments$constraints$bd$maxout[1], n), ncol = 1)
  mhf$arguments$constraints$bd$maxin <- matrix(rep(n, n), ncol = 1)
  mhf$arguments$constraints$bd$minout <-
    mhf$arguments$constraints$bd$minin <- matrix(rep(0, n), ncol = 1)


  ## 4. Update MHproposal.diss ##
  mhd$arguments$constraints$bd <- mhf$arguments$constraints$bd


  ## 5. Output ##
  p <- list(model.form = mf, model.diss = md,
            MHproposal.form = mhf, MHproposal.diss = mhd)

  dat$p <- p
  return(dat)
}


#' @title Adjustment for the Edges Coefficient with Changing Network Size
#'
#' @description Adjusts the edges coefficients in a dynamic network model
#'              to preserve the mean degree.
#'
#' @inheritParams aging_het
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
#' The network model parameters stored in \code{dat$nwparam} are updated.
#'
#' @references
#' Krivitsky PN, Handcock MS, and Morris M. "Adjusting for network size and
#' composition effects in exponential-family random graph models." Statistical
#' Methodology. 2011; 8.4: 319-339.
#'
#' @keywords module het
#'
#' @export
#'
edges_correct_het <- function(dat, at) {

  # Popsize
  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)

  # New Coefs
  coef.form <- get_nwparam(dat)$coef.form
  coef.form[1] <- coef.form[1] + log(old.num) - log(new.num)
  dat$nwparam[[1]]$coef.form <- coef.form

  return(dat)
}

