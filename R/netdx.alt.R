
#' @title netdx - shamp alternative
#'
#' @description Runs dynamic diagnostics on an ERGM/STERGM estimated through ergm.ego and altered for EpiModelHIV.
#'
#' @param x is a fit object.
#'        
#'
#' @return
#' This function returns \code{netdx}.
#'
#' @keywords module netwrok diagnostics
#' @export
#'


netdx.alt <- function (x, nsims = 1, dynamic = TRUE, nsteps, nwstats.formula = "formation", 
          set.control.ergm, set.control.stergm, keep.tedgelist = FALSE, 
          verbose = TRUE, ncores = 1) 
{
#  if (class(x) != "netest") {
#    stop("x must be an object of class netest", call. = FALSE)
#  }
  #if (class(x$fit) == "network") {
  #  nw <- x$fit
  #}
  #else {
    nw <- x$fit$newnetwork
    fit <- x$fit
  #}
  formation <- x$formation
  coef.form <- x$coef.form
  dissolution <- x$coef.diss$dissolution
  coef.diss <- x$coef.diss
  constraints <- x$constraints
  target.stats <- x$target.stats
  edapprox <- x$edapprox
  if (dynamic == TRUE && missing(nsteps)) {
    stop("Specify number of time steps with nsteps", call. = FALSE)
  }
  if (dynamic == FALSE && nwstats.formula == "formation") {
    nwstats.formula <- x$formation
  }
  if (verbose == TRUE) {
    cat("\nNetwork Diagnostics")
    cat("\n-----------------------\n")
  }
  if (verbose == TRUE) {
    cat("- Simulating", nsims, "networks")
  }
  if (edapprox == FALSE) {
    if (missing(set.control.stergm)) {
      set.control.stergm <- control.simulate.stergm()
    }
    if (nsims == 1 || ncores == 1) {
      diag.sim <- list()
      if (verbose == TRUE & nsims > 1) {
        cat("\n  |")
      }
      for (i in 1:nsims) {
        diag.sim[[i]] <- simulate(fit, time.slices = nsteps, 
                                  monitor = nwstats.formula, nsim = 1, control = set.control.stergm)
        if (verbose == TRUE & nsims > 1) {
          cat("*")
        }
      }
      if (verbose == TRUE & nsims > 1) {
        cat("|")
      }
    }
    else {
      cluster.size <- min(nsims, ncores)
      registerDoParallel(cluster.size)
      diag.sim <- foreach(i = 1:nsims) %dopar% {
        #simulate(fit, time.slices = nsteps, monitor = nwstats.formula, 
        #         nsim = 1, control = set.control.stergm)
        simulate(fit, time.slices = nsteps, monitor = nwstats.formula, nsim = 1)
      }
    }
  }
  if (edapprox == TRUE) {
    if (missing(set.control.ergm)) {
      set.control.ergm <- control.simulate.ergm()
    }
    if (missing(set.control.stergm)) {
      set.control.stergm <- control.simulate.network()
    }
    if (dynamic == TRUE) {
      if (nsims == 1 || ncores == 1) {
        diag.sim <- list()
        if (verbose == TRUE & nsims > 1) {
          cat("\n  |")
        }
        for (i in 1:nsims) {
          #if (class(x$fit) == "network") {
          #  fit.sim <- simulate(formation, basis = nw, 
          #                      coef = x$coef.form.crude, constraints = constraints)
          #}
          #else {
            fit.sim <- simulate(fit, control = set.control.ergm)
          #}
          diag.sim[[i]] <- simulate(fit.sim, formation = formation, 
                                    dissolution = dissolution, coef.form = coef.form, 
                                    coef.diss = coef.diss$coef.crude, time.slices = nsteps, 
                                    constraints = constraints, monitor = nwstats.formula, 
                                    nsim = 1, control = set.control.stergm)
          if (verbose == TRUE & nsims > 1) {
            cat("*")
          }
        }
        if (verbose == TRUE & nsims > 1) {
          cat("|")
        }
      }
      else {
        cluster.size <- min(nsims, ncores)
        registerDoParallel(cluster.size)
        diag.sim <- foreach(i = 1:nsims) %dopar% {
          #if (class(x$fit) == "network") {
          #  fit.sim <- simulate(formation, basis = nw, 
          #                      coef = x$coef.form.crude, constraints = constraints)
          #}
          #else {
            fit.sim <- simulate(fit)
          #}
          simulate(fit.sim, formation = formation, dissolution = dissolution, 
                   coef.form = coef.form, coef.diss = coef.diss$coef.crude, 
                   time.slices = nsteps, constraints = constraints, 
                   monitor = nwstats.formula, nsim = 1, control = set.control.stergm)
        }
      }
    }
    if (dynamic == FALSE) {
      #if (class(x$fit) == "network") {
      #  diag.sim <- simulate(formation, basis = nw, coef = x$coef.form.crude, 
      #                       constraints = constraints, nsim = nsims, statsonly = TRUE, 
      #                       monitor = nwstats.formula)
      #}
      #else {
        #diag.sim <- simulate(fit, nsim = nsims, statsonly = TRUE, 
        #                     control = set.control.ergm, monitor = nwstats.formula)
        
        diag.sim <- simulate(fit, nsim = nsims, statsonly = TRUE, 
                             monitor = nwstats.formula)
      #}
    }
  }
  if (verbose == TRUE) {
    cat("\n- Calculating formation statistics")
  }
  if (dynamic == TRUE) {
    stats <- list()
    for (i in 1:length(diag.sim)) {
      stats[[i]] <- as.matrix(attributes(diag.sim[[i]])$stats)[1:nsteps, 
                                                               , drop = FALSE]
    }
    if (nsims > 1) {
      merged.stats <- matrix(NA, nrow = nrow(stats[[1]]) * 
                               nsims, ncol = ncol(stats[[1]]))
      for (i in 1:ncol(stats[[1]])) {
        merged.stats[, i] <- as.numeric(sapply(stats, 
                                               function(x) c(x[, i])))
      }
      colnames(merged.stats) <- colnames(stats[[1]])
    }
    else {
      merged.stats <- stats[[1]]
    }
  }
  else {
    stats <- list(diag.sim[, !duplicated(colnames(diag.sim)), 
                           drop = FALSE])
    merged.stats <- diag.sim[, !duplicated(colnames(diag.sim)), 
                             drop = FALSE]
  }
  stats.means <- colMeans(merged.stats)
  stats.sd <- apply(merged.stats, 2, sd)
  stats.table <- data.frame(sorder = 1:length(names(stats.means)), 
                            names = names(stats.means), stats.means, stats.sd)
  ts.attr.names <- x$target.stats.names
  if (length(ts.attr.names) != length(target.stats)) {
    target.stats <- target.stats[which(target.stats > 0)]
  }
  ts.out <- data.frame(names = ts.attr.names, targets = target.stats)
  stats.table <- merge(ts.out, stats.table, all = TRUE)
  stats.table <- stats.table[do.call("order", stats.table[, 
                                                          "sorder", drop = FALSE]), , drop = FALSE]
  rownames(stats.table) <- stats.table$names
  stats.table$reldiff <- (stats.table$stats.means - stats.table$targets)/stats.table$targets
  stats.table.formation <- stats.table[, c(2, 4, 6, 5)]
  colnames(stats.table.formation) <- c("Target", "Sim Mean", 
                                       "Pct Diff", "Sim SD")
  if (dynamic == TRUE) {
    if (verbose == TRUE) {
      cat("\n- Calculating duration statistics")
    }
    sim.df <- list()
    for (i in 1:length(diag.sim)) {
      sim.df[[i]] <- as.data.frame(diag.sim[[i]])
    }
    ncens <- which(sim.df[[1]]$onset.censored == FALSE & 
                     sim.df[[1]]$terminus.censored == FALSE)
    durVec <- sim.df[[1]]$duration[ncens]
    if (nsims > 1) {
      for (i in 2:length(diag.sim)) {
        ncens <- which(sim.df[[i]]$onset.censored == 
                         FALSE & sim.df[[i]]$terminus.censored == FALSE)
        durVec <- c(durVec, sim.df[[i]]$duration[ncens])
      }
    }
    if (nsims == 1 || ncores == 1) {
      pages <- list()
      if (verbose == TRUE & nsims > 1) {
        cat("\n  |")
      }
      for (i in 1:length(diag.sim)) {
        pages[[i]] <- edgelist_meanage(el = sim.df[[i]])
        if (verbose == TRUE & nsims > 1) {
          cat("*")
        }
      }
      if (verbose == TRUE & nsims > 1) {
        cat("|")
      }
    }
    else {
      cluster.size <- min(nsims, ncores)
      registerDoParallel(cluster.size)
      pages <- foreach(i = 1:nsims) %dopar% {
        edgelist_meanage(el = sim.df[[i]])
      }
    }
    if (verbose == TRUE) {
      cat("\n- Calculating dissolution statistics")
    }
    if (nsims == 1 || ncores == 1) {
      if (verbose == TRUE & nsims > 1) {
        cat("\n  |")
      }
      prop.diss <- list()
      for (i in 1:length(diag.sim)) {
        prop.diss[[i]] <- sapply(1:nsteps, function(x) sum(sim.df[[i]]$terminus == 
                                                             x)/sum(sim.df[[i]]$onset < x & sim.df[[i]]$terminus >= 
                                                                      x))
        if (verbose == TRUE & nsims > 1) {
          cat("*")
        }
      }
      if (verbose == TRUE & nsims > 1) {
        cat("|")
      }
    }
    else {
      cluster.size <- min(nsims, ncores)
      registerDoParallel(cluster.size)
      prop.diss <- foreach(i = 1:nsims) %dopar% {
        sapply(1:nsteps, function(x) sum(sim.df[[i]]$terminus == 
                                           x)/sum(sim.df[[i]]$onset < x & sim.df[[i]]$terminus >= 
                                                    x))
      }
    }
    if (verbose == TRUE) {
      cat("\n ")
    }
    duration.mean <- mean(durVec)
    duration.sd <- sd(durVec)
    duration.expected <- exp(coef.diss$coef.crude[1]) + 1
    duration.pctdiff <- (duration.mean - duration.expected)/duration.expected
    dissolution.mean <- mean(unlist(prop.diss))
    dissolution.sd <- sd(unlist(prop.diss))
    dissolution.expected <- 1/(exp(coef.diss$coef.crude[1]) + 
                                 1)
    dissolution.pctdiff <- (dissolution.mean - dissolution.expected)/dissolution.expected
    stats.table.dissolution <- data.frame(Targets = c(duration.expected, 
                                                      dissolution.expected), Sim_Means = c(duration.mean, 
                                                                                           dissolution.mean), Pct_Diff = c(duration.pctdiff, 
                                                                                                                           dissolution.pctdiff), Sim_SD = c(duration.sd, dissolution.sd))
    colnames(stats.table.dissolution) <- c("Target", "Sim Mean", 
                                           "Pct Diff", "Sim SD")
    rownames(stats.table.dissolution) <- c("Edge Duration", 
                                           "Pct Edges Diss")
  }
  out <- list()
  out$nw <- nw
  out$formation <- formation
  out$coef.form <- coef.form
  out$dissolution <- dissolution
  out$coef.diss <- coef.diss
  out$constraints <- constraints
  out$edapprox <- edapprox
  out$target.stats <- ts.out
  out$nsims <- nsims
  out$dynamic <- dynamic
  out$stats <- stats
  out$stats.table.formation <- stats.table.formation
  if (dynamic == TRUE) {
    out$nsteps <- nsteps
    out$stats.table.dissolution <- stats.table.dissolution
    out$edgelist <- sim.df
    out$pages <- pages
    out$prop.diss <- prop.diss
    if (keep.tedgelist == TRUE) {
      out$tedgelist <- sim.df
    }
  }
  class(out) <- "netdx"
  return(out)
}
