
#' @export
updatenwp.mard <- function(dat, at) {

  n <- attributes(dat$el[[1]])$n
  maxdyads <- choose(n, 2)

  for (i in 1:3) {

    p <- dat$p[[i]]

    mf <- p$model.form
    md <- p$model.diss
    mhf <- p$MHproposal.form
    mhd <- p$MHproposal.diss

    if (i == 1) {

      ## Update model.form ##

      # edges
      mf$terms[[1]]$maxval <- maxdyads

      # nodefactor("deg.pers")
      nodecov <- dat$attr$deg.pers
      u <- sort(unique(nodecov))
      u <- u[-1] # remove base values here
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      ui <- seq(along = u)
      inputs <- c(ui, nodecov)
      mf$terms[[2]]$inputs <- c(length(ui), length(mf$terms[[2]]$coef.names), length(inputs), inputs)

      # absdiff("sqrt.age")
      nodecov <- dat$attr$sqrt.age
      power <- 1
      inputs <- c(power, nodecov)
      mf$terms[[3]]$inputs <- c(0, length(mf$terms[[3]]$coef.names), length(inputs), inputs)

      # nodematch("role.class)
      nodecov <- dat$attr$role.class
      u <- sort(unique(nodecov))
      u <- u[1:2] # keep = 1:2
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      dontmatch <- nodecov == (length(u) + 1)
      nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
      ui <- seq(along = u)
      inputs <- c(ui, nodecov)
      mf$terms[[4]]$inputs <- c(0, length(mf$terms[[4]]$coef.names), length(inputs), inputs)

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

    }
    if (i == 2) {

      ## Update model.form ##

      # edges
      mf$terms[[1]]$maxval <- maxdyads

      # nodefactor("deg.main")
      nodecov <- dat$attr$deg.main
      u <- sort(unique(nodecov))
      u <- u[-1] # remove base values here
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      ui <- seq(along = u)
      inputs <- c(ui, nodecov)
      mf$terms[[2]]$inputs <- c(length(ui), length(mf$terms[[2]]$coef.names), length(inputs), inputs)

      # concurrent
      mf$terms[[3]]$maxval <- n

      # absdiff("sqrt.age")
      nodecov <- dat$attr$sqrt.age
      power <- 1
      inputs <- c(power, nodecov)
      mf$terms[[4]]$inputs <- c(0, length(mf$terms[[4]]$coef.names), length(inputs), inputs)

      # nodematch("role.class)
      nodecov <- dat$attr$role.class
      u <- sort(unique(nodecov))
      u <- u[1:2] # keep = 1:2
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      dontmatch <- nodecov == (length(u) + 1)
      nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
      ui <- seq(along = u)
      inputs <- c(ui, nodecov)
      mf$terms[[5]]$inputs <- c(0, length(mf$terms[[5]]$coef.names), length(inputs), inputs)

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

    }
    if (i == 3) {

      ## Update model.form ##

      # edges
      mf$terms[[1]]$maxval <- maxdyads

      # nodefactor(c("deg.main", "deg.pers"))
      nodecov <- do.call(paste, c(sapply(c("deg.main", "deg.pers"),
                                         function(oneattr) dat$attr[[oneattr]],
                                         simplify = FALSE), sep = "."))
      u <- sort(unique(nodecov))
      u <- u[-1] # remove base values here
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      ui <- seq(along = u)
      inputs <- c(ui, nodecov)
      mf$terms[[2]]$inputs <- c(length(ui), length(mf$terms[[2]]$coef.names), length(inputs), inputs)


      # nodefactor("riskg", base = 3)
      nodecov <- dat$attr$riskg
      u <- sort(unique(nodecov))
      u <- u[-3] # remove base values here
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      ui <- seq(along = u)
      inputs <- c(ui, nodecov)
      mf$terms[[3]]$inputs <- c(length(ui), length(mf$terms[[3]]$coef.names), length(inputs), inputs)

      # absdiff("sqrt.age")
      nodecov <- dat$attr$sqrt.age
      power <- 1
      inputs <- c(power, nodecov)
      mf$terms[[4]]$inputs <- c(0, length(mf$terms[[4]]$coef.names), length(inputs), inputs)

      # nodematch("role.class)
      nodecov <- dat$attr$role.class
      u <- sort(unique(nodecov))
      u <- u[1:2] # keep = 1:2
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      dontmatch <- nodecov == (length(u) + 1)
      nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
      ui <- seq(along = u)
      inputs <- c(ui, nodecov)
      mf$terms[[5]]$inputs <- c(0, length(mf$terms[[5]]$coef.names), length(inputs), inputs)

      ## combined maxval ##
      mf$maxval[1] <- maxdyads

      ## Update model.diss ##
      md$terms[[1]]$maxval <- maxdyads
      md$maxval <- maxdyads

      ## Update MHproposal.form ##
      # no changes

      ## Update MHproposal.diss ##
      # no changes

    }

    dat$p[[i]] <- list(model.form = mf, model.diss = md,
                       MHproposal.form = mhf, MHproposal.diss = mhd)
  }


  return(dat)
}
