
#' @title Apportion Least-Remainder Method
#'
#' @description Apportions a vector of values given a specified frequency
#'              distribution of those values such that the length of the output
#'              is robust to rounding and other instabilities.
#'
#' @param vector.length Length for the output vector.
#' @param values Values for the output vector.
#' @param proportions Proportion distribution with one number for each value. This
#'        must sum to 1.
#' @param shuffled If \code{TRUE}, randomly shuffle the order of the vector.
#'
#' @export
#'
apportion_lr <- function(vector.length, values,
                         proportions, shuffled = FALSE) {

  if (vector.length != round(vector.length)) {
    stop("argument vector.length must be a positive integer")
  }
  if (vector.length <= 0) {
    stop("argument vector.length must be a positive integer")
  }
  if (is.vector(values) == FALSE) {
    stop("argument values must be a vector")
  }
  if (!(length(proportions) == length(values) && round(sum(proportions), 10) == 1) &&
     (!(length(proportions) == length(values) - 1 && round(sum(proportions), 10) <= 1 &&
        round(sum(proportions), 10) >= 0))) {
    stop("error in proportions length or proportions sum")
  }

  if (length(proportions) == length(values) - 1) {
    proportions <- c(proportions, 1 - round(sum(proportions), 10))
  }
  result <- rep(NA, vector.length)
  exp.nums <- proportions * vector.length
  counts <- floor(exp.nums)
  remainders <- exp.nums - counts
  leftovers <- vector.length - sum(counts)
  if (leftovers > 0) {
    additions <- order(remainders, decreasing = TRUE)[1:leftovers]
    counts[additions]   <- counts[additions] + 1
  }
  result <- rep(values, counts)
  if (shuffled == TRUE) {
    result <- sample(result, length(result))
  }

  return(result)
}


#' @title Get Arguments from EpiModel Parameterization Functions
#'
#' @description Returns a list of argument names and values for use for parameter
#'              processing functions.
#'
#' @param formal.args The output of \code{formals(sys.function())}.
#' @param dot.args The output of \code{list(...)}.
#'
#' @export
#'
get_args <- function(formal.args, dot.args){
  p <- list()
  formal.args[["..."]] <- NULL
  for (arg in names(formal.args)) {
    p[arg] <- list(get(arg, pos = parent.frame()))
  }

  names.dot.args <- names(dot.args)
  if (length(dot.args) > 0) {
    for (i in 1:length(dot.args)) {
      p[[names.dot.args[i]]] <- dot.args[[i]]
    }
  }
  return(p)
}


#' @title Proportionally Reallocate PrEP Adherence Class Probability
#'
#' @description Shifts probabilities from the high-adherence category to the lower
#'              three adherence categories while maintaining the proportional
#'              distribution of those lower categories.
#'
#' @param in.pcp Input vector of length four for the \code{prep.class.prob}
#'        parameter.
#' @param reall The pure percentage points to shift from the high adherence
#'        group to the lower three groups.
#'
#' @export
#'
reallocate_pcp <- function(in.pcp = c(0.211, 0.07, 0.1, 0.619), reall = 0) {

  dist <- in.pcp[1]/sum(in.pcp[1:3])
  dist[2] <- in.pcp[2]/sum(in.pcp[1:3])
  dist[3] <- in.pcp[3]/sum(in.pcp[1:3])

  out.pcp <- rep(NA, 4)
  out.pcp[1:3] <- in.pcp[1:3] - (dist * reall)
  out.pcp[4] <- 1 - sum(out.pcp[1:3])

  return(out.pcp)
}


#' @title Truncate Simulation Time Series
#'
#' @description Left-truncates a simulation epidemiological summary statistics and
#'              network statistics at a specified time step.
#'
#' @param x Object of class \code{netsim}.
#' @param at Time step at which to left-truncate the time series.
#'
#' @details
#' This function would be used when running a follow-up simulation from time steps
#' \code{b} to \code{c} after a burnin period from time \code{a} to \code{b},
#' where the final time window of interest for data analysis is \code{b} to \code{c}
#' only.
#'
#' @export
#'
truncate_sim <- function(x, at) {

  rows <- at:(x$control$nsteps)

  # epi
  x$epi <- lapply(x$epi, function(r) r[rows, ])

  # control settings
  x$control$start <- 1
  x$control$nsteps <- max(seq_along(rows))

  return(x)
}



#' @title Source All Files in a Directory
#'
#' @description Loops over all files in a directory to source them to the
#'              Global Environment.
#'
#' @param path Directory of files to source.
#' @param verbose Print names of sourced files to console.
#' @param ... Additional arguments passed to \code{source}.
#'
#' @export
#'
sourceDir <- function(path, verbose = TRUE, ...) {

  fn <- list.files(path, pattern = "\\.[Rr]$")
  if (length(fn) == 0) {
    stop("No R files in that path", call. = FALSE)
  }

  if (verbose == TRUE) {
    cat("\n Sourced Files in", path,
        "\n -----------------")
  }

  for (nm in fn) {
    if (verbose == TRUE) cat("\n", nm)
    source(file.path(path, nm), ...)
  }
}


nbsdtosize <- function(mu, sd) {
  mu ^ 2 / (sd ^ 2 - mu)
}

get_attr <- function(x, sim = 1) {
  if (is.null(x$attr)) {
    stop("No attr on x")
  } else {
    x$attr[[1]]
  }
}

cut_age <- function(age, breaks = c(0, 29, 39, Inf)) {
  cut(age, breaks = breaks, labels = FALSE)
}

keep.attr <- function(attrList, keep) {
  lapply(attrList, function(x) x[keep])
}



#######################    tergmLite funtions that are no longer being exported.

#' @title Prep stergm
#'
#' @description set up the structure for the stergm function. No longer exported in tegmLite
#'
#' @param path Directory of files to source.
#' @param verbose Print names of sourced files to console.
#' @param ... Additional arguments passed to \code{source}.
#'
#' @export
#'
stergm_prep <- function(nw,
                        formation,
                        dissolution,
                        coef.form,
                        coef.diss,
                        constraints,
                        control = control.simulate.network()) {
  
  if (inherits(nw, "network") == FALSE) {
    stop("A network object must be given")
  }
  
  formation <- statnet.common::nonsimp_update.formula(formation, nw ~ ., from.new = "nw")
  dissolution <- statnet.common::nonsimp_update.formula(dissolution, nw ~ ., from.new = "nw")
  
  model.form <- ergm::ergm_model(formation, nw, role = "formation")
  if (!missing(coef.form) && ergm::nparam(model.form) != length(coef.form)) {
    stop("coef.form has ", length(coef.form), " elements, while the model requires ",
         ergm::nparam(model.form), " parameters.")
  }
  model.diss <- ergm::ergm_model(dissolution, nw, role = "dissolution")
  if (!missing(coef.diss) && ergm::nparam(model.diss) != length(coef.diss)) {
    stop("coef.diss has ", length(coef.diss), " elements, while the model requires ",
         ergm::nparam(model.diss), " parameters.")
  }
  
  MHproposal.form <- ergm::ergm_proposal(constraints, control$MCMC.prop.args.form,
                                         nw, weights = control$MCMC.prop.weights.form,
                                         class = "f")
  MHproposal.diss <- ergm::ergm_proposal(constraints, control$MCMC.prop.args.diss,
                                         nw, weights = control$MCMC.prop.weights.diss,
                                         class = "d")
  
  MHproposal.form$arguments$constraints$.attributes <- NULL
  MHproposal.diss$arguments$constraints$.attributes <- NULL
  
  out <- list(model.form = model.form, model.diss = model.diss,
              MHproposal.form = MHproposal.form, MHproposal.diss = MHproposal.diss)
  
  return(out)
}


#' @title PrEP ergm
#'
#' @description set up the structure for the ergm function. No longer exported in tegmLite
#'
#' @param path Directory of files to source.
#' @param verbose Print names of sourced files to console.
#' @param ... Additional arguments passed to \code{source}.
#'
#' @export
#'
ergm_prep <- function(nw,
                      formation,
                      coef,
                      constraints,
                      control = ergm::control.simulate.ergm()) {
  
  form <- statnet.common::nonsimp_update.formula(formation, nw ~ ., from.new = "nw")
  m <- ergm::ergm_model(form, nw, response = NULL, role = "static")
  
  MHproposal <- ergm_proposal(constraints, arguments = control$MCMC.prop.args,
                              nw = nw, weights = control$MCMC.prop.weights, class = "c",
                              reference = ~Bernoulli, response = NULL)
  
  MHproposal$arguments$constraints$.attributes <- NULL
  
  out <- list(model.form = m, MHproposal = MHproposal)
  return(out)
}



# Note eefit is the ergm.ego fit object, disscoef is the dissolution coef object

#' @title Ego to est for KenyaTM
#'
#' @description temporary version of ego.netest for imposing changes in the Kenya formation models
#'
#' @param path Directory of files to source.
#' @param verbose Print names of sourced files to console.
#' @param ... Additional arguments passed to \code{source}.
#'
#' @export
#'
ego.netest.KTM <- function (eefit,
                        disscoef, #set.formula.nw=FALSE, start.net=NULL, 
                        obs.adjust=FALSE, obs.offset = 52) {
  
  ## the netest object contains
  ##[1] "fit"                "formation"          "target.stats"       "target.stats.names" "coef.form"          "coef.form.crude"    "coef.diss"
  ##[8] "constraints"        "edapprox"
  
  # Preliminaries
  
  ## Initialize output list
  
  out <- list()
  
  ## fit is an ergm object that contains.
  ##[1] "coef"          "sample"        "iterations"    "MCMCtheta"     "loglikelihood" "gradient"      "hessian"       "covar"
  ##[9] "failure"       "network"       "newnetworks"   "newnetwork"    "coef.init"     "est.cov"       "coef.hist"     "stats.hist"
  ##[17] "steplen.hist"  "control"       "etamap"        "formula"       "target.stats"  "target.esteq"  "constraints"   "reference"
  ##[25] "estimate"      "offset"        "drop"          "estimable"     "null.lik"
  
  ## Delete ergm.ego fit object components we will not need
  ## Keeping ppopsize and popsize just in case
  eefit$egodata <- eefit$sample.obs <- eefit$constrained <- NULL
  eefit$v <- eefit$m <- eefit$ergm.formula <- eefit$offset.coef <- NULL
  eefit$ergm.offset.coef <- eefit$ergm.covar <- eefit$DtDe <- NULL
  
  
  # First: The fit component of the netest object
  #        Assign the ergm.ego fit object, and modify as needed to remove
  #        the network size adjustments unique to ergm.ego.
  
  out$fit <- eefit
  class(out$fit) <- "ergm"
  
  
  ## edges coef for ergm = (edges + network size adjustment) for ergm.ego
  out$fit$coef <- eefit$coef[2:length(eefit$coef)]
  out$fit$coef[1] <- out$fit$coef[1] + eefit$coef[1]
  t<-c(out$fit$coef[1:5],-Inf,out$fit$coef[6:length(out$fit$coef)]) 
  names(t)<-c(names(out$fit$coef[1:5]),"nodefactor.age_grp.6",names(out$fit$coef[6:length(out$fit$coef)]))
  out$fit$coef <-t
  
  ## Sample (modified)
  ncol<-length(eefit$sample[1,])
  out$fit$sample<-eefit$sample[,-1]
  t<-cbind(out$fit$sample[,1:5],rep(0,length(out$fit$sample[,1])),out$fit$sample[,6:length(out$fit$sample[1,])]) 
  colnames(t)<-c(names(out$fit$sample[1,1:5]),"nodefactor.age_grp.6",names(out$fit$sample[1,6:length(out$fit$sample[1,])]))
  out$fit$sample <-t
  class(out$fit$sample) <- "mcmc"
  
  # Iterations
  out$fit$iterations <- eefit$iterations
  
  ## MCMCtheta coeficient (modified)
  out$fit$MCMCtheta<-eefit$MCMCtheta[2:length(eefit$MCMCtheta)]
  out$fit$MCMCtheta[1]<-out$fit$MCMCtheta[1] + eefit$MCMCtheta[1]
  t<-c(out$fit$MCMCtheta[1:5],-Inf,out$fit$MCMCtheta[6:length(out$fit$MCMCtheta)]) 
  names(t)<-c(names(out$fit$MCMCtheta[1:5]),"nodefactor.age_grp.6",names(out$fit$MCMCtheta[6:length(out$fit$MCMCtheta)]))
  out$fit$MCMCtheta <-t
  
  # Loglikelihood
  out$fit$loglikelihood <- eefit$loglikelihood
  
  ## Gradient (modified).
  out$fit$gradient<-eefit$gradient[-1]
  out$fit$gradient<-c(out$fit$gradient[1:5],NA,out$fit$gradient[6:length(out$fit$gradient)])
  
  ## Hessian matrix (modified) (r1,c1).
  ncol <- length(eefit$hessian[1,])
  out$fit$hessian <- eefit$hessian[2:ncol,2:ncol]
  t<-cbind(out$fit$hessian[,1:5],rep(NA,length(out$fit$hessian[,1])), out$fit$hessian[,6:length(out$fit$hessian[1,])])
  t<-rbind(t[1:5,],rep(NA,length(t[1,])), t[6:length(t[,1]),])
  rownames(t)<-c(rownames(t[1:5,]),"nodefactor.age_grp.6",rownames(t[7:length(t[1,]),]))
  colnames(t)<-c(colnames(t[,1:5]),"nodefactor.age_grp.6",colnames(t[,7:length(t[,1])]))
  out$fit$hessian <- t
  
  ## Covar matrix (modified) (r1,c1).
  ncol <- length(eefit$covar[1,])
  out$fit$covar <- eefit$covar[2:ncol,2:ncol]
  t<-cbind(out$fit$covar[,1:5],rep(NA,length(out$fit$covar[,1])), out$fit$covar[,6:length(out$fit$covar[1,])])
  t<-rbind(t[1:5,],rep(NA,length(t[1,])), t[6:length(t[,1]),])
  rownames(t)<-c(rownames(t[1:5,]),"nodefactor.age_grp.6",rownames(t[7:length(t[1,]),]))
  colnames(t)<-c(colnames(t[,1:5]),"nodefactor.age_grp.6",colnames(t[,7:length(t[,1])]))
  out$fit$covar <- t
  
  # Failure
  out$fit$failure <- eefit$failure
  
  # network, newnetworks and newnetwork
  out$fit$network <- eefit$network
  out$fit$newnetworks <- eefit$newnetworks
  out$fit$newnetwork <- eefit$newnetwork
  
  ## coef.init (modified)
  out$fit$coef.init<-eefit$coef.init[2:length(eefit$coef.init)]
  out$fit$coef.init[1]<-out$fit$coef.init[1] + eefit$coef.init[1]
  t<-c(out$fit$coef.init[1:5],-Inf,out$fit$coef.init[6:length(out$fit$coef.init)]) 
  names(t)<-c(names(out$fit$coef.init[1:5]),"nodefactor.age_grp.6",names(out$fit$coef.init[6:length(out$fit$coef.init)]))
  out$fit$coef.init <-t
  
  ## est.cov matrix adjusted to remove network size adjustment (r1,c1).
  ncol <- length(eefit$est.cov[1,])
  out$fit$est.cov <- eefit$est.cov[2:ncol,2:ncol]
  t<-cbind(out$fit$est.cov[,1:5],rep(0,length(out$fit$est.cov[,1])), out$fit$est.cov[,6:length(out$fit$est.cov[1,])])
  t<-rbind(t[1:5,],rep(0,length(t[1,])), t[6:length(t[,1]),])
  rownames(t)<-c(rownames(t[1:5,]),"nodefactor.age_grp.6",rownames(t[7:length(t[1,]),]))
  colnames(t)<-c(colnames(t[,1:5]),"nodefactor.age_grp.6",colnames(t[,7:length(t[,1])]))
  out$fit$est.cov <- t
  
 
  ## coef.hist (modified)
  ncol <- length(eefit$coef.hist[1,])
  eefit$coef.hist[,2] <- eefit$coef.hist[,1] + eefit$coef.hist[,2]
  out$fit$coef.hist <- eefit$coef.hist[,-1]
  t<-cbind(out$fit$coef.hist[,1:5],rep(-Inf,length(out$fit$coef.hist[,1])), out$fit$coef.hist[,6:length(out$fit$coef.hist[1,])])
  colnames(t)<-c(colnames(out$fit$coef.hist[,1:5]),"nodefactor.age_grp.6",colnames(out$fit$coef.hist[,6:length(out$fit$coef.hist[1,])]))
  out$fit$coef.hist <- t
  
  ## stats.hist (modified)
  out$fit$stats.hist <- eefit$stats.hist[,-1]
  t<-cbind(out$fit$stats.hist [,1:5],rep(0,length(out$fit$stats.hist [,1])),out$fit$stats.hist [,6:length(out$fit$stats.hist [1,])]) 
  colnames(t)<-c(names(out$fit$stats.hist [1,1:5]),"nodefactor.age_grp.6",names(out$fit$stats.hist[1,6:length(out$fit$stats.hist[1,])]))
  out$fit$stats.hist  <-t
  
  # steplen.hist
  out$fit$steplen.hist <- eefit$steplen.hist
  
  ## control (modified)
  out$fit$control <- eefit$control
  out$fit$control$init <- eefit$control$init[-1]
  out$fit$control$init[1] <- eefit$control$init[1] + eefit$control$init[2]
  t<-c(out$fit$control$init[1:5],-Inf,out$fit$control$init[6:length(out$fit$control$init)]) 
  names(t)<-c(names(out$fit$control$init[1:5]),"nodefactor.age_grp.6",names(out$fit$control$init[6:length(out$fit$control$init)]))
  out$fit$control$init <-t
  
  ## etamap (modified) (first element of each attribute dropped)
  out$fit$etamap <- eefit$etamap
  #out$fit$etamap$canonical <- eefit$etamap$canonical[-(length(eefit$etamap$canonical))]
  out$fit$etamap$canonical <- eefit$etamap$canonical
  #out$fit$etamap$offsetmap <- eefit$etamap$offsetmap[-1]
  out$fit$etamap$offsetmap <- c(eefit$etamap$offsetmap[2:5], FALSE, eefit$etamap$offsetmap[6:length(eefit$etamap$offsetmap)])
  out$fit$etamap$offset <- eefit$etamap$offset[-1]
  #out$fit$etamap$offsettheta <- eefit$etamap$offsettheta[-1]
  out$fit$etamap$offsettheta <- c(eefit$etamap$offsettheta[2:5],FALSE,eefit$etamap$offsettheta[6:length(eefit$etamap$offsettheta)])
  out$fit$etamap$curved <- eefit$etamap$curved
  #out$fit$etamap$etalength <- eefit$etamap$etalength - 1
  out$fit$etamap$etalength <- eefit$etamap$etalength
  # formula, 
  out$fit$formula <- eefit$formula
  
  # Commenting out all the code that set the network -- that
  # should be done in another initialization file, not here.
  
  #  z<-(as.formula("nw ~ o"))
  #  out$fit$formula[2] <- z[2]
  
  # #JKB 4/24/18 - nw object needs to be added to the environment of the 
  # #formation formula; patterned after netest() lines 182-183
  # #DTH 4/46/2018 - testing cause of Males-Male ties
  # # deleting all ties to creat an empty network
  # 
  # #    ego.netest <- function (eefit, disscoef, 
  # #                            set.formula.nw=FALSE, start.net=NULL, 
  # #                            obs.adjust=FALSE, obs.offset = 52)
  # 
  # #JKB 6/5/18: allow for an alternate starting network if using ppop.wt='round'
  # if (is.null(start.net)) {
  #   nw <- eefit$network
  # } else nw <- start.net
  # 
  # count <- network::network.edgecount(nw)
  # network::delete.edges(nw,1:count)
  # if (set.formula.nw) environment(out$fit$formula) <- environment()
  
  # target stats
  
  if(is.na(eefit$target.stats[1])){
    out$fit$target.stats <- eefit$target.stats[-1]
  } else out$fit$target.stats <- eefit$target.stats
  
  t<-c(out$fit$target.stats[1:5],0,out$fit$target.stats[6:length(out$fit$target.stats)])
  names(t)<-c(names(out$fit$target.stats[1:5]),"nodefactor.age_grp.6",names(out$fit$target.stats[6:length(out$fit$target.stats)]))
  out$fit$target.stats<-t
  
  # target.esteq
  out$fit$target.esteq <- eefit$target.esteq
  t<-matrix(c(out$fit$target.esteq[,1:5],0,out$fit$target.esteq[,6:length(out$fit$target.esteq)]),nrow=1)
  colnames(t)<-c(names(out$fit$target.esteq[,1:5]),"nodefactor.age_grp.6",names(out$fit$target.esteq[,6:length(out$fit$target.esteq)]))
  out$fit$target.esteq<-t
  
  # constraints
  out$fit$constraints <- eefit$constraints
  
  # reference
  out$fit$reference <- eefit$reference
  
  # estimate
  out$fit$estimate <- eefit$estimate
  
  ## offset (modified)
  #out$fit$offset <- eefit$offset[-1]
  out$fit$offset <- eefit$offset
  out$fit$offset[1] <- FALSE
  
  ## drop (modified)
  #out$fit$drop <- eefit$drop[-1]
  out$fit$drop <- eefit$drop
  
  ## estimable (modified)
  #out$fit$estimable <- eefit$estimable[-1]
  out$fit$estimable <- eefit$estimable
  
  # null.lik
  out$fit$null.lik <- eefit$null.lik[-1]
  
  
  # Next: The other netest components (other than fit).
  
  ## the netest object contains these components:
  ##[1] "fit"                "formation"          "target.stats"       "target.stats.names" "coef.form"          "coef.form.crude"    "coef.diss"
  ##[8] "constraints"        "edapprox"
  
  formation <- eefit$formula[-2]
  target.stats <- as.numeric(out$fit$target.stats)

  # Remove the offset term(s) from the target stats components
  
  target.stats <- target.stats[!out$fit$offset]
  target.stats.names <- names(out$fit$target.stats)[!out$fit$offset]
  
  # coef.form and coef.form.crude
  
  coef.form <- coef.form.crude <- out$fit$coef
  #  coef.form.crude <- coef.form
  
  ##  Adjust edges coef if partner count obs period != sim time step.
  ##  Edges coefficient must be adjusted to reflect formation 
  ##  rates per simulated timestep (in both fit$coef  and coef.form).
  ##  The default is set to 52 (obs period = 1 yr, sim step = 1 wk)
  
  if (obs.adjust){
    coef.form[1] <- out$fit$coef[1] <- coef.form[1] -log(obs.offset)
    #    out$fit$coef[1] <- out$fit$coef[1] -log(ot.offset)
  }
  
  # for edapprox=T:  Adjust all formation coefs by their dissolution coefs
  # This has been updated with Sam's code from 
  # https://github.com/statnet/EpiModel/blob/master/R/netest.R
  # but this requires careful model specification to ensure
  # order and position of formation coefs = dissolution coefs
  
  edapprox <- TRUE
  
  if (disscoef$coef.crude[1] > -Inf) {
    len.cd <- length(disscoef$coef.crude)
    coef.form[1:len.cd] <- coef.form[1:len.cd] - disscoef$coef.crude
  }
  
  constraints <- eefit$constraints
  
  # Finally, assign the new components to the netest object
  # NOTE: this doesn't include the weights from the egodata object
  
  out$formation <- formation
  out$target.stats <- target.stats
  out$target.stats.names <- target.stats.names
  out$coef.form <- coef.form
  out$coef.form.crude <- coef.form.crude
  out$coef.diss <- disscoef
  out$constraints <- constraints
  out$edapprox <- edapprox
  
  class(out) <- "netest"
  return(out)
}

#' @title Convert ergm.ego fit to netest object
#'
#' @description Converts an ergm.ego output object and a dissolution
#' coef (output from Epimodel function \code{dissolution_coefs} 
#' to a netest object for use in EpiModel.  Started from \code{ee.netest}.
#'
#' @param eefit An object of class \code{ergm.ego} output from the 
#'              \code{ergm.ego} function in the ergm.ego package.
#' @param disscoef An object of class \code{disscoef} output from the 
#'                 \code{dissolution_coefs} function in the EpiModel package.
#' @param obs.adjust Are the partner counts taken over a period 
#'                    different than the simulation step length? Typically
#'                    used for one-time partner models.
#' @param obs.offset A value representing the partner count observation 
#'                   period in simulation time steps.
#'                   The edges coefficient in the formation model is adjusted 
#'                   by \code{-log(obs.offset)}
#'                   The default (52) assumes the partner count is over the 
#'                   last year, and the simulation step is one week.  
#'        
#'
#' @return This function returns a \code{netest} for use with EpiModel functions. The output is organized as a list of 9 objects, 
#'         including \code{fit}, \code{formation}, \code{target.stats}, \code{target.stats.names}, 
#'         \code{coef.form}, \code{dissolution}, \code{coef.diss}, \code{constraints}, \code{edapprox}. 
#'         \item{"fit"}{Fitting object from \code{ergm.ego}.}
#'         \item{"fit$coef"}{Estimated coefficents from \code{ergm.ego}.}
#'         \item{"fit$sample"}{A list of matrices storing MCMC samples for each coefficient.}
#'         \item{"fit$iterations"}{The number of Newton-Raphson iterations required before convergence.}
#'         \item{"fit$MCMCtheta"}{The value of \eqn{\theta} used to produce the Markov Chain Monte Carlo samples.}
#'         \item{"fit$loglikelihood"}{The approximate change in log-likelihood in the last iteration. The value is only approximate because it is estimated based on the MCMC random sample.}
#'         \item{"fit$gradiant"}{The value of the gradient vector of the approximated loglikelihood function, evaluated at the maximizer. This vector should be very close to zero.}
#'         \item{"fit$hessian"}{}
#'         \item{"fit$covar"}{Approximate covariance matrix for the MLE, based on the inverse Hessian of the approximated loglikelihood evaluated at the maximizer.}
#'         \item{"fit$failure"}{Logical: Did the MCMC estimation fail?}
#'         \item{"fit$network"}{Original network.}
#'         \item{"fit$newnetwork"}{The final network at the end of the MCMC simulation.}
#'         \item{"fit$newnetworks"}{}
#'         \item{"fit$coef.init"}{The initial value \eqn{\theta}.}
#'         \item{"fit$est.cov"}{The covariance matrix of the model statistics in the final MCMC sample.}
#'         \item{"fit$coef.hist"&"$stats.hist"&"$steplen.hist"}{For the MCMLE method, the history of coefficients, Hummel step lengths, and average model statistics for each iteration.}
#'         \item{"fit$control"}{A list of control parameters for algorithm tuning.}
#'         \item{"fit$etamap"}{The set of functions mapping the true parameter theta to the canonical parameter eta (irrelevant except in a curved exponential family model).}
#'         \item{"fit$formula"}{The original formula used in \code{ergm.ego} function}
#'         \item{"fit$target.stats"}{The target stats used for estimation.}
#'         \item{"fit$target.esteq"}{Used for curved models to preserve the target mean values of the curved terms. It is identical to target.stats for non-curved models.}
#'         \item{"fit$constraints"}{Constraints used during estimation (passed through from the Arguments).}
#'         \item{"fit$reference"}{The reference measure used during estimation (passed through from the Arguments).}
#'         \item{"fit$estimate"}{The estimation method used (passed through from the Arguments).}
#'         \item{"fit$offset"}{Vector of logical telling which model parameters are to be set at a fixed value (i.e., not estimated).}
#'         \item{"fit$drop"}{}
#'         \item{"fit$estimable"}{A logical vector indicating which terms could not be estimated due to a constraints constraint fixing that term at a constant value.}
#'         \item{"fit$null.lik"}{Log-likelihood of the null model. Valid only for unconstrained models.}
#'         \item{"formation"}{The original formation formula}
#'         \item{"target.stats"}{The target stats used for estimation}
#'         \item{"target.stats.names"}{The names of the target stats.}
#'         \item{"coef.form"}{The coefficients of the formation model.}
#'         \item{"coef.diss"}{The coefficients of the dissolution model.}
#'         \item{"constraints"}{}
#'         \item{"edapprox"}{Logical}
#' @examples
#'
#' \dontrun{
#'  ergm.ego model
#'   library(ergm.ego)
#'  data(faux.mesa.high)
#'  fmh.ego <- as.egodata(faux.mesa.high)
#'  ego_fit <- ergm.ego(fmh.ego~edges+degree(0:3)
#'                     +nodefactor("Race")+nodematch("Race")
#'                     +nodefactor("Sex")+nodematch("Sex")
#'                     +absdiff("Grade"),
#'                     ppopsize="samp") # uses size of fmh for ppop
#'  
#'  # Also works - use to test having an offset, this one prevents same sex ties
#'  ego_fit2 <- ergm.ego(fmh.ego~edges+offset(nodematch("Sex", diff=FALSE)),
#'                     offset.coef=c(-Inf),
#'                     ppopsize="samp") # as above
#'  
#'  #----------------------------
#'  # Homogeneous dissolution model (note, model specified now with EpiModel)
#'  library(EpiModel)
#'  diss = ~offset(edges)
#'  ## Need to set two times:  mean partnership duration, mean lifespan
#'  ## here we use lifespan = 100*partnershp duration
#'  coef.diss <- dissolution_coefs(dissolution = diss,
#'                                 duration = 40, 
#'                                 d.rate = 1/(100*40))
#'  
#'  ## Translate the ergm.ego fit object into a netest object
#'  netest_fit <- egonet::ego.netest(ego_fit, coef.diss)
#'  
#'  # netdx, static, note monitoring of non-model terms
#'  dx1 <- netdx(netest_fit, nsims = 1e4, dynamic = FALSE,
#'               nwstats.formula = ~edges + meandeg + concurrent)
#'  dx1
#'  plot(dx1, method = "b", stats = c("edges", "concurrent"))
#'  
#'  # netdx, dynamic
#'  dx2 <- netdx(netest_fit, nsims = 5, nsteps = 1000,
#'               nwstats.formula = ~edges + meandeg + concurrent,
#'               set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))
#'  dx2
#'  plot(dx2, stats = c("edges", "meandeg"), plots.joined = FALSE)
#'  plot(dx2, type = "duration")
#'  plot(dx2, type = "dissolution", qnts.col = "orange2")
#'  plot(dx2, type = "dissolution", method = "b", col = "bisque")
#'  
#'  
#'  #----------------------------
#'  # Dissolution varies by Race
#'  ## Note for this to work, all terms in the dissolution model
#'  ## must also be in the formation model, and in the same order in
#'  ## both models, with the common terms at the front of the formation
#'  ## model specification
#'  
#'  dissR <- ~offset(edges) + offset(nodematch("Race", diff=TRUE))
#'  coef.dissR <- dissolution_coefs(dissolution = dissR,
#'                                  duration = c(30, 10, 20, 30, 10, 40))
#'  netestR_fit <- egonet::ego.netest(ego_fit, coef.dissR)
#'  
#'  dxR1 <- netdx(netestR_fit, nsims=1e4, dynamic=FALSE) # static draws
#'  plot(dxR1, method = "b")
#'  
#'  # see what happens when sequential=F, static draws from same initial net
#'  dxR2 <- netdx(netestR_fit, nsims=5, nsteps=1000, 
#'                dynamic=FALSE, sequential=FALSE, 
#'                set.control.ergm=control.simulate.ergm(MCMC.burnin=1e6))
#'  dxR2
#'  plot(dxR2, plots.joined = FALSE)
#'  
#'  # Now with dynamic=T, note it doesn't look great, may need larger interval
#'  dxR3 <- netdx(netestR_fit, nsims=5, nsteps=1000, 
#'                set.control.ergm=control.simulate.ergm(MCMC.burnin=1e6))
#'  plot(dxR3, plots.joined = FALSE)
#'  plot(dxR3, type = "duration")
#'  plot(dxR3, type = "dissolution", qnts.col = "orange2")
#'  plot(dxR3, type = "dissolution", method = "b", col = "bisque")
#'
#' }
#'
#' @keywords module convert object class
#' @export

# Note eefit is the ergm.ego fit object, disscoef is the dissolution coef object

ego.netest <- function (eefit,
                        disscoef, #set.formula.nw=FALSE, start.net=NULL, 
                        obs.adjust=FALSE, obs.offset = 52) {
  
  ## the netest object contains
  ##[1] "fit"                "formation"          "target.stats"       "target.stats.names" "coef.form"          "coef.form.crude"    "coef.diss"
  ##[8] "constraints"        "edapprox"
  
  # Preliminaries
  
  ## Initialize output list
  
  out <- list()
  
  ## fit is an ergm object that contains.
  ##[1] "coef"          "sample"        "iterations"    "MCMCtheta"     "loglikelihood" "gradient"      "hessian"       "covar"
  ##[9] "failure"       "network"       "newnetworks"   "newnetwork"    "coef.init"     "est.cov"       "coef.hist"     "stats.hist"
  ##[17] "steplen.hist"  "control"       "etamap"        "formula"       "target.stats"  "target.esteq"  "constraints"   "reference"
  ##[25] "estimate"      "offset"        "drop"          "estimable"     "null.lik"
  
  ## Delete ergm.ego fit object components we will not need
  ## Keeping ppopsize and popsize just in case
  eefit$egodata <- eefit$sample.obs <- eefit$constrained <- NULL
  eefit$v <- eefit$m <- eefit$ergm.formula <- eefit$offset.coef <- NULL
  eefit$ergm.offset.coef <- eefit$ergm.covar <- eefit$DtDe <- NULL
  
  
  # First: The fit component of the netest object
  #        Assign the ergm.ego fit object, and modify as needed to remove
  #        the network size adjustments unique to ergm.ego.
  
  out$fit <- eefit
  class(out$fit) <- "ergm"
  
  
  ## edges coef for ergm = (edges + network size adjustment) for ergm.ego
  out$fit$coef <- eefit$coef[2:length(eefit$coef)]
  out$fit$coef[1] <- out$fit$coef[1] + eefit$coef[1]
  
  ## Sample (modified)
  ncol<-length(eefit$sample[1,])
  out$fit$sample<-eefit$sample[,-1]
  
  # Iterations
  out$fit$iterations <- eefit$iterations
  
  ## MCMCtheta coeficient (modified)
  out$fit$MCMCtheta<-eefit$MCMCtheta[2:length(eefit$MCMCtheta)]
  out$fit$MCMCtheta[1]<-out$fit$MCMCtheta[1] + eefit$MCMCtheta[1]
  
  # Loglikelihood
  out$fit$loglikelihood <- eefit$loglikelihood
  
  ## Gradient (modified).
  out$fit$gradient<-eefit$gradient[-1]
  
  ## Hessian matrix (modified) (r1,c1).
  ncol <- length(eefit$hessian[1,])
  out$fit$hessian <- eefit$hessian[2:ncol,2:ncol]
  
  ## Covar matrix (modified) (r1,c1).
  ncol <- length(eefit$covar[1,])
  out$fit$covar <- eefit$covar[2:ncol,2:ncol]
  
  # Failure
  out$fit$failure <- eefit$failure
  
  # network, newnetworks and newnetwork
  out$fit$network <- eefit$network
  out$fit$newnetworks <- eefit$newnetworks
  out$fit$newnetwork <- eefit$newnetwork
  
  ## coef.init (modified)
  out$fit$coef.init<-eefit$coef.init[2:length(eefit$coef.init)]
  out$fit$coef.init[1]<-out$fit$coef.init[1] + eefit$coef.init[1]
  
  ## est.cov matrix adjusted to remove network size adjustment (r1,c1).
  ncol <- length(eefit$est.cov[1,])
  out$fit$est.cov <- eefit$est.cov[2:ncol,2:ncol]
  
  ## coef.hist (modified)
  ncol <- length(eefit$coef.hist[1,])
  eefit$coef.hist[,2] <- eefit$coef.hist[,1] + eefit$coef.hist[,2]
  out$fit$coef.hist <- eefit$coef.hist[,-1]
  
  ## stats.hist (modified)
  out$fit$stats.hist <- eefit$stats.hist[,-1]
  
  # steplen.hist
  out$fit$steplen.hist <- eefit$steplen.hist
  
  ## control (modified)
  out$fit$control <- eefit$control
  out$fit$control$init <- eefit$control$init[-1]
  out$fit$control$init[1] <- eefit$control$init[1] + eefit$control$init[2]
  
  ## etamap (modified) (first element of each attribute dropped)
  out$fit$etamap <- eefit$etamap
  out$fit$etamap$canonical <- eefit$etamap$canonical[-(length(eefit$etamap$canonical))]
  out$fit$etamap$offsetmap <- eefit$etamap$offsetmap[-1]
  out$fit$etamap$offset <- eefit$etamap$offset[-1]
  out$fit$etamap$offsettheta <- eefit$etamap$offsettheta-1
  out$fit$etamap$curved <- eefit$etamap$curved
  out$fit$etamap$etalength <- eefit$etamap$etalength - 1
  
  # formula, 
  out$fit$formula <- eefit$formula
  
  # Commenting out all the code that set the network -- that
  # should be done in another initialization file, not here.
  
  #  z<-(as.formula("nw ~ o"))
  #  out$fit$formula[2] <- z[2]
  
  # #JKB 4/24/18 - nw object needs to be added to the environment of the 
  # #formation formula; patterned after netest() lines 182-183
  # #DTH 4/46/2018 - testing cause of Males-Male ties
  # # deleting all ties to creat an empty network
  # 
  # #    ego.netest <- function (eefit, disscoef, 
  # #                            set.formula.nw=FALSE, start.net=NULL, 
  # #                            obs.adjust=FALSE, obs.offset = 52)
  # 
  # #JKB 6/5/18: allow for an alternate starting network if using ppop.wt='round'
  # if (is.null(start.net)) {
  #   nw <- eefit$network
  # } else nw <- start.net
  # 
  # count <- network::network.edgecount(nw)
  # network::delete.edges(nw,1:count)
  # if (set.formula.nw) environment(out$fit$formula) <- environment()
  
  # target stats
  
  if(is.na(eefit$target.stats[1])){
    out$fit$target.stats <- eefit$target.stats[-1]
  } else out$fit$target.stats <- eefit$target.stats
  
  # target.esteq
  out$fit$target.esteq <- eefit$target.esteq
  
  # constraints
  out$fit$constraints <- eefit$constraints
  
  # reference
  out$fit$reference <- eefit$reference
  
  # estimate
  out$fit$estimate <- eefit$estimate
  
  ## offset (modified)
  out$fit$offset <- eefit$offset[-1]
  
  ## drop (modified)
  out$fit$drop <- eefit$drop[-1]
  
  ## estimable (modified)
  out$fit$estimable <- eefit$estimable[-1]
  
  # null.lik
  out$fit$null.lik <- eefit$null.lik[-1]
  
  
  # Next: The other netest components (other than fit).
  
  ## the netest object contains these components:
  ##[1] "fit"                "formation"          "target.stats"       "target.stats.names" "coef.form"          "coef.form.crude"    "coef.diss"
  ##[8] "constraints"        "edapprox"
  
  formation <- eefit$formula[-2]
  target.stats <- as.numeric(out$fit$target.stats)
  
  # Remove the offset term(s) from the target stats components
  
  target.stats <- target.stats[!out$fit$offset]
  target.stats.names <- names(out$fit$target.stats)[!out$fit$offset]
  
  # coef.form and coef.form.crude
  
  coef.form <- coef.form.crude <- out$fit$coef
  #  coef.form.crude <- coef.form
  
  ##  Adjust edges coef if partner count obs period != sim time step.
  ##  Edges coefficient must be adjusted to reflect formation 
  ##  rates per simulated timestep (in both fit$coef  and coef.form).
  ##  The default is set to 52 (obs period = 1 yr, sim step = 1 wk)
  
  if (obs.adjust){
    coef.form[1] <- out$fit$coef[1] <- coef.form[1] -log(obs.offset)
    #    out$fit$coef[1] <- out$fit$coef[1] -log(ot.offset)
  }
  
  # for edapprox=T:  Adjust all formation coefs by their dissolution coefs
  # This has been updated with Sam's code from 
  # https://github.com/statnet/EpiModel/blob/master/R/netest.R
  # but this requires careful model specification to ensure
  # order and position of formation coefs = dissolution coefs
  
  edapprox <- TRUE
  
  if (disscoef$coef.crude[1] > -Inf) {
    len.cd <- length(disscoef$coef.crude)
    coef.form[1:len.cd] <- coef.form[1:len.cd] - disscoef$coef.crude
  }
  
  constraints <- eefit$constraints
  
  # Finally, assign the new components to the netest object
  # NOTE: this doesn't include the weights from the egodata object
  
  out$formation <- formation
  out$target.stats <- target.stats
  out$target.stats.names <- target.stats.names
  out$coef.form <- coef.form
  out$coef.form.crude <- coef.form.crude
  out$coef.diss <- disscoef
  out$constraints <- constraints
  out$edapprox <- edapprox
  
  class(out) <- "netest"
  return(out)
}

check.attr <- function(x){
  for(i in 1 : length(dat$attr)){
    print(length(dat$attr[[i]]))
  }
}




