
#' @title deprecated do not use
#' @details This function has been deprecated. Use `param.net` instead.
#' See https://github.com/EpiModel/EpiModelHIV-Template for details.
#'
#' @export
param_msm <- function() {
  stop(
    "This function has been deprecated. Use `param.net` instead.",
    "See https://github.com/EpiModel/EpiModelHIV-Template for details."
  )
}


#' @title Epidemic Model Initial Conditions
#'
#' @description Sets the initial conditions for a stochastic epidemic models
#'              simulated with [`netsim`].
#'
#' @param prev.gono.uret Initial prevalence of urethral gonorrhea.
#' @param prev.gono.rect Initial prevalence of rectal gonorrhea.
#' @param prev.chla.uret Initial prevalence of urethral chlamydia.
#' @param prev.chla.rect Initial prevalence of rectal chlamydia.
#' @param prev.syph Initial prevalence of syphilis.
#' @param ... Additional arguments passed to function.
#'
#' @return
#' A list object of class `init.net`, which can be passed to EpiModel function [`netsim`].
#'
#' @export
#'
init_msm <- function(prev.gono.rect,
                     prev.gono.uret,
                     prev.chla.rect,
                     prev.chla.uret,
                     prev.syph,
                     ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  class(p) <- "init.net"
  return(p)
}


#' @title Epidemic Model Control Settings
#'
#' @description Sets the controls for stochastic network models simulated with [`netsim`].
#'
#' @param simno Unique ID for the simulation run, used for file naming purposes if used
#'        with the `EpiModelHPC` package.
#' @param nsims Number of simulations to run overall.
#' @param ncores Number of parallel cores per run, if parallel processing is used.
#' @param nsteps Number of time steps per simulation.
#' @param start Starting time step for simulation, with default to `1` to run new simulation. This
#'        may also be set to 1 greater than the final time step of a previous simulation to resume
#'        the simulation with different parameters.
#' @param cumulative.edgelist If `TRUE`, models tracks a cumulative edgelist (a historical list of
#'        all partners within each network over time), useful for partner services models.
#' @param truncate.el.cuml Number of time steps of the cumulative edgelist to retain. See help file
#'        for `EpiModel::update_cumulative_edgelist` for details
#' @param initialize.FUN Module function to use for initialization of the epidemic model.
#' @param aging.FUN Module function for aging.
#' @param departure.FUN Module function for general and disease-related departures.
#' @param arrival.FUN Module function for entries into the sexually active population.
#' @param hivtest.FUN Module function for HIV diagnostic disease testing.
#' @param hivtx.FUN Module function for ART initiation and adherence.
#' @param prep.FUN Module function for PrEP initiation and utilization.
#' @param hivprogress.FUN Module function for HIV disease progression.
#' @param hivvl.FUN Module function for HIV viral load evolution.
#' @param resim_nets.FUN Module function for network resimulation at each time
#'        step.
#' @param summary_nets.FUN Module function for network statistic extraction at
#'        each time step.
#' @param acts.FUN Module function to simulate the number of sexual acts within
#'        partnerships.
#' @param condoms.FUN Module function to simulate condom use within acts.
#' @param position.FUN Module function to simulate sexual position within acts.
#' @param hivtrans.FUN Module function to stochastically simulate HIV
#'        transmission over acts given individual and dyadic attributes.
#' @param sti.FUN Module function calling the function responsibles for
#'        transmission, progression, screening and treatment of Gonorrhea,
#'        Chlamydia and Syphilis.
#' @param prev.FUN Module function to calculate prevalence summary statistics.
#' @param verbose.FUN Module function to print model progress to the console or
#'        external text files.
#' @param cleanup.FUN Module function to tidy the `netsim_dat` object at the end
#'        of each step.
#' @param save.nwstats Calculate and save network statistics as defined in the
#'        `simnet` modules.
#' @param nwstats.formula Right-hand side formula for network statistics to
#'        monitor, with default `"formation"` equal to the model formula.
#'        Supports \code{\link{multilayer}} specification.
#' @param tergmLite Logical indicating usage of `tergmLite` (non-modifiable
#'        currently for `EpiModelHIV`).
#' @param save.network If `TRUE`, networkLite object is saved at simulation end.
#' @param tergmLite.track.duration If `TRUE`, track duration information
#'        for models in `tergmLite` simulations. Supports
#'        \code{\link{multilayer}} specification.
#' @param save.other A character vector of elements on the main `netsim_dat`
#'        data object to save out after each simulation.
#' @param verbose If `TRUE`, print out simulation progress to the console
#'        if in interactive mode or text files if in batch mode.
#' @param set.control.ergm Control arguments passed to `ergm`'s
#'        `simulate_formula.network`. Supports \code{\link{multilayer}}
#'        specification.
#' @param set.control.tergm Control arguments passed to `tergm`'s
#'        `simulate_formula.network`. Supports \code{\link{multilayer}}
#'        specification.
#' @param dat.updates See the general documentation for `dat.updates` at [`EpiModel::control.net`].
#'        The standard implementation for the current MSM model is found in
#'        [`standard_resimnet_updates`].
#' @param debug Enables the debugging Epi Trackers (prefixed with `dbg_`)
#' @param renal.FUN Module responsible of handling the renal function
#' @param hepatic.FUN Module responsible of handling the hepatic function
#' @param prep_otc.FUN Module for the Over The Counter (OTC) PrEP
#' @param tx_resist.FUN Module for handling the acquisition of HIV cART
#'        resistance mutations
#' @param ... Additional arguments passed to the function.
#'
#' @inheritParams EpiModel::control.net
#'
#' @return
#' A list object of class `control.net`, which can be passed to the EpiModel function [`netsim`].
#'
#' @export
#'
control_msm <- function(simno = 1,
                        nsims = 1,
                        ncores = 1,
                        nsteps = 100,
                        start = 1,
                        cumulative.edgelist = TRUE,
                        truncate.el.cuml = 0,
                        debug = FALSE, # enables the debugging trackers
                        #### Modules List ####
                        initialize.FUN = initialize_msm,
                        # demography
                        aging.FUN = aging_msm,
                        departure.FUN = departure_msm,
                        arrival.FUN = arrival_msm,
                        renal.FUN = renal_msm,
                        hepatic.FUN = hepatic_msm,
                        # network
                        resim_nets.FUN = EpiModel::resim_nets,
                        summary_nets.FUN = EpiModel::summary_nets,
                        # sexual activity - make the act list
                        acts.FUN = acts_msm,
                        condoms.FUN = condoms_msm,
                        position.FUN = position_msm,
                        # HIV epidemic
                        prep.FUN = prep_msm,
                        prep_otc.FUN = prep_otc_msm,
                        hivtest.FUN = hivtest_msm,
                        hivtx.FUN = hivtx_msm,
                        hivprogress.FUN = hivprogress_msm,
                        hivvl.FUN = hivvl_msm,
                        hivtrans.FUN = hivtrans_msm,
                        tx_resist.FUN = tx_resist_msm,
                        # STI epidemics
                        sti.FUN = sti_msm,
                        # end of step
                        prev.FUN = prevalence_msm,
                        verbose.FUN = EpiModel::verbose.net,
                        cleanup.FUN = cleanup_msm,
                        ######################
                        save.nwstats = FALSE,
                        nwstats.formula = "formation",
                        tergmLite = TRUE,
                        save.network = FALSE,
                        tergmLite.track.duration = FALSE,
                        set.control.ergm = control.simulate.formula(
                          MCMC.burnin = 2e5
                        ),
                        set.control.tergm = control.simulate.formula.tergm(
                          MCMC.burnin.min = 5000
                        ),
                        save.run = TRUE,
                        save.other = c(),
                        save.cumulative.edgelist = FALSE,
                        verbose = TRUE,
                        dat.updates = standard_resimnet_updates,
                        ...) {
  formal.args <- formals(sys.function())
  dot.args <- list(...)
  p <- get_args(formal.args, dot.args)

  p$skip.check <- TRUE
  p$save.transmat <- FALSE

  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
                                  USE.NAMES = FALSE) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)
  p[["f.names"]] <- c(p[["bi.mods"]], p[["user.mods"]])

  if (is.null(p$verbose.int)) {
    p$verbose.int <- 1
  }

  p$resimulate.network <- TRUE
  p$save.diss.stats <- FALSE

  p <- set.control.class("control.net", p)
  return(p)
}

#' @title Standard Updates to Network Attributes During Resimulation
#'
#' @description This function is updates the network-related degree attributes on the three-layer
#' MSM sexual network that account for the resimulated network structure.
#'
#' @inheritParams aging_msm
#' @param network Integer for network number (values 1, 2, or 3 for main, casual, and one-time
#'                networks).
#'
#' @details
#' This function is called between network resimulations in [`EpiModel::resim_nets`], passed into
#' [`control_msm`] through the `dat.updates` argument. This implementation updates degree attributes
#' calculated from the current network snapshot for use as ERGM terms in the other network layers
#' (e.g., degree in the casual network is a function of the degree in the main network). See the
#' general documentation for `dat.updates` at [`EpiModel::control.net`].
#'
#'
#' @export
#'
standard_resimnet_updates <- function(dat, at, network) {
  if (network == 0L) {
    dat <- set_attr(dat, "deg.casl", EpiModel::get_degree(dat$run$el[[2]]))
  } else if (network == 1L) {
    dat <- set_attr(dat, "deg.main", EpiModel::get_degree(dat$run$el[[1]]))
  } else if (network == 2L) {
    dat <- set_attr(dat, "deg.casl", EpiModel::get_degree(dat$run$el[[2]]))
    dat <- set_attr(dat, "deg.tot", pmin(get_attr(dat, "deg.main") +
                                           EpiModel::get_degree(dat$run$el[[2]]),
                                         3))
  }
  return(dat)
}
