#' Network-Based Epidemic Modeling of HIV/STI Transmission
#'
#' \tabular{ll}{
#'    Package: \tab EpiModelHIV\cr
#'    Type: \tab Package\cr
#'    Version: \tab 2.0.5\cr
#'    Date: \tab 2020-10-08\cr
#'    License: \tab GPL-3\cr
#'    LazyLoad: \tab yes\cr
#' }
#'
#' @details
#' EpiModelHIV provides extensions to our general EpiModel package to allow
#' for simulating HIV transmission dynamics among two populations: men who
#' have sex with men (MSM) in the United States and heterosexual adults in
#' sub-Saharan Africa.
#'
#' @name EpiModelHIV-package
#' @aliases EpiModelHIV
#'
#' @import EpiModel EpiModelHPC network networkDynamic tergmLite tergm ergm bindata
#' @importFrom stats rbinom rgeom rmultinom rpois runif simulate rnbinom plogis predict
#' @importFrom dplyr group_by summarise
#' @importFrom utils tail
#'
#' @docType package
#' @keywords package msm het
#'
NULL
