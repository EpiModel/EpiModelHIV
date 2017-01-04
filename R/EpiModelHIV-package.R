#' HIV Transmission Dynamics among MSM and Heterosexuals
#'
#' \tabular{ll}{
#'    Package: \tab EpiModelHIV\cr
#'    Type: \tab Package\cr
#'    Version: \tab 1.0.0\cr
#'    Date: \tab 2016-06-25\cr
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
#' @importFrom stats rbinom rgeom rmultinom rpois runif simulate rnbinom plogis
#' @importFrom Rcpp sourceCpp
#' @importFrom dplyr group_by summarise
#'
#' @useDynLib EpiModelHIV
#'
#' @docType package
#' @keywords package msm het
#'
NULL
