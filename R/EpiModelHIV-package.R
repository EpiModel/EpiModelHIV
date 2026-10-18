#' Network-Based Epidemic Modeling of HIV/STI Transmission
#'
#' \tabular{ll}{
#'    Package: \tab EpiModelHIV\cr
#'    Type: \tab Package\cr
#'    Version: \tab 2.6.0\cr
#'    Date: \tab 2022-12-21\cr
#'    License: \tab GPL-3\cr
#'    LazyLoad: \tab yes\cr
#' }
#'
#' @details
#' EpiModelHIV provides extensions to our general EpiModel package to allow
#' for simulating HIV/STI transmission dynamics.
#'
#' @name EpiModelHIV-package
#' @aliases EpiModelHIV
#'
#' @import EpiModel network networkDynamic tergm ergm
#' @importFrom stats rbinom rgeom rmultinom rpois runif simulate rnbinom plogis
#'                   predict
#' @importFrom utils tail
#' @importFrom statnet.common set.control.class
#' @importFrom dplyr select filter
#' @importFrom networkLite networkLite
#' @importFrom rlang .data
#'
#' @keywords package msm het
#'
"_PACKAGE"

# Avoids NOTE about Rglpk import
f1 <- function() {
  Rglpk::Rglpk_read_file(
    system.file(file.path("examples", "plan.lp"), package = "Rglpk"),
    type = "CPLEX_LP"
  )
}
