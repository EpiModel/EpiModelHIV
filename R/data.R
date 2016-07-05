
#' Fitted Network Models for 10000 Person Network
#'
#' A dataset containing a list of three outputs from \code{netest}, one for each
#' network. This would be used primarily for more efficient package testing.
#'
#' @docType data
#' @keywords datasets
#' @format A list with three elements.
#' @name est
NULL

#' Network Statistics Network Models for 10000 Person Network
#'
#' A dataset containing the necessary network statistics for efficient
#' package testing.
#'
#' @docType data
#' @keywords datasets
#' @format A list with 18 elements.
#' @name st
NULL

#' Life Table for Ghana
#'
#' A dataset containing age-sex-period specific mortality rates.
#'
#' \itemize{
#'   \item \strong{year}: calendar year for rates, with options of 1990, 2000,
#'         and 2011.
#'   \item \strong{male}: binary indicator for male (vs female) rates.
#'   \item \strong{agStart}: lower age interval at which the rate applies.
#'   \item \strong{agEnd}: upper age interval at which the rate applies.
#'   \item \strong{mrate}: the mortality rate.
#' }
#'
#' @docType data
#' @keywords datasets
#' @format A data frame with 132 rows and 5 variables.
#' @name ltGhana
NULL
