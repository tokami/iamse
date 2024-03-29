# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @name simpop
#'
#' @title Simulate an age-based population
#'
#' @details Simulate a pop
#'
#' @param logF fishing mortality
#' @param dat List with species data
#' @param set List with MSE settings
#' @param tvy year index for all time-variant (tv) processes (so far: Msel, Ms, sel)
#' @param opt If 1 the function returns the yield in the last year,
#' if 2 the function returns a list with yield, TSB, SSB, and ESB over
#' the whole simulation period.
#'
#' @export
NULL

initdist <- function(MAA, FAA, R0, spawning, indage0) {
    .Call('_iamse_initdist', PACKAGE = 'iamse', MAA, FAA, R0, spawning, indage0)
}

simpop <- function(logFM, dat, set, out) {
    .Call('_iamse_simpop', PACKAGE = 'iamse', logFM, dat, set, out)
}

