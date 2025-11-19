#' iamse: Intra-Annual Management Strategy Evaluation (IAMSE)
#'
#' Provides a framework for intra-annual age-based operating models used in
#' management strategy evaluation (MSE) of exploited fish stocks.
#'
#' The package includes tools to:
#' \itemize{
#'   \item define and simulate age-based operating models at intra-annual
#'         time steps,
#'   \item specify harvest control rules (HCRs) and management procedures,
#'   \item run management strategy evaluations (MSEs), and
#'   \item summarise and visualise performance metrics and trade-offs.
#' }
#'
#' @name iamse
#' @aliases iamse iamse-package
#' @author T.K. Mildenberger
#'
#' @useDynLib iamse, .registration = TRUE
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"
