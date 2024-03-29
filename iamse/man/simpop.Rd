% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{simpop}
\alias{simpop}
\title{Simulate an age-based population}
\arguments{
\item{logF}{fishing mortality}

\item{dat}{List with species data}

\item{set}{List with MSE settings}

\item{tvy}{year index for all time-variant (tv) processes (so far: Msel, Ms, sel)}

\item{opt}{If 1 the function returns the yield in the last year,
if 2 the function returns a list with yield, TSB, SSB, and ESB over
the whole simulation period.}
}
\description{
Simulate an age-based population
}
\details{
Simulate a pop
}
