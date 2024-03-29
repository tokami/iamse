% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{recfunc}
\alias{recfunc}
\title{recfunc}
\usage{
recfunc(
  h,
  SPR0,
  SSB,
  R0 = 1e+06,
  method = "bevholt",
  bp = 0,
  beta = 0,
  gamma = 0
)
}
\arguments{
\item{h}{- steepness}

\item{SPR0}{- spawning biomass produced by one recrut in its lifetime}

\item{SSB}{- spawning biomass}

\item{R0}{- recruitment in unfished population}

\item{method}{- SR type}

\item{bp}{- breakpoint for hockey-stick SR}
}
\description{
Function to calculate recruitment (Beverton - Holt)
}
