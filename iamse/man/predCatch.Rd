% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{predCatch}
\alias{predCatch}
\title{predCatch}
\usage{
predCatch(
  logFM,
  NAA,
  MAAy,
  sel,
  weight,
  seasons,
  ns,
  y,
  h2,
  asmax,
  mat,
  pzbm,
  spawning,
  R0,
  SR,
  bp,
  recBeta,
  recGamma,
  eR,
  indage0,
  iaFM,
  TAC = NULL,
  out = 0
)
}
\arguments{
\item{seasons}{vector with season indices}

\item{ns}{number of seasons}

\item{h}{steepness}
}
\description{
predCatch
}
\details{
get predicted catch for TAC period or difference between provided
    and predicted catch
}
