# Predict catch over a TAC period

\`predCatch()\` computes the predicted catch over a TAC period, given
fishing mortality, numbers at age, natural mortality, selectivity, and
other biological and stock–recruitment parameters. If a target TAC is
supplied, the function can also return the difference between the
predicted catch and the specified TAC.

## Usage

``` r
predCatch(
  logFM,
  NAA,
  MAA,
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
  TAC = NULL,
  out = 0
)
```

## Arguments

- logFM:

  Array, matrix, or vector of log fishing mortality values (e.g. by age
  and season) for the TAC period.

- NAA:

  Numbers-at-age at the start of the TAC period (typically a vector or
  matrix indexed by age and possibly season).

- MAA:

  Natural mortality-at-age (and, if applicable, season) used in the
  projection.

- sel:

  Selectivity-at-age (and possibly season) applied to fishing mortality
  when computing catch.

- weight:

  Mean weight-at-age (and possibly season) used to convert numbers
  caught into catch biomass.

- seasons:

  Integer vector with the season indices included in the TAC period
  (e.g. quarters or months within a year).

- ns:

  Integer giving the total number of seasons per year.

- y:

  Integer index of the current year in the projection.

- h2:

  Steepness (or related transformed parameter) of the stock–recruitment
  relationship, used when updating recruitment.

- asmax:

  Integer giving the maximum spawning age (or related age limit) used
  when computing spawning biomass and recruitment.

- mat:

  Maturity-at-age (and possibly season), used to calculate spawning
  biomass.

- pzbm:

  Additional parameter controlling the timing or proportion of biomass
  contributing to spawning (exact role depends on the internal
  implementation).

- spawning:

  Indicator (e.g. logical or integer vector) identifying the spawning
  season(s) within the year.

- R0:

  Unfished recruitment level used in the stock–recruitment relationship.

- SR:

  Character string or code specifying the type of stock–recruitment
  function (e.g. Beverton–Holt, Ricker), as used in the operating model.

- bp:

  Additional parameter (or vector of parameters) for the
  stock–recruitment relationship (e.g. breakpoints or shape parameters).

- recBeta:

  Parameter (or vector) controlling the slope or shape of the
  recruitment function.

- recGamma:

  Parameter (or vector) controlling the curvature or shape of the
  recruitment function.

- eR:

  Stochastic recruitment deviation(s) applied to the expected
  recruitment (e.g. lognormal errors).

- indage0:

  Integer index specifying the recruitment age (age at which new
  recruits enter the age structure).

- TAC:

  Optional numeric value giving the target total allowable catch for the
  TAC period. If \`NULL\` (default), only the predicted catch is
  returned. If non-\`NULL\`, the function can return the difference
  between predicted catch and \`TAC\`, depending on \`out\`.

- out:

  Integer flag controlling the type of output returned. The precise
  meaning depends on the internal implementation but typically
  distinguishes between returning the predicted catch, the difference
  between catch and \`TAC\`, or possibly additional diagnostics.

## Value

A numeric value or vector containing the predicted catch over the TAC
period, or the difference between predicted and target catch if \`TAC\`
is supplied and \`out\` is set accordingly. The exact structure of the
return value depends on the internal implementation. This function is
intended for internal use within IAMSE.

## Details

This is an internal helper used by higher-level IAMSE functions when
searching for fishing mortality levels that match a given TAC or when
projecting catches under a given management rule.
