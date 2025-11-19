# Solve for fishing mortality corresponding to a target TAC

\`get.f()\` finds the fishing mortality (FM) that produces a given total
allowable catch (TAC) over a TAC period, accounting for seasonal
structure, natural mortality, selectivity, and recruitment. Internally,
it repeatedly calls the catch-prediction logic (e.g. \[predCatch()\])
and uses a numerical root-finding or optimisation routine to match
predicted catch to the target \`TAC\`.

## Usage

``` r
get.f(
  TAC,
  NAA,
  MAA,
  sel,
  weight,
  seasons,
  ns,
  y,
  h,
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
  lastFM = 0.01,
  upper = log(100)
)
```

## Arguments

- TAC:

  Numeric value giving the target total allowable catch for the TAC
  period (e.g. over the seasons specified in \`seasons\`).

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

- h:

  Steepness (or related parameter) of the stock–recruitment
  relationship, used when updating recruitment along the projection.

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

- lastFM:

  Numeric value giving the starting (or typical) fishing mortality used
  as an initial guess in the numerical search. Default is \`0.01\`.

- upper:

  Numeric value giving the upper bound (on the log scale) for the search
  over fishing mortality. Default is \`log(100)\`, corresponding to a
  relatively high maximum F.

## Value

A numeric value giving the fishing mortality (or F-multiplier) that
produces a predicted catch matching \`TAC\` as closely as possible under
the given biological and seasonal assumptions.

## Details

This function is typically called internally by higher-level IAMSE
functions when implementing HCRs that are specified in terms of TAC
(rather than directly in terms of fishing mortality).

The function typically:

1.  Defines an objective function measuring the difference between the
    predicted catch (for a given FM) and the target \`TAC\`.

2.  Searches over FM (on a log scale) between a lower bound and
    \`upper\`, using \`lastFM\` as an initial value.

3.  Returns the FM that minimises this difference (or sets the root to
    zero), within numerical tolerance.

The exact optimisation method depends on the internal implementation
(e.g. use of \[stats::optimize()\] or \[stats::uniroot()\]).
