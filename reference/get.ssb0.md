# Calculate spawning stock biomass at F = 0 (SSB0)

\`get.ssb0()\` computes spawning stock biomass under zero fishing
mortality (often referred to as SSB0), based on natural mortality,
maturity, weight, and fecundity-at-age, accounting for seasonal
structure and the timing of spawning. Optionally, a non-zero fishing
mortality (\`FM\`) can be supplied to obtain spawning biomass under that
F instead.

## Usage

``` r
get.ssb0(
  M,
  mat,
  weight,
  fecun = 1,
  asmax,
  ns,
  spawning,
  R0,
  indage0,
  season,
  FM = NULL
)
```

## Arguments

- M:

  Natural mortality-at-age (and, if applicable, season), typically given
  as a vector or matrix of instantaneous natural mortality rates.

- mat:

  Maturity ogive, giving the proportion mature at each age (and season,
  if relevant).

- weight:

  Mean weight-at-age (and possibly season), used to convert numbers of
  mature fish into spawning biomass.

- fecun:

  Fecundity-at-age (and possibly season). Defaults to \`1\`, in which
  case fecundity is assumed to be proportional to weight and maturity
  only.

- asmax:

  Integer giving the maximum age class (or the number of age classes)
  used in the calculation.

- ns:

  Integer giving the number of seasons per year.

- spawning:

  Indicator (e.g. logical or integer vector of length \`ns\`)
  identifying the spawning season(s) within the year.

- R0:

  Unfished recruitment level used to scale the numbers-at-age under F
  = 0. Depending on the implementation, this may be used to convert
  per-recruit quantities into absolute biomass.

- indage0:

  Integer index specifying the recruitment age (age at which new
  recruits enter the age structure).

- season:

  Integer or vector specifying the season(s) for which SSB is evaluated;
  typically aligned with \`spawning\` and the seasonal structure.

- FM:

  Optional fishing mortality-at-age (and season). If \`NULL\` (default),
  fishing mortality is set to zero and SSB0 is computed. If
  non-\`NULL\`, the function returns spawning biomass under the combined
  mortality \`M + FM\`.

## Value

A numeric value giving spawning stock biomass under F = 0 (SSB0), or
under the supplied \`FM\` if not \`NULL\`. The exact scaling (per
recruit vs. absolute biomass) depends on how \`R0\` is used internally.
