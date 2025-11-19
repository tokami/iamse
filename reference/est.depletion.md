# Estimate depletion under alternative fishing mortalities

\`est.depletion()\` evaluates the depletion level (typically biomass
relative to unfished biomass, \\B/B_0\\, or a related quantity) for a
range of fishing mortalities. The function simulates the IAMSE operating
model over a grid of fishing mortality values between \`fmin\` and
\`fmax\`, using \`nrep\` stochastic replicates, and summarises the
resulting depletion distribution according to the chosen \`method\`.

## Usage

``` r
est.depletion(
  dat,
  set = NULL,
  fmin = 1e-04,
  fmax = 10,
  nrep = 100,
  verbose = TRUE,
  method = "percentile",
  B0K = FALSE,
  tol = 1e-04,
  do.opt = TRUE
)
```

## Arguments

- dat:

  Data object as returned by \[check.dat()\], containing life-history
  and stock information for the operating model.

- set:

  Optional settings list as returned by \[check.set()\]. If \`NULL\`
  (default), internal defaults are used where possible. When supplied,
  \`set\` controls aspects such as the number of years, projection
  horizon, and noise structure.

- fmin:

  Numeric value giving the minimum fishing mortality (or F-multiplier)
  considered in the depletion calculation. Default is \`0.0001\`.

- fmax:

  Numeric value giving the maximum fishing mortality considered in the
  depletion calculation. Default is \`10\`.

- nrep:

  Integer specifying the number of stochastic replicates used for each
  fishing mortality level. Larger values give more stable depletion
  estimates at the cost of increased computation time. Default is
  \`100\`.

- verbose:

  Logical; if \`TRUE\`, print progress messages and basic diagnostics
  during the calculation. Default is \`TRUE\`.

- method:

  Character string indicating how the depletion distribution is
  summarised across replicates. The default \`"percentile"\` typically
  uses quantiles of the simulated depletion distribution. Other methods
  may be available depending on the internal implementation.

- B0K:

  Logical; if \`TRUE\`, treat the carrying capacity \\K\\ as equivalent
  to unfished biomass \\B_0\\ when computing depletion. If \`FALSE\`
  (default), \\B_0\\ is used explicitly if available.

- tol:

  Numeric tolerance used in any optimisation or root-finding procedures
  (e.g. when solving for the fishing mortality associated with a target
  depletion level). Default is \`0.0001\`.

- do.opt:

  Logical; if \`TRUE\` (default), perform optimisation to identify
  fishing mortality values corresponding to target depletion levels. If
  \`FALSE\`, the function only evaluates depletion over the grid defined
  by \`fmin\` and \`fmax\` without attempting to optimise.

## Value

An object (typically a list or data frame) summarising depletion as a
function of fishing mortality, including summary statistics across
replicates. The exact structure depends on the internal implementation.
The object is returned invisibly if the function is mainly used for its
side effects (e.g. storing results in \`dat\`).

## Details

Optionally, an optimisation step can be performed to identify the
fishing mortality that leads to a target depletion level, controlled by
the convergence tolerance \`tol\` and the flag \`do.opt\`.
