# Estimate stochastic reference points

\`est.ref.levels.stochastic()\` estimates reference points such as
\\F\_{MSY}\\, \\B\_{MSY}\\, MSY, and unfished biomass \\B_0\\ using
stochastic simulations. For each candidate fishing mortality level (up
to \`fmax\`), the operating model is simulated over many replicates,
taking into account process, observation, and parameter uncertainty as
specified in \`set\`. Reference points are then derived from the
distribution of simulated long-term outcomes.

## Usage

``` r
est.ref.levels.stochastic(
  dat,
  set = NULL,
  fmax = 10,
  ncores = parallel::detectCores() - 1,
  ref = c("Fmsy", "Bmsy", "MSY", "B0", "ESBmsy", "SSBmsy", "ESB0"),
  plot = FALSE,
  get.final = FALSE
)
```

## Arguments

- dat:

  Data object as returned by \[check.dat()\], containing life-history
  and stock information for one or more stocks. The list element
  \`dat\$ref\` is created or updated with the estimated stochastic
  reference levels.

- set:

  Settings list as returned by \[check.set()\]. This should define the
  number of stochastic replicates used for the reference-point
  calculation (e.g. \`set\$refN\`) and the noise structure (e.g.
  \`set\$noise\`). If \`NULL\` (default), internal defaults are used
  where possible.

- fmax:

  Numeric value giving the maximum fishing mortality (or F-multiplier)
  considered when searching for MSY-based reference points. The internal
  grid of F values typically runs from 0 to \`fmax\`. Default is \`10\`.

- ncores:

  Integer giving the number of CPU cores to use for parallel
  computation. Default is \`parallel::detectCores() - 1\`. On non-Unix
  systems (e.g. Windows), \`mclapply()\` falls back to serial execution
  regardless of \`ncores\`.

- ref:

  Character vector specifying which reference points to estimate.
  Typical entries include:

  - \`"Fmsy"\`: Fishing mortality at MSY.

  - \`"Bmsy"\`: Biomass at MSY.

  - \`"MSY"\`: Maximum sustainable yield.

  - \`"B0"\`: Unfished (virgin) biomass.

  - \`"ESBmsy"\`: Equilibrium spawning biomass at F\\\_{MSY}\\.

  - \`"SSBmsy"\`: Spawning stock biomass at F\\\_{MSY}\\.

  - \`"ESB0"\`: Equilibrium spawning biomass at F = 0.

  The default is \`c("Fmsy", "Bmsy", "MSY", "B0", "ESBmsy", "SSBmsy",
  "ESB0")\`.

- plot:

  Logical; if \`TRUE\`, produce diagnostic plots showing, for example,
  simulated yield and biomass as a function of fishing mortality and the
  resulting distribution of reference-point estimates. Default is
  \`FALSE\`.

- get.final:

  Logical; if \`TRUE\`, return additional information on the final
  simulation results used to derive the reference points (e.g.
  replicate-level outcomes at the estimated \\F\_{MSY}\\). If \`FALSE\`
  (default), only summary reference levels are stored in \`dat\$ref\`
  and returned.

## Value

The updated \`dat\` object, with \`dat\$ref\` containing stochastic
reference levels for each stock. If \`get.final = TRUE\`, the return
value may include additional elements with the underlying simulation
results used to construct the reference points. The object is also
returned invisibly.

## Details

Compared to \[est.ref.levels()\], which is deterministic, this function
provides stochastic, simulation-based estimates that are more consistent
with the full uncertainty structure used in IAMSE MSE runs. Parallel
computation across stocks or fishing mortality levels is supported via
\[parallel::mclapply()\] on Unix-like systems.

The function typically:

1.  For each stock, simulates the operating model over a range of F
    values from 0 to \`fmax\`, using the uncertainty specification in
    \`set\`.

2.  Aggregates long-term outcomes across replicates for each F.

3.  Identifies the F that maximises yield (and associated biomass
    levels) in a stochastic sense (e.g. using medians across
    replicates).

4.  Stores the requested stochastic reference points in \`dat\$ref\`.

Parallelisation is implemented via \[parallel::mclapply()\], which uses
forking on Unix-like systems. On Windows, the computation is executed in
serial even if \`ncores \> 1\`.
