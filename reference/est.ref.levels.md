# Estimate deterministic reference points

\`est.ref.levels()\` estimates deterministic reference points for each
stock in an IAMSE data object over a grid of fishing mortality values.
Typical reference points include \\F\_{MSY}\\, \\B\_{MSY}\\, MSY, and
unfished biomass \\B_0\\. The function evaluates the operating model at
a range of fishing mortalities (\`fvec\`), identifies the values that
maximise long-term yield, and stores the resulting reference levels in
\`dat\$ref\`.

## Usage

``` r
est.ref.levels(
  dat,
  set = NULL,
  fvec = seq(0, 5, 0.1),
  ncores = parallel::detectCores() - 1,
  ref = c("Fmsy", "Bmsy", "MSY", "ESBmsy", "SSBmsy", "B0"),
  plot = FALSE
)
```

## Arguments

- dat:

  Data object as returned by \[check.dat()\], containing life-history
  and stock information for one or more stocks. The list element
  \`dat\$ref\` is created or updated with the estimated reference
  levels.

- set:

  Optional settings list as returned by \[check.set()\]. If \`NULL\`
  (default), internal defaults are used where needed. When provided,
  \`set\` can control aspects of the operating model and reference-point
  calculation.

- fvec:

  Numeric vector giving the grid of fishing mortality values (e.g. \\F\\
  or \\F/F\_{MSY}\\) over which yield and biomass are evaluated. Default
  is \`seq(0, 5, 0.1)\`.

- ncores:

  Integer giving the number of CPU cores to use for parallel computation
  across stocks. Default is \`parallel::detectCores() - 1\`. On non-Unix
  systems (e.g. Windows), \`mclapply()\` falls back to serial execution
  regardless of \`ncores\`.

- ref:

  Character vector specifying which reference points to estimate.
  Typical entries include:

  - \`"Fmsy"\`: Fishing mortality at MSY.

  - \`"Bmsy"\`: Biomass at MSY.

  - \`"MSY"\`: Maximum sustainable yield.

  - \`"ESBmsy"\`: Equilibrium spawning biomass at F\\\_{MSY}\\.

  - \`"SSBmsy"\`: Spawning stock biomass at F\\\_{MSY}\\.

  - \`"B0"\`: Unfished (virgin) biomass.

  The default is \`c("Fmsy", "Bmsy", "MSY", "ESBmsy", "SSBmsy", "B0")\`.

- plot:

  Logical; if \`TRUE\`, produce diagnostic plots of the relationship
  between fishing mortality and yield/biomass (e.g. equilibrium yield
  curves, biomass vs. F). Default is \`FALSE\`.

## Value

The updated \`dat\` object, with the component \`dat\$ref\` containing a
table or list of estimated reference levels for each stock. The object
is also returned invisibly.

## Details

This function provides a fast, deterministic alternative to
\[est.ref.levels.stochastic()\], which uses stochastic simulations. It
can use multiple cores via \`parallel::mclapply()\` on Unix-like
systems.

The function typically:

1.  Loops over stocks defined in \`dat\`.

2.  For each stock, evaluates equilibrium quantities over \`fvec\`.

3.  Identifies \\F\_{MSY}\\ and associated biomass and yield metrics.

4.  Stores the requested reference points in \`dat\$ref\`.

Parallelisation is implemented via \[parallel::mclapply()\], which uses
forking on Unix-like systems. On Windows, the computation is executed in
serial even if \`ncores \> 1\`.
