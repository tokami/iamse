# Estimate productivity curve

\`est.productivity()\` estimates a productivity relationship for the
IAMSE operating model, typically linking biomass (or depletion) to
surplus production or long-term yield. The operating model is simulated
over a range of fishing mortalities up to \`fmax\` for \`ny\` years, and
summary quantities are extracted to characterise the stock's
productivity.

## Usage

``` r
est.productivity(
  dat,
  set = NULL,
  ny = 100,
  fmax = 10,
  nf = 1000,
  tsSplit = 8,
  plot = TRUE
)
```

## Arguments

- dat:

  Data object as returned by \[check.dat()\], containing life-history
  and stock information for one or more stocks.

- set:

  Optional settings list as returned by \[check.set()\]. If \`NULL\`
  (default), internal defaults are used where possible. When supplied,
  \`set\` controls aspects such as the number of replicates, noise
  structure, and other operating-model options.

- ny:

  Integer giving the number of years to simulate when estimating the
  productivity relationship. Larger values allow the system to approach
  equilibrium or typical dynamic behaviour under each fishing level.
  Default is \`100\`.

- fmax:

  Numeric value giving the maximum fishing mortality (or F-multiplier)
  considered when constructing the productivity curve. Default is
  \`10\`.

- nf:

  Integer giving the number of fishing-mortality levels (or points along
  the F gradient) to consider between 0 and \`fmax\`. Default is \`1e3\`
  (1000 points), providing a relatively fine grid.

- tsSplit:

  Integer controlling how the simulated time series are split into
  segments when summarising productivity (for example, discarding an
  initial burn-in and using the remaining periods for analysis). The
  precise interpretation depends on the internal implementation. Default
  is \`8\`.

- plot:

  Logical; if \`TRUE\` (default), produce diagnostic plots of the
  productivity relationship (e.g. yield or surplus production vs.
  biomass or F). If \`FALSE\`, no plots are produced.

## Value

An object (typically a list or data frame) containing productivity
summaries as a function of fishing mortality and/or biomass for each
stock. The exact structure depends on the internal implementation. The
object is usually returned invisibly if the function is called primarily
for its side effects (e.g. plotting).

## Details

The function can be used to construct productivity curves or clouds that
relate state variables (e.g. biomass, depletion) to production, which
can then be analysed or plotted for diagnostics and comparison across
stocks.
