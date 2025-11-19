# Plot trade-offs between IAMSE performance metrics

\`plotiamse.tradeoff()\` produces a two-dimensional trade-off plot for
performance metrics summarising IAMSE management strategy evaluation
(MSE) results. Each point usually represents a harvest control rule
(HCR), positioned according to two selected metrics (e.g. average catch
vs. risk of falling below a biomass limit).

## Usage

``` r
plotiamse.tradeoff(
  mets,
  metrics = c("avRelCatch", "PBBlim"),
  hcrs = NA,
  plot.legend = TRUE,
  xlim = NULL,
  ylim = NULL
)
```

## Arguments

- mets:

  Data frame or list returned by \[est.metrics()\], containing
  performance metrics for each HCR (and possibly for each stock or
  scenario). The object must include columns with the names specified in
  \`metrics\`.

- metrics:

  Character vector of length two specifying the names of the metrics to
  plot on the x- and y-axes, respectively. Default is \`c("avRelCatch",
  "PBBlim")\`. See \[est.metrics()\] for available metric names.

- hcrs:

  Optional character vector or numeric indices specifying which HCRs to
  include in the trade-off plot. If \`NA\` (default), all HCRs in
  \`mets\` are plotted.

- plot.legend:

  Logical; if \`TRUE\` (default), draw a legend identifying the plotted
  HCRs.

- xlim:

  Optional numeric vector of length two giving the x-axis limits. If
  \`NULL\` (default), limits are chosen automatically based on the
  selected metric.

- ylim:

  Optional numeric vector of length two giving the y-axis limits. If
  \`NULL\` (default), limits are chosen automatically based on the
  selected metric.

## Value

This function is called for its side effect of producing a plot.
Invisibly returns \`NULL\`.

## Details

The function takes the output from \[est.metrics()\] and plots one
metric on the x-axis and another on the y-axis, optionally with a legend
identifying the different HCRs.
