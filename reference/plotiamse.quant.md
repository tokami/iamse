# Plot IAMSE diagnostic / uncertainty quantities

\`plotiamse.quant()\` visualises selected diagnostic or
uncertainty-related quantities computed from an IAMSE management
strategy evaluation. Typical examples include the variability of biomass
and fishing mortality relative to reference points, or the variability
of catch over time.

## Usage

``` r
plotiamse.quant(
  dat,
  set,
  resMSE,
  hcrs = NULL,
  quants = c("bpbmsy.sd", "fmfmsy.sd", "cp.sd"),
  hline = NULL,
  plot.legend = TRUE
)
```

## Arguments

- dat:

  Data object as returned by \[check.dat()\], containing life-history
  and reference-point information used in the MSE.

- set:

  Settings list as returned by \[check.set()\] and modified for the MSE
  run (e.g. containing \`nysim\`, \`nrep\`, and \`hcr\`).

- resMSE:

  Result object returned by \[run.mse()\], containing simulated time
  series for each HCR and replicate.

- hcrs:

  Optional character vector or numeric indices specifying which HCRs to
  include in the plot. If \`NULL\` (default), all HCRs present in
  \`set\$hcr\` (and \`resMSE\`) are used.

- quants:

  Character vector specifying the names of the quantities to plot.
  Defaults to \`c("bpbmsy.sd", "fmfmsy.sd", "cp.sd")\`, which typically
  correspond to the standard deviation of biomass relative to
  B\\\_{MSY}\\, fishing mortality relative to F\\\_{MSY}\\, and catch
  performance (or a related quantity). The available names depend on how
  summary statistics were computed upstream.

- hline:

  Optional numeric value or vector giving horizontal reference line(s)
  to be added to the plot (e.g. a threshold or target level). If
  \`NULL\` (default), no reference lines are drawn.

- plot.legend:

  Logical; if \`TRUE\` (default), add a legend identifying the plotted
  HCRs.

## Value

This function is called for its side effect of producing a plot.
Invisibly returns \`NULL\`.

## Details

The function extracts the specified quantities for each harvest control
rule (HCR) and produces a comparative plot (e.g. one panel per quantity
or grouped bars by HCR, depending on the internal implementation).
Optional horizontal reference lines can be added to facilitate
interpretation.
