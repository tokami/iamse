# Plot fishing mortality time series from an IAMSE MSE

\`plotiamse.f()\` produces time-series plots of fishing mortality for
each harvest control rule (HCR) evaluated in an IAMSE management
strategy evaluation. The function aggregates results across replicates
and can display median trajectories, uncertainty envelopes, and optional
smoothed trend lines.

## Usage

``` r
plotiamse.f(
  dat,
  set,
  resMSE,
  trendline = TRUE,
  uncert = TRUE,
  med = TRUE,
  hcrs = NA,
  ylim = NULL,
  plot.legend = TRUE,
  ylab = "Fishing mortality",
  xlab = ""
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

- trendline:

  Logical; if \`TRUE\`, draw a smoothed trend or additional line
  highlighting the general trajectory of fishing mortality over time
  (e.g. using a smoother or running mean). Default is \`TRUE\`.

- uncert:

  Logical; if \`TRUE\`, display uncertainty across replicates, typically
  in the form of quantile bands or envelopes around the median. Default
  is \`TRUE\`.

- med:

  Logical; if \`TRUE\`, draw the median fishing mortality trajectory
  across replicates for each HCR. If \`FALSE\`, only uncertainty or
  other summaries are shown, depending on the internal implementation.
  Default is \`TRUE\`.

- hcrs:

  Optional character vector or indices specifying which HCRs to plot. If
  \`NA\` (default), all HCRs present in \`set\$hcr\` (and \`resMSE\`)
  are plotted.

- ylim:

  Optional numeric vector of length two giving the y-axis limits. If
  \`NULL\` (default), limits are chosen automatically to encompass the
  plotted fishing mortality values.

- plot.legend:

  Logical; if \`TRUE\` (default), add a legend identifying the plotted
  HCRs.

- ylab:

  Character string specifying the y-axis label. Default is \`"Fishing
  mortality"\`.

- xlab:

  Character string specifying the x-axis label. Default is an empty
  string \`""\`. Often set to \`"Year"\` or similar by the user.

## Value

This function is called for its side effect of producing a plot.
Invisibly returns \`NULL\`.

## Details

By default, one line per HCR is drawn, and uncertainty across stochastic
replicates is indicated (for example, by plotting quantile bands; the
exact representation depends on the internal implementation). The
function uses base R graphics.
