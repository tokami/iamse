# Define an index-based harvest control rule

\`def.hcr.index()\` defines one or more harvest control rules (HCRs)
that use changes in a relative biomass index to set catch advice (TAC)
or fishing mortality. The function constructs HCR functions and assigns
them to the specified environment, typically the global environment, so
they can be used within an MSE framework.

## Usage

``` r
def.hcr.index(
  id = "r23",
  x = 2,
  y = 3,
  stab = TRUE,
  lower = 0.8,
  upper = 1.2,
  clType = "TAC",
  clyears = 1,
  red = NA,
  redyears = 3,
  redAlways = FALSE,
  ffmsySD = 0,
  bbtriggerSD = 0,
  rightRef = 1,
  assessmentInterval = 1,
  env = globalenv()
)
```

## Arguments

- id:

  Character vector giving the identifier(s) for the HCR(s) (e.g.
  \`"r23"\`). These names are used when referring to the HCR in the MSE
  setup and outputs.

- x:

  Numeric scalar or vector used to parameterise the index-based rule.
  Together with \`y\`, it controls how changes in the biomass index
  translate into changes in TAC (e.g. slope, threshold, or tuning
  factor).

- y:

  Numeric scalar or vector used alongside \`x\` to tune the response of
  the HCR to the biomass index. If \`x\` and \`y\` are vectors, they
  must have the same length.

- stab:

  Logical; if \`TRUE\` (default), apply stabilisation bounds to limit
  year-to-year changes in the TAC or fishing mortality.

- lower:

  Numeric value giving the lower stabilisation bound, typically
  interpreted as the minimum allowed ratio \\\text{TAC}\_{y+1} /
  \text{TAC}\_y\\. The default is \`0.8\` (maximum 20% decrease).

- upper:

  Numeric value giving the upper stabilisation bound, typically
  interpreted as the maximum allowed ratio \\\text{TAC}\_{y+1} /
  \text{TAC}\_y\\. The default is \`1.2\` (maximum 20% increase).

- clType:

  Character string specifying the type of control variable returned by
  the HCR. Common options include \`"TAC"\` for catch limits or \`"F"\`
  for fishing mortality targets. The default is \`"TAC"\`.

- clyears:

  Integer specifying the number of recent years used to calculate the
  reference catch level when deriving TAC advice (e.g. mean catch over
  the last \`clyears\` years).

- red:

  Optional numeric reduction factor applied to the catch when biomass is
  low (e.g. \`0.8\` for a 20% reduction). If \`NA\` (default), no
  automatic biomass-based reduction is applied.

- redyears:

  Integer giving the number of years over which the biomass condition
  for reduction is evaluated (e.g. how many years below a threshold
  before the reduction is triggered).

- redAlways:

  Logical; if \`TRUE\`, the reduction rule is evaluated and applied
  whenever the criteria are met. If \`FALSE\`, the reduction may only be
  applied under more restrictive conditions, depending on the internal
  implementation.

- ffmsySD:

  Non-negative numeric specifying the standard deviation of the (log)
  assessment error in \\F/F\_{MSY}\\ associated with this HCR, when
  simulating perceived stock status.

- bbtriggerSD:

  Non-negative numeric specifying the standard deviation of the (log)
  assessment error in biomass-related indicators (e.g.
  \\B/B\_{trigger}\\) associated with this HCR.

- rightRef:

  Integer index specifying which set of reference points should be used
  when evaluating this HCR (e.g. an index into a reference table).

- assessmentInterval:

  Positive integer giving the interval (in years) between full
  assessments and updates of the advice. The default \`1\` corresponds
  to annual updates.

- env:

  Environment in which the generated HCR function(s) are created and
  stored. Defaults to \[globalenv()\].

## Value

A character vector with the name(s) of the HCR function(s) created. The
main purpose of the function is its side effect of defining these HCRs
in \`env\`.

## Details

The generated HCR functions usually take arguments of the form
\`function(x, Data, reps, ...)\`, where \`Data\` contains the catch and
index time series and \`reps\` is the number of stochastic TAC draws (if
used). Internally, the rule uses recent values of a biomass index,
applies stabilisation bounds on year-to-year changes, and can optionally
reduce catch when biomass is low.

Several arguments (e.g. \`id\`, \`x\`, \`y\`) can be provided as vectors
of equal length, in which case multiple HCRs are created in a single
call.

## Examples

``` r
if (FALSE) { # \dontrun{
  ## Define a single index-based HCR with default tuning
  hcr_names <- def.hcr.index(id = "index_HCR")
  hcr_names

  ## Define several index-based HCRs at once with different tuning
  hcr_names2 <- def.hcr.index(
    id = c("index_low", "index_high"),
    x  = c(0.8, 1.0),
    y  = c(1.2, 1.5)
  )
  hcr_names2

  ## The created functions are available in the chosen environment
  ls(pattern = "index_", envir = globalenv())
} # }
```
