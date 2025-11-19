# Define a SPiCT-based harvest control rule

\`def.hcr.spict()\` defines one or more harvest control rules (HCRs)
that use the SPiCT model (stochastic surplus production model in
continuous time) to generate catch advice. For each call, a function is
created that:

## Usage

``` r
def.hcr.spict(
  id = "spict-msy",
  fractiles = list(catch = 0.5, ffmsy = 0.5, bbmsy = 0.5, bmsy = 0.5, fmsy = 0.5),
  breakpointB = 0,
  evalBreakpointB = 0,
  safeguard = list(limitB = 0, prob = 0.95),
  dteuler = 1/4,
  reportmode = 1,
  stabilise = 0,
  priorlogn = c(log(2), 2, 1),
  priorlogsdf = c(5, 2, 0),
  priorlogsdc = c(log(0.2), 2, 0),
  priorlogalpha = c(log(1), 2, 1),
  priorlogbeta = c(log(1), 2, 1),
  priorlogbkfrac = c(log(0.5), 2, 0),
  priorlogr = c(log(0.4), 2, 0),
  fixn = NA,
  bfac = NA,
  bref = "current",
  brefType = "target",
  nyBref = 5,
  btar = "bmsy",
  probtar = 0.4,
  brule = 0,
  red = NA,
  redyears = 3,
  redAlways = FALSE,
  rai = 0.2,
  manstartdY = 0,
  assessmentInterval = 1,
  cpvec = FALSE,
  halfAssessInt = FALSE,
  intC = NA,
  nonconvHCR = "conscat",
  clType = "TAC",
  clyears = 1,
  stab = FALSE,
  lower = 0.8,
  upper = 1.2,
  bm = FALSE,
  env = globalenv(),
  scenario = NA
)
```

## Arguments

- id:

  Character string giving the name/ID of the HCR. Default is
  \`"spict-msy"\`. If vectors are supplied for some arguments, multiple
  HCRs can be created in one call, each with its own \`id\`.

- fractiles:

  List of fractiles used to derive advice from SPiCT output. Typical
  elements include quantiles for predicted catch and reference points,
  e.g. \`list(catch = 0.5, ffmsy = 0.5, bbmsy = 0.5, bmsy = 0.5, fmsy =
  0.5)\`. Values are probabilities in \\0, 1\].

- breakpointB:

  Numeric biomass breakpoint (on the scale of \\B/B\_{MSY}\\ or related
  indicator) that can modify the rule below a given biomass level (e.g.
  more conservative advice when biomass is low).

- evalBreakpointB:

  Numeric or integer specifying how the breakpoint should be evaluated
  (e.g. over which period or using which biomass indicator). The exact
  interpretation depends on the internal implementation.

- safeguard:

  List describing an additional biomass safeguard, typically of the form
  \`list(limitB = \<threshold\>, prob = \<probability\>)\`. For example,
  the rule may reduce or cap TAC if the probability that biomass is
  below \`limitB\` exceeds \`1 - prob\`.

- dteuler:

  Time step of the forward Euler scheme used by SPiCT to solve the
  underlying continuous-time model. Default is \`1/4\` (quarterly).

- reportmode:

  Integer flag passed to SPiCT controlling the level of detail in
  reporting and output.

- stabilise:

  Numeric or logical argument passed to SPiCT's \`stabilise\` option,
  controlling how strongly the SPiCT fit is stabilised (e.g. in the
  presence of limited data). Default is \`0\` (no additional
  stabilisation).

- priorlogn:

  Numeric vector giving the prior for the initial depletion (log scale),
  typically of the form \`c(mean, sd, p)\` where \`p\` is an indicator
  controlling whether the prior is active.

- priorlogsdf:

  Numeric vector giving the prior for the process-noise standard
  deviation on fishing mortality (log scale).

- priorlogsdc:

  Numeric vector giving the prior for the process-noise standard
  deviation on catches (log scale).

- priorlogalpha:

  Numeric vector giving the prior for the production function exponent
  \\\alpha\\ (log scale).

- priorlogbeta:

  Numeric vector giving the prior for the production function exponent
  \\\beta\\ (log scale).

- priorlogbkfrac:

  Numeric vector giving the prior for the fraction of carrying capacity
  at which biomass starts (log scale).

- priorlogr:

  Numeric vector giving the prior for the intrinsic growth rate \\r\\
  (log scale).

- fixn:

  Optional value indicating whether the depletion parameter (\`n\`)
  should be fixed in the SPiCT fit. Default \`NA\` means it is
  estimated.

- bfac:

  Optional numeric factor used to scale biomass-based reference points
  when computing advice. Default \`NA\` means no additional scaling.

- bref:

  Character string specifying how the biomass reference level is
  defined, e.g. \`"current"\`, \`"lowest"\`, \`"highest"\`,
  \`"average"\` or \`"last"\`. This reference is used together with
  \`brefType\`, \`nyBref\`, and \`btar\` to construct biomass-based
  triggers.

- brefType:

  Character string specifying the type of biomass reference, typically
  \`"target"\` or \`"limit"\`. This affects how conservative the rule
  becomes when biomass is near or below the reference.

- nyBref:

  Integer giving the number of years used to define the biomass
  reference (e.g. average biomass over the last \`nyBref\` years).

- btar:

  Character string specifying the target biomass level used in the rule,
  e.g. \`"bmsy"\` for biomass at MSY.

- probtar:

  Numeric probability (in \\0, 1\]) indicating how conservative the
  target rule should be (e.g. requiring that the probability of being
  above \`btar\` exceeds \`probtar\`).

- brule:

  Numeric code or scaling factor controlling how the biomass rule is
  applied below the target (e.g. linear reduction, more aggressive
  reductions, or no reduction).

- red:

  Optional numeric reduction factor applied to the advised catch when
  additional criteria are met (e.g. low biomass, repeated assessment
  failure). For example, \`0.8\` corresponds to a 20% reduction. If
  \`NA\`, no such reduction is applied.

- redyears:

  Integer indicating over how many years the reduction criteria are
  evaluated (e.g. number of years below a threshold before the reduction
  is triggered).

- redAlways:

  Logical; if \`TRUE\`, the reduction rule is evaluated and applied
  every time the criteria are met. If \`FALSE\`, the reduction may be
  applied more sparingly, depending on the internal implementation.

- rai:

  Numeric "raising factor" that controls how quickly the rule responds
  when moving back towards higher catches (e.g. after a reduction).
  Default is \`0.2\`.

- manstartdY:

  Numeric management start year offset relative to the assessment time
  series (e.g. \`0\` for the last year, negative values for earlier
  years). Default is \`0\`.

- assessmentInterval:

  Positive integer giving the interval (in years) between full SPiCT
  assessments and advice updates. Default is \`1\`, i.e. annual
  assessments.

- cpvec:

  Logical; if \`TRUE\`, use a vector of catch profiles or control points
  when computing advice. Default is \`FALSE\`.

- halfAssessInt:

  Logical; if \`TRUE\`, the assessment is effectively shifted by half an
  assessment interval (e.g. mid-year updates). Default is \`FALSE\`.

- intC:

  Numeric or \`NA\`. If non-\`NA\`, this can be used to specify a
  predefined catch level in an initial period or over an interval.
  Default is \`NA\`.

- nonconvHCR:

  Character string giving the name of an alternative HCR to be used if
  the SPiCT assessment does not converge. Default is \`"conscat"\`,
  which typically corresponds to a constant-catch rule defined by
  \[def.hcr.conscat()\].

- clType:

  Character string specifying the type of control variable returned by
  the HCR, e.g. \`"TAC"\` for catch advice. Used in combination with the
  stabilisation cap.

- clyears:

  Integer giving the number of years used to compute the reference level
  for the stabilisation cap (e.g. mean catch over the last \`clyears\`
  years). Default is \`1\`.

- stab:

  Logical; if \`TRUE\`, apply an uncertainty or stabilisation cap that
  limits year-to-year changes in the advice (e.g. TAC). Default is
  \`FALSE\`.

- lower:

  Numeric lower bound for the stabilisation cap, typically interpreted
  as the minimum allowed ratio \\\text{Advice}\_{y+1} /
  \text{Advice}\_y\\. Default is \`0.8\` (maximum 20% decrease).

- upper:

  Numeric upper bound for the stabilisation cap, typically interpreted
  as the maximum allowed ratio \\\text{Advice}\_{y+1} /
  \text{Advice}\_y\\. Default is \`1.2\` (maximum 20% increase).

- bm:

  Logical or numeric. If \`FALSE\` (default), no explicit benchmarks are
  applied. If a positive integer is supplied (e.g. \`5\`), the biomass
  reference \`Bref\` may be redefined every \`bm\` years (benchmarking).

- env:

  Environment in which the generated HCR function(s) are created and
  stored. Defaults to \[globalenv()\].

- scenario:

  Optional label (character or \`NA\`) used to identify the scenario
  associated with this HCR in subsequent summaries and plots.

## Value

A character vector containing the name(s) of the HCR function(s)
created. The main purpose of the function is its side effect of defining
SPiCT-based HCRs in \`env\`.

## Details

1\. Fits or updates a SPiCT assessment using catch and index data, 2.
Extracts stock status and predicted catch distributions, and 3.
Translates these into TAC recommendations according to specified
fractiles, biomass breakpoints, safeguards, and stabilisation rules.

The generated HCR function(s) are assigned to the environment \`env\`
(typically the global environment) and can then be used within a
management strategy evaluation (MSE) framework.

This function assumes that the \*\*spict\*\* package is available and
that the underlying catch and index data can be passed to SPiCT in the
expected format. Many arguments (e.g. priors and Euler time step) are
passed directly to SPiCT's input list and control how the assessment
behaves. For details on SPiCT itself, see the \*\*spict\*\*
documentation.
