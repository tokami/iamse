# Define a constant-catch harvest control rule

\`def.hcr.conscat()\` defines a simple harvest control rule (HCR) that
recommends a constant catch (TAC) over time. The baseline TAC can be
supplied directly via \`constantC\` or estimated from recent catches,
and optional reduction rules can be used to step down the TAC when
performance is poor (e.g. low biomass or repeated assessment failure).

## Usage

``` r
def.hcr.conscat(
  id = "conscat",
  constantC = NULL,
  clyears = 1,
  red = NA,
  redyears = 2,
  redAlways = FALSE,
  assessmentInterval = 1,
  ffmsySD = 0,
  bbtriggerSD = 0,
  rightRef = 1,
  set. = check.set(),
  env = globalenv()
)
```

## Arguments

- id:

  Character string giving the identifier for this HCR (for example
  \`"conscat"\`). This label is used when selecting HCRs in \`set\$hcr\`
  and in summaries and plots.

- constantC:

  Optional numeric giving the constant catch level (e.g. in tonnes). If
  \`NULL\` (default), the constant catch is derived internally from
  recent catches, typically using the last \`clyears\` years.

- clyears:

  Integer specifying the number of historical years used to calculate
  the baseline catch level when \`constantC\` is \`NULL\`.

- red:

  Optional numeric reduction factor applied to the constant catch (e.g.
  \`0.8\` corresponds to a 20 automatic reduction is applied and the TAC
  remains constant.

- redyears:

  Integer giving the number of years over which the reduction criteria
  are evaluated. The exact interpretation depends on the internal
  implementation (e.g. how many years of poor performance trigger a
  reduction).

- redAlways:

  Logical; if \`TRUE\`, the reduction rule is evaluated in all years
  where the criteria are met. If \`FALSE\`, the reduction may be applied
  only once or under more restrictive conditions, depending on the
  internal implementation.

- assessmentInterval:

  Positive integer giving the interval (in years) between assessments
  and TAC updates. The default \`1\` corresponds to an annual advice
  cycle.

- ffmsySD:

  Non-negative numeric giving the standard deviation of the (log)
  assessment error in \\F/F\_{MSY}\\ used when simulating the assessed
  status for this HCR.

- bbtriggerSD:

  Non-negative numeric giving the standard deviation of the (log)
  assessment error in biomass-related indicators (e.g.
  \\B/B\_{trigger}\\) used for this HCR.

- rightRef:

  Integer index specifying which reference-point set should be used when
  evaluating this HCR (e.g. an index into \`dat\$ref\`).

- set.:

  An \`iamse\` settings list as returned by \[check.set()\]. This object
  is updated to include the newly defined HCR.

- env:

  Environment in which the HCR function is created and stored. Defaults
  to \[globalenv()\]. In most cases the default should be used.

## Value

Invisibly returns the updated \`set.\` object. The function is called
for its side effects (defining and registering the HCR).

## Details

The function is usually called for its side effects: it registers a new
HCR with identifier \`id\` in the \`iamse\` settings object \`set.\`
and, by default, also in the supplied environment \`env\`. The HCR can
then be included in \`set\$hcr\` and evaluated with \[run.mse()\].
