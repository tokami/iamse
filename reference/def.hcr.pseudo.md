# Define harvest control rule with pseudo assessment

def.hcr.pseudo

## Usage

``` r
def.hcr.pseudo(
  id = "pseudo-msy",
  fractiles = list(catch = 0.5, ffmsy = 0.5, bbmsy = 0.5, bmsy = 0.5, fmsy = 0.5),
  breakpointB = 0,
  clType = "TAC",
  clyears = 1,
  stab = FALSE,
  lower = 0.8,
  upper = 1.2,
  env = globalenv()
)
```

## Arguments

- id:

  Name/ID of HCR. Default: "pseudo-msy"

- fractiles:

  Fractiles. List

- breakpointB:

  breakpointb

- clType:

  Catch type for uncertainty cap. Default: "TAC".

- clyears:

  Number of years for uncertainty cap. Default: 1.

- stab:

  Uncertainty cap. Default: FALSE.

- lower:

  Upper bound of uncertainty cap. Default: 0.8.

- upper:

  Upper bound of uncertainty cap. Default: 1.2.

- env:

  Environment. Default: globalenv()
