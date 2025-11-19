# Estimate stochastic productivity relationship

Estimate stochastic productivity relationship

## Usage

``` r
est.productivity.stochastic(
  dat,
  set = NULL,
  fmax = 10,
  nf = 1000,
  prob = c(0.1, 0.9),
  ncores = parallel::detectCores() - 1,
  plot = TRUE
)
```

## Arguments

- dat:

  info

- set:

  info

- fmax:

  info

- nf:

  info

- prob:

  info

- ncores:

  info

- plot:

  info

## Value

An object (typically a list or data frame) containing stochastic
productivity summaries as a function of fishing mortality (and possibly
biomass or depletion) for each stock, including the quantiles specified
in `prob`. The exact structure depends on the internal implementation.
The object is usually returned invisibly if the function is called
primarily for its side effects (for example plotting).
