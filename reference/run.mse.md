# run.mse

Run MSE

## Usage

``` r
run.mse(
  dat,
  set,
  fest = NULL,
  ncores = parallel::detectCores() - 1,
  globals = NULL,
  verbose = TRUE
)
```

## Arguments

- dat:

  data

- set:

  settings

- fest:

  F

- ncores:

  for parallel

- globals:

  list with extra objects to be passed on to cluseters (only if ncores
  \> 1)

- verbose:

  logical; print info

## Value

a list with MSE results
