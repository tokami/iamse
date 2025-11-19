# Intra-annual Management Strategy Evaluation (iamse)

`iamse` is an R package for **intra-annual, age-based management strategy evaluation (MSE)**.
It is designed to explore how different harvest control rules (HCRs) perform when key processes
(e.g. growth, recruitment, fishing mortality) vary within the year and through time.

---

## Highlights

1. **Sub-annual time steps**
   - Explicit intra-annual dynamics (e.g. seasons or quarters) in the operating model.

2. **Optional length-dependent processes**
   - Length-based growth, selectivity, and mortality can be included where data allow.

3. **Stochastic age–length key**
   - Flexible, stochastic age–length conversion to link length-based observations to age-based dynamics.

4. **Stochastic reference points**
   - Simulation-based estimation of reference points (e.g. MSY, F_MSY, B_MSY, B_0) using
     `est.ref.levels.stochastic()`.

5. **Dynamic (time-varying) processes and reference points**
   - Time-varying life-history parameters and reference points to reflect changing productivity.

---

## Getting started

Once the package is installed and loaded, a typical workflow is:

1. Prepare life-history and stock data with `check.dat()`.
2. Set up simulation and noise settings with `check.set()`.
3. Estimate reference points (deterministic or stochastic).
4. Define one or more HCRs (e.g. `def.hcr.conscat()`, `def.hcr.ref()`, `def.hcr.spict()`).
5. Run the MSE with `run.mse()`.
6. Summarise and visualise results using the `plotiamse.*()` and `est.metrics()` functions.

A basic example (mirroring the demo vignette):

```r
library(iamse)

## Example life-history data
data("stocklist")
dat0 <- stocklist$mlh
dat0$CVlen <- 0.1
dat  <- check.dat(dat0)

## Simulation settings
set <- check.set()
set$nysim   <- 20   # projection years
set$nrep    <- 10   # number of replicates (low for demo)
set$refMethod <- "median"

## Noise structure (simplified)
set$noise$rep$H  <- c(0.1, 0, 1)
set$noise$rep$M  <- c(0.1, 0, 1)
set$noise$rep$R0 <- c(0.1, 0, 1)
set$noise$time$R <- c(0.7, 0.4, 1)
set$noise$rep$F  <- c(0.15, 0, 1)
set$noise$time$I <- c(0.3, 0, 1)
set$noise$time$C <- c(0.1, 0, 1)

## Stochastic reference points
set$refN <- 100
dat <- est.ref.levels.stochastic(dat, set, fmax = 5,
                                 ncores = 2, plot = TRUE)

dat$ref$Blim <- 0.2 * dat$ref$B0

## Define some HCRs
iamse::def.hcr.ref(set. = set)                      # "noF"
iamse::def.hcr.ref(consF = "fmsy", set. = set)      # "refFmsy"
iamse::def.hcr.conscat(set. = set)                  # "conscat"

set$hcr <- c("noF", "refFmsy", "conscat")

## Run MSE (single core for portability)
resMSE <- run.mse(dat, set, ncores = 1, verbose = TRUE)

## Plot time series
par(mfrow = c(3, 1), mar = c(4, 4, 1, 1))
plotiamse.f(dat, set, resMSE)
plotiamse.cw(dat, set, resMSE)
plotiamse.b(dat, set, resMSE)
```

---

## Documentation

A **vignette** will be available to provide a more complete introduction to:

- The structure of the operating model,
- How to define and combine HCRs,
- Options for stochastic reference points and productivity analysis,
- Interpretation of key performance metrics and trade-offs.

Until then, the main help pages are a good starting point:

- `?check.dat` – prepare and validate life-history / stock input.
- `?check.set` – define simulation and noise settings.
- `?run.mse` – run the MSE given data, settings, and HCRs.

Function-specific documentation is available via the usual R help system, e.g.:

```r
?check.dat
?run.mse
?def.hcr.conscat
?plotiamse.tradeoff
```

---

## Package requirements

The package depends on several other R packages. Key ones include:

- [`Rcpp`](https://www.rcpp.org) (and `LinkingTo: Rcpp`),
- `doBy`,
- `future`, `future.apply`,
- `parallel`,
- (optionally) [`spict`](https://github.com/DTUAqua/spict) for SPiCT-based HCRs.

For the full list of dependencies, see the `DESCRIPTION` file of the package.

---

## Installing the `iamse` package

You can install the development version of `iamse` from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("tokami/iamse")
```

Then load it with:

```r
library(iamse)
```

---

## Citation

Please use:

```r
citation("iamse")
```

to get the recommended citations for the package and associated methodology.

---

## Funding

The development of **iamse** was co-funded by the European Union.

---

<figure>
  <img src="man/figures/EN_Co-fundedbytheEU_RGB_POS.png" width="250"
       alt="Co-funded by the EU logo" />
  <figcaption aria-hidden="true">Co-funded by the EU logo</figcaption>
</figure>
