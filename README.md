Management strategy evaluation (MSE)
=====

An R-package for mse

Highlights:

1. Sub-annual time steps.
2. Stochastic estimation of reference points (MSY, Fmsy, Bmsy).
3. Stochastic age-length key.
4. Length-dependent processes.


## Help files

A vignette for the package is available [`here`](https://github.com/mawp/spict/raw/master/spict/vignettes/vignette.pdf), and serves as an introduction to the package and its functionality. The vignette also contains description of the more advanced features of the package.

The package also contains reasonable documentation in the form of help texts associated with each function (some may not be fully up-to-date). These can be accessed in the usual R manner by writing e.g. ```?checkDat```. A good place to start (in addition to reading the vignette) is to read ```?spict::check.inp``` and ```?runMSE```.

## Citation

```citation("mse")```

## Package requirements

The package requires [`Rcpp`](http://www.tmb-project.org) to be installed. TMB is now a part of CRAN and can therefore be installed using the install.packages() command. For more information about TMB click [`here`](https://github.com/kaskr/adcomp).

## Installing the spict package

To install spict from GitHub use

```
library(devtools)
install_github("tokami/mse")            # master branch
```
