# Summarise an IAMSE object

This is the \`summary()\` method for objects of class \`"iamse"\`. It
provides a compact overview of the main components of an IAMSE object,
such as the number of stocks, harvest control rules (HCRs), simulation
years and replicates, and possibly key performance metrics or reference
points, depending on what is stored in \`object\`.

## Usage

``` r
# S3 method for class 'iamse'
summary(object, ...)
```

## Arguments

- object:

  An object of class \`"iamse"\`, usually returned by IAMSE helper
  functions (e.g. a wrapper around data, settings, and/or results).

- ...:

  Further arguments passed to or from other methods. Currently not used,
  but included for method compatibility.

## Value

An object of class \`"summary.iamse"\` (or similar), which is typically
printed to the console. The function is usually called for its side
effects (printing a readable summary), and the return value is returned
invisibly.

## Details

Typically, users will not call \`summary.iamse()\` directly, but will
instead use \[summary()\] on an \`"iamse"\` object:

      summary(x)
