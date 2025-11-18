# surveysd

[![Travis build
status](https://img.shields.io/travis/statistikat/surveysd.svg?logo=travis)](https://app.travis-ci.com/statistikat/surveysd)
[![Coverage
Status](https://img.shields.io/coveralls/github/statistikat/surveysd.svg?colorB=red&logo=codecov)](https://coveralls.io/github/statistikat/surveysd?branch=master)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-maturing-blue.svg?logo=github)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![GitHub last
commit](https://img.shields.io/github/last-commit/statistikat/surveysd.svg?logo=github)](https://github.com/statistikat/surveysd/commits/master)
[![GitHub code size in
bytes](https://img.shields.io/github/languages/code-size/statistikat/surveysd.svg?logo=github)](https://github.com/statistikat/surveysd)
[![Downloads](http://cranlogs.r-pkg.org/badges/surveysd)](https://CRAN.R-project.org/package=surveysd)
[![CRAN](https://img.shields.io/cran/v/surveysd.svg?colorB=green&logo=R&logoColor=blue&label=CRAN)](https://CRAN.R-project.org/package=surveysd)
[![R-CMD-check](https://github.com/statistikat/surveysd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/statistikat/surveysd/actions/workflows/R-CMD-check.yaml)

[![Mentioned in Awesome Official
Statistics](https://awesome.re/mentioned-badge.svg)](https://github.com/SNStatComp/awesome-official-statistics-software)

This is the development place for the R-package `surveysd`. The package
can be used to estimate the standard deviation of estimates in complex
surveys using bootstrap weights.

## Installation

``` r
# Install release version from CRAN
install.packages("surveysd")

# Install development version from GitHub
devtools::install_github("statistikat/surveysd")
```

## Concept

Bootstrapping has long been around and used widely to estimate
confidence intervals and standard errors of point estimates. This
package aims to combine all necessary steps for applying a calibrated
bootstrapping procedure with custom estimating functions.

## Workflow

A typical workflow with this package consists of three steps. To see
these concepts in practice, please refer to the [getting started
vignette](https://statistikat.github.io/surveysd/articles/surveysd.html).

- Calibrated weights can be generated with the function
  [`ipf()`](https://statistikat.github.io/surveysd/reference/ipf.md)
  using an iterative proportional updating algorithm.
- Bootstrap samples are drawn with rescaled bootstrapping in the
  function
  [`draw.bootstrap()`](https://statistikat.github.io/surveysd/reference/draw.bootstrap.md).
- These samples can then be calibrated with an iterative proportional
  updating algorithm using
  [`recalib()`](https://statistikat.github.io/surveysd/reference/recalib.md).
- Finally, estimation functions can be applied over all bootstrap
  replicates with
  [`calc.stError()`](https://statistikat.github.io/surveysd/reference/calc.stError.md).

## Further reading

More information can be found on the [github-pages
site](https://statistikat.github.io/surveysd/) for surveysd.

- The methodology is covered in the [methodology
  vignette](https://statistikat.github.io/surveysd/articles/methodology.html).
- A more comprehensive documentation of
  [`calc.stError()`](https://statistikat.github.io/surveysd/reference/calc.stError.md)
  is available in the [error estimation
  vignette](https://statistikat.github.io/surveysd/articles/error_estimation.html).
