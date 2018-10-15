
surveysd <img src="man/figures/logo.png" align="right" alt="" />
================================================================

[![Build Status](https://travis-ci.org/statistikat/surveysd.svg?branch=master)](https://travis-ci.org/statistikat/surveysd) [![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![GitHub last commit](https://img.shields.io/github/last-commit/statistikat/surveysd.svg)](https://github.com/statistikat/surveysd/commits/master) [![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/statistikat/surveysd.svg)](https://github.com/statistikat/surveysd)

<!--[![Coverage Status](https://coveralls.io/repos/github/statistikat/surveysd/badge.svg?branch=master)](https://coveralls.io/github/statistikat/surveysd?branch=master)-->
<!--[![CRAN](http://www.r-pkg.org/badges/version/surveysd)](https://CRAN.R-project.org/package=surveysd)-->
<!--[![Downloads](http://cranlogs.r-pkg.org/badges/surveysd)](https://CRAN.R-project.org/package=surveysd)-->
<!--[![Mentioned in Awesome Official Statistics ](https://awesome.re/mentioned-badge.svg)](http://www.awesomeofficialstatistics.org)-->
This is the development place for the R-package `surveysd`. This package can be used to estimate the standard deviation of estimates in complex surveys using bootstrap weights.

Installation
------------

This package can be installed like any other R package on github via `install_github`

``` r
devtools::install_github("statistikat/surveysd")
```

Concept
-------

Bootstrapping has long been around and used widely to estimate confidence intervals and standard errors of point estimates. This package aims to combine all necessary steps for applying a calibrated bootstrapping procedure with custom estimating functions.

Workflow
--------

-   Bootstrap samples are drawn with rescaled bootstrapping as described in `bla` in the function `draw.bootstrap()`.
-   These samples can then be calibrated with an iterative proportional updating algorithm using `recalib()`.
-   Finally, estimation functions can be applied over all bootstrap replicates with `calc.stError()`.

Further reading
---------------

-   Simple syntax examples can be found in the [getting started vignette](https://statistikat.github.io/surveysd/articles/surveysd.html).
-   The methodology is covered in the [methodology vignette](https://statistikat.github.io/surveysd/articles/Methodology.html).
