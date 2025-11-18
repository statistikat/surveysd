# Package index

## data import

- [`demo.eusilc()`](https://statistikat.github.io/surveysd/reference/demo.eusilc.md)
  : Generate multiple years of EU-SILC data

## bootstrapping and calibration

Functions to draw bootstrap samples and calibrate each sample accoring
to population totals. The methods used in these function relate to a
rescaled bootstrap as described in
[`vignette("methodology")`](https://statistikat.github.io/surveysd/articles/methodology.md).

- [`ipf()`](https://statistikat.github.io/surveysd/reference/ipf.md) :
  Iterative Proportional Fitting
- [`draw.bootstrap()`](https://statistikat.github.io/surveysd/reference/draw.bootstrap.md)
  : Draw bootstrap replicates
- [`rescaled.bootstrap()`](https://statistikat.github.io/surveysd/reference/rescaled.bootstrap.md)
  : Draw bootstrap replicates
- [`generate.HHID()`](https://statistikat.github.io/surveysd/reference/generate.HHID.md)
  : Generate new houshold ID for survey data with rotating panel design
  taking into account split households
- [`get.selection()`](https://statistikat.github.io/surveysd/reference/get.selection.md)
  : Get sample selection (~deltas) from drawn bootstrap replicates
- [`recalib()`](https://statistikat.github.io/surveysd/reference/recalib.md)
  : Calibrate weights

## estimation of standard errors

Apply estimators to each sample to generate standard errors as well as
confidence intervals. See
[`vignette("error_estimation")`](https://statistikat.github.io/surveysd/articles/error_estimation.md)
for more details.

- [`weightedRatio()`](https://statistikat.github.io/surveysd/reference/PointEstimates.md)
  [`weightedSum()`](https://statistikat.github.io/surveysd/reference/PointEstimates.md)
  : Weighted Point Estimates
- [`calc.stError()`](https://statistikat.github.io/surveysd/reference/calc.stError.md)
  : Calcualte point estimates and their standard errors using bootstrap
  weights.
- [`plot(`*`<surveysd>`*`)`](https://statistikat.github.io/surveysd/reference/plot.surveysd.md)
  : Plot surveysd-Objects
- [`print(`*`<surveysd>`*`)`](https://statistikat.github.io/surveysd/reference/print.surveysd.md)
  : Print function for surveysd objects

## utility functions

Misc helper functions that are used in or related to
[`ipf()`](https://statistikat.github.io/surveysd/reference/ipf.md).

- [`ipf_step_ref()`](https://statistikat.github.io/surveysd/reference/ipf_step.md)
  [`ipf_step()`](https://statistikat.github.io/surveysd/reference/ipf_step.md)
  [`ipf_step_f()`](https://statistikat.github.io/surveysd/reference/ipf_step.md)
  [`combine_factors()`](https://statistikat.github.io/surveysd/reference/ipf_step.md)
  : Perform one step of iterative proportional updating
- [`kishFactor()`](https://statistikat.github.io/surveysd/reference/kishFactor.md)
  : Kish Factor
- [`geometric_mean_reference()`](https://statistikat.github.io/surveysd/reference/cpp_mean.md)
  : Calculate mean by factors
- [`computeLinear()`](https://statistikat.github.io/surveysd/reference/computeFrac.md)
  [`computeLinearG1_old()`](https://statistikat.github.io/surveysd/reference/computeFrac.md)
  [`computeLinearG1()`](https://statistikat.github.io/surveysd/reference/computeFrac.md)
  [`computeFrac()`](https://statistikat.github.io/surveysd/reference/computeFrac.md)
  : Numerical weighting functions
- [`summary(`*`<ipf>`*`)`](https://statistikat.github.io/surveysd/reference/summary.ipf.md)
  : Generate Summary Output for IPF Calibration
- [`print(`*`<summary.ipf>`*`)`](https://statistikat.github.io/surveysd/reference/print.summary.ipf.md)
  : Print method for IPF calibration summary
