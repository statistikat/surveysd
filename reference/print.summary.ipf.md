# Print method for IPF calibration summary

Provides a concise summary of an IPF (Iterative Proportional Fitting)
calibration summary object. It extracts the calibration weight from
`all_formulas`, computes the Kish factor for the weights, and prints the
first 10 rows of any `calib_results_` tables. Useful for a quick
overview of calibration results. Additional details can be explored with
[`str()`](https://rdrr.io/r/utils/str.html) or
[`names()`](https://rdrr.io/r/base/names.html).

## Usage

``` r
# S3 method for class 'summary.ipf'
print(x, ...)
```

## Arguments

- x:

  An object of class `summary.ipf`, as returned by
  [`summary.ipf`](https://statistikat.github.io/surveysd/reference/summary.ipf.md).

- ...:

  Additional arguments (currently ignored).

## Value

The input object `x`, invisibly (for chaining).
