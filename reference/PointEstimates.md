# Weighted Point Estimates

Predefined functions for weighted point estimates in package `surveysd`.

## Usage

``` r
weightedRatio(x, w)

weightedSum(x, w)
```

## Arguments

- x:

  numeric vector

- w:

  weight vector

## Value

Each of the functions return a single numeric value

## Details

Predefined functions are weighted ratio and weighted sum.

## Examples

``` r
x <- 1:10
w <- 10:1
weightedRatio(x,w)
#> [1] 18.18182
x <- 1:10
w <- 10:1
weightedSum(x,w)
#> [1] 220
```
