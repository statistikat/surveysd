# Kish Factor

Compute the design effect due to unequal weighting.

## Usage

``` r
kishFactor(w, na.rm = FALSE)
```

## Arguments

- w:

  a numeric vector with weights

- na.rm:

  a logical value indicating whether NA values should be stripped before
  the computation proceeds.

## Value

The function will return the the kish factor

## Details

The factor is computed acording to 'Weighting for Unequal P_i', Leslie
Kish, Journal of Official Statistics, Vol. 8. No. 2, 1992 \$\$ deff =
\sqrt n \sum_j w_j^2 / (\sum_j w_j)^2\$\$

## Author

Alexander Kowarik

## Examples

``` r
kishFactor(rep(1,10))
#> [1] 1
kishFactor(rlnorm(10))
#> [1] 1.1811
```
