# Calculate mean by factors

These functions calculate the arithmetic and geometric mean of the
weight for each class. `geometric_mean` and `arithmetic_mean` return a
`numeric` vector of the same length as `w` which stores the averaged
weight for each observation. `geometric_mean_reference` returns the same
value by reference, i.e. the input value `w` gets overwritten by the
updated weights. See examples.

## Usage

``` r
geometric_mean_reference(w, classes)
```

## Arguments

- w:

  An numeric vector. All entries should be positive.

- classes:

  A factor variable. Must have the same length as `w`.

## Examples

``` r
if (FALSE) { # \dontrun{

## create random data
nobs <- 10
classLabels <- letters[1:3]
dat = data.frame(
  weight = exp(rnorm(nobs)),
  household = factor(sample(classLabels, nobs, replace = TRUE))
)
dat

## calculate weights with geometric_mean
geom_weight <- geometric_mean(dat$weight, dat$household)
cbind(dat, geom_weight)

## calculate weights with arithmetic_mean
arith_weight <- arithmetic_mean(dat$weight, dat$household)
cbind(dat, arith_weight)

## calculate weights "by reference"
geometric_mean_reference(dat$weight, dat$household)
dat
} # }
```
