# Perform one step of iterative proportional updating

C++ routines to invoke a single iteration of the Iterative proportional
updating (IPU) scheme. Targets and classes are assumed to be one
dimensional in the `ipf_step` functions. `combine_factors` aggregates
several vectors of type factor into a single one to allow
multidimensional ipu-steps. See examples.

## Usage

``` r
ipf_step_ref(w, classes, targets)

ipf_step(w, classes, targets)

ipf_step_f(w, classes, targets)

combine_factors(dat, targets)
```

## Arguments

- w:

  a numeric vector of weights. All entries should be positive.

- classes:

  a factor variable. Must have the same length as `w`.

- targets:

  key figure to target with the ipu scheme. A numeric verctor of the
  same length as `levels(classes)`. This can also be a `table` produced
  by `xtabs`. See examples.

- dat:

  a `data.frame` containing the factor variables to be combined.

## Details

`ipf_step` returns the adjusted weights. `ipf_step_ref` does the same,
but updates `w` by reference rather than returning. `ipf_step_f` returns
a multiplicator: adjusted weights divided by unadjusted weights.
`combine_factors` is designed to make `ipf_step` work with contingency
tables produced by [xtabs](https://rdrr.io/r/stats/xtabs.html).

## Examples

``` r

############# one-dimensional ipu ##############

## create random data
nobs <- 10
classLabels <- letters[1:3]
dat = data.frame(
  weight = exp(rnorm(nobs)),
  household = factor(sample(classLabels, nobs, replace = TRUE))
)
dat
#>       weight household
#> 1  3.5593264         a
#> 2  0.3192965         b
#> 3  4.9140467         a
#> 4  0.9747069         c
#> 5  2.4586285         b
#> 6  1.2461800         b
#> 7  2.9306996         b
#> 8  4.7771154         c
#> 9  0.8780809         c
#> 10 7.0632484         b

## create targets (same lenght as classLabels!)
targets <- 3:5

## calculate weights
new_weight <- ipf_step(dat$weight, dat$household, targets)
cbind(dat, new_weight)
#>       weight household new_weight
#> 1  3.5593264         a 1.26018046
#> 2  0.3192965         b 0.09111008
#> 3  4.9140467         a 1.73981954
#> 4  0.9747069         c 0.73508382
#> 5  2.4586285         b 0.70156063
#> 6  1.2461800         b 0.35559289
#> 7  2.9306996         b 0.83626437
#> 8  4.7771154         c 3.60270374
#> 9  0.8780809         c 0.66221244
#> 10 7.0632484         b 2.01547203

## check solution
xtabs(new_weight ~ dat$household)
#> dat$household
#> a b c 
#> 3 4 5 

## calculate weights "by reference"
ipf_step_ref(dat$weight, dat$household, targets)
dat
#>        weight household
#> 1  1.26018046         a
#> 2  0.09111008         b
#> 3  1.73981954         a
#> 4  0.73508382         c
#> 5  0.70156063         b
#> 6  0.35559289         b
#> 7  0.83626437         b
#> 8  3.60270374         c
#> 9  0.66221244         c
#> 10 2.01547203         b

############# multidimensional ipu ##############

## load data
factors <- c("time", "sex", "smoker", "day")
tips <- data.frame(sex=c("Female","Male","Male"), day=c("Sun","Mon","Tue"),
time=c("Dinner","Lunch","Lunch"), smoker=c("No","Yes","No"))
tips <- tips[factors]

## combine factors
con <- xtabs(~., tips)
cf <- combine_factors(tips, con)
cbind(tips, cf)[sample(nrow(tips), 10, replace = TRUE),]
#>       time    sex smoker day cf
#> 3    Lunch   Male     No Tue 20
#> 2    Lunch   Male    Yes Mon  8
#> 3.1  Lunch   Male     No Tue 20
#> 2.1  Lunch   Male    Yes Mon  8
#> 2.2  Lunch   Male    Yes Mon  8
#> 2.3  Lunch   Male    Yes Mon  8
#> 3.2  Lunch   Male     No Tue 20
#> 2.4  Lunch   Male    Yes Mon  8
#> 1   Dinner Female     No Sun  9
#> 2.5  Lunch   Male    Yes Mon  8

## adjust weights
weight <- rnorm(nrow(tips)) + 5
adjusted_weight <- ipf_step(weight, cf, con)

## check outputs
con2 <- xtabs(adjusted_weight ~ ., data = tips)
sum((con - con2)^2)
#> [1] 0
```
