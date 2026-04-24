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
#> 1  0.5231661         c
#> 2  1.6672206         b
#> 3  2.0636623         c
#> 4  0.7314631         a
#> 5  2.1868641         c
#> 6  2.8496370         a
#> 7  0.7967923         b
#> 8  0.9417180         b
#> 9  3.0598142         b
#> 10 0.9139348         b

## create targets (same lenght as classLabels!)
targets <- 3:5

## calculate weights
new_weight <- ipf_step(dat$weight, dat$household, targets)
cbind(dat, new_weight)
#>       weight household new_weight
#> 1  0.5231661         c  0.5479680
#> 2  1.6672206         b  0.9037063
#> 3  2.0636623         c  2.1614948
#> 4  0.7314631         a  0.6127696
#> 5  2.1868641         c  2.2905372
#> 6  2.8496370         a  2.3872304
#> 7  0.7967923         b  0.4318962
#> 8  0.9417180         b  0.5104522
#> 9  3.0598142         b  1.6585528
#> 10 0.9139348         b  0.4953925

## check solution
xtabs(new_weight ~ dat$household)
#> dat$household
#> a b c 
#> 3 4 5 

## calculate weights "by reference"
ipf_step_ref(dat$weight, dat$household, targets)
dat
#>       weight household
#> 1  0.5479680         c
#> 2  0.9037063         b
#> 3  2.1614948         c
#> 4  0.6127696         a
#> 5  2.2905372         c
#> 6  2.3872304         a
#> 7  0.4318962         b
#> 8  0.5104522         b
#> 9  1.6585528         b
#> 10 0.4953925         b

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
#> 2    Lunch   Male    Yes Mon  8
#> 3    Lunch   Male     No Tue 20
#> 3.1  Lunch   Male     No Tue 20
#> 1   Dinner Female     No Sun  9
#> 2.1  Lunch   Male    Yes Mon  8
#> 2.2  Lunch   Male    Yes Mon  8
#> 3.2  Lunch   Male     No Tue 20
#> 2.3  Lunch   Male    Yes Mon  8
#> 1.1 Dinner Female     No Sun  9
#> 1.2 Dinner Female     No Sun  9

## adjust weights
weight <- rnorm(nrow(tips)) + 5
adjusted_weight <- ipf_step(weight, cf, con)

## check outputs
con2 <- xtabs(adjusted_weight ~ ., data = tips)
sum((con - con2)^2)
#> [1] 0
```
