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
#> 1  1.4389507         a
#> 2  0.7604319         b
#> 3  1.2640370         b
#> 4  0.5975856         a
#> 5  0.3327697         c
#> 6  0.3952720         b
#> 7  1.6274510         a
#> 8  1.8242116         c
#> 9  0.8413260         c
#> 10 0.8239812         c

## create targets (same lenght as classLabels!)
targets <- 3:5

## calculate weights
new_weight <- ipf_step(dat$weight, dat$household, targets)
cbind(dat, new_weight)
#>       weight household new_weight
#> 1  1.4389507         a  1.1781843
#> 2  0.7604319         b  1.2570469
#> 3  1.2640370         b  2.0895410
#> 4  0.5975856         a  0.4892912
#> 5  0.3327697         c  0.4353016
#> 6  0.3952720         b  0.6534121
#> 7  1.6274510         a  1.3325245
#> 8  1.8242116         c  2.3862820
#> 9  0.8413260         c  1.1005527
#> 10 0.8239812         c  1.0778637

## check solution
xtabs(new_weight ~ dat$household)
#> dat$household
#> a b c 
#> 3 4 5 

## calculate weights "by reference"
ipf_step_ref(dat$weight, dat$household, targets)
dat
#>       weight household
#> 1  1.1781843         a
#> 2  1.2570469         b
#> 3  2.0895410         b
#> 4  0.4892912         a
#> 5  0.4353016         c
#> 6  0.6534121         b
#> 7  1.3325245         a
#> 8  2.3862820         c
#> 9  1.1005527         c
#> 10 1.0778637         c

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
#> 1   Dinner Female     No Sun  9
#> 2    Lunch   Male    Yes Mon  8
#> 2.1  Lunch   Male    Yes Mon  8
#> 1.1 Dinner Female     No Sun  9
#> 1.2 Dinner Female     No Sun  9
#> 1.3 Dinner Female     No Sun  9
#> 1.4 Dinner Female     No Sun  9
#> 1.5 Dinner Female     No Sun  9
#> 2.2  Lunch   Male    Yes Mon  8
#> 2.3  Lunch   Male    Yes Mon  8

## adjust weights
weight <- rnorm(nrow(tips)) + 5
adjusted_weight <- ipf_step(weight, cf, con)

## check outputs
con2 <- xtabs(adjusted_weight ~ ., data = tips)
sum((con - con2)^2)
#> [1] 0
```
