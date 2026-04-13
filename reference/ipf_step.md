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
#> 1  1.1850357         b
#> 2  1.2409946         c
#> 3  0.7089801         c
#> 4  2.1531990         c
#> 5  1.1281514         a
#> 6  0.5278370         c
#> 7  0.8531578         a
#> 8  1.4710711         c
#> 9  1.0081016         b
#> 10 2.9852986         b

## create targets (same lenght as classLabels!)
targets <- 3:5

## calculate weights
new_weight <- ipf_step(dat$weight, dat$household, targets)
cbind(dat, new_weight)
#>       weight household new_weight
#> 1  1.1850357         b  0.9153619
#> 2  1.2409946         c  1.0168616
#> 3  0.7089801         c  0.5809330
#> 4  2.1531990         c  1.7643151
#> 5  1.1281514         a  1.7081908
#> 6  0.5278370         c  0.4325056
#> 7  0.8531578         a  1.2918092
#> 8  1.4710711         c  1.2053846
#> 9  1.0081016         b  0.7786920
#> 10 2.9852986         b  2.3059462

## check solution
xtabs(new_weight ~ dat$household)
#> dat$household
#> a b c 
#> 3 4 5 

## calculate weights "by reference"
ipf_step_ref(dat$weight, dat$household, targets)
dat
#>       weight household
#> 1  0.9153619         b
#> 2  1.0168616         c
#> 3  0.5809330         c
#> 4  1.7643151         c
#> 5  1.7081908         a
#> 6  0.4325056         c
#> 7  1.2918092         a
#> 8  1.2053846         c
#> 9  0.7786920         b
#> 10 2.3059462         b

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
#> 1   Dinner Female     No Sun  9
#> 2    Lunch   Male    Yes Mon  8
#> 3.1  Lunch   Male     No Tue 20
#> 2.1  Lunch   Male    Yes Mon  8
#> 2.2  Lunch   Male    Yes Mon  8
#> 2.3  Lunch   Male    Yes Mon  8
#> 1.1 Dinner Female     No Sun  9
#> 1.2 Dinner Female     No Sun  9
#> 2.4  Lunch   Male    Yes Mon  8

## adjust weights
weight <- rnorm(nrow(tips)) + 5
adjusted_weight <- ipf_step(weight, cf, con)

## check outputs
con2 <- xtabs(adjusted_weight ~ ., data = tips)
sum((con - con2)^2)
#> [1] 0
```
