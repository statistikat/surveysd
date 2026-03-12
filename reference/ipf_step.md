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
#>        weight household
#> 1  0.79662907         b
#> 2  0.69127271         b
#> 3  0.83058150         b
#> 4  1.23208006         b
#> 5  0.09165295         b
#> 6  1.78321648         a
#> 7  0.15028812         c
#> 8  2.32318788         b
#> 9  6.91417013         a
#> 10 0.51374262         a

## create targets (same lenght as classLabels!)
targets <- 3:5

## calculate weights
new_weight <- ipf_step(dat$weight, dat$household, targets)
cbind(dat, new_weight)
#>        weight household new_weight
#> 1  0.79662907         b 0.53416603
#> 2  0.69127271         b 0.46352112
#> 3  0.83058150         b 0.55693225
#> 4  1.23208006         b 0.82615027
#> 5  0.09165295         b 0.06145632
#> 6  1.78321648         a 0.58078106
#> 7  0.15028812         c 5.00000000
#> 8  2.32318788         b 1.55777400
#> 9  6.91417013         a 2.25189658
#> 10 0.51374262         a 0.16732236

## check solution
xtabs(new_weight ~ dat$household)
#> dat$household
#> a b c 
#> 3 4 5 

## calculate weights "by reference"
ipf_step_ref(dat$weight, dat$household, targets)
dat
#>        weight household
#> 1  0.53416603         b
#> 2  0.46352112         b
#> 3  0.55693225         b
#> 4  0.82615027         b
#> 5  0.06145632         b
#> 6  0.58078106         a
#> 7  5.00000000         c
#> 8  1.55777400         b
#> 9  2.25189658         a
#> 10 0.16732236         a

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
#> 3    Lunch   Male     No Tue 20
#> 3.1  Lunch   Male     No Tue 20
#> 3.2  Lunch   Male     No Tue 20
#> 2.1  Lunch   Male    Yes Mon  8
#> 1.1 Dinner Female     No Sun  9
#> 1.2 Dinner Female     No Sun  9
#> 3.3  Lunch   Male     No Tue 20
#> 3.4  Lunch   Male     No Tue 20

## adjust weights
weight <- rnorm(nrow(tips)) + 5
adjusted_weight <- ipf_step(weight, cf, con)

## check outputs
con2 <- xtabs(adjusted_weight ~ ., data = tips)
sum((con - con2)^2)
#> [1] 0
```
