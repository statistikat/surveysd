# Error estimation

For the most part, this document will present the functionalities of the
function
[`surveysd::calc.stError()`](https://statistikat.github.io/surveysd/reference/calc.stError.md)
which generates point estimates and standard errors for user-supplied
estimation functions.

## Prerequisites

In order to use a dataset with
[`calc.stError()`](https://statistikat.github.io/surveysd/reference/calc.stError.md),
several weight columns have to be present. Each weight column
corresponds to a bootstrap sample. In the following examples, we will
use the data from
[`demo.eusilc()`](https://statistikat.github.io/surveysd/reference/demo.eusilc.md)
and attach the bootstrap weights using
[`draw.bootstrap()`](https://statistikat.github.io/surveysd/reference/draw.bootstrap.md)
and
[`recalib()`](https://statistikat.github.io/surveysd/reference/recalib.md).
Please refer to the documentation of those functions for more detail.

``` r
library(surveysd)

set.seed(1234)
eusilc <- demo.eusilc(prettyNames = TRUE)
dat_boot <- draw.bootstrap(eusilc, REP = 10, hid = "hid", weights = "pWeight",
                           strata = "region", period = "year")
dat_boot_calib <- recalib(dat_boot, conP.var = "gender", conH.var = "region",
                          epsP = 1e-2, epsH = 2.5e-2, verbose = FALSE)
dat_boot_calib[, onePerson := nrow(.SD) == 1, by = .(year, hid)]

## print part of the dataset
dat_boot_calib[1:5, .(year, povertyRisk, eqIncome, onePerson, pWeight, w1, w2, w3, w4, w5)]
```

    ##     year povertyRisk eqIncome onePerson  pWeight          w1           w2
    ##    <num>      <lgcl>    <num>    <lgcl>    <num>       <num>        <num>
    ## 1:  2010       FALSE 16090.69     FALSE 504.5696   0.4518189 1007.6108352
    ## 2:  2010       FALSE 16090.69     FALSE 504.5696   0.4518189 1007.6108352
    ## 3:  2010       FALSE 16090.69     FALSE 504.5696   0.4518189 1007.6108352
    ## 4:  2010       FALSE 27076.24     FALSE 493.3824 993.8351804    0.4385411
    ## 5:  2010       FALSE 27076.24     FALSE 493.3824 993.8351804    0.4385411
    ##           w3        w4          w5
    ##        <num>     <num>       <num>
    ## 1: 1007.6765 0.4535474 997.8806647
    ## 2: 1007.6765 0.4535474 997.8806647
    ## 3: 1007.6765 0.4535474 997.8806647
    ## 4:  987.3836 0.4433913   0.4336001
    ## 5:  987.3836 0.4433913   0.4336001

## Estimator functions

The parameters `fun` and `var` in
[`calc.stError()`](https://statistikat.github.io/surveysd/reference/calc.stError.md)
define the estimator to be used in the error analysis. There are two
built-in estimator functions
[`weightedSum()`](https://statistikat.github.io/surveysd/reference/PointEstimates.md)
and
[`weightedRatio()`](https://statistikat.github.io/surveysd/reference/PointEstimates.md)
which can be used as follows.

``` r
povertyRate <- calc.stError(dat_boot_calib, var = "povertyRisk", fun = weightedRatio)
totalIncome <- calc.stError(dat_boot_calib, var = "eqIncome", fun = weightedSum)
```

Those functions calculate the ratio of persons at risk of poverty (in
percent) and the total income. By default, the results are calculated
separately for each reference period.

``` r
povertyRate$Estimates
```

    ## Key: <year, n, N, estimate_type>
    ##     year     n       N estimate_type val_povertyRisk stE_povertyRisk
    ##    <num> <int>   <num>        <char>           <num>           <num>
    ## 1:  2010 14827 8182222        direct        14.44422       0.4124800
    ## 2:  2011 14827 8182222        direct        14.77393       0.3000355
    ## 3:  2012 14827 8182222        direct        15.04515       0.3160549
    ## 4:  2013 14827 8182222        direct        14.89013       0.4722579
    ## 5:  2014 14827 8182222        direct        15.14556       0.2692546
    ## 6:  2015 14827 8182222        direct        15.53640       0.2061388
    ## 7:  2016 14827 8182222        direct        15.08315       0.3997422
    ## 8:  2017 14827 8182222        direct        15.42019       0.4979087

``` r
totalIncome$Estimates
```

    ## Key: <year, n, N, estimate_type>
    ##     year     n       N estimate_type val_eqIncome stE_eqIncome
    ##    <num> <int>   <num>        <char>        <num>        <num>
    ## 1:  2010 14827 8182222        direct 162750998071    939434272
    ## 2:  2011 14827 8182222        direct 161926931417   1009565325
    ## 3:  2012 14827 8182222        direct 162576509628   1176900540
    ## 4:  2013 14827 8182222        direct 163199507862   1236142040
    ## 5:  2014 14827 8182222        direct 163986275009   1100511806
    ## 6:  2015 14827 8182222        direct 163416275447    898033567
    ## 7:  2016 14827 8182222        direct 162706205137   1056587828
    ## 8:  2017 14827 8182222        direct 164314959107   1880465101

Columns that use the `val_` prefix denote the point estimate belonging
to the “main weight” of the dataset, which is `pWeight` in case of the
dataset used here.

Columns with the `stE_` prefix denote standard errors calculated with
bootstrap replicates. The replicates result in using `w1`, `w2`, …,
`w10` instead of `pWeight` when applying the estimator.

`n` denotes the number of observations for the year and `N` denotes the
total weight of those persons.

### Custom estimators

In order to define a custom estimator function to be used in `fun`, the
function needs to have at least two arguments like the example below.

``` r
## define custom estimator
myWeightedSum <- function(x, w) {
  sum(x*w)
}

## check if results are equal to the one using `surveysd::weightedSum()`
totalIncome2 <- calc.stError(dat_boot_calib, var = "eqIncome", fun = myWeightedSum)
all.equal(totalIncome$Estimates, totalIncome2$Estimates)
```

    ## [1] TRUE

The parameters `x` and `w` can be assumed to be vectors with equal
length with `w` being numeric weight vector and `x` being the column
defined in the `var` argument. It will be called once for each period
(in this case `year`) and for each weight column (in this case
`pWeight`, `w1`, `w2`, …, `w10`).

Custom estimators using additional parameters can also be supplied and
parameter `add.arg` can be used to set the additional arguments for the
custom estimator.

``` r
## use add.arg-argument
fun <- function(x, w, b) {
  sum(x*w*b)
}
add.arg = list(b="onePerson")

err.est <- calc.stError(dat_boot_calib, var = "povertyRisk", fun = fun,
                        period.mean = 0, add.arg=add.arg)
err.est$Estimates
```

    ## Key: <year, n, N, estimate_type>
    ##     year     n       N estimate_type val_povertyRisk stE_povertyRisk
    ##    <num> <int>   <num>        <char>           <num>           <num>
    ## 1:  2010 14827 8182222        direct        273683.9        9905.294
    ## 2:  2011 14827 8182222        direct        261883.6        7593.573
    ## 3:  2012 14827 8182222        direct        243083.9       13394.882
    ## 4:  2013 14827 8182222        direct        238004.4       16860.720
    ## 5:  2014 14827 8182222        direct        218572.1       14619.394
    ## 6:  2015 14827 8182222        direct        219984.1       11200.726
    ## 7:  2016 14827 8182222        direct        201753.9       10828.578
    ## 8:  2017 14827 8182222        direct        196881.2       15022.007

``` r
# compare with direct computation
compare.value <- dat_boot_calib[,fun(povertyRisk,pWeight,b=onePerson),
                                 by=c("year")]
all((compare.value$V1-err.est$Estimates$val_povertyRisk)==0)
```

    ## [1] TRUE

The above chunk computes the weighted poverty ratio for single person
households.

### Adjust variable depending on bootstrap weights

In our example the variable `povertyRisk` is a boolean and is `TRUE` if
the income is less than 60% of the weighted median income. Thus it
directly depends on the original weight vector `pWeight`. To further
reduce the estimated error one should calculate for each bootstrap
replicate weight $w$ the weighted median income $medIncome_{w}$ and then
define $povertyRisk_{w}$ as

\$\$ povertyRisk_w = \cases{1 \quad\text{if Income}\<0.6\cdot
medIncome\_{w}\\ 0 \quad\text{else}} \$\$

The estimator can then be applied to the new variable $povertyRisk_{w}$.
This can be realized using a custom estimator function.

``` r
# custom estimator to first derive poverty threshold 
# and then estimate a weighted ratio
povmd <- function(x, w) {
 md <- laeken::weightedMedian(x, w)*0.6
 pmd60 <- x < md
 # weighted ratio is directly estimated inside the function
 return(sum(w[pmd60])/sum(w)*100)
}

err.est <- calc.stError(
  dat_boot_calib, var = "povertyRisk", fun = weightedRatio,
  fun.adjust.var = povmd, adjust.var = "eqIncome")
err.est$Estimates
```

    ## Key: <year, n, N, estimate_type>
    ##     year     n       N estimate_type val_povertyRisk stE_povertyRisk
    ##    <num> <int>   <num>        <char>           <num>           <num>
    ## 1:  2010 14827 8182222        direct        14.44422               0
    ## 2:  2011 14827 8182222        direct        14.77393               0
    ## 3:  2012 14827 8182222        direct        15.04515               0
    ## 4:  2013 14827 8182222        direct        14.89013               0
    ## 5:  2014 14827 8182222        direct        15.14556               0
    ## 6:  2015 14827 8182222        direct        15.53640               0
    ## 7:  2016 14827 8182222        direct        15.08315               0
    ## 8:  2017 14827 8182222        direct        15.42019               0

The approach shown above is only valid if no grouping variables are
supplied (parameter `group = NULL`). If grouping variables are supplied
one should use parameters `fun.adjust.var` and `adjust.var` such that
the $povertyRisk_{w}$ is first calculated for each `period` and then
used for each grouping in `group`.

``` r
# using fun.adjust.var and adjust.var to estimate povmd60 indicator
# for each period and bootstrap weight before applying the weightedRatio
povmd2 <- function(x, w) {
 md <- laeken::weightedMedian(x, w)*0.6
 pmd60 <- x < md
 return(as.integer(pmd60))
}

# set adjust.var="eqIncome" so the income vector is used to estimate
# the povmd60 indicator for each bootstrap weight
# and the resulting indicators are passed to function weightedRatio
group <- "gender"
err.est <- calc.stError(
  dat_boot_calib, var = "povertyRisk", fun = weightedRatio, group = "gender",
  fun.adjust.var = povmd2, adjust.var = "eqIncome")
err.est$Estimates
```

    ## Key: <year, n, N, gender, estimate_type>
    ##      year     n       N gender estimate_type val_povertyRisk stE_povertyRisk
    ##     <num> <int>   <num> <fctr>        <char>           <num>           <num>
    ##  1:  2010  7267 3979572   male        direct        12.02660       0.4443128
    ##  2:  2010  7560 4202650 female        direct        16.73351       0.3484212
    ##  3:  2010 14827 8182222   <NA>        direct        14.44422       0.3157465
    ##  4:  2011  7267 3979572   male        direct        12.81921       0.3881769
    ##  5:  2011  7560 4202650 female        direct        16.62488       0.5057501
    ##  6:  2011 14827 8182222   <NA>        direct        14.77393       0.3594573
    ##  7:  2012  7267 3979572   male        direct        13.76065       0.4948190
    ##  8:  2012  7560 4202650 female        direct        16.26147       0.3970795
    ##  9:  2012 14827 8182222   <NA>        direct        15.04515       0.3443189
    ## 10:  2013  7267 3979572   male        direct        13.88962       0.5389362
    ## 11:  2013  7560 4202650 female        direct        15.83754       0.5849565
    ## 12:  2013 14827 8182222   <NA>        direct        14.89013       0.4842279
    ## 13:  2014  7267 3979572   male        direct        14.50351       0.3684713
    ## 14:  2014  7560 4202650 female        direct        15.75353       0.4739899
    ## 15:  2014 14827 8182222   <NA>        direct        15.14556       0.3772557
    ## 16:  2015  7267 3979572   male        direct        15.12289       0.3950068
    ## 17:  2015  7560 4202650 female        direct        15.92796       0.3616357
    ## 18:  2015 14827 8182222   <NA>        direct        15.53640       0.3435780
    ## 19:  2016  7267 3979572   male        direct        14.57968       0.4080958
    ## 20:  2016  7560 4202650 female        direct        15.55989       0.5020602
    ## 21:  2016 14827 8182222   <NA>        direct        15.08315       0.4061863
    ## 22:  2017  7267 3979572   male        direct        14.94816       0.4273152
    ## 23:  2017  7560 4202650 female        direct        15.86717       0.5340769
    ## 24:  2017 14827 8182222   <NA>        direct        15.42019       0.4265635
    ##      year     n       N gender estimate_type val_povertyRisk stE_povertyRisk

### Multiple estimators

In case an estimator should be applied to several columns of the
dataset, `var` can be set to a vector containing all necessary columns.

``` r
multipleRates <- calc.stError(dat_boot_calib, var = c("povertyRisk", "onePerson"), fun = weightedRatio)
multipleRates$Estimates
```

    ## Key: <year, n, N, estimate_type>
    ##     year     n       N estimate_type val_povertyRisk stE_povertyRisk
    ##    <num> <int>   <num>        <char>           <num>           <num>
    ## 1:  2010 14827 8182222        direct        14.44422       0.4313838
    ## 2:  2011 14827 8182222        direct        14.77393       0.3243568
    ## 3:  2012 14827 8182222        direct        15.04515       0.4026730
    ## 4:  2013 14827 8182222        direct        14.89013       0.4760682
    ## 5:  2014 14827 8182222        direct        15.14556       0.4339296
    ## 6:  2015 14827 8182222        direct        15.53640       0.5777787
    ## 7:  2016 14827 8182222        direct        15.08315       0.4080630
    ## 8:  2017 14827 8182222        direct        15.42019       0.5485506
    ##    val_onePerson stE_onePerson
    ##            <num>         <num>
    ## 1:      14.85737     0.4313838
    ## 2:      14.85737     0.3243568
    ## 3:      14.85737     0.4026730
    ## 4:      14.85737     0.4760682
    ## 5:      14.85737     0.4339296
    ## 6:      14.85737     0.5777787
    ## 7:      14.85737     0.4080630
    ## 8:      14.85737     0.5485506

Here we see the relative number of persons at risk of poverty and the
relative number of one-person households.

## Grouping

The `groups` argument can be used to calculate estimators for different
subsets of the data. This argument can take the grouping variable as a
string that refers to a column name (usually a factor) in `dat`. If set,
all estimators are not only split by the reference period but also by
the grouping variable. For simplicity, only one reference period of the
above data is used.

``` r
dat2 <- subset(dat_boot_calib, year == 2010)
for (att  in c("period", "weights", "b.rep"))
  attr(dat2, att) <- attr(dat_boot_calib, att)
```

To calculate the ratio of persons at risk of poverty for each federal
state of Austria, `group = "region"` can be used.

``` r
povertyRates <- calc.stError(dat2, var = "povertyRisk", fun = weightedRatio, group = "region")
povertyRates$Estimates
```

    ## Key: <year, n, N, region, estimate_type>
    ##      year     n       N        region estimate_type val_povertyRisk
    ##     <num> <int>   <num>        <fctr>        <char>           <num>
    ##  1:  2010   549  260564    Burgenland        direct        19.53984
    ##  2:  2010   733  377355    Vorarlberg        direct        16.53731
    ##  3:  2010   924  535451      Salzburg        direct        13.78734
    ##  4:  2010  1078  563648     Carinthia        direct        13.08627
    ##  5:  2010  1317  701899         Tyrol        direct        15.30819
    ##  6:  2010  2295 1167045        Styria        direct        14.37464
    ##  7:  2010  2322 1598931        Vienna        direct        17.23468
    ##  8:  2010  2804 1555709 Lower Austria        direct        13.84362
    ##  9:  2010  2805 1421620 Upper Austria        direct        10.88977
    ## 10:  2010 14827 8182222          <NA>        direct        14.44422
    ##     stE_povertyRisk
    ##               <num>
    ##  1:       3.2383552
    ##  2:       2.3306155
    ##  3:       2.2953986
    ##  4:       1.2913181
    ##  5:       1.8957772
    ##  6:       0.8869132
    ##  7:       1.3046408
    ##  8:       0.7826762
    ##  9:       0.9339947
    ## 10:       0.4124800

The last row with `region = NA` denotes the aggregate over all regions.
Note that the columns `N` and `n` now show the weighted and unweighted
number of persons in each region.

### Several grouping variables

In case more than one grouping variable is used, there are several
options of calling
[`calc.stError()`](https://statistikat.github.io/surveysd/reference/calc.stError.md)
depending on whether combinations of grouping levels should be regarded
or not. We will consider the variables `gender` and `region` as our
grouping variables and show three options on how
[`calc.stError()`](https://statistikat.github.io/surveysd/reference/calc.stError.md)
can be called.

#### Option 1: All regions and all genders

Calculate the point estimate and standard error for each region and each
gender. The number of rows in the output is therefore

$$n_{\text{periods}} \cdot \left( n_{\text{regions}} + n_{\text{genders}} + 1 \right) = 1 \cdot (9 + 2 + 1) = 12.$$

The last row is again the estimate for the whole period.

``` r
povertyRates <- calc.stError(dat2, var = "povertyRisk", fun = weightedRatio, 
                             group = c("gender", "region"))
povertyRates$Estimates
```

    ## Key: <year, n, N, gender, region, estimate_type>
    ##      year     n       N gender        region estimate_type val_povertyRisk
    ##     <num> <int>   <num> <fctr>        <fctr>        <char>           <num>
    ##  1:  2010   549  260564   <NA>    Burgenland        direct        19.53984
    ##  2:  2010   733  377355   <NA>    Vorarlberg        direct        16.53731
    ##  3:  2010   924  535451   <NA>      Salzburg        direct        13.78734
    ##  4:  2010  1078  563648   <NA>     Carinthia        direct        13.08627
    ##  5:  2010  1317  701899   <NA>         Tyrol        direct        15.30819
    ##  6:  2010  2295 1167045   <NA>        Styria        direct        14.37464
    ##  7:  2010  2322 1598931   <NA>        Vienna        direct        17.23468
    ##  8:  2010  2804 1555709   <NA> Lower Austria        direct        13.84362
    ##  9:  2010  2805 1421620   <NA> Upper Austria        direct        10.88977
    ## 10:  2010  7267 3979572   male          <NA>        direct        12.02660
    ## 11:  2010  7560 4202650 female          <NA>        direct        16.73351
    ## 12:  2010 14827 8182222   <NA>          <NA>        direct        14.44422
    ##     stE_povertyRisk
    ##               <num>
    ##  1:       3.2383552
    ##  2:       2.3306155
    ##  3:       2.2953986
    ##  4:       1.2913181
    ##  5:       1.8957772
    ##  6:       0.8869132
    ##  7:       1.3046408
    ##  8:       0.7826762
    ##  9:       0.9339947
    ## 10:       0.4599945
    ## 11:       0.5072902
    ## 12:       0.4124800

#### Option 2: All combinations of `region` and `gender`

Split the data by all combinations of the two grouping variables. This
will result in a larger output-table of the size

$$n_{\text{periods}} \cdot \left( n_{\text{regions}} \cdot n_{\text{genders}} + 1 \right) = 1 \cdot (9 \cdot 2 + 1) = 19.$$

``` r
povertyRates <- calc.stError(dat2, var = "povertyRisk", fun = weightedRatio, 
                             group = list(c("gender", "region")))
povertyRates$Estimates
```

    ## Key: <year, n, N, gender, region, estimate_type>
    ##      year     n         N gender        region estimate_type val_povertyRisk
    ##     <num> <int>     <num> <fctr>        <fctr>        <char>           <num>
    ##  1:  2010   261  122741.8   male    Burgenland        direct       17.414524
    ##  2:  2010   288  137822.2 female    Burgenland        direct       21.432598
    ##  3:  2010   359  182732.9   male    Vorarlberg        direct       12.973259
    ##  4:  2010   374  194622.1 female    Vorarlberg        direct       19.883637
    ##  5:  2010   440  253143.7   male      Salzburg        direct        9.156964
    ##  6:  2010   484  282307.3 female      Salzburg        direct       17.939382
    ##  7:  2010   517  268581.4   male     Carinthia        direct       10.552149
    ##  8:  2010   561  295066.6 female     Carinthia        direct       15.392924
    ##  9:  2010   650  339566.5   male         Tyrol        direct       12.857542
    ## 10:  2010   667  362332.5 female         Tyrol        direct       17.604861
    ## 11:  2010  1128  571011.7   male        Styria        direct       11.671247
    ## 12:  2010  1132  774405.4   male        Vienna        direct       15.590616
    ## 13:  2010  1167  596033.3 female        Styria        direct       16.964539
    ## 14:  2010  1190  824525.6 female        Vienna        direct       18.778813
    ## 15:  2010  1363  684272.5   male Upper Austria        direct        9.074690
    ## 16:  2010  1387  772593.2 female Lower Austria        direct       16.372949
    ## 17:  2010  1417  783115.8   male Lower Austria        direct       11.348283
    ## 18:  2010  1442  737347.5 female Upper Austria        direct       12.574206
    ## 19:  2010 14827 8182222.0   <NA>          <NA>        direct       14.444218
    ##     stE_povertyRisk
    ##               <num>
    ##  1:       3.4678020
    ##  2:       3.5372448
    ##  3:       2.0822938
    ##  4:       2.8106479
    ##  5:       1.8033118
    ##  6:       2.9546983
    ##  7:       1.6941596
    ##  8:       1.3452298
    ##  9:       2.3588456
    ## 10:       1.9160929
    ## 11:       0.9317810
    ## 12:       1.7825863
    ## 13:       1.0493806
    ## 14:       1.2113023
    ## 15:       0.6731709
    ## 16:       1.2390922
    ## 17:       0.6078340
    ## 18:       1.3523566
    ## 19:       0.4124800

#### Option 3: Cobination of Option 1 and Option 2

In this case, the estimates and standard errors are calculated for

- every gender,
- every region and
- every combination of region and gender.

The number of rows in the output is therefore

$$n_{\text{periods}} \cdot \left( n_{\text{regions}} \cdot n_{\text{genders}} + n_{\text{regions}} + n_{\text{genders}} + 1 \right) = 1 \cdot (9 \cdot 2 + 9 + 2 + 1) = 30.$$

``` r
povertyRates <- calc.stError(dat2, var = "povertyRisk", fun = weightedRatio, 
                             group = list("gender", "region", c("gender", "region")))
povertyRates$Estimates
```

    ## Key: <year, n, N, gender, region, estimate_type>
    ##      year     n         N gender        region estimate_type val_povertyRisk
    ##     <num> <int>     <num> <fctr>        <fctr>        <char>           <num>
    ##  1:  2010   261  122741.8   male    Burgenland        direct       17.414524
    ##  2:  2010   288  137822.2 female    Burgenland        direct       21.432598
    ##  3:  2010   359  182732.9   male    Vorarlberg        direct       12.973259
    ##  4:  2010   374  194622.1 female    Vorarlberg        direct       19.883637
    ##  5:  2010   440  253143.7   male      Salzburg        direct        9.156964
    ##  6:  2010   484  282307.3 female      Salzburg        direct       17.939382
    ##  7:  2010   517  268581.4   male     Carinthia        direct       10.552149
    ##  8:  2010   549  260564.0   <NA>    Burgenland        direct       19.539837
    ##  9:  2010   561  295066.6 female     Carinthia        direct       15.392924
    ## 10:  2010   650  339566.5   male         Tyrol        direct       12.857542
    ## 11:  2010   667  362332.5 female         Tyrol        direct       17.604861
    ## 12:  2010   733  377355.0   <NA>    Vorarlberg        direct       16.537310
    ## 13:  2010   924  535451.0   <NA>      Salzburg        direct       13.787343
    ## 14:  2010  1078  563648.0   <NA>     Carinthia        direct       13.086268
    ## 15:  2010  1128  571011.7   male        Styria        direct       11.671247
    ## 16:  2010  1132  774405.4   male        Vienna        direct       15.590616
    ## 17:  2010  1167  596033.3 female        Styria        direct       16.964539
    ## 18:  2010  1190  824525.6 female        Vienna        direct       18.778813
    ## 19:  2010  1317  701899.0   <NA>         Tyrol        direct       15.308190
    ## 20:  2010  1363  684272.5   male Upper Austria        direct        9.074690
    ## 21:  2010  1387  772593.2 female Lower Austria        direct       16.372949
    ## 22:  2010  1417  783115.8   male Lower Austria        direct       11.348283
    ## 23:  2010  1442  737347.5 female Upper Austria        direct       12.574206
    ## 24:  2010  2295 1167045.0   <NA>        Styria        direct       14.374637
    ## 25:  2010  2322 1598931.0   <NA>        Vienna        direct       17.234683
    ## 26:  2010  2804 1555709.0   <NA> Lower Austria        direct       13.843623
    ## 27:  2010  2805 1421620.0   <NA> Upper Austria        direct       10.889773
    ## 28:  2010  7267 3979571.7   male          <NA>        direct       12.026600
    ## 29:  2010  7560 4202650.3 female          <NA>        direct       16.733508
    ## 30:  2010 14827 8182222.0   <NA>          <NA>        direct       14.444218
    ##      year     n         N gender        region estimate_type val_povertyRisk
    ##     stE_povertyRisk
    ##               <num>
    ##  1:       3.4678020
    ##  2:       3.5372448
    ##  3:       2.0822938
    ##  4:       2.8106479
    ##  5:       1.8033118
    ##  6:       2.9546983
    ##  7:       1.6941596
    ##  8:       3.2383552
    ##  9:       1.3452298
    ## 10:       2.3588456
    ## 11:       1.9160929
    ## 12:       2.3306155
    ## 13:       2.2953986
    ## 14:       1.2913181
    ## 15:       0.9317810
    ## 16:       1.7825863
    ## 17:       1.0493806
    ## 18:       1.2113023
    ## 19:       1.8957772
    ## 20:       0.6731709
    ## 21:       1.2390922
    ## 22:       0.6078340
    ## 23:       1.3523566
    ## 24:       0.8869132
    ## 25:       1.3046408
    ## 26:       0.7826762
    ## 27:       0.9339947
    ## 28:       0.4599945
    ## 29:       0.5072902
    ## 30:       0.4124800
    ##     stE_povertyRisk

## Group differences

If differences between groups need to be calculated, e.g difference of
poverty rates between `gender = "male"` and `gender = "female"`,
parameter `group.diff` can be utilised. Setting `group.diff = TRUE` the
differences and the standard error of these differences for all
variables defined in `groups` will be calculated.

``` r
povertyRates <- calc.stError(dat2, var = "povertyRisk", fun = weightedRatio, 
                             group = c("gender", "region"),
                             group.diff = TRUE)
povertyRates$Estimates
```

    ## Key: <year, n, N, gender, region, estimate_type>
    ##      year       n         N        gender                        region
    ##     <num>   <num>     <num>        <fctr>                        <fctr>
    ##  1:  2010   549.0  260564.0          <NA>                    Burgenland
    ##  2:  2010   641.0  318959.5          <NA>       Burgenland - Vorarlberg
    ##  3:  2010   733.0  377355.0          <NA>                    Vorarlberg
    ##  4:  2010   736.5  398007.5          <NA>         Burgenland - Salzburg
    ##  5:  2010   813.5  412106.0          <NA>        Burgenland - Carinthia
    ##  6:  2010   828.5  456403.0          <NA>         Salzburg - Vorarlberg
    ##  7:  2010   905.5  470501.5          <NA>        Carinthia - Vorarlberg
    ##  8:  2010   924.0  535451.0          <NA>                      Salzburg
    ##  9:  2010   933.0  481231.5          <NA>            Burgenland - Tyrol
    ## 10:  2010  1001.0  549549.5          <NA>          Carinthia - Salzburg
    ## 11:  2010  1025.0  539627.0          <NA>            Tyrol - Vorarlberg
    ## 12:  2010  1078.0  563648.0          <NA>                     Carinthia
    ## 13:  2010  1120.5  618675.0          <NA>              Salzburg - Tyrol
    ## 14:  2010  1197.5  632773.5          <NA>             Carinthia - Tyrol
    ## 15:  2010  1317.0  701899.0          <NA>                         Tyrol
    ## 16:  2010  1422.0  713804.5          <NA>           Burgenland - Styria
    ## 17:  2010  1435.5  929747.5          <NA>           Burgenland - Vienna
    ## 18:  2010  1514.0  772200.0          <NA>           Styria - Vorarlberg
    ## 19:  2010  1527.5  988143.0          <NA>           Vienna - Vorarlberg
    ## 20:  2010  1609.5  851248.0          <NA>             Salzburg - Styria
    ## 21:  2010  1623.0 1067191.0          <NA>             Salzburg - Vienna
    ## 22:  2010  1676.5  908136.5          <NA>    Burgenland - Lower Austria
    ## 23:  2010  1677.0  841092.0          <NA>    Burgenland - Upper Austria
    ## 24:  2010  1686.5  865346.5          <NA>            Carinthia - Styria
    ## 25:  2010  1700.0 1081289.5          <NA>            Carinthia - Vienna
    ## 26:  2010  1768.5  966532.0          <NA>    Lower Austria - Vorarlberg
    ## 27:  2010  1769.0  899487.5          <NA>    Upper Austria - Vorarlberg
    ## 28:  2010  1806.0  934472.0          <NA>                Styria - Tyrol
    ## 29:  2010  1819.5 1150415.0          <NA>                Tyrol - Vienna
    ## 30:  2010  1864.0 1045580.0          <NA>      Lower Austria - Salzburg
    ## 31:  2010  1864.5  978535.5          <NA>      Salzburg - Upper Austria
    ## 32:  2010  1941.0 1059678.5          <NA>     Carinthia - Lower Austria
    ## 33:  2010  1941.5  992634.0          <NA>     Carinthia - Upper Austria
    ## 34:  2010  2060.5 1128804.0          <NA>         Lower Austria - Tyrol
    ## 35:  2010  2061.0 1061759.5          <NA>         Tyrol - Upper Austria
    ## 36:  2010  2295.0 1167045.0          <NA>                        Styria
    ## 37:  2010  2308.5 1382988.0          <NA>               Styria - Vienna
    ## 38:  2010  2322.0 1598931.0          <NA>                        Vienna
    ## 39:  2010  2549.5 1361377.0          <NA>        Lower Austria - Styria
    ## 40:  2010  2550.0 1294332.5          <NA>        Styria - Upper Austria
    ## 41:  2010  2563.0 1577320.0          <NA>        Lower Austria - Vienna
    ## 42:  2010  2563.5 1510275.5          <NA>        Upper Austria - Vienna
    ## 43:  2010  2804.0 1555709.0          <NA>                 Lower Austria
    ## 44:  2010  2804.5 1488664.5          <NA> Lower Austria - Upper Austria
    ## 45:  2010  2805.0 1421620.0          <NA>                 Upper Austria
    ## 46:  2010  7267.0 3979571.7          male                          <NA>
    ## 47:  2010  7413.5 4091111.0 male - female                          <NA>
    ## 48:  2010  7560.0 4202650.3        female                          <NA>
    ## 49:  2010 14827.0 8182222.0          <NA>                          <NA>
    ##      year       n         N        gender                        region
    ##        estimate_type val_povertyRisk stE_povertyRisk
    ##               <char>           <num>           <num>
    ##  1:           direct     19.53983651       3.2383552
    ##  2: group difference      3.00252634       2.8689577
    ##  3:           direct     16.53731017       2.3306155
    ##  4: group difference      5.75249330       3.4112935
    ##  5: group difference      6.45356876       3.2367403
    ##  6: group difference     -2.74996696       3.5145999
    ##  7: group difference     -3.45104242       2.4970123
    ##  8:           direct     13.78734321       2.2953986
    ##  9: group difference      4.23164602       4.1855567
    ## 10: group difference     -0.70107546       2.9800719
    ## 11: group difference     -1.22911968       3.6953361
    ## 12:           direct     13.08626775       1.2913181
    ## 13: group difference     -1.52084728       3.4672708
    ## 14: group difference     -2.22192274       1.5913467
    ## 15:           direct     15.30819049       1.8957772
    ## 16: group difference      5.16519923       3.1382729
    ## 17: group difference      2.30515330       2.6806674
    ## 18: group difference     -2.16267289       2.3894575
    ## 19: group difference      0.69737304       1.6564706
    ## 20: group difference     -0.58729407       2.4648529
    ## 21: group difference     -3.44734000       2.8265825
    ## 22: group difference      5.69621369       3.4740716
    ## 23: group difference      8.65006312       3.4877400
    ## 24: group difference     -1.28836953       1.9976524
    ## 25: group difference     -4.14841546       1.9739810
    ## 26: group difference     -2.69368735       2.8144534
    ## 27: group difference     -5.64753678       2.4036757
    ## 28: group difference     -0.93355321       2.3374554
    ## 29: group difference     -1.92649272       2.7594162
    ## 30: group difference      0.05627961       2.0378119
    ## 31: group difference      2.89756982       2.8799229
    ## 32: group difference     -0.75735506       1.5516228
    ## 33: group difference      2.19649436       1.3016042
    ## 34: group difference     -1.46456768       1.9234989
    ## 35: group difference      4.41841710       1.9847285
    ## 36:           direct     14.37463728       0.8869132
    ## 37: group difference     -2.86004593       1.0714784
    ## 38:           direct     17.23468321       1.3046408
    ## 39: group difference     -0.53101447       1.2999287
    ## 40: group difference      3.48486389       1.5345470
    ## 41: group difference     -3.39106040       1.7357561
    ## 42: group difference     -6.34490982       1.6128502
    ## 43:           direct     13.84362281       0.7826762
    ## 44: group difference      2.95384943       1.1561239
    ## 45:           direct     10.88977339       0.9339947
    ## 46:           direct     12.02659998       0.4599945
    ## 47: group difference     -4.70690810       0.4978371
    ## 48:           direct     16.73350808       0.5072902
    ## 49:           direct     14.44421817       0.4124800
    ##        estimate_type val_povertyRisk stE_povertyRisk

The resulting output table contains 49 rows. 12 rows for all the direct
estimators

$$n_{\text{periods}} \cdot \left( n_{\text{regions}} + n_{\text{genders}} + 1 \right) = 1 \cdot (9 + 2 + 1) = 12,$$

and another 37 for all the differences within the variable `"gender"`
and `"region"` seperately. Variable `"gender"` has 2 unique values
(`unique(dat2$gender)`) resulting in 1 difference, ~ `gender = "male"` -
`gender = "female"` and variable `"region"` has 9 unique values
(`unique(dat2$region)`) resulting in

$$8 + 7 + 6 + 5 + 4 + 3 + 2 + 1 = \sum\limits_{1 = 1}^{9 - 1}i = 36$$

estimates. Thus the output contains 1 + 36 = 37 estimates with respect
to group differences.

If a combintaion of grouping variables is used in `group` and
`group.diff = TRUE` then differences between combinations will only be
calculated if one of the grouping variables differs. For example the
difference between the following groups would be calculated

- `gender = "female" & region = "Vienna"` -
  `gender = "male" & region = "Vienna"`
- `gender = "female" & region = "Vienna"` -
  `gender = "female" & region = "Salzburg"`
- `gender = "male" & region = "Salzburg"` -
  `gender = "female" & region = "Salzburg"`

The difference between `gender = "female" & region = "Vienna"` and
`gender = "male" & region = "Salzburg"` however would not be calculated.

Thus this leads to

$$2 \cdot \left( \sum\limits_{1 = 1}^{9 - 1}i \right) + 9 \cdot 1 = 81$$

results with respect to the differences. The Output contains an
additional column `estimate_type` and

``` r
povertyRates <- calc.stError(dat2, var = "povertyRisk", fun = weightedRatio, 
                             group = list(c("gender", "region")),
                             group.diff = TRUE)
povertyRates$Estimates[,.N,by=.(estimate_type)]
```

    ##       estimate_type     N
    ##              <char> <int>
    ## 1:           direct    19
    ## 2: group difference    81

## Differences between survey periods

Differences of estimates between `period`s can be calculated using
parameter `period.diff`. `period.diff` expects a character vector (if
not `NULL`) specifying for which `period`s the differences should be
calcualed for. The inputs should be specified in the form
`"period2" - "period1"`.

``` r
povertyRates <- calc.stError(dat_boot_calib[year>2013], var = "povertyRisk", fun = weightedRatio, 
                             period.diff = c("2017 - 2016", "2016 - 2015", "2015 - 2014"))
povertyRates$Estimates
```

    ## Key: <year, n, N, estimate_type>
    ##         year     n       N     estimate_type val_povertyRisk stE_povertyRisk
    ##       <char> <num>   <num>            <char>           <num>           <num>
    ## 1:      2014 14827 8182222            direct      15.1455601       0.2692546
    ## 2:      2015 14827 8182222            direct      15.5364014       0.2061388
    ## 3: 2015-2014 14827 8182222 period difference       0.3908413       0.3075315
    ## 4:      2016 14827 8182222            direct      15.0831502       0.3997422
    ## 5: 2016-2015 14827 8182222 period difference      -0.4532512       0.4047544
    ## 6:      2017 14827 8182222            direct      15.4201916       0.4979087
    ## 7: 2017-2016 14827 8182222 period difference       0.3370414       0.3858008

If additional grouping variables are supplied to
[`calc.stError()`](https://statistikat.github.io/surveysd/reference/calc.stError.md)
die differences across `period`s are also carried out for all variables
in `group`.

``` r
povertyRates <- calc.stError(dat_boot_calib[year>2013], var = "povertyRisk", fun = weightedRatio, 
                             group = "gender",
                             period.diff = c("2017 - 2016", "2016 - 2015", "2015 - 2014"))
povertyRates$Estimates
```

    ## Key: <year, n, N, gender, estimate_type>
    ##          year     n       N gender     estimate_type val_povertyRisk
    ##        <char> <num>   <num> <fctr>            <char>           <num>
    ##  1:      2014  7267 3979572   male            direct      14.5035068
    ##  2:      2014  7560 4202650 female            direct      15.7535328
    ##  3:      2014 14827 8182222   <NA>            direct      15.1455601
    ##  4:      2015  7267 3979572   male            direct      15.1228904
    ##  5:      2015  7560 4202650 female            direct      15.9279630
    ##  6:      2015 14827 8182222   <NA>            direct      15.5364014
    ##  7: 2015-2014  7267 3979572   male period difference       0.6193836
    ##  8: 2015-2014  7560 4202650 female period difference       0.1744301
    ##  9: 2015-2014 14827 8182222   <NA> period difference       0.3908413
    ## 10:      2016  7267 3979572   male            direct      14.5796824
    ## 11:      2016  7560 4202650 female            direct      15.5598937
    ## 12:      2016 14827 8182222   <NA>            direct      15.0831502
    ## 13: 2016-2015  7267 3979572   male period difference      -0.5432080
    ## 14: 2016-2015  7560 4202650 female period difference      -0.3680693
    ## 15: 2016-2015 14827 8182222   <NA> period difference      -0.4532512
    ## 16:      2017  7267 3979572   male            direct      14.9481591
    ## 17:      2017  7560 4202650 female            direct      15.8671684
    ## 18:      2017 14827 8182222   <NA>            direct      15.4201916
    ## 19: 2017-2016  7267 3979572   male period difference       0.3684767
    ## 20: 2017-2016  7560 4202650 female period difference       0.3072748
    ## 21: 2017-2016 14827 8182222   <NA> period difference       0.3370414
    ##          year     n       N gender     estimate_type val_povertyRisk
    ##     stE_povertyRisk
    ##               <num>
    ##  1:       0.3237634
    ##  2:       0.3342519
    ##  3:       0.2692546
    ##  4:       0.2876297
    ##  5:       0.2347633
    ##  6:       0.2061388
    ##  7:       0.3403725
    ##  8:       0.3885999
    ##  9:       0.3075315
    ## 10:       0.4320988
    ## 11:       0.4523708
    ## 12:       0.3997422
    ## 13:       0.4227987
    ## 14:       0.4192416
    ## 15:       0.4047544
    ## 16:       0.5528365
    ## 17:       0.5335792
    ## 18:       0.4979087
    ## 19:       0.4242154
    ## 20:       0.3706255
    ## 21:       0.3858008
    ##     stE_povertyRisk

## Averages across periods

With parameter `period.mean` averages across `period`s are calculated
additional. The parameter accepts only odd integer values. The resulting
table will contain the direct estimates as well as rolling averages of
length `period.mean`.

``` r
povertyRates <- calc.stError(dat_boot_calib[year>2013], var = "povertyRisk", fun = weightedRatio, 
                             period.mean = 3)
povertyRates$Estimates
```

    ## Key: <year, n, N, estimate_type>
    ##              year     n       N  estimate_type val_povertyRisk stE_povertyRisk
    ##            <char> <num>   <num>         <char>           <num>           <num>
    ## 1:           2014 14827 8182222         direct        15.14556       0.2692546
    ## 2: 2014_2015_2016 14827 8182222 period average        15.25504       0.2067857
    ## 3:           2015 14827 8182222         direct        15.53640       0.2061388
    ## 4: 2015_2016_2017 14827 8182222 period average        15.34658       0.2830507
    ## 5:           2016 14827 8182222         direct        15.08315       0.3997422
    ## 6:           2017 14827 8182222         direct        15.42019       0.4979087

if in addition the parameters `group` and/or `period.diff` are specified
then differences and groupings of averages will be calculated.

``` r
povertyRates <- calc.stError(dat_boot_calib[year>2013], var = "povertyRisk", fun = weightedRatio, 
                             period.mean = 3, period.diff = "2016 - 2015",
                             group = "gender")
povertyRates$Estimates
```

    ## Key: <year, n, N, gender, estimate_type>
    ##               year     n       N gender                      estimate_type
    ##             <char> <num>   <num> <fctr>                             <char>
    ##  1:           2014  7267 3979572   male                             direct
    ##  2:           2014  7560 4202650 female                             direct
    ##  3:           2014 14827 8182222   <NA>                             direct
    ##  4: 2014_2015_2016  7267 3979572   male                     period average
    ##  5: 2014_2015_2016  7560 4202650 female                     period average
    ##  6: 2014_2015_2016 14827 8182222   <NA>                     period average
    ##  7:           2015  7267 3979572   male                             direct
    ##  8:           2015  7560 4202650 female                             direct
    ##  9:           2015 14827 8182222   <NA>                             direct
    ## 10: 2015_2016_2017  7267 3979572   male                     period average
    ## 11: 2015_2016_2017  7560 4202650 female                     period average
    ## 12: 2015_2016_2017 14827 8182222   <NA>                     period average
    ## 13:           2016  7267 3979572   male                             direct
    ## 14:           2016  7560 4202650 female                             direct
    ## 15:           2016 14827 8182222   <NA>                             direct
    ## 16:      2016-2015  7267 3979572   male                  period difference
    ## 17:      2016-2015  7560 4202650 female                  period difference
    ## 18:      2016-2015 14827 8182222   <NA>                  period difference
    ## 19: 2016-2015_mean  7267 3979572   male difference between period averages
    ## 20: 2016-2015_mean  7560 4202650 female difference between period averages
    ## 21: 2016-2015_mean 14827 8182222   <NA> difference between period averages
    ## 22:           2017  7267 3979572   male                             direct
    ## 23:           2017  7560 4202650 female                             direct
    ## 24:           2017 14827 8182222   <NA>                             direct
    ##               year     n       N gender                      estimate_type
    ##     val_povertyRisk stE_povertyRisk
    ##               <num>           <num>
    ##  1:     14.50350682       0.3237634
    ##  2:     15.75353283       0.3342519
    ##  3:     15.14556006       0.2692546
    ##  4:     14.73535987       0.2635717
    ##  5:     15.74712982       0.2507021
    ##  6:     15.25503720       0.2067857
    ##  7:     15.12289042       0.2876297
    ##  8:     15.92796296       0.2347633
    ##  9:     15.53640136       0.2061388
    ## 10:     14.88357729       0.3291229
    ## 11:     15.78500836       0.3396680
    ## 12:     15.34658105       0.2830507
    ## 13:     14.57968239       0.4320988
    ## 14:     15.55989368       0.4523708
    ## 15:     15.08315018       0.3997422
    ## 16:     -0.54320803       0.4227987
    ## 17:     -0.36806928       0.4192416
    ## 18:     -0.45325118       0.4047544
    ## 19:      0.14821741       0.2049635
    ## 20:      0.03787854       0.1815158
    ## 21:      0.09154385       0.1809570
    ## 22:     14.94815906       0.5528365
    ## 23:     15.86716845       0.5335792
    ## 24:     15.42019160       0.4979087
    ##     val_povertyRisk stE_povertyRisk
