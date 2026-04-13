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

    ##     year povertyRisk eqIncome onePerson  pWeight           w1        w2
    ##    <num>      <lgcl>    <num>    <lgcl>    <num>        <num>     <num>
    ## 1:  2010       FALSE 16090.69     FALSE 504.5696 1017.1567980 0.4462212
    ## 2:  2010       FALSE 16090.69     FALSE 504.5696 1017.1567980 0.4462212
    ## 3:  2010       FALSE 16090.69     FALSE 504.5696 1017.1567980 0.4462212
    ## 4:  2010       FALSE 27076.24     FALSE 493.3824    0.4424102 0.4362065
    ## 5:  2010       FALSE 27076.24     FALSE 493.3824    0.4424102 0.4362065
    ##           w3           w4          w5
    ##        <num>        <num>       <num>
    ## 1: 1002.6493 1014.5680291   0.4498068
    ## 2: 1002.6493 1014.5680291   0.4498068
    ## 3: 1002.6493 1014.5680291   0.4498068
    ## 4:  980.2386    0.4416282 989.1489374
    ## 5:  980.2386    0.4416282 989.1489374

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
    ## 1:  2010 14827 8182222        direct        14.44422       0.7763632
    ## 2:  2011 14827 8182222        direct        14.77393       0.6689374
    ## 3:  2012 14827 8182222        direct        15.04515       0.5660591
    ## 4:  2013 14827 8182222        direct        14.89013       0.5233140
    ## 5:  2014 14827 8182222        direct        15.14556       0.4499388
    ## 6:  2015 14827 8182222        direct        15.53640       0.4542393
    ## 7:  2016 14827 8182222        direct        15.08315       0.4510617
    ## 8:  2017 14827 8182222        direct        15.42019       0.4195721

``` r
totalIncome$Estimates
```

    ## Key: <year, n, N, estimate_type>
    ##     year     n       N estimate_type val_eqIncome stE_eqIncome
    ##    <num> <int>   <num>        <char>        <num>        <num>
    ## 1:  2010 14827 8182222        direct 162750998071    984494436
    ## 2:  2011 14827 8182222        direct 161926931417   1168201212
    ## 3:  2012 14827 8182222        direct 162576509628   1510078186
    ## 4:  2013 14827 8182222        direct 163199507862   1756061906
    ## 5:  2014 14827 8182222        direct 163986275009   1544364971
    ## 6:  2015 14827 8182222        direct 163416275447   1407931979
    ## 7:  2016 14827 8182222        direct 162706205137   1488015338
    ## 8:  2017 14827 8182222        direct 164314959107   1123272745

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
    ## 1:  2010 14827 8182222        direct        273683.9       14342.813
    ## 2:  2011 14827 8182222        direct        261883.6       18596.933
    ## 3:  2012 14827 8182222        direct        243083.9       10306.230
    ## 4:  2013 14827 8182222        direct        238004.4        8017.457
    ## 5:  2014 14827 8182222        direct        218572.1       10412.411
    ## 6:  2015 14827 8182222        direct        219984.1        9711.024
    ## 7:  2016 14827 8182222        direct        201753.9        7588.736
    ## 8:  2017 14827 8182222        direct        196881.2       10699.702

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
    ##  1:  2010  7267 3979572   male        direct        12.02660       0.7369632
    ##  2:  2010  7560 4202650 female        direct        16.73351       0.8571299
    ##  3:  2010 14827 8182222   <NA>        direct        14.44422       0.7485161
    ##  4:  2011  7267 3979572   male        direct        12.81921       0.6261007
    ##  5:  2011  7560 4202650 female        direct        16.62488       0.8417165
    ##  6:  2011 14827 8182222   <NA>        direct        14.77393       0.6958582
    ##  7:  2012  7267 3979572   male        direct        13.76065       0.5475589
    ##  8:  2012  7560 4202650 female        direct        16.26147       0.7557210
    ##  9:  2012 14827 8182222   <NA>        direct        15.04515       0.6036640
    ## 10:  2013  7267 3979572   male        direct        13.88962       0.5558253
    ## 11:  2013  7560 4202650 female        direct        15.83754       0.6751942
    ## 12:  2013 14827 8182222   <NA>        direct        14.89013       0.5193916
    ## 13:  2014  7267 3979572   male        direct        14.50351       0.4936863
    ## 14:  2014  7560 4202650 female        direct        15.75353       0.6099725
    ## 15:  2014 14827 8182222   <NA>        direct        15.14556       0.4592074
    ## 16:  2015  7267 3979572   male        direct        15.12289       0.5957205
    ## 17:  2015  7560 4202650 female        direct        15.92796       0.8113203
    ## 18:  2015 14827 8182222   <NA>        direct        15.53640       0.6156644
    ## 19:  2016  7267 3979572   male        direct        14.57968       0.5945621
    ## 20:  2016  7560 4202650 female        direct        15.55989       0.8794678
    ## 21:  2016 14827 8182222   <NA>        direct        15.08315       0.6896737
    ## 22:  2017  7267 3979572   male        direct        14.94816       0.6078366
    ## 23:  2017  7560 4202650 female        direct        15.86717       0.7523452
    ## 24:  2017 14827 8182222   <NA>        direct        15.42019       0.6523009
    ##      year     n       N gender estimate_type val_povertyRisk stE_povertyRisk
    ##     <num> <int>   <num> <fctr>        <char>           <num>           <num>

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
    ## 1:  2010 14827 8182222        direct        14.44422       0.6062829
    ## 2:  2011 14827 8182222        direct        14.77393       0.5076536
    ## 3:  2012 14827 8182222        direct        15.04515       0.4806795
    ## 4:  2013 14827 8182222        direct        14.89013       0.4389647
    ## 5:  2014 14827 8182222        direct        15.14556       0.4193012
    ## 6:  2015 14827 8182222        direct        15.53640       0.5987976
    ## 7:  2016 14827 8182222        direct        15.08315       0.4462862
    ## 8:  2017 14827 8182222        direct        15.42019       0.5115568
    ##    val_onePerson stE_onePerson
    ##            <num>         <num>
    ## 1:      14.85737     0.6062829
    ## 2:      14.85737     0.5076536
    ## 3:      14.85737     0.4806795
    ## 4:      14.85737     0.4389647
    ## 5:      14.85737     0.4193012
    ## 6:      14.85737     0.5987976
    ## 7:      14.85737     0.4462862
    ## 8:      14.85737     0.5115568

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
    ##  1:       2.7134079
    ##  2:       2.0227058
    ##  3:       2.0110792
    ##  4:       1.9021917
    ##  5:       1.6437091
    ##  6:       1.4364691
    ##  7:       1.5386759
    ##  8:       1.5650570
    ##  9:       1.0629040
    ## 10:       0.7763632

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
    ##  1:       2.7134079
    ##  2:       2.0227058
    ##  3:       2.0110792
    ##  4:       1.9021917
    ##  5:       1.6437091
    ##  6:       1.4364691
    ##  7:       1.5386759
    ##  8:       1.5650570
    ##  9:       1.0629040
    ## 10:       0.7377296
    ## 11:       0.9086065
    ## 12:       0.7763632

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
    ##  1:       2.3891776
    ##  2:       3.1874147
    ##  3:       2.3250528
    ##  4:       1.9360523
    ##  5:       1.8686773
    ##  6:       2.3761722
    ##  7:       2.0631231
    ##  8:       2.1470231
    ##  9:       1.7976068
    ## 10:       2.0485364
    ## 11:       1.3844393
    ## 12:       1.7976311
    ## 13:       1.8868485
    ## 14:       1.6801057
    ## 15:       0.9273187
    ## 16:       2.0394249
    ## 17:       1.2115344
    ## 18:       1.2783744
    ## 19:       0.7763632

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
    ##     <num> <int>     <num> <fctr>        <fctr>        <char>           <num>
    ##     stE_povertyRisk
    ##               <num>
    ##  1:       2.3891776
    ##  2:       3.1874147
    ##  3:       2.3250528
    ##  4:       1.9360523
    ##  5:       1.8686773
    ##  6:       2.3761722
    ##  7:       2.0631231
    ##  8:       2.7134079
    ##  9:       2.1470231
    ## 10:       1.7976068
    ## 11:       2.0485364
    ## 12:       2.0227058
    ## 13:       2.0110792
    ## 14:       1.9021917
    ## 15:       1.3844393
    ## 16:       1.7976311
    ## 17:       1.8868485
    ## 18:       1.6801057
    ## 19:       1.6437091
    ## 20:       0.9273187
    ## 21:       2.0394249
    ## 22:       1.2115344
    ## 23:       1.2783744
    ## 24:       1.4364691
    ## 25:       1.5386759
    ## 26:       1.5650570
    ## 27:       1.0629040
    ## 28:       0.7377296
    ## 29:       0.9086065
    ## 30:       0.7763632
    ##     stE_povertyRisk
    ##               <num>

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
    ##     <num>   <num>     <num>        <fctr>                        <fctr>
    ##        estimate_type val_povertyRisk stE_povertyRisk
    ##               <char>           <num>           <num>
    ##  1:           direct     19.53983651       2.7134079
    ##  2: group difference      3.00252634       3.4390815
    ##  3:           direct     16.53731017       2.0227058
    ##  4: group difference      5.75249330       3.5736779
    ##  5: group difference      6.45356876       2.1797188
    ##  6: group difference     -2.74996696       3.0666142
    ##  7: group difference     -3.45104242       2.3819563
    ##  8:           direct     13.78734321       2.0110792
    ##  9: group difference      4.23164602       3.0873390
    ## 10: group difference     -0.70107546       2.5538239
    ## 11: group difference     -1.22911968       3.0625120
    ## 12:           direct     13.08626775       1.9021917
    ## 13: group difference     -1.52084728       2.3681476
    ## 14: group difference     -2.22192274       2.5106483
    ## 15:           direct     15.30819049       1.6437091
    ## 16: group difference      5.16519923       2.4570250
    ## 17: group difference      2.30515330       3.5841203
    ## 18: group difference     -2.16267289       2.0553855
    ## 19: group difference      0.69737304       2.7067472
    ## 20: group difference     -0.58729407       1.9984498
    ## 21: group difference     -3.44734000       1.1935650
    ## 22: group difference      5.69621369       3.0410389
    ## 23: group difference      8.65006312       2.3070680
    ## 24: group difference     -1.28836953       1.6542625
    ## 25: group difference     -4.14841546       2.3869868
    ## 26: group difference     -2.69368735       2.6926491
    ## 27: group difference     -5.64753678       2.7756289
    ## 28: group difference     -0.93355321       2.3499792
    ## 29: group difference     -1.92649272       1.8773306
    ## 30: group difference      0.05627961       2.1712574
    ## 31: group difference      2.89756982       2.0900060
    ## 32: group difference     -0.75735506       2.6364478
    ## 33: group difference      2.19649436       2.1039008
    ## 34: group difference     -1.46456768       2.7974845
    ## 35: group difference      4.41841710       1.5999008
    ## 36:           direct     14.37463728       1.4364691
    ## 37: group difference     -2.86004593       1.8253605
    ## 38:           direct     17.23468321       1.5386759
    ## 39: group difference     -0.53101447       1.4466614
    ## 40: group difference      3.48486389       1.7418955
    ## 41: group difference     -3.39106040       1.8921525
    ## 42: group difference     -6.34490982       1.9569850
    ## 43:           direct     13.84362281       1.5650570
    ## 44: group difference      2.95384943       1.8647130
    ## 45:           direct     10.88977339       1.0629040
    ## 46:           direct     12.02659998       0.7377296
    ## 47: group difference     -4.70690810       0.6013164
    ## 48:           direct     16.73350808       0.9086065
    ## 49:           direct     14.44421817       0.7763632
    ##        estimate_type val_povertyRisk stE_povertyRisk
    ##               <char>           <num>           <num>

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
    ## 1:      2014 14827 8182222            direct      15.1455601       0.4499388
    ## 2:      2015 14827 8182222            direct      15.5364014       0.4542393
    ## 3: 2015-2014 14827 8182222 period difference       0.3908413       0.2964377
    ## 4:      2016 14827 8182222            direct      15.0831502       0.4510617
    ## 5: 2016-2015 14827 8182222 period difference      -0.4532512       0.4690850
    ## 6:      2017 14827 8182222            direct      15.4201916       0.4195721
    ## 7: 2017-2016 14827 8182222 period difference       0.3370414       0.2325975

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
    ##        <char> <num>   <num> <fctr>            <char>           <num>
    ##     stE_povertyRisk
    ##               <num>
    ##  1:       0.6477092
    ##  2:       0.4554481
    ##  3:       0.4499388
    ##  4:       0.6076495
    ##  5:       0.5404325
    ##  6:       0.4542393
    ##  7:       0.4432807
    ##  8:       0.2778951
    ##  9:       0.2964377
    ## 10:       0.5296535
    ## 11:       0.5305984
    ## 12:       0.4510617
    ## 13:       0.5771921
    ## 14:       0.4136583
    ## 15:       0.4690850
    ## 16:       0.4528624
    ## 17:       0.4730819
    ## 18:       0.4195721
    ## 19:       0.1974976
    ## 20:       0.3768628
    ## 21:       0.2325975
    ##     stE_povertyRisk
    ##               <num>

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
    ## 1:           2014 14827 8182222         direct        15.14556       0.4499388
    ## 2: 2014_2015_2016 14827 8182222 period average        15.25504       0.4007557
    ## 3:           2015 14827 8182222         direct        15.53640       0.4542393
    ## 4: 2015_2016_2017 14827 8182222 period average        15.34658       0.3672687
    ## 5:           2016 14827 8182222         direct        15.08315       0.4510617
    ## 6:           2017 14827 8182222         direct        15.42019       0.4195721

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
    ##             <char> <num>   <num> <fctr>                             <char>
    ##     val_povertyRisk stE_povertyRisk
    ##               <num>           <num>
    ##  1:     14.50350682       0.6477092
    ##  2:     15.75353283       0.4554481
    ##  3:     15.14556006       0.4499388
    ##  4:     14.73535987       0.5181154
    ##  5:     15.74712982       0.4775339
    ##  6:     15.25503720       0.4007557
    ##  7:     15.12289042       0.6076495
    ##  8:     15.92796296       0.5404325
    ##  9:     15.53640136       0.4542393
    ## 10:     14.88357729       0.4607284
    ## 11:     15.78500836       0.4313383
    ## 12:     15.34658105       0.3672687
    ## 13:     14.57968239       0.5296535
    ## 14:     15.55989368       0.5305984
    ## 15:     15.08315018       0.4510617
    ## 16:     -0.54320803       0.5771921
    ## 17:     -0.36806928       0.4136583
    ## 18:     -0.45325118       0.4690850
    ## 19:      0.14821741       0.1479679
    ## 20:      0.03787854       0.1428562
    ## 21:      0.09154385       0.1157052
    ## 22:     14.94815906       0.4528624
    ## 23:     15.86716845       0.4730819
    ## 24:     15.42019160       0.4195721
    ##     val_povertyRisk stE_povertyRisk
    ##               <num>           <num>
