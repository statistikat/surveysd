# Calcualte point estimates and their standard errors using bootstrap weights.

Calculate point estimates as well as standard errors of variables in
surveys. Standard errors are estimated using bootstrap weights (see
[draw.bootstrap](https://statistikat.github.io/surveysd/reference/draw.bootstrap.md)
and
[recalib](https://statistikat.github.io/surveysd/reference/recalib.md)).
In addition the standard error of an estimate can be calcualted using
the survey data for 3 or more consecutive periods, which results in a
reduction of the standard error.

## Usage

``` r
calc.stError(
  dat,
  weights = attr(dat, "weights"),
  b.weights = attr(dat, "b.rep"),
  period = attr(dat, "period"),
  var = NULL,
  fun = weightedRatio,
  relative.share = FALSE,
  group = NULL,
  group.diff = FALSE,
  fun.adjust.var = NULL,
  adjust.var = NULL,
  period.diff = NULL,
  period.mean = NULL,
  bias = FALSE,
  size.limit = 20,
  cv.limit = 10,
  p = NULL,
  add.arg = NULL,
  national = FALSE
)
```

## Arguments

- dat:

  either data.frame or data.table containing the survey data. Surveys
  can be a panel survey or rotating panel survey, but does not need to
  be. For rotating panel survey bootstrap weights can be created using
  [draw.bootstrap](https://statistikat.github.io/surveysd/reference/draw.bootstrap.md)
  and
  [recalib](https://statistikat.github.io/surveysd/reference/recalib.md).

- weights:

  character specifying the name of the column in `dat` containing the
  original sample weights. Used to calculate point estimates.

- b.weights:

  character vector specifying the names of the columns in `dat`
  containing bootstrap weights. Used to calculate standard errors.

- period:

  character specifying the name of the column in `dat` containing the
  sample periods.

- var:

  character vector containing variable names in `dat` on which `fun`
  shall be applied for each sample period. If `var = NULL` the results
  will reflect the sum of `weights`.

- fun:

  function which will be applied on `var` for each sample period.
  Predefined functions are
  [weightedRatio](https://statistikat.github.io/surveysd/reference/PointEstimates.md),
  [weightedSum](https://statistikat.github.io/surveysd/reference/PointEstimates.md),
  but can also take any other function which returns a double or integer
  and uses weights as its second argument.

- relative.share:

  boolean, if `TRUE` point estimates resulting from `fun` will be
  divided by the point estimate at population level per `period`.

- group:

  character vectors or list of character vectors containig variables in
  `dat`. For each list entry `dat` will be split in subgroups according
  to the containing variables as well as `period`. The pointestimates
  are then estimated for each subgroup seperately. If `group=NULL` the
  data will split into sample periods by default.

- group.diff:

  boolen, if `TRUE` differences and the standard error between groups
  defined in `group` are calculated. See details for more explanations.

- fun.adjust.var:

  can be either `NULL` or a function. This argument can be used to apply
  a function for each `period` and bootstrap weight to the data. The
  resulting estimates will be passed down to `fun`. See details for more
  explanations.

- adjust.var:

  can be either `NULL` or a character specifying the first argument in
  `fun.adjust.var`.

- period.diff:

  character vectors, defining periods for which the differences in the
  point estimate as well it's standard error is calculated. Each entry
  must have the form of `"period1 - period2"`. Can be NULL

- period.mean:

  odd integer, defining the range of periods over which the sample mean
  of point estimates is additionally calcualted.

- bias:

  boolean, if `TRUE` the sample mean over the point estimates of the
  bootstrap weights is returned.

- size.limit:

  integer defining a lower bound on the number of observations on `dat`
  in each group defined by `period` and the entries in `group`. Warnings
  are returned if the number of observations in a subgroup falls below
  `size.limit`. In addition the concerned groups are available in the
  function output.

- cv.limit:

  non-negativ value defining a upper bound for the standard error in
  relation to the point estimate. If this relation exceed `cv.limit`,
  for a point estimate, they are flagged and available in the function
  output.

- p:

  numeric vector containing values between 0 and 1. Defines which
  quantiles for the distribution of `var` are additionally estimated.

- add.arg:

  additional arguments which will be passed to fun. Can be either a
  named list or vector. The names of the object correspond to the
  function arguments and the values to column names in dat, see also
  examples.

- national:

  DEPRECATED use `relative.share` instead! boolean, if TRUE point
  estimates resulting from fun will be divided by the point estimate at
  the national level.

## Value

Returns a list containing:

- `Estimates`: data.table containing period differences and/or k period
  averages for estimates of `fun` applied to `var` as well as the
  corresponding standard errors, which are calculated using the
  bootstrap weights. In addition the sample size, `n`, and poplutaion
  size for each group is added to the output.

- `smallGroups`: data.table containing groups for which the number of
  observation falls below `size.limit`.

- `cvHigh`: data.table containing a boolean variable which indicates for
  each estimate if the estimated standard error exceeds `cv.limit`.

- `stEDecrease`: data.table indicating for each estimate the theoretical
  increase in sample size which is gained when averaging over k periods.
  Only returned if `period.mean` is not `NULL`.

## Details

`calc.stError` takes survey data (`dat`) and returns point estimates as
well as their standard Errors defined by `fun` and `var` for each sample
period in `dat`. `dat` must be household data where household members
correspond to multiple rows with the same household identifier. The data
should at least contain the following columns:

- Column indicating the sample period;

- Column indicating the household ID;

- Column containing the household sample weights;

- Columns which contain the bootstrap weights (see output of
  [recalib](https://statistikat.github.io/surveysd/reference/recalib.md));

- Columns listed in `var` as well as in `group`

For each variable in `var` as well as sample period the function `fun`
is applied using the original as well as the bootstrap sample weights.  
The point estimate is then selected as the result of `fun` when using
the original sample weights and it's standard error is estimated with
the result of `fun` using the bootstrap sample weights.  
  
`fun` can be any function which returns a double or integer and uses
sample weights as it's second argument. The predifined options are
`weightedRatio` and `weightedSum`.  
  
For the option `weightedRatio` a weighted ratio (in \\ calculated for
`var` equal to 1, e.g
`sum(weight[var==1])/sum(weight[!is.na(var)])*100`.  
Additionally using the option `national=TRUE` the weighted ratio (in \\
divided by the weighted ratio at the national level for each `period`.  
If `group` is not `NULL` but a vector of variables from `dat` then `fun`
is applied on each subset of `dat` defined by all combinations of values
in `group`.  
For instance if `group = "sex"` with "sex" having the values "Male" and
"Female" in `dat` the point estimate and standard error is calculated on
the subsets of `dat` with only "Male" or "Female" value for "sex". This
is done for each value of `period`. For variables in `group` which have
`NA`s in `dat` the rows containing the missings will be discarded.  
When `group` is a list of character vectors, subsets of `dat` and the
following estimation of the point estimate, including the estimate for
the standard error, are calculated for each list entry.  
  
If `group.diff = TRUE` difference between groups definded by `group` are
calculated. Differences are only calculated within each variables of
`group`, e.g `group = c("gender", "region")` will calcualate estimates
of each group and also differences within `"gender"` and `"region"`
seperately. If grouping is done with multiple variables e.g
`group = list(c("gender","region")`) (~ grouping by `"gender"` x
`"region"`) differences are calculated only between groups where one of
the grouping variables is different. For instance the difference between
`gender = "female" & region = "Vienna"` and
`gender = "male" & region = "Vienna"` OR
`gender = "female" & region = "Vienna"` and
`gender = "female" & region = "Salzburg"` will be calculated. The
difference between `gender = "female" & region = "Vienna"` and
`gender = "male" & region = "Salzburg"` will not be calculated. The
order of difference is determined by order of value (alpha-numerical
order) or if grouping contains factor variables the factor levels
determin the order.  
The optional parameters `fun.adjust.var` and `adjust.var` can be used if
the values in `var` are dependent on the `weights`. As is for instance
the case for the poverty thershhold calculated from EU-SILC. In such a
case an additional function can be supplied using `fun.adjust.var` as
well as its first argument `adjust.var`, which needs to be part of the
data set `dat`. Then, before applying `fun` on variable `var` for all
`period` and groups, the function `fun.adjust.var` is applied to
`adjust.var` using each of the bootstrap weights seperately (NOTE:
weight is used as the second argument of `fun.adjust.var`). Thus
creating i=1,...,`length(b.weights)` additional variables. For applying
`fun` on `var` the estimates for the bootstrap replicate will now use
each of the corresponding new additional variables. So instead of
\$\$fun(var,weights,...),fun(var,b.weights\[1\],...),
fun(var,b.weights\[2\],...),...\$\$ the function `fun` will be applied
in the way
\$\$fun(var,weights,...),fun(var.1,b.weights\[1\],...),fun(var.2,
b.weights\[2\],...),...\$\$

where `var.1`, `var.2`, `...` correspond to the estimates resulting from
`fun.adjust.var` and `adjust.var`. NOTE: This procedure is especially
usefull if the `var` is dependent on `weights` and `fun` is applied on
subgroups of the data set. Then it is not possible to capture this
procedure with `fun` and `var`, see examples for a more hands on
explanation.  
When defining `period.diff` the difference of point estimates between
periods as well their standard errors are calculated.  
The entries in `period.diff` must have the form of `"period1 - period2"`
which means that the results of the point estimates for `period2` will
be substracted from the results of the point estimates for `period1`.  
  
Specifying `period.mean` leads to an improvement in standard error by
averaging the results for the point estimates, using the bootstrap
weights, over `period.mean` periods. Setting, for instance,
`period.mean = 3` the results in averaging these results over each
consecutive set of 3 periods.  
Estimating the standard error over these averages gives an improved
estimate of the standard error for the central period, which was used
for averaging.  
The averaging of the results is also applied in differences of point
estimates. For instance defining `period.diff = "2015-2009"` and
`period.mean = 3` the differences in point estimates of 2015 and 2009,
2016 and 2010 as well as 2014 and 2008 are calcualated and finally the
average over these 3 differences is calculated. The periods set in
`period.diff` are always used as the middle periods around which the
mean over `period.mean` years is build.  
Setting `bias` to `TRUE` returns the calculation of a mean over the
results from the bootstrap replicates. In the output the corresponding
columns is labeled *\_mean* at the end.  
  
If `fun` needs more arguments they can be supplied in `add.arg`. This
can either be a named list or vector.  
  
The parameter `size.limit` indicates a lower bound of the sample size
for subsets in `dat` created by `group`. If the sample size of a subset
falls below `size.limit` a warning will be displayed.  
In addition all subsets for which this is the case can be selected from
the output of `calc.stError` with `$smallGroups`.  
With the parameter `cv.limit` one can set an upper bound on the
coefficient of variantion. Estimates which exceed this bound are flagged
with `TRUE` and are available in the function output with `$cvHigh`.
`cv.limit` must be a positive integer and is treated internally as \\
for `cv.limit=1` the estimate will be flagged if the coefficient of
variantion exceeds 1\\  
When specifying `period.mean`, the decrease in standard error for
choosing this method is internally calcualted and a rough estimate for
an implied increase in sample size is available in the output with
`$stEDecrease`. The rough estimate for the increase in sample size uses
the fact that for a sample of size \\n\\ the sample estimate for the
standard error of most point estimates converges with a factor
\\1/\sqrt{n}\\ against the true standard error \\\sigma\\.

## See also

[draw.bootstrap](https://statistikat.github.io/surveysd/reference/draw.bootstrap.md)  
[recalib](https://statistikat.github.io/surveysd/reference/recalib.md)

## Author

Johannes Gussenbauer, Alexander Kowarik, Statistics Austria

## Examples

``` r
# Import data and calibrate

library(surveysd)
library(data.table)
setDTthreads(1)
set.seed(1234)
eusilc <- demo.eusilc(n = 4,prettyNames = TRUE)
dat_boot <- draw.bootstrap(eusilc, REP = 3, hid = "hid", weights = "pWeight",
                           strata = "region", period = "year")
dat_boot_calib <- recalib(dat_boot, conP.var = "gender", conH.var = "region")
#> Iteration stopped after 2 steps
#> Convergence reached
#> Iteration stopped after 1 steps
#> Convergence reached
#> Iteration stopped after 4 steps
#> Convergence reached

# estimate weightedRatio for povertyRisk per period

err.est <- calc.stError(dat_boot_calib, var = "povertyRisk",
                        fun = weightedRatio)
err.est$Estimates
#> Key: <year, n, N, estimate_type>
#>     year     n       N estimate_type val_povertyRisk stE_povertyRisk
#>    <num> <int>   <num>        <char>           <num>           <num>
#> 1:  2010 14827 8182222        direct        14.44422       0.5270298
#> 2:  2011 14827 8182222        direct        14.77393       0.3942901
#> 3:  2012 14827 8182222        direct        15.04515       0.4114000
#> 4:  2013 14827 8182222        direct        14.89013       0.3463014

# calculate weightedRatio for povertyRisk and fraction of one-person
# households per period

dat_boot_calib[, onePerson := .N == 1, by = .(year, hid)]
#>          hid  hsize        region     pid       age gender   ecoStat
#>        <int> <fctr>        <fctr>   <int>    <fctr> <fctr>    <fctr>
#>     1:     1      3         Tyrol     101   (25,45] female part time
#>     2:     1      3         Tyrol     102   (25,45]   male full time
#>     3:     1      3         Tyrol     103 (-Inf,16]   male      <NA>
#>     4:     2      4         Tyrol     201   (25,45] female  domestic
#>     5:     2      4         Tyrol     202   (25,45]   male full time
#>    ---                                                              
#> 59304:  8999      4 Lower Austria  899904 (-Inf,16] female education
#> 59305: 10500      1 Upper Austria 1050001   (25,45] female full time
#> 59306:  5999      1         Tyrol  599901   (25,45]   male full time
#> 59307:  7500      2         Tyrol  750001   (45,65]   male full time
#> 59308:  7500      2         Tyrol  750002   (45,65] female  disabled
#>        citizenship   py010n py050n  py090n py100n py110n py120n py130n py140n
#>             <fctr>    <num>  <num>   <num>  <num>  <num>  <num>  <num>  <num>
#>     1:          AT  9756.25      0    0.00      0      0      0      0      0
#>     2:       Other 12471.60      0    0.00      0      0      0      0      0
#>     3:        <NA>       NA     NA      NA     NA     NA     NA     NA     NA
#>     4:          AT 12487.03      0    0.00      0      0      0      0      0
#>     5:          AT 42821.23      0    0.00      0      0      0      0      0
#>    ---                                                                       
#> 59304:          AT     0.00      0    0.00      0      0      0      0      0
#> 59305:          AT 13962.56      0    0.00      0      0      0      0      0
#> 59306:          AT 14685.18      0    0.00      0      0      0      0      0
#> 59307:          AT 20606.82      0    0.00      0      0      0      0      0
#> 59308:          AT     0.00      0 3825.63      0      0      0      0      0
#>        hy040n  hy050n hy070n hy080n hy090n hy110n hy130n hy145n  eqSS  eqIncome
#>         <num>   <num>  <num>  <num>  <num>  <num>  <num>  <num> <num>     <num>
#>     1: 4273.9 2428.11      0      0  33.39      0      0      0   1.8 16090.694
#>     2: 4273.9 2428.11      0      0  33.39      0      0      0   1.8 16090.694
#>     3: 4273.9 2428.11      0      0  33.39      0      0      0   1.8 16090.694
#>     4:    0.0 1549.72      0      0   2.13      0      0      0   2.1 27076.243
#>     5:    0.0 1549.72      0      0   2.13      0      0      0   2.1 27076.243
#>    ---                                                                         
#> 59304:    0.0 1955.19      0      0   0.00      0      0      0   2.5 26970.023
#> 59305:    0.0    0.00      0      0 424.85      0      0      0   1.0  6923.625
#> 59306:    0.0    0.00      0      0 120.65      0      0      0   1.0 14805.830
#> 59307:    0.0    0.00      0      0   0.00      0      0      0   1.5 24680.877
#> 59308:    0.0    0.00      0      0   0.00      0      0      0   1.5 24680.877
#>           db090  pWeight  year povertyRisk           w1           w2
#>           <num>    <num> <num>      <lgcl>        <num>        <num>
#>     1: 504.5696 504.5696  2010       FALSE    0.4470694    0.4486785
#>     2: 504.5696 504.5696  2010       FALSE    0.4470694    0.4486785
#>     3: 504.5696 504.5696  2010       FALSE    0.4470694    0.4486785
#>     4: 493.3824 493.3824  2010       FALSE  982.8594451  986.3259754
#>     5: 493.3824 493.3824  2010       FALSE  982.8594451  986.3259754
#>    ---                                                              
#> 59304: 556.4260 556.4260  2013       FALSE 1107.0247723 1112.8575590
#> 59305: 643.2557 643.2557  2013        TRUE 1289.8469103    0.6060913
#> 59306: 679.7288 679.7288  2013       FALSE    0.5962466    0.6044353
#> 59307: 567.1544 567.1544  2013       FALSE 1127.8618193 1133.8043941
#> 59308: 567.1544 567.1544  2013       FALSE 1127.8618193 1133.8043941
#>                  w3 onePerson
#>               <num>    <lgcl>
#>     1:    0.4547508     FALSE
#>     2:    0.4547508     FALSE
#>     3:    0.4547508     FALSE
#>     4:    0.4465770     FALSE
#>     5:    0.4465770     FALSE
#>    ---                       
#> 59304: 1123.0762927     FALSE
#> 59305:    0.5885571      TRUE
#> 59306:    0.6336001      TRUE
#> 59307: 1144.2154706     FALSE
#> 59308: 1144.2154706     FALSE
err.est <- calc.stError(dat_boot_calib, var = c("povertyRisk", "onePerson"),
                        fun = weightedRatio)
err.est$Estimates
#> Key: <year, n, N, estimate_type>
#>     year     n       N estimate_type val_povertyRisk stE_povertyRisk
#>    <num> <int>   <num>        <char>           <num>           <num>
#> 1:  2010 14827 8182222        direct        14.44422       0.4935789
#> 2:  2011 14827 8182222        direct        14.77393       0.3697396
#> 3:  2012 14827 8182222        direct        15.04515       0.3370720
#> 4:  2013 14827 8182222        direct        14.89013       0.3024828
#>    val_onePerson stE_onePerson
#>            <num>         <num>
#> 1:      14.85737     0.4935789
#> 2:      14.85737     0.3697396
#> 3:      14.85737     0.3370720
#> 4:      14.85737     0.3024828


# estimate weightedRatio for povertyRisk per period and gender and
# period x region x gender 

group <- list("gender", c("gender", "region"))
err.est <- calc.stError(dat_boot_calib, var = "povertyRisk",
                        fun = weightedRatio, group = group)
err.est$Estimates
#> Key: <year, n, N, gender, region, estimate_type>
#>      year     n         N gender        region estimate_type val_povertyRisk
#>     <num> <int>     <num> <fctr>        <fctr>        <char>           <num>
#>  1:  2010   261  122741.8   male    Burgenland        direct       17.414524
#>  2:  2010   288  137822.2 female    Burgenland        direct       21.432598
#>  3:  2010   359  182732.9   male    Vorarlberg        direct       12.973259
#>  4:  2010   374  194622.1 female    Vorarlberg        direct       19.883637
#>  5:  2010   440  253143.7   male      Salzburg        direct        9.156964
#>  6:  2010   484  282307.3 female      Salzburg        direct       17.939382
#>  7:  2010   517  268581.4   male     Carinthia        direct       10.552149
#>  8:  2010   561  295066.6 female     Carinthia        direct       15.392924
#>  9:  2010   650  339566.5   male         Tyrol        direct       12.857542
#> 10:  2010   667  362332.5 female         Tyrol        direct       17.604861
#> 11:  2010  1128  571011.7   male        Styria        direct       11.671247
#> 12:  2010  1132  774405.4   male        Vienna        direct       15.590616
#> 13:  2010  1167  596033.3 female        Styria        direct       16.964539
#> 14:  2010  1190  824525.6 female        Vienna        direct       18.778813
#> 15:  2010  1363  684272.5   male Upper Austria        direct        9.074690
#> 16:  2010  1387  772593.2 female Lower Austria        direct       16.372949
#> 17:  2010  1417  783115.8   male Lower Austria        direct       11.348283
#> 18:  2010  1442  737347.5 female Upper Austria        direct       12.574206
#> 19:  2010  7267 3979571.7   male          <NA>        direct       12.026600
#> 20:  2010  7560 4202650.3 female          <NA>        direct       16.733508
#> 21:  2010 14827 8182222.0   <NA>          <NA>        direct       14.444218
#> 22:  2011   261  122741.8   male    Burgenland        direct       16.465679
#> 23:  2011   288  137822.2 female    Burgenland        direct       21.245704
#> 24:  2011   359  182732.9   male    Vorarlberg        direct       12.590109
#> 25:  2011   374  194622.1 female    Vorarlberg        direct       17.125259
#> 26:  2011   440  253143.7   male      Salzburg        direct       11.189644
#> 27:  2011   484  282307.3 female      Salzburg        direct       18.213248
#> 28:  2011   517  268581.4   male     Carinthia        direct       12.789630
#> 29:  2011   561  295066.6 female     Carinthia        direct       16.745958
#> 30:  2011   650  339566.5   male         Tyrol        direct       12.531181
#> 31:  2011   667  362332.5 female         Tyrol        direct       16.522466
#> 32:  2011  1128  571011.7   male        Styria        direct       10.163611
#> 33:  2011  1132  774405.4   male        Vienna        direct       16.616484
#> 34:  2011  1167  596033.3 female        Styria        direct       15.322461
#> 35:  2011  1190  824525.6 female        Vienna        direct       18.582643
#> 36:  2011  1363  684272.5   male Upper Austria        direct       11.360200
#> 37:  2011  1387  772593.2 female Lower Austria        direct       16.109340
#> 38:  2011  1417  783115.8   male Lower Austria        direct       12.419104
#> 39:  2011  1442  737347.5 female Upper Austria        direct       14.426605
#> 40:  2011  7267 3979571.7   male          <NA>        direct       12.819212
#> 41:  2011  7560 4202650.3 female          <NA>        direct       16.624882
#> 42:  2011 14827 8182222.0   <NA>          <NA>        direct       14.773925
#> 43:  2012   261  122741.8   male    Burgenland        direct       16.005251
#> 44:  2012   288  137822.2 female    Burgenland        direct       19.646025
#> 45:  2012   359  182732.9   male    Vorarlberg        direct       12.650134
#> 46:  2012   374  194622.1 female    Vorarlberg        direct       16.089942
#> 47:  2012   440  253143.7   male      Salzburg        direct       11.305057
#> 48:  2012   484  282307.3 female      Salzburg        direct       18.150559
#> 49:  2012   517  268581.4   male     Carinthia        direct       14.157852
#> 50:  2012   561  295066.6 female     Carinthia        direct       15.926182
#> 51:  2012   650  339566.5   male         Tyrol        direct       15.104689
#> 52:  2012   667  362332.5 female         Tyrol        direct       17.449297
#> 53:  2012  1128  571011.7   male        Styria        direct       11.019424
#> 54:  2012  1132  774405.4   male        Vienna        direct       16.099733
#> 55:  2012  1167  596033.3 female        Styria        direct       13.902894
#> 56:  2012  1190  824525.6 female        Vienna        direct       16.796126
#> 57:  2012  1363  684272.5   male Upper Austria        direct       13.557275
#> 58:  2012  1387  772593.2 female Lower Austria        direct       16.563011
#> 59:  2012  1417  783115.8   male Lower Austria        direct       13.606128
#> 60:  2012  1442  737347.5 female Upper Austria        direct       15.494037
#> 61:  2012  7267 3979571.7   male          <NA>        direct       13.760647
#> 62:  2012  7560 4202650.3 female          <NA>        direct       16.261469
#> 63:  2012 14827 8182222.0   <NA>          <NA>        direct       15.045149
#> 64:  2013   261  122741.8   male    Burgenland        direct       20.636575
#> 65:  2013   288  137822.2 female    Burgenland        direct       21.223485
#> 66:  2013   359  182732.9   male    Vorarlberg        direct       12.311667
#> 67:  2013   374  194622.1 female    Vorarlberg        direct       15.716327
#> 68:  2013   440  253143.7   male      Salzburg        direct       11.427862
#> 69:  2013   484  282307.3 female      Salzburg        direct       17.234804
#> 70:  2013   517  268581.4   male     Carinthia        direct       12.594461
#> 71:  2013   561  295066.6 female     Carinthia        direct       13.368498
#> 72:  2013   650  339566.5   male         Tyrol        direct       14.564544
#> 73:  2013   667  362332.5 female         Tyrol        direct       16.559763
#> 74:  2013  1128  571011.7   male        Styria        direct       11.942973
#> 75:  2013  1132  774405.4   male        Vienna        direct       16.144697
#> 76:  2013  1167  596033.3 female        Styria        direct       14.377224
#> 77:  2013  1190  824525.6 female        Vienna        direct       17.419665
#> 78:  2013  1363  684272.5   male Upper Austria        direct       13.461762
#> 79:  2013  1387  772593.2 female Lower Austria        direct       14.978617
#> 80:  2013  1417  783115.8   male Lower Austria        direct       13.710926
#> 81:  2013  1442  737347.5 female Upper Austria        direct       15.272207
#> 82:  2013  7267 3979571.7   male          <NA>        direct       13.889623
#> 83:  2013  7560 4202650.3 female          <NA>        direct       15.837536
#> 84:  2013 14827 8182222.0   <NA>          <NA>        direct       14.890134
#>      year     n         N gender        region estimate_type val_povertyRisk
#>     stE_povertyRisk
#>               <num>
#>  1:      3.91582094
#>  2:      4.61118944
#>  3:      1.97092591
#>  4:      3.28637313
#>  5:      2.92829528
#>  6:      3.81125582
#>  7:      2.05173509
#>  8:      0.52378948
#>  9:      1.51156359
#> 10:      0.73178455
#> 11:      1.18582696
#> 12:      1.91045605
#> 13:      1.16785450
#> 14:      1.87783396
#> 15:      1.49386669
#> 16:      1.27805725
#> 17:      0.30524904
#> 18:      2.14297092
#> 19:      0.60534991
#> 20:      0.52757574
#> 21:      0.52702985
#> 22:      2.46601320
#> 23:      3.58191387
#> 24:      0.52690363
#> 25:      0.74189709
#> 26:      3.35317800
#> 27:      3.37077372
#> 28:      1.23382184
#> 29:      0.63584451
#> 30:      0.64657212
#> 31:      1.49536068
#> 32:      0.61875626
#> 33:      0.42434709
#> 34:      0.82878872
#> 35:      0.44278814
#> 36:      1.30875690
#> 37:      1.93272371
#> 38:      0.78188947
#> 39:      1.74422460
#> 40:      0.18133964
#> 41:      0.61597086
#> 42:      0.39429009
#> 43:      3.30564720
#> 44:      4.62630849
#> 45:      2.51209290
#> 46:      1.02468931
#> 47:      2.43407212
#> 48:      2.31629871
#> 49:      1.34654132
#> 50:      2.55905021
#> 51:      1.28095658
#> 52:      1.18312497
#> 53:      0.75816641
#> 54:      0.84608517
#> 55:      0.16500396
#> 56:      1.07938417
#> 57:      1.70355764
#> 58:      1.56487271
#> 59:      1.39075566
#> 60:      2.16016623
#> 61:      0.15255888
#> 62:      0.68130895
#> 63:      0.41140002
#> 64:      2.28114079
#> 65:      0.76807268
#> 66:      0.67249468
#> 67:      1.41269402
#> 68:      1.67223967
#> 69:      2.90634765
#> 70:      0.72893737
#> 71:      2.84323993
#> 72:      1.02899397
#> 73:      1.74306998
#> 74:      0.86664536
#> 75:      0.12330786
#> 76:      1.20865559
#> 77:      0.91592787
#> 78:      0.31718765
#> 79:      1.22834999
#> 80:      1.13432978
#> 81:      1.25790729
#> 82:      0.08370471
#> 83:      0.59567303
#> 84:      0.34630135
#>     stE_povertyRisk

# use average over 3 periods for standard error estimation
# and calculate estimate for difference of
# period 2011 and 2012 inclulding standard errors
period.diff <- c("2012-2011")
err.est <- calc.stError(
  dat_boot_calib, var = "povertyRisk", fun = weightedRatio,
  period.diff = period.diff,  # <- take difference of periods 2012 and 2011
  period.mean = 3)  # <- average over 3 periods
err.est$Estimates
#> Key: <year, n, N, estimate_type>
#>              year     n       N                      estimate_type
#>            <char> <num>   <num>                             <char>
#> 1:           2010 14827 8182222                             direct
#> 2: 2010_2011_2012 14827 8182222                     period average
#> 3:           2011 14827 8182222                             direct
#> 4: 2011_2012_2013 14827 8182222                     period average
#> 5:           2012 14827 8182222                             direct
#> 6:      2012-2011 14827 8182222                  period difference
#> 7: 2012-2011_mean 14827 8182222 difference between period averages
#> 8:           2013 14827 8182222                             direct
#>    val_povertyRisk stE_povertyRisk
#>              <num>           <num>
#> 1:      14.4442182      0.52702985
#> 2:      14.7544308      0.44413839
#> 3:      14.7739255      0.39429009
#> 4:      14.9030692      0.38279602
#> 5:      15.0451487      0.41140002
#> 6:       0.2712233      0.02777690
#> 7:       0.1486385      0.06431972
#> 8:      14.8901335      0.34630135

# for more examples see https://statistikat.github.io/surveysd/articles/error_estimation.html
```
