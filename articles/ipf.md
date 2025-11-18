# Iterative Proportional Fitting

This vignette explains the usage of the
[`ipf()`](https://statistikat.github.io/surveysd/reference/ipf.md)
function, which has been used for calibrating the labour force survey of
Austria for several years. It is based on the Iterative Proportional
Fitting algorithm and gives some flexibility about the details of the
implementation. See (Meraner, Gumprecht, and Kowarik 2016) or
[`vignette("methodology")`](https://statistikat.github.io/surveysd/articles/methodology.md)
for more details.

## Setup

We will assume the output of
[`demo.eusilc()`](https://statistikat.github.io/surveysd/reference/demo.eusilc.md)
is our population. From this population, a sample without replacement is
drawn. The sample covers 10 percent of the population. We assign a
weight of one for all observations of the population and a weight of ten
for all observations of the sample.

``` r
library(surveysd)
population <- demo.eusilc(1, prettyNames = TRUE)
population[, pWeight := 1]
pop_sample <- population[sample(1:.N, floor(.N*0.10)), ]
pop_sample[, pWeight := 10]
```

## One constraint, one variable

We will start with an example where we want to adapt the weights of
`pop_sample` such that the weighted number of males and females matches
the ones of `population`. We can see that this is currently not the
case.

``` r
(gender_distribution <- xtabs(pWeight ~ gender, population))
#> gender
#>   male female 
#>   7267   7560
xtabs(pWeight ~ gender, pop_sample)
#> gender
#>   male female 
#>   7380   7440
```

Due to random sampling (rather than stratified sampling), there are
differences between the gender distributions. We can pass
`gender_distribution` as a parameter to
[`ipf()`](https://statistikat.github.io/surveysd/reference/ipf.md) to
obtain modified weights.

``` r
pop_sample_c <- ipf(pop_sample, conP = list(gender_distribution), w = "pWeight")
```

The resulting dataset, `pop_sample_c` is similar to `pop_sample` but has
an additional column with the adjusted weights.

``` r
dim(pop_sample)
#> [1] 1482   30
dim(pop_sample_c)
#> [1] 1482   31
setdiff(names(pop_sample_c), names(pop_sample))
#> [1] "calibWeight"
```

We can now calculate the weighted number of males and females according
to this new weight. This will result in a match for the constraints.

``` r
xtabs(calibWeight ~ gender, pop_sample_c)
#> gender
#>   male female 
#>   7267   7560
xtabs(pWeight ~ gender, population)
#> gender
#>   male female 
#>   7267   7560
```

In this simple case, `ipf` just performs a post stratification step.
This means, that all males and all females have the same weight.

``` r
xtabs(~ calibWeight + gender, pop_sample_c)
#>                   gender
#> calibWeight        male female
#>   9.84688346883469  738      0
#>   10.1612903225806    0    744
```

All males have been weighted down (`calibWeight < 10`) to compensate for
the overrepresentation in the sample.

## One constraint, two variables

Let’s now assume that we want to put constraints on the number of males
and females for each age group. The numbers from the original population
can be obtained with [`xtabs()`](https://rdrr.io/r/stats/xtabs.html).

``` r
(con_ga <- xtabs(pWeight ~ gender + age, population))
#>         age
#> gender   (-Inf,16] (16,25] (25,45] (45,65] (65, Inf]
#>   male        1528     855    2165    1822       897
#>   female      1375     848    2255    1845      1237
xtabs(pWeight ~ gender + age, pop_sample)
#>         age
#> gender   (-Inf,16] (16,25] (25,45] (45,65] (65, Inf]
#>   male        1480     870    2310    1830       890
#>   female      1620     780    2180    1690      1170
```

Again, we can see that those constraints are not met. Supplying the
contingency table `con_ga` to
[`ipf()`](https://statistikat.github.io/surveysd/reference/ipf.md) will
again resolve this.

``` r
pop_sample_c2 <- ipf(pop_sample, conP = list(con_ga), w = "pWeight")
xtabs(pWeight ~ gender + age, population)
#>         age
#> gender   (-Inf,16] (16,25] (25,45] (45,65] (65, Inf]
#>   male        1528     855    2165    1822       897
#>   female      1375     848    2255    1845      1237
xtabs(calibWeight ~ gender + age, pop_sample_c2)
#>         age
#> gender   (-Inf,16] (16,25] (25,45] (45,65] (65, Inf]
#>   male        1528     855    2165    1822       897
#>   female      1375     848    2255    1845      1237
```

## Two constraints

Now we assume that we know the number of persons living in each nuts2
region from registry data.

``` r
registry_table <- xtabs(pWeight ~ region, population)
```

However, those registry data does not provide any information about age
or `gender`. Therefore, the two contingency tables (`con_ga` and
`registry_table`) have to be specified independently. This can be done
by supplying a list to `conP`

``` r
pop_sample_c2 <- ipf(pop_sample, conP = list(con_ga, registry_table), w = "pWeight")
xtabs(pWeight ~ gender + age, population)
#>         age
#> gender   (-Inf,16] (16,25] (25,45] (45,65] (65, Inf]
#>   male        1528     855    2165    1822       897
#>   female      1375     848    2255    1845      1237
xtabs(calibWeight ~ gender + age, pop_sample_c2)
#>         age
#> gender   (-Inf,16]   (16,25]   (25,45]   (45,65] (65, Inf]
#>   male   1528.0003  854.9994 2164.9996 1821.9998  897.0004
#>   female 1374.9999  847.9999 2255.0012 1845.0003 1236.9992
xtabs(pWeight ~ region, population)
#> region
#>    Burgenland     Carinthia Lower Austria      Salzburg        Styria 
#>           549          1078          2804           924          2295 
#>         Tyrol Upper Austria        Vienna    Vorarlberg 
#>          1317          2805          2322           733
xtabs(calibWeight ~ region, pop_sample_c2)
#> region
#>    Burgenland     Carinthia Lower Austria      Salzburg        Styria 
#>           549          1078          2804           924          2295 
#>         Tyrol Upper Austria        Vienna    Vorarlberg 
#>          1317          2805          2322           733
```

this time, the constraints are not matched perfectly. That is, because
we provided more than one constraint. therefore, the
[`ipf()`](https://statistikat.github.io/surveysd/reference/ipf.md)
algorithm had to work iteratively.

## Household Constraints

If the dataset has a household structure, household constraints can be
passed via the parameter `conH`. If this parameter is used, it is also
necessary to supply `hid`, which defines the column names that contains
household ids.

``` r
(conH1 <- xtabs(pWeight ~ hsize + region, data = population[!duplicated(hid)]))
#>      region
#> hsize Burgenland Carinthia Lower Austria Salzburg Styria Tyrol Upper Austria
#>     1         58       117           325      103    264   118           262
#>     2         82       126           345      102    260   149           321
#>     3         37        80           189       55    187    79           203
#>     4         33        63           169       71    122   102           168
#>     5         16        39           103       30     83    48           114
#>      region
#> hsize Vienna Vorarlberg
#>     1    431         67
#>     2    355         72
#>     3    175         44
#>     4     96         53
#>     5     50         34
pop_sample_hh <- ipf(pop_sample, hid = "hid", conH = list(conH1), w = "pWeight",
                     bound = 10)
xtabs(calibWeight ~ hsize + region, data = pop_sample_hh[!duplicated(hid)])
#>      region
#> hsize Burgenland Carinthia Lower Austria Salzburg Styria Tyrol Upper Austria
#>     1         58       117           325      103    264   118           262
#>     2         82       126           345      102    260   149           321
#>     3         37        80           189       55    187    79           203
#>     4         33        63           169       71    122   102           168
#>     5         16        39           103       30     83    48           114
#>      region
#> hsize Vienna Vorarlberg
#>     1    431         67
#>     2    355         72
#>     3    175         44
#>     4     96         53
#>     5     50         34
```

## Tolerances

If `conP` or `conH` contain several contingency tables or if `conP` and
`conH` are used at the same time, the ipf algorithm will operate
iteratively. This means that the calibrated dataset will satisfy the
constraints only approximately. The default tolerances of the
approximation can be overwritten using the parameters `conP` and `conH`.

Lowering the tolerances will improve the match between the constraints
and the contingency tables according to the calibrated weights. However,
lower tolerances will also make it so more iterations are necessary
until a convergence is met. If the constraints are too small, ipf will
return with a warning that indicates that a convergence could not be
reached.

``` r
ipf(pop_sample, conP = list(con_ga, registry_table), w = "pWeight",
    verbose = TRUE, epsP = 0.01)
#> Iteration stopped after 3 steps
#> Convergence reached
ipf(pop_sample, conP = list(con_ga, registry_table), w = "pWeight",
    verbose = TRUE, epsP = 0.0001)
#> Iteration stopped after 3 steps
#> Convergence reached
```

We see that changing the tolerances from `0.01` (one percent) to
`0.0001` increases the number of required iterations.

## References

Meraner, Angelika, Daniela Gumprecht, and Alexander Kowarik. 2016.
“Weighting Procedure of the Austrian Microcensus Using Administrative
Data.” *Austrian Journal of Statistics* 45 (June): 3.
<https://doi.org/10.17713/ajs.v45i3.120>.
