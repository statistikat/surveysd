# Introduction to surveysd

The goal of surveysd is to combine all necessary steps to use calibrated
bootstrapping with custom estimation functions. This vignette will cover
the usage of the most important functions. For insights in the theory
used in this package, refer to
[`vignette("methodology")`](https://statistikat.github.io/surveysd/articles/methodology.md).

### Load dummy data

A test data set based on `data(eusilc, package = "laeken")` can be
created with
[`demo.eusilc()`](https://statistikat.github.io/surveysd/reference/demo.eusilc.md)

``` r
library(surveysd)

set.seed(1234)
eusilc <- demo.eusilc(n = 2, prettyNames = TRUE)

eusilc[1:5, .(year, povertyRisk, gender, pWeight)]
```

    ##     year povertyRisk gender  pWeight
    ##    <num>      <lgcl> <fctr>    <num>
    ## 1:  2010       FALSE female 504.5696
    ## 2:  2010       FALSE   male 504.5696
    ## 3:  2010       FALSE   male 504.5696
    ## 4:  2010       FALSE female 493.3824
    ## 5:  2010       FALSE   male 493.3824

### Draw bootstrap replicates

Use stratified resampling without replacement to generate 10 samples.
Those samples are consistent with respect to the reference periods.

``` r
dat_boot <- draw.bootstrap(eusilc, REP = 10, hid = "hid", weights = "pWeight", 
                           strata = "region", period = "year")
```

### Calibrate bootstrap replicates

Calibrate each sample according to the distribution of `gender` (on a
personal level) and `region` (on a household level).

``` r
dat_boot_calib <- recalib(dat_boot, conP.var = "gender", conH.var = "region",
                          epsP = 1e-2, epsH = 2.5e-2, verbose = FALSE)
dat_boot_calib[1:5, .(year, povertyRisk, gender, pWeight, w1, w2, w3, w4)]
```

    ##     year povertyRisk gender  pWeight           w1          w2        w3
    ##    <num>      <lgcl> <fctr>    <num>        <num>       <num>     <num>
    ## 1:  2010       FALSE female 504.5696 1012.0710607   0.4450574 0.4502766
    ## 2:  2010       FALSE   male 504.5696 1012.0710607   0.4450574 0.4502766
    ## 3:  2010       FALSE   male 504.5696 1012.0710607   0.4450574 0.4502766
    ## 4:  2010       FALSE female 493.3824    0.4406142 977.7492539 0.4409384
    ## 5:  2010       FALSE   male 493.3824    0.4406142 977.7492539 0.4409384
    ##             w4
    ##          <num>
    ## 1:   0.4442893
    ## 2:   0.4442893
    ## 3:   0.4442893
    ## 4: 977.0283364
    ## 5: 977.0283364

### Estimate with respect to a grouping variable

Estimate relative amount of persons at risk of poverty per period and
`gender`.

``` r
err.est <- calc.stError(dat_boot_calib, var = "povertyRisk", fun = weightedRatio, group = "gender")
err.est$Estimates
```

    ## Key: <year, n, N, gender, estimate_type>
    ##     year     n       N gender estimate_type val_povertyRisk stE_povertyRisk
    ##    <num> <int>   <num> <fctr>        <char>           <num>           <num>
    ## 1:  2010  7267 3979572   male        direct        12.02660       0.5890890
    ## 2:  2010  7560 4202650 female        direct        16.73351       0.7699564
    ## 3:  2010 14827 8182222   <NA>        direct        14.44422       0.6007732
    ## 4:  2011  7267 3979572   male        direct        12.81921       0.5928896
    ## 5:  2011  7560 4202650 female        direct        16.62488       0.6689113
    ## 6:  2011 14827 8182222   <NA>        direct        14.77393       0.5108341

The output contains estimates (`val_povertyRisk`) as well as standard
errors (`stE_povertyRisk`) measured in percent. The rows with
`gender = NA` denotes the aggregate over all genders for the
corresponding year.

### Estimate with respect to several variables

Estimate relative amount of persons at risk of poverty per period for
each `region`, `gender`, and combination of both.

``` r
group <- list("gender", "region", c("gender", "region"))
err.est <- calc.stError(dat_boot_calib, var = "povertyRisk", fun = weightedRatio, group = group)
head(err.est$Estimates)
```

    ## Key: <year, n, N, gender, region, estimate_type>
    ##     year     n        N gender     region estimate_type val_povertyRisk
    ##    <num> <int>    <num> <fctr>     <fctr>        <char>           <num>
    ## 1:  2010   261 122741.8   male Burgenland        direct       17.414524
    ## 2:  2010   288 137822.2 female Burgenland        direct       21.432598
    ## 3:  2010   359 182732.9   male Vorarlberg        direct       12.973259
    ## 4:  2010   374 194622.1 female Vorarlberg        direct       19.883637
    ## 5:  2010   440 253143.7   male   Salzburg        direct        9.156964
    ## 6:  2010   484 282307.3 female   Salzburg        direct       17.939382
    ##    stE_povertyRisk
    ##              <num>
    ## 1:        4.106583
    ## 2:        4.699786
    ## 3:        1.775553
    ## 4:        2.286644
    ## 5:        1.342583
    ## 6:        2.475953

``` r
## skipping 54 more rows
```
