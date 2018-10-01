[![Build Status](https://travis-ci.org/statistikat/surveysd.svg?branch=master)](https://travis-ci.org/statistikat/surveysd)
<!--[![Coverage Status](https://coveralls.io/repos/github/statistikat/surveysd/badge.svg?branch=master)](https://coveralls.io/github/statistikat/surveysd?branch=master)-->
<!--[![CRAN](http://www.r-pkg.org/badges/version/surveysd)](https://CRAN.R-project.org/package=surveysd)-->
<!--[![Downloads](http://cranlogs.r-pkg.org/badges/surveysd)](https://CRAN.R-project.org/package=surveysd)-->
<!--[![Mentioned in Awesome Official Statistics ](https://awesome.re/mentioned-badge.svg)](http://www.awesomeofficialstatistics.org)-->

## surveysd
This is the development place for R-package surveysd

Estimate standard deviation of estimates in complex surveys using bootstrap weights.


#### Load dummy data
```{r}
library(surveysd)
library(laeken)
library(data.table)

eusilc <- surveysd:::demo.eusilc()
``` 

#### Draw bootstrap replicates

```{r}
dat_boot <- draw.bootstrap(eusilc, REP = 10, hid = "db030", weights = "rb050", 
                           strata = c("db040"), period = "year")
```
   
#### Calibrate bootstrap replicates
```{r}
dat_boot_calib <- recalib(copy(dat_boot), conP.var = c("rb090"), conH.var = c("db040"))
```

#### Estimate weighted ratio for variable `povmd60` per period and `rb090`

```{r}
err.est <- calc.stError(dat_boot_calib, var = "povmd60", fun = weightedRatio,
                        group = "rb090", period.diff = NULL, period.mean = NULL)
```

#### Estimate weighted ratio for povmd60 per period and `rb090`, `db040` and combination of both

```{r}
group <- list("rb090", "db040", c("rb090","db040"))
err.est <- calc.stError(dat_boot_calib, var = "povmd60", fun = weightedRatio,
                        group = group, period.diff = NULL, period.mean = NULL)
```
