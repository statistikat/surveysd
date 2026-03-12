# Plot surveysd-Objects

Plot results of
[`calc.stError()`](https://statistikat.github.io/surveysd/reference/calc.stError.md)

## Usage

``` r
# S3 method for class 'surveysd'
plot(
  x,
  variable = x$param$var[1],
  type = c("summary", "grouping"),
  groups = NULL,
  sd.type = c("dot", "ribbon"),
  ...
)
```

## Arguments

- x:

  object of class 'surveysd' output of function
  [calc.stError](https://statistikat.github.io/surveysd/reference/calc.stError.md)

- variable:

  Name of the variable for which standard errors have been calcualated
  in `dat`

- type:

  can bei either `"summary"` or `"grouping"`, default value is
  `"summary"`. For `"summary"` a barplot is created giving an overview
  of the number of estimates having the flag `smallGroup`, `cvHigh`,
  both or none of them. For 'grouping' results for point estimate and
  standard error are plotted for pre defined groups.

- groups:

  If `type='grouping'` variables must be defined by which the data is
  grouped. Only 2 levels are supported as of right now. If only one
  group is defined the higher group will be the estimate over the whole
  period. Results are plotted for the first argument in `groups` as well
  as for the combination of `groups[1]` and `groups[2]`.

- sd.type:

  can bei either `'ribbon'` or `'dot'` and is only used if
  `type='grouping'`. Default is `"dot"` For `sd.type='dot'` point
  estimates are plotted and flagged if the corresponding standard error
  and/or the standard error using the mean over k-periods exceeded the
  value `cv.limit` (see
  [calc.stError](https://statistikat.github.io/surveysd/reference/calc.stError.md)).
  For `sd.type='ribbon'` the point estimates including ribbons, defined
  by point estimate +- estimated standard error are plotted. The
  calculated standard errors using the mean over k periods are plotted
  using less transparency. Results for the higher level (~`groups[1]`)
  are coloured grey.

- ...:

  additional arguments supplied to plot.

## Examples

``` r
library(surveysd)

set.seed(1234)
eusilc <- demo.eusilc(n = 3, prettyNames = TRUE)

dat_boot <- draw.bootstrap(eusilc, REP = 3, hid = "hid", weights = "pWeight",
                           strata = "region", period = "year")

# calibrate weight for bootstrap replicates
dat_boot_calib <- recalib(dat_boot, conP.var = "gender", conH.var = "region")
#> Iteration stopped after 1 steps
#> Convergence reached
#> Iteration stopped after 1 steps
#> Convergence reached
#> 10:Not yet converged for P-Constraint1
#> 10:Not yet converged for H-Constraint1
#> 20:Not yet converged for P-Constraint1
#> 20:Not yet converged for H-Constraint1
#> 30:Not yet converged for P-Constraint1
#> 30:Not yet converged for H-Constraint1
#> 40:Not yet converged for P-Constraint1
#> 40:Not yet converged for H-Constraint1
#> 50:Not yet converged for P-Constraint1
#> 50:Not yet converged for H-Constraint1
#> 60:Not yet converged for P-Constraint1
#> 60:Not yet converged for H-Constraint1
#> 70:Not yet converged for P-Constraint1
#> 70:Not yet converged for H-Constraint1
#> 80:Not yet converged for P-Constraint1
#> 80:Not yet converged for H-Constraint1
#> 90:Not yet converged for P-Constraint1
#> 90:Not yet converged for H-Constraint1
#> 100:Not yet converged for P-Constraint1
#>      year gender     maxFac     N  epsP CalibMargin PopMargin
#>    <fctr> <fctr>      <num> <int> <num>       <num>     <num>
#> 1:   2012   male 0.02215357  7267  0.01     4067733   3979572
#> 2:   2012 female 0.02215357  7560  0.01     4295754   4202650
#> 3:   2011 female 0.01689261  7560  0.01     4273644   4202650
#> 4:   2011   male 0.01689261  7267  0.01     4046797   3979572
#> 5:   2010 female 0.01340012  7560  0.01     4258966   4202650
#> 6:   2010   male 0.01340012  7267  0.01     4032898   3979572
#> -----------------------------------------
#> 100:Not yet converged for H-Constraint1
#>      year        region     maxFac     N  epsH sumCalibWeight PopMargin
#>    <fctr>        <fctr>      <num> <int> <num>          <num>     <num>
#> 1:   2012      Salzburg 0.02167343   924  0.02       214917.8    219679
#> 2:   2012    Vorarlberg 0.02167343   733  0.02       141852.5    144995
#> 3:   2012     Carinthia 0.02167343  1078  0.02       228679.9    233746
#> 4:   2012         Tyrol 0.02167343  1317  0.02       272969.7    279017
#> 5:   2012    Burgenland 0.02167343   549  0.02       107465.3    109846
#> 6:   2012        Styria 0.02167343  2295  0.02       479738.1    490366
#> 7:   2012 Lower Austria 0.02167343  2804  0.02       633330.5    647361
#> 8:   2012 Upper Austria 0.02167343  2805  0.02       554721.9    567011
#> 9:   2012        Vienna 0.02167343  2322  0.02       795500.8    813124
#> -----------------------------------------
#> 110:Not yet converged for P-Constraint1
#> 110:Not yet converged for H-Constraint1
#> 120:Not yet converged for P-Constraint1
#> 120:Not yet converged for H-Constraint1
#> 130:Not yet converged for P-Constraint1
#> 130:Not yet converged for H-Constraint1
#> 140:Not yet converged for P-Constraint1
#> 140:Not yet converged for H-Constraint1
#> 150:Not yet converged for P-Constraint1
#> 150:Not yet converged for H-Constraint1
#> 160:Not yet converged for P-Constraint1
#> 160:Not yet converged for H-Constraint1
#> 170:Not yet converged for P-Constraint1
#> 170:Not yet converged for H-Constraint1
#> 180:Not yet converged for P-Constraint1
#> 180:Not yet converged for H-Constraint1
#> 190:Not yet converged for P-Constraint1
#> 190:Not yet converged for H-Constraint1
#> 200:Not yet converged for P-Constraint1
#>      year gender     maxFac     N  epsP CalibMargin PopMargin
#>    <fctr> <fctr>      <num> <int> <num>       <num>     <num>
#> 1:   2012   male 0.02215357  7267  0.01     4067733   3979572
#> 2:   2012 female 0.02215357  7560  0.01     4295754   4202650
#> 3:   2011   male 0.01689261  7267  0.01     4046797   3979572
#> 4:   2011 female 0.01689261  7560  0.01     4273644   4202650
#> 5:   2010   male 0.01340012  7267  0.01     4032898   3979572
#> 6:   2010 female 0.01340012  7560  0.01     4258966   4202650
#> -----------------------------------------
#> 200:Not yet converged for H-Constraint1
#>      year        region     maxFac     N  epsH sumCalibWeight PopMargin
#>    <fctr>        <fctr>      <num> <int> <num>          <num>     <num>
#> 1:   2012        Vienna 0.02167343  2322  0.02       795500.8    813124
#> 2:   2012         Tyrol 0.02167343  1317  0.02       272969.7    279017
#> 3:   2012    Burgenland 0.02167343   549  0.02       107465.3    109846
#> 4:   2012 Upper Austria 0.02167343  2805  0.02       554721.9    567011
#> 5:   2012        Styria 0.02167343  2295  0.02       479738.1    490366
#> 6:   2012    Vorarlberg 0.02167343   733  0.02       141852.5    144995
#> 7:   2012      Salzburg 0.02167343   924  0.02       214917.8    219679
#> 8:   2012     Carinthia 0.02167343  1078  0.02       228679.9    233746
#> 9:   2012 Lower Austria 0.02167343  2804  0.02       633330.5    647361
#> -----------------------------------------
#> Warning: Not converged in 200 steps
#> No convergence reached
#> Calibration failed for bootstrap replicates w3 
#> Corresponding bootstrap replicates will be discarded
#> Returning 2 calibrated bootstrap weights

# estimate weightedRatio for povmd60 per period
group <- list("gender", "region", c("gender", "region"))
err.est <- calc.stError(dat_boot_calib, var = "povertyRisk",
                        fun = weightedRatio,
                        group = group , period.mean = NULL)


plot(err.est)


# plot results for gender
# dotted line is the result on the national level
plot(err.est, type = "grouping", groups = "gender")
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's shape values.



# plot results for rb090 in each db040
# with standard errors as ribbons
plot(err.est, type = "grouping", groups = c("gender", "region"), sd.type = "ribbon")


```
