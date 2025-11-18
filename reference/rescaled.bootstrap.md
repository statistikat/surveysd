# Draw bootstrap replicates

Draw bootstrap replicates from survey data using either the rescaled
bootstrap for stratified multistage sampling, presented by J. Preston
(2009) or the Rao-Wu boostrap by J. N. K. Rao and C. F. J. Wu (1988)

## Usage

``` r
rescaled.bootstrap(
  dat,
  method = c("Preston", "Rao-Wu"),
  REP = 1000,
  strata = "DB050>1",
  cluster = "DB060>DB030",
  fpc = "N.cluster>N.households",
  single.PSU = c("merge", "mean"),
  return.value = c("data", "replicates"),
  run.input.checks = TRUE,
  already.selected = NULL,
  seed = NULL
)
```

## Arguments

- dat:

  either data frame or data table containing the survey sample

- method:

  for bootstrap replicates, either `"Preston"` or `"Rao-Wu"`

- REP:

  integer indicating the number of bootstraps to be drawn

- strata:

  string specifying the column name in `dat` that is used for
  stratification. For multistage sampling multiple column names can be
  specified by `strata=c("strata1","strata2","strata3")` or
  `strata=c("strata1>strata2>strata3")`. See Details for more
  information.

- cluster:

  string specifying the column name in `dat` that is used for
  clustering. For instance given a household sample the column
  containing the household ID should be supplied. For multistage
  sampling multiple column names can be specified by
  `cluster=c("cluster1","cluster2","cluster3")` or
  `cluster=c("cluster1>cluster2>cluster3")`. See Details for more
  information.

- fpc:

  string specifying the column name in `dat` that contains the number of
  PSUs at the first stage. For multistage sampling the number of PSUs at
  each stage must be specified by `strata=c("fpc1","fpc2","fpc3")` or
  `strata=c("fpc1>fpc2>fpc3")`.

- single.PSU:

  either "merge" or "mean" defining how single PSUs need to be dealt
  with. For `single.PSU="merge"` single PSUs at each stage are merged
  with the strata or cluster with the next least number of PSUs. If
  multiple of those exist one will be select via random draw. For
  `single.PSU="mean"` single PSUs will get the mean over all bootstrap
  replicates at the stage which did not contain single PSUs.

- return.value:

  either "data", "replicates" and/or "selection" specifying the return
  value of the function. For "data" the survey data is returned as class
  `data.table`, for "replicates" only the bootstrap replicates are
  returned as `data.table`. For "selection" list of data.tables with
  length of `length(strata)` is returned containing 1:`REP` 0-1 columns
  indicating if a PSU was selected for each sampling stage.

- run.input.checks:

  logical, if TRUE the input will be checked before applying the
  bootstrap procedure

- already.selected:

  list of data.tables or `NULL` where each data.table contains columns
  in `cluster`, `strata` and additionally 1:`REP` columns containing 0-1
  values which indicate if a PSU was selected for each bootstrap
  replicate. Each of the data.tables corresponds to one of the sampling
  stages. First entry in the list corresponds to the first sampling
  stage and so on.

- seed:

  integer specifying the seed for the random number generator.

## Value

returns the complete data set including the bootstrap replicates or just
the bootstrap replicates, depending on `return.value="data"` or
`return.value="replicates"` respectively.

## Details

For specifying multistage sampling designs the column names in
`strata`,`cluster` and `fpc` need to be seperated by "\>".  
For multistage sampling the strings are read from left to right meaning
that the first vector entry or column name before the first "\>" is
taken as the column for stratification/clustering/number of PSUs at the
first and the last vector entry or column after the last "\>" is taken
as the column for stratification/clustering/number of PSUs at the last
stage. If for some stages the sample was not stratified or clustered one
must specify this by "1" or "I", e.g.
`strata=c("strata1","I","strata3")` or `strata=c("strata1>I>strata3")`
if there was no stratification at the second stage or
`cluster=c("cluster1","cluster2","I")` respectively
`cluster=c("cluster1>cluster2>I")` if there were no clusters at the last
stage.  
The number of PSUs at each stage is not calculated internally and must
be specified for any sampling design. For single stage sampling using
stratification this can usually be done by adding over all sample
weights of each PSU by each strata-code.  
Spaces in each of the strings will be removed, so if column names
contain spaces they should be renamed before calling this procedure!  
If `already.selected` is supplied the sampling of bootstrap replicates
considers if speficif PSUs have already been selected by a previous
survey wave. For a specific `strata` and `cluster` this could lead to
more than `floor(n/2)` records selected. In that case records will be
de-selected such that `floor(n/2)` records, with `n` as the total number
of records, are selected for each `strata` and `cluster`. This parameter
ist mostly used by
[draw.bootstrap](https://statistikat.github.io/surveysd/reference/draw.bootstrap.md)
in order to consider the rotation of the sampling units over time.

## References

Preston, J. (2009). Rescaled bootstrap for stratified multistage
sampling. Survey Methodology. 35. 227-234. Rao, J. N. K., and C. F. J.
Wu. (1988). Resampling Inference with Complex Survey Data. Journal of
the American Statistical Association 83 (401): 231–41.

## Author

Johannes Gussenbauer, Eileen Vattheuer, Statistics Austria

## Examples

``` r
library(surveysd)
library(data.table)
setDTthreads(1)
set.seed(1234)
eusilc <- demo.eusilc(n = 1,prettyNames = TRUE)

eusilc[,N.households:=uniqueN(hid),by=region]
#>          hid  hsize        region    pid       age gender   ecoStat citizenship
#>        <int> <fctr>        <fctr>  <int>    <fctr> <fctr>    <fctr>      <fctr>
#>     1:     1      3         Tyrol    101   (25,45] female part time          AT
#>     2:     1      3         Tyrol    102   (25,45]   male full time       Other
#>     3:     1      3         Tyrol    103 (-Inf,16]   male      <NA>        <NA>
#>     4:     2      4         Tyrol    201   (25,45] female  domestic          AT
#>     5:     2      4         Tyrol    202   (25,45]   male full time          AT
#>    ---                                                                         
#> 14823:  5997      4 Lower Austria 599704 (-Inf,16] female education          AT
#> 14824:  5998      1 Upper Austria 599801   (25,45] female full time          AT
#> 14825:  5999      1         Tyrol 599901   (25,45]   male full time          AT
#> 14826:  6000      2         Tyrol 600001   (45,65]   male full time          AT
#> 14827:  6000      2         Tyrol 600002   (45,65] female  disabled          AT
#>          py010n py050n  py090n py100n py110n py120n py130n py140n hy040n
#>           <num>  <num>   <num>  <num>  <num>  <num>  <num>  <num>  <num>
#>     1:  9756.25      0    0.00      0      0      0      0      0 4273.9
#>     2: 12471.60      0    0.00      0      0      0      0      0 4273.9
#>     3:       NA     NA      NA     NA     NA     NA     NA     NA 4273.9
#>     4: 12487.03      0    0.00      0      0      0      0      0    0.0
#>     5: 42821.23      0    0.00      0      0      0      0      0    0.0
#>    ---                                                                  
#> 14823:     0.00      0    0.00      0      0      0      0      0    0.0
#> 14824: 13962.56      0    0.00      0      0      0      0      0    0.0
#> 14825: 14685.18      0    0.00      0      0      0      0      0    0.0
#> 14826: 20606.82      0    0.00      0      0      0      0      0    0.0
#> 14827:     0.00      0 3825.63      0      0      0      0      0    0.0
#>         hy050n hy070n hy080n hy090n hy110n hy130n hy145n  eqSS eqIncome
#>          <num>  <num>  <num>  <num>  <num>  <num>  <num> <num>    <num>
#>     1: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     2: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     3: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     4: 1549.72      0      0   2.13      0      0      0   2.1 27076.24
#>     5: 1549.72      0      0   2.13      0      0      0   2.1 27076.24
#>    ---                                                                 
#> 14823: 1955.19      0      0   0.00      0      0      0   2.5 26508.20
#> 14824:    0.00      0      0 424.85      0      0      0   1.0 14387.41
#> 14825:    0.00      0      0 120.65      0      0      0   1.0 14805.83
#> 14826:    0.00      0      0   0.00      0      0      0   1.5 16288.30
#> 14827:    0.00      0      0   0.00      0      0      0   1.5 16288.30
#>           db090  pWeight  year povertyRisk N.households
#>           <num>    <num> <num>      <lgcl>        <int>
#>     1: 504.5696 504.5696  2010       FALSE          496
#>     2: 504.5696 504.5696  2010       FALSE          496
#>     3: 504.5696 504.5696  2010       FALSE          496
#>     4: 493.3824 493.3824  2010       FALSE          496
#>     5: 493.3824 493.3824  2010       FALSE          496
#>    ---                                                 
#> 14823: 556.4260 556.4260  2010       FALSE         1131
#> 14824: 643.2557 643.2557  2010       FALSE         1068
#> 14825: 679.7288 679.7288  2010       FALSE          496
#> 14826: 567.1544 567.1544  2010       FALSE          496
#> 14827: 567.1544 567.1544  2010       FALSE          496
eusilc.bootstrap <- rescaled.bootstrap(eusilc,REP=10,strata="region",
                                       cluster="hid",fpc="N.households")

eusilc[,new_strata:=paste(region,hsize,sep="_")]
#>          hid  hsize        region    pid       age gender   ecoStat citizenship
#>        <int> <fctr>        <fctr>  <int>    <fctr> <fctr>    <fctr>      <fctr>
#>     1:     1      3         Tyrol    101   (25,45] female part time          AT
#>     2:     1      3         Tyrol    102   (25,45]   male full time       Other
#>     3:     1      3         Tyrol    103 (-Inf,16]   male      <NA>        <NA>
#>     4:     2      4         Tyrol    201   (25,45] female  domestic          AT
#>     5:     2      4         Tyrol    202   (25,45]   male full time          AT
#>    ---                                                                         
#> 14823:  5997      4 Lower Austria 599704 (-Inf,16] female education          AT
#> 14824:  5998      1 Upper Austria 599801   (25,45] female full time          AT
#> 14825:  5999      1         Tyrol 599901   (25,45]   male full time          AT
#> 14826:  6000      2         Tyrol 600001   (45,65]   male full time          AT
#> 14827:  6000      2         Tyrol 600002   (45,65] female  disabled          AT
#>          py010n py050n  py090n py100n py110n py120n py130n py140n hy040n
#>           <num>  <num>   <num>  <num>  <num>  <num>  <num>  <num>  <num>
#>     1:  9756.25      0    0.00      0      0      0      0      0 4273.9
#>     2: 12471.60      0    0.00      0      0      0      0      0 4273.9
#>     3:       NA     NA      NA     NA     NA     NA     NA     NA 4273.9
#>     4: 12487.03      0    0.00      0      0      0      0      0    0.0
#>     5: 42821.23      0    0.00      0      0      0      0      0    0.0
#>    ---                                                                  
#> 14823:     0.00      0    0.00      0      0      0      0      0    0.0
#> 14824: 13962.56      0    0.00      0      0      0      0      0    0.0
#> 14825: 14685.18      0    0.00      0      0      0      0      0    0.0
#> 14826: 20606.82      0    0.00      0      0      0      0      0    0.0
#> 14827:     0.00      0 3825.63      0      0      0      0      0    0.0
#>         hy050n hy070n hy080n hy090n hy110n hy130n hy145n  eqSS eqIncome
#>          <num>  <num>  <num>  <num>  <num>  <num>  <num> <num>    <num>
#>     1: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     2: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     3: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     4: 1549.72      0      0   2.13      0      0      0   2.1 27076.24
#>     5: 1549.72      0      0   2.13      0      0      0   2.1 27076.24
#>    ---                                                                 
#> 14823: 1955.19      0      0   0.00      0      0      0   2.5 26508.20
#> 14824:    0.00      0      0 424.85      0      0      0   1.0 14387.41
#> 14825:    0.00      0      0 120.65      0      0      0   1.0 14805.83
#> 14826:    0.00      0      0   0.00      0      0      0   1.5 16288.30
#> 14827:    0.00      0      0   0.00      0      0      0   1.5 16288.30
#>           db090  pWeight  year povertyRisk N.households      new_strata
#>           <num>    <num> <num>      <lgcl>        <int>          <char>
#>     1: 504.5696 504.5696  2010       FALSE          496         Tyrol_3
#>     2: 504.5696 504.5696  2010       FALSE          496         Tyrol_3
#>     3: 504.5696 504.5696  2010       FALSE          496         Tyrol_3
#>     4: 493.3824 493.3824  2010       FALSE          496         Tyrol_4
#>     5: 493.3824 493.3824  2010       FALSE          496         Tyrol_4
#>    ---                                                                 
#> 14823: 556.4260 556.4260  2010       FALSE         1131 Lower Austria_4
#> 14824: 643.2557 643.2557  2010       FALSE         1068 Upper Austria_1
#> 14825: 679.7288 679.7288  2010       FALSE          496         Tyrol_1
#> 14826: 567.1544 567.1544  2010       FALSE          496         Tyrol_2
#> 14827: 567.1544 567.1544  2010       FALSE          496         Tyrol_2
eusilc[,N.housholds:=uniqueN(hid),by=new_strata]
#>          hid  hsize        region    pid       age gender   ecoStat citizenship
#>        <int> <fctr>        <fctr>  <int>    <fctr> <fctr>    <fctr>      <fctr>
#>     1:     1      3         Tyrol    101   (25,45] female part time          AT
#>     2:     1      3         Tyrol    102   (25,45]   male full time       Other
#>     3:     1      3         Tyrol    103 (-Inf,16]   male      <NA>        <NA>
#>     4:     2      4         Tyrol    201   (25,45] female  domestic          AT
#>     5:     2      4         Tyrol    202   (25,45]   male full time          AT
#>    ---                                                                         
#> 14823:  5997      4 Lower Austria 599704 (-Inf,16] female education          AT
#> 14824:  5998      1 Upper Austria 599801   (25,45] female full time          AT
#> 14825:  5999      1         Tyrol 599901   (25,45]   male full time          AT
#> 14826:  6000      2         Tyrol 600001   (45,65]   male full time          AT
#> 14827:  6000      2         Tyrol 600002   (45,65] female  disabled          AT
#>          py010n py050n  py090n py100n py110n py120n py130n py140n hy040n
#>           <num>  <num>   <num>  <num>  <num>  <num>  <num>  <num>  <num>
#>     1:  9756.25      0    0.00      0      0      0      0      0 4273.9
#>     2: 12471.60      0    0.00      0      0      0      0      0 4273.9
#>     3:       NA     NA      NA     NA     NA     NA     NA     NA 4273.9
#>     4: 12487.03      0    0.00      0      0      0      0      0    0.0
#>     5: 42821.23      0    0.00      0      0      0      0      0    0.0
#>    ---                                                                  
#> 14823:     0.00      0    0.00      0      0      0      0      0    0.0
#> 14824: 13962.56      0    0.00      0      0      0      0      0    0.0
#> 14825: 14685.18      0    0.00      0      0      0      0      0    0.0
#> 14826: 20606.82      0    0.00      0      0      0      0      0    0.0
#> 14827:     0.00      0 3825.63      0      0      0      0      0    0.0
#>         hy050n hy070n hy080n hy090n hy110n hy130n hy145n  eqSS eqIncome
#>          <num>  <num>  <num>  <num>  <num>  <num>  <num> <num>    <num>
#>     1: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     2: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     3: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     4: 1549.72      0      0   2.13      0      0      0   2.1 27076.24
#>     5: 1549.72      0      0   2.13      0      0      0   2.1 27076.24
#>    ---                                                                 
#> 14823: 1955.19      0      0   0.00      0      0      0   2.5 26508.20
#> 14824:    0.00      0      0 424.85      0      0      0   1.0 14387.41
#> 14825:    0.00      0      0 120.65      0      0      0   1.0 14805.83
#> 14826:    0.00      0      0   0.00      0      0      0   1.5 16288.30
#> 14827:    0.00      0      0   0.00      0      0      0   1.5 16288.30
#>           db090  pWeight  year povertyRisk N.households      new_strata
#>           <num>    <num> <num>      <lgcl>        <int>          <char>
#>     1: 504.5696 504.5696  2010       FALSE          496         Tyrol_3
#>     2: 504.5696 504.5696  2010       FALSE          496         Tyrol_3
#>     3: 504.5696 504.5696  2010       FALSE          496         Tyrol_3
#>     4: 493.3824 493.3824  2010       FALSE          496         Tyrol_4
#>     5: 493.3824 493.3824  2010       FALSE          496         Tyrol_4
#>    ---                                                                 
#> 14823: 556.4260 556.4260  2010       FALSE         1131 Lower Austria_4
#> 14824: 643.2557 643.2557  2010       FALSE         1068 Upper Austria_1
#> 14825: 679.7288 679.7288  2010       FALSE          496         Tyrol_1
#> 14826: 567.1544 567.1544  2010       FALSE          496         Tyrol_2
#> 14827: 567.1544 567.1544  2010       FALSE          496         Tyrol_2
#>        N.housholds
#>              <int>
#>     1:          79
#>     2:          79
#>     3:          79
#>     4:         102
#>     5:         102
#>    ---            
#> 14823:         169
#> 14824:         262
#> 14825:         118
#> 14826:         149
#> 14827:         149
eusilc.bootstrap <- rescaled.bootstrap(eusilc,REP=10,strata=c("new_strata"),
                                       cluster="hid",fpc="N.households")

```
