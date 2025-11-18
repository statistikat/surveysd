# Draw bootstrap replicates

Draw bootstrap replicates from survey data with rotating panel design.
Survey information, like ID, sample weights, strata and population
totals per strata, should be specified to ensure meaningfull survey
bootstraping.

## Usage

``` r
draw.bootstrap(
  dat,
  method = "Preston",
  REP = 1000,
  hid = NULL,
  weights,
  period = NULL,
  strata = NULL,
  cluster = NULL,
  totals = NULL,
  single.PSU = c("merge", "mean"),
  boot.names = NULL,
  split = FALSE,
  pid = NULL,
  seed = NULL,
  already.selected = NULL
)
```

## Arguments

- dat:

  either data.frame or data.table containing the survey data with
  rotating panel design.

- method:

  for bootstrap replicates, either `"Preston"` or `"Rao-Wu"`

- REP:

  integer indicating the number of bootstrap replicates.

- hid:

  character specifying the name of the column in `dat` containing the
  household id. If `NULL` (the default), the household structure is not
  regarded.

- weights:

  character specifying the name of the column in `dat` containing the
  sample weights.

- period:

  character specifying the name of the column in `dat` containing the
  sample periods. If `NULL` (the default), it is assumed that all
  observations belong to the same period.

- strata:

  character vector specifying the name(s) of the column in `dat` by
  which the population was stratified. If `strata` is a vector
  stratification will be assumed as the combination of column names
  contained in `strata`. Setting in addition `cluster` not NULL
  stratification will be assumed on multiple stages, where each
  additional entry in `strata` specifies the stratification variable for
  the next lower stage. see Details for more information.

- cluster:

  character vector specifying cluster in the data. If not already
  specified in `cluster` household ID is taken es the lowest level
  cluster.

- totals:

  character specifying the name of the column in `dat` containing the
  the totals per strata and/or cluster. Is ONLY optional if `cluster` is
  `NULL` or equal `hid` and `strata` contains one columnname! Then the
  households per strata will be calcualted using the `weights` argument.
  If clusters and strata for multiple stages are specified `totals`
  needs to be a vector of `length(strata)` specifying the column on
  `dat` that contain the total number of PSUs at each stage. `totals` is
  interpreted from left the right, meaning that the first argument
  corresponds to the number of PSUs at the first and the last argument
  to the number of PSUs at the last stage.

- single.PSU:

  either "merge" or "mean" defining how single PSUs need to be dealt
  with. For `single.PSU="merge"` single PSUs at each stage are merged
  with the strata or cluster with the next least number of PSUs. If
  multiple of those exist one will be select via random draw. For
  `single.PSU="mean"` single PSUs will get the mean over all bootstrap
  replicates at the stage which did not contain single PSUs.

- boot.names:

  character indicating the leading string of the column names for each
  bootstrap replica. If NULL defaults to "w".

- split:

  logical, if TRUE split households are considered using `pid`, for more
  information see Details.

- pid:

  column in `dat` specifying the personal identifier. This identifier
  needs to be unique for each person throught the whole data set. Can be
  `NULL`.

- seed:

  integer specifying the seed for the random number generator.

- already.selected:

  list of data.tables indicating if record was drawn in previous
  `period`. `length(already.selected)` must be equal to the number of
  sampling stages specified. See
  [`get.selection()`](https://statistikat.github.io/surveysd/reference/get.selection.md)
  for an example.

## Value

the survey data with the number of REP bootstrap replicates added as
columns.

Returns a data.table containing the original data as well as the number
of `REP` columns containing the bootstrap replicates for each
repetition.  
The columns of the bootstrap replicates are by default labeled
"w*Number*" where *Number* goes from 1 to `REP`. If the column names of
the bootstrap replicates should start with a different character or
string the parameter `boot.names` can be used.

## Details

`draw.bootstrap` takes `dat` and draws `REP` bootstrap replicates from
it. `dat` must be household data where household members correspond to
multiple rows with the same household identifier. For most practical
applications, the following columns should be available in the dataset
and passed via the corresponding parameters:

- Column indicating the sample period (parameter `period`).

- Column indicating the household ID (parameter `hid`)

- Column containing the household sample weights (parameter `weights`);

- Columns by which population was stratified during the sampling process
  (parameter: `strata`).

As methods either the either the rescaled bootstrap for stratified
multistage sampling, presented by J. Preston (2009)
(`method = "Preston"`) or the Rao-Wu boostrap by J. N. K. Rao and C. F.
J. Wu (1988) (`method = "Rao-Wu"`) are supported.  
For single stage sampling design a column the argument `totals` is
optional, meaning that a column of the number of PSUs at the first stage
does not need to be supplied. The number of PSUs is calculated and added
to `dat` using `strata` and `weights`. By setting `cluster` to NULL
single stage sampling design is always assumed and if `strata` contains
of multiple column names the combination of all those column names will
be used for stratification.  
In the case of multi stage sampling design the argument `totals` needs
to be specified and needs to have the same number of arguments as
`strata`.  
If `cluster` is `NULL` or does not contain `hid` at the last stage,
`hid` will automatically be used as the final cluster. If, besides
`hid`, clustering in additional stages is specified the number of column
names in `strata` and `cluster` (including `hid`) must be the same. If
for any stage there was no clustering or stratification one can set "1"
or "I" for this stage.  
For example `strata=c("REGION","I"),cluster=c("MUNICIPALITY","HID")`
would speficy a 2 stage sampling design where at the first stage the
municipalities where drawn stratified by regions and at the 2nd stage
housholds are drawn in each municipality without stratification.  
Bootstrap replicates are drawn for each survey period consecutively
(`period`) using the function
[rescaled.bootstrap](https://statistikat.github.io/surveysd/reference/rescaled.bootstrap.md).
Bootstrap replicates are drawn consistently in the way that in each
`period` and sampling stage always \\\lfloor n/2 \rfloor\\ clusters are
selected in each strata.  
This ensures that the bootstrap replicates follow the same logic as the
sampled households, making the bootstrap replicates more comparable to
the actual sample units.  
If `split` ist set to `TRUE` and `pid` is specified, the bootstrap
replicates are carried forward using the personal identifiers instead of
the household identifier. This takes into account the issue of a
houshold splitting up. Any person in this new split household will get
the same bootstrap replicate as the person that has come from an other
household in the survey. People who enter already existing households
will also get the same bootstrap replicate as the other households
members had in the previous periods.  
`already.selected` can be specified in order to construct bootstrap
replicates considering that already drawn bootstrap replicates from a
previous `period` exist for the same survey. See
[`get.selection()`](https://statistikat.github.io/surveysd/reference/get.selection.md)
for more explanations and examples.

## References

Preston, J. (2009). Rescaled bootstrap for stratified multistage
sampling. Survey Methodology. 35. 227-234. Rao, J. N. K., and C. F. J.
Wu. (1988). Resampling Inference with Complex Survey Data. Journal of
the American Statistical Association 83 (401): 231–41.

## See also

[`data.table`](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
for more information on data.table objects.

## Author

Johannes Gussenbauer, Alexander Kowarik, Statistics Austria

## Examples

``` r
library(surveysd)
library(data.table)
setDTthreads(1)
set.seed(1234)
eusilc <- demo.eusilc(n = 3, prettyNames = TRUE)

## draw replicates without stratification or clustering
dat_boot <- draw.bootstrap(eusilc, REP = 1, weights = "pWeight",
                           period = "year")

## use stratification w.r.t. region and clustering w.r.t. households
dat_boot <- draw.bootstrap(
  eusilc, REP = 1, hid = "hid", weights = "pWeight",
  strata = "region", period = "year")

## use multi-level clustering
dat_boot <- draw.bootstrap(
  eusilc, REP = 1, hid = "hid", weights = "pWeight",
  strata = c("region", "hsize"), period = "year")


# create spit households
eusilc[, pidsplit := pid]
#>          hid  hsize        region    pid       age gender   ecoStat citizenship
#>        <int> <fctr>        <fctr>  <int>    <fctr> <fctr>    <fctr>      <fctr>
#>     1:     1      3         Tyrol    101   (25,45] female part time          AT
#>     2:     1      3         Tyrol    102   (25,45]   male full time       Other
#>     3:     1      3         Tyrol    103 (-Inf,16]   male      <NA>        <NA>
#>     4:     2      4         Tyrol    201   (25,45] female  domestic          AT
#>     5:     2      4         Tyrol    202   (25,45]   male full time          AT
#>    ---                                                                         
#> 44477:  8999      4 Lower Austria 899904 (-Inf,16] female education          AT
#> 44478:  9000      1 Upper Austria 900001   (25,45] female full time          AT
#> 44479:  5999      1         Tyrol 599901   (25,45]   male full time          AT
#> 44480:  7500      2         Tyrol 750001   (45,65]   male full time          AT
#> 44481:  7500      2         Tyrol 750002   (45,65] female  disabled          AT
#>          py010n py050n  py090n py100n py110n py120n py130n py140n hy040n
#>           <num>  <num>   <num>  <num>  <num>  <num>  <num>  <num>  <num>
#>     1:  9756.25      0    0.00      0      0      0      0      0 4273.9
#>     2: 12471.60      0    0.00      0      0      0      0      0 4273.9
#>     3:       NA     NA      NA     NA     NA     NA     NA     NA 4273.9
#>     4: 12487.03      0    0.00      0      0      0      0      0    0.0
#>     5: 42821.23      0    0.00      0      0      0      0      0    0.0
#>    ---                                                                  
#> 44477:     0.00      0    0.00      0      0      0      0      0    0.0
#> 44478: 13962.56      0    0.00      0      0      0      0      0    0.0
#> 44479: 14685.18      0    0.00      0      0      0      0      0    0.0
#> 44480: 20606.82      0    0.00      0      0      0      0      0    0.0
#> 44481:     0.00      0 3825.63      0      0      0      0      0    0.0
#>         hy050n hy070n hy080n hy090n hy110n hy130n hy145n  eqSS eqIncome
#>          <num>  <num>  <num>  <num>  <num>  <num>  <num> <num>    <num>
#>     1: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     2: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     3: 2428.11      0      0  33.39      0      0      0   1.8 16090.69
#>     4: 1549.72      0      0   2.13      0      0      0   2.1 27076.24
#>     5: 1549.72      0      0   2.13      0      0      0   2.1 27076.24
#>    ---                                                                 
#> 44477: 1955.19      0      0   0.00      0      0      0   2.5 26970.02
#> 44478:    0.00      0      0 424.85      0      0      0   1.0 20242.09
#> 44479:    0.00      0      0 120.65      0      0      0   1.0 14805.83
#> 44480:    0.00      0      0   0.00      0      0      0   1.5 24680.88
#> 44481:    0.00      0      0   0.00      0      0      0   1.5 24680.88
#>           db090  pWeight  year povertyRisk pidsplit
#>           <num>    <num> <num>      <lgcl>    <int>
#>     1: 504.5696 504.5696  2010       FALSE      101
#>     2: 504.5696 504.5696  2010       FALSE      102
#>     3: 504.5696 504.5696  2010       FALSE      103
#>     4: 493.3824 493.3824  2010       FALSE      201
#>     5: 493.3824 493.3824  2010       FALSE      202
#>    ---                                             
#> 44477: 556.4260 556.4260  2012       FALSE   899904
#> 44478: 643.2557 643.2557  2012       FALSE   900001
#> 44479: 679.7288 679.7288  2012       FALSE   599901
#> 44480: 567.1544 567.1544  2012       FALSE   750001
#> 44481: 567.1544 567.1544  2012       FALSE   750002
year <- eusilc[, unique(year)]
year <- year[-1]
leaf_out <- c()
for(y in year) {
  split.person <- eusilc[
    year == (y-1) & !duplicated(hid) & !(hid %in% leaf_out),
    sample(pid, 20)
  ]
  overwrite.person <- eusilc[
    (year == (y)) & !duplicated(hid) & !(hid %in% leaf_out),
    .(pid = sample(pid, 20))
  ]
  overwrite.person[, c("pidsplit", "year_curr") := .(split.person, y)]

  eusilc[overwrite.person, pidsplit := i.pidsplit,
         on = .(pid, year >= year_curr)]
  leaf_out <- c(leaf_out,
                eusilc[pid %in% c(overwrite.person$pid,
                                  overwrite.person$pidsplit),
                unique(hid)])
}

dat_boot <- draw.bootstrap(
  eusilc, REP = 1, hid = "hid", weights = "pWeight",
  strata = c("region", "hsize"), period = "year", split = TRUE,
  pid = "pidsplit")
# split households were considered e.g. household and
# split household were both selected or not selected
dat_boot[, data.table::uniqueN(w1), by = pidsplit][V1 > 1]
#> Empty data.table (0 rows and 2 cols): pidsplit,V1

```
