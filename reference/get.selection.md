# Get sample selection (~deltas) from drawn bootstrap replicates

Reconstruct sample selection, e.g. record was drawn or not drawn (delta
= 0/1) in each sampling stage from bootstrap replicates.
`get.selection()` needs the `cluster`, `strata` and `hid`/`pid`
information (if not `NULL`) to correctly reconstruct if a record was
drawn in each sampling stage for each bootstrap replicate. Is only
needed if bootstrap replicates are drawn for a survey with existing
bootstrap replicates from a previous period, see parameter
`already.selected` in function
[`draw.bootstrap()`](https://statistikat.github.io/surveysd/reference/draw.bootstrap.md).

## Usage

``` r
get.selection(
  dat,
  b.rep = attr(dat, "b.rep"),
  strata = attr(dat, "strata"),
  cluster = attr(dat, "cluster"),
  hid = attr(dat, "hid"),
  pid = attr(dat, "pid")
)
```

## Arguments

- dat:

  either data.frame or data.table containing the survey data with
  rotating panel design. Should contain only survey data from a single
  time period.

- b.rep:

  character specifying the names of the columns in `dat` containing
  bootstrap replicates.

- strata:

  character vector specifying the name(s) of the column in `dat` by
  which the population was stratified.

- cluster:

  character vector specifying cluster in the data.

- hid:

  character specifying the name of the column in `dat` containing the
  household id. If `NULL` (the default), the household structure is not
  regarded. `hid` and `pid` cannot both be `NULL`.

- pid:

  pid column in `dat` specifying the personal identifier. This
  identifier needs to be unique for each person throught the whole data
  set. `hid` and `pid` cannot both be `NULL`.

## Value

Returns a list of data.tables. The length of the list equals the number
of sampling stages specified. Each list entry contains a `data.table`
with variables for sampling stage and/or `hid`/`pid` as well as
`length(attr(dat,"b.rep"))` columns each indicating if record/cluster
was drawn in the respective sampling stage for the i-th boostrap
replicate.

## Examples

``` r
library(surveysd)
library(data.table)
setDTthreads(1)
set.seed(1234)
eusilc <- demo.eusilc(n = 3, prettyNames = TRUE)

## draw replicates with stratification
dat_boot <- draw.bootstrap(eusilc[year<2012], REP = 3, weights = "pWeight",
                           strata = "region", hid = "hid",
                           period = "year")

## get selection matrix for year 2011 
dat_selection <- get.selection(dat_boot[year==2011])
print(dat_selection)
#> $SamplingStage1
#> Key: <region, hid>
#>           region   hid delta_1_1 delta_1_2 delta_1_3
#>           <fctr> <int>    <lgcl>    <lgcl>    <lgcl>
#>    1: Burgenland    12     FALSE     FALSE      TRUE
#>    2: Burgenland    59     FALSE     FALSE      TRUE
#>    3: Burgenland   112     FALSE      TRUE      TRUE
#>    4: Burgenland   135     FALSE      TRUE      TRUE
#>    5: Burgenland   170     FALSE      TRUE     FALSE
#>   ---                                               
#> 5996: Vorarlberg  7384     FALSE     FALSE      TRUE
#> 5997: Vorarlberg  7396      TRUE     FALSE     FALSE
#> 5998: Vorarlberg  7437     FALSE     FALSE      TRUE
#> 5999: Vorarlberg  7445     FALSE     FALSE     FALSE
#> 6000: Vorarlberg  7488      TRUE     FALSE      TRUE
#> 

## draw bootstrap replicates for year 2012
## respecting already selected units for year 2011 ~ dat_selection
## in order to mimic rotating panel design
dat_boot_2012 <- draw.bootstrap(eusilc[year==2012], REP = 3, weights = "pWeight",
                                strata = "region", hid = "hid",
                                period = "year", 
                                already.selected = dat_selection)


```
