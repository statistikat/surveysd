# Generate multiple years of EU-SILC data

Create a dummy dataset to be used for demonstrating the functionalities
of the `surveysd` package based on
[laeken::eusilc](https://rdrr.io/pkg/laeken/man/eusilc.html). Please
refer to the documentation page of the original data for details about
the variables.

## Usage

``` r
demo.eusilc(n = 8, prettyNames = FALSE)
```

## Arguments

- n:

  Number of years to generate. Should be at least 1

- prettyNames:

  Create easy-to-read names for certain variables. Recommended for
  demonstration purposes. Otherwise, use the original codes documented
  in [laeken::eusilc](https://rdrr.io/pkg/laeken/man/eusilc.html).

## Details

If `prettyNames` is `TRUE`, the following variables will be available in
an easy-to-read manner.

- `hid` Household id. Consistent with respect to the reference period
  (`year`)

- `hsize` Size of the household. derived from `hid` and `period`

- `region` Federal state of austria where the household is located

- `pid` Personal id. Consistent with respect to the reference period
  (`year`)

- `age` Age-class of the respondent

- `gender` A persons gender (`"male"`, `"Female"`)

- `ecoStat` Ecnomic status (`"part time"`, `"full time"`,
  `"unemployed"`, ...)

- `citizenship` Citizenship (`"AT"`, `"EU"`, `"other"`)

- `pWeight` Personal sample weight inside the reference period

- `year`. Simulated reference period

- `povertyRisk`. Logical variable determining whether a respondent is at
  risk of poverty

## Examples

``` r
demo.eusilc(n = 1, prettyNames = TRUE)[, c(1:8, 26, 28:30)]
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
#>        eqIncome  pWeight  year povertyRisk
#>           <num>    <num> <num>      <lgcl>
#>     1: 16090.69 504.5696  2010       FALSE
#>     2: 16090.69 504.5696  2010       FALSE
#>     3: 16090.69 504.5696  2010       FALSE
#>     4: 27076.24 493.3824  2010       FALSE
#>     5: 27076.24 493.3824  2010       FALSE
#>    ---                                    
#> 14823: 26508.20 556.4260  2010       FALSE
#> 14824: 14387.41 643.2557  2010       FALSE
#> 14825: 14805.83 679.7288  2010       FALSE
#> 14826: 16288.30 567.1544  2010       FALSE
#> 14827: 16288.30 567.1544  2010       FALSE
```
