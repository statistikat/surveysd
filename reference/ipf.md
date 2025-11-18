# Iterative Proportional Fitting

Adjust sampling weights to given totals based on household-level and/or
individual level constraints.

## Usage

``` r
ipf(
  dat,
  hid = NULL,
  conP = NULL,
  conH = NULL,
  epsP = 1e-06,
  epsH = 0.01,
  verbose = FALSE,
  w = NULL,
  bound = 4,
  maxIter = 200,
  meanHH = TRUE,
  allPthenH = TRUE,
  returnNA = TRUE,
  looseH = FALSE,
  numericalWeighting = computeLinear,
  check_hh_vars = TRUE,
  conversion_messages = FALSE,
  nameCalibWeight = "calibWeight",
  minMaxTrim = NULL,
  print_every_n = 100
)
```

## Arguments

- dat:

  a `data.table` containing household ids (optionally), base weights
  (optionally), household and/or personal level variables (numerical or
  categorical) that should be fitted.

- hid:

  name of the column containing the household-ids within `dat` or NULL
  if such a variable does not exist.

- conP:

  list or (partly) named list defining the constraints on person level.
  The list elements are contingency tables in array representation with
  dimnames corresponding to the names of the relevant calibration
  variables in `dat`. If a numerical variable is to be calibrated, the
  respective list element has to be named with the name of that
  numerical variable. Otherwise the list element shoud NOT be named.

- conH:

  list or (partly) named list defining the constraints on household
  level. The list elements are contingency tables in array
  representation with dimnames corresponding to the names of the
  relevant calibration variables in `dat`. If a numerical variable is to
  be calibrated, the respective list element has to be named with the
  name of that numerical variable. Otherwise the list element shoud NOT
  be named.

- epsP:

  numeric value or list (of numeric values and/or arrays) specifying the
  convergence limit(s) for `conP`. The list can contain numeric values
  and/or arrays which must appear in the same order as the corresponding
  constraints in `conP`. Also, an array must have the same dimensions
  and dimnames as the corresponding constraint in `conP`.

- epsH:

  numeric value or list (of numeric values and/or arrays) specifying the
  convergence limit(s) for `conH`. The list can contain numeric values
  and/or arrays which must appear in the same order as the corresponding
  constraints in `conH`. Also, an array must have the same dimensions
  and dimnames as the corresponding constraint in `conH`.

- verbose:

  if TRUE, some progress information will be printed.

- w:

  name if the column containing the base weights within `dat` or NULL if
  such a variable does not exist. In the latter case, every observation
  in `dat` is assigned a starting weight of 1.

- bound:

  numeric value specifying the multiplier for determining the weight
  trimming boundary if the change of the base weights should be
  restricted, i.e. if the weights should stay between 1/`bound`\*`w` and
  `bound`\*`w`.

- maxIter:

  numeric value specifying the maximum number of iterations that should
  be performed.

- meanHH:

  if TRUE, every person in a household is assigned the mean of the
  person weights corresponding to the household. If `"geometric"`, the
  geometric mean is used rather than the arithmetic mean.

- allPthenH:

  if TRUE, all the person level calibration steps are performed before
  the houshold level calibration steps (and `meanHH`, if specified). If
  FALSE, the houshold level calibration steps (and `meanHH`, if
  specified) are performed after everey person level calibration step.
  This can lead to better convergence properties in certain cases but
  also means that the total number of calibration steps is increased.

- returnNA:

  if TRUE, the calibrated weight will be set to NA in case of no
  convergence.

- looseH:

  if FALSE, the actual constraints `conH` are used for calibrating all
  the hh weights. If TRUE, only the weights for which the lower and
  upper thresholds defined by `conH` and `epsH` are exceeded are
  calibrated. They are however not calibrated against the actual
  constraints `conH` but against these lower and upper thresholds, i.e.
  `conH`-`conH`\*`epsH` and `conH`+`conH`\*`epsH`.

- numericalWeighting:

  See
  [numericalWeighting](https://statistikat.github.io/surveysd/reference/computeFrac.md)

- check_hh_vars:

  If `TRUE` check for non-unique values inside of a household for
  variables in household constraints

- conversion_messages:

  show a message, if inputs need to be reformatted. This can be useful
  for speed optimizations if ipf is called several times with similar
  inputs (for example bootstrapping)

- nameCalibWeight:

  character defining the name of the variable for the newly generated
  calibrated weight.

- minMaxTrim:

  numeric vector of length2, first element a minimum value for weights
  to be trimmed to, second element a maximum value for weights to be
  trimmed to.

- print_every_n:

  number of interation steps after which a summary table is printed. The
  summary table shows all constraints which are not yet reached
  according to `epsP` and `epsH`

## Value

The function will return the input data `dat` with the calibrated
weights `calibWeight` as an additional column as well as attributes. If
no convergence has been reached in `maxIter` steps, and `returnNA` is
`TRUE` (the default), the column `calibWeights` will only consist of
`NA`s. The attributes of the table are attributes derived from the
`data.table` class as well as the following.

|                                |                                                                                           |
|--------------------------------|-------------------------------------------------------------------------------------------|
| `converged`                    | Did the algorithm converge in `maxIter` steps?                                            |
| `iterations`                   | The number of iterations performed.                                                       |
| `conP`, `conH`, `epsP`, `epsH` | See Arguments.                                                                            |
| `conP_adj`, `conH_adj`         | Adjusted versions of `conP` and `conH`                                                    |
| `formP`, `formH`               | Formulas that were used to calculate `conP_adj` and `conH_adj` based on the output table. |

## Details

This function implements the weighting procedure described here:
[doi:10.17713/ajs.v45i3.120](https://doi.org/10.17713/ajs.v45i3.120) .
Usage examples can be found in the corresponding vignette
([`vignette("ipf")`](https://statistikat.github.io/surveysd/articles/ipf.md)).

`conP` and `conH` are contingency tables, which can be created with
`xtabs`. The `dimnames` of those tables should match the names and
levels of the corresponding columns in `dat`.

`maxIter`, `epsP` and `epsH` are the stopping criteria. `epsP` and
`epsH` describe relative tolerances in the sense that \$\$1-epsP \<
\frac{w\_{i+1}}{w_i} \< 1+epsP\$\$ will be used as convergence
criterium. Here i is the iteration step and wi is the weight of a
specific person at step i.

The algorithm performs best if all varables occuring in the constraints
(`conP` and `conH`) as well as the household variable are coded as
`factor`-columns in `dat`. Otherwise, conversions will be necessary
which can be monitored with the `conversion_messages` argument. Setting
`check_hh_vars` to `FALSE` can also incease the performance of the
scheme.

## Author

Alexander Kowarik, Gregor de Cillia

## Examples

``` r
if (FALSE) { # \dontrun{

# load data
eusilc <- demo.eusilc(n = 1, prettyNames = TRUE)

# personal constraints
conP1 <- xtabs(pWeight ~ age, data = eusilc)
conP2 <- xtabs(pWeight ~ gender + region, data = eusilc)
conP3 <- xtabs(pWeight*eqIncome ~ gender, data = eusilc)

# household constraints
conH1 <- xtabs(pWeight ~ hsize + region, data = eusilc[!duplicated(hid)])

# simple usage ------------------------------------------

calibweights1 <- ipf(
  eusilc,
  conP = list(conP1, conP2, eqIncome = conP3),
  bound = NULL,
  verbose = TRUE
)

# compare personal weight with the calibweigth
calibweights1[, .(hid, pWeight, calibWeight)]

# advanced usage ----------------------------------------

# use an array of tolerances
epsH1 <- conH1
epsH1[1:4, ] <- 0.005
epsH1[5, ] <- 0.2

# create an initial weight for the calibration
eusilc[, regSamp := .N, by = region]
eusilc[, regPop := sum(pWeight), by = region]
eusilc[, baseWeight := regPop/regSamp]

calibweights2 <- ipf(
  eusilc,
  conP = list(conP1, conP2),
  conH = list(conH1),
  epsP = 1e-6,
  epsH = list(epsH1),
  bound = 4,
  w = "baseWeight",
  verbose = TRUE
)

# show an adjusted version of conP and the original
attr(calibweights2, "conP_adj")
attr(calibweights2, "conP")
} # }
```
