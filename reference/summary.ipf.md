# Generate Summary Output for IPF Calibration

Generates a detailed summary of an Iterative Proportional Fitting (IPF)
calibration, providing a complete tool for evaluating the calibration's
success and the validity of the resulting weights.

The output is a list of data.tables for a comprehensive evaluation,
including:

**Calibration Results:**

- `calib_results_conP_*` and `calib_results_conH_*`: Key diagnostic
  tables that compare calibrated margins to population targets and
  assess the goodness of fit via metrics like `maxFac`.

**Data and Diagnostics:**

- `weighted data`: An excerpt of the final dataset with the calculated
  calibration weights.

- `distribution of the weights`: A statistical overview of the weight
  distribution (min, max, CV).

**Detailed Margin Comparisons:**

- `conP_*`, `conH_*`, `*_adjusted`, `*_original`, `*_rel_diff_*`: Tables
  that compare original sample margins, calibrated margins, and
  population targets, along with their relative differences.

## Usage

``` r
# S3 method for class 'ipf'
summary(object, ...)
```

## Arguments

- object:

  object of class ipf

- ...:

  additional arguments

## Value

a list of the following outputs

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
conH1 <- xtabs(pWeight ~ hsize + region, data = eusilc)

# simple usage ------------------------------------------

calibweights1 <- ipf(
 eusilc,
 conP = list(conP1, conP2, eqIncome = conP3),
 bound = NULL,
 verbose = TRUE
)
output <- summary(calibweights1)
# the output can easily be exported to an Excel file, e.g. with
# library(openxlsx)
# write.xlsx(output, "SummaryIPF.xlsx")
} # }
```
