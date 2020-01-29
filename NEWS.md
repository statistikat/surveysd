# surveysd 1.2.0

* Add vignette about `ipf()` (#23)
* Add additional check to `ipf()` so the procedure does not fail silently ([33fba8d](https://github.com/statistikat/surveysd/commit/33fba8d))
* Resolve bug that appeared if constraints contained only one household or only
  one person (#17)
* Improve code style according to suggestions from `lintr::lint_package()`
* `ipf()` now throws a warning if convergence can't be reached

# surveysd 1.1.1

* `recalib()` accepts `conP` and `conH` in the same way as `ipf` (#21)
* `recalib()` has arguments `epsP` and `epsH`, to make convergence limits more
  transparent ([3def439](https://github.com/statistikat/surveysd/commit/3def439))
* fixed bug in `ipf()` when supplying `hid="hid"` (#20). Thanks [@asiripanich](https://github.com/asiripanich)

# surveysd 1.1.0

* uncouple surveysd from simPop (#13)
    * use `demo.eusilc()` for all examples and unit tests in this package
    * remove simPop from the "Suggests" field (DESCRIPTION)
* remove hardcoded variable names in default parameters (#15)
* allow calls to `recalib()` without specifying `conP.vars` or `conH.vars`

# surveysd 1.0.2

* fix bug with NSE
* update documentation for `calc.stError()`

# surveysd 1.0.1

## bugfixes

* fix bug with fpc calculation ([0824834](https://github.com/statistikat/surveysd/commit/0824834))
* fix issue when household column is not a factor (#12)
* fix bug with dummy data ([5adb6b7](https://github.com/statistikat/surveysd/commit/5adb6b7))
* resolve issues when datasets have columns with certain names (#7, #10)

## linters

* improve readability of code and resolve linters
* set up automated linter-checking

## documentation

* automatize the gh-pages builds
    * css updates
* new badges (code coverage, cran)
* update setup for vignettes

# surveysd 1.0.0

* copy `simPop::ipu2()` and some related functions to surveysd (#6)
    * add bugfixes for ported functions
* simplify tests
* update references for the methodology vignette
* fix sampling with full population in some strata

# surveysd 0.2.3

* Fix citation
* Resolve issue that occured when columns had multiple classes

# surveysd 0.2.2

* Update some default parameters (#4, #5)
* Split README into README + Get Started
* Bugfixes
* Update documentation titles, function categories and vignette titles

# surveysd 0.2.1

* Create vignette for error estimation
* Bugfixes for plotting and the `...` argument in `calc.stError()`
* Documentation updates. Make example plots available in gh-pages

# surveysd 0.2.0

* Automatically pass variable names for household id, weights and
  reference period from `draw.bootstrap()` to downstream functions
* Reorganize names of R source files

## documentation

* Use markown pre-processing for documentation
* Add github pages, more badges and a logo
* export `draw.bootstrap()` and include `prettyNames` argument. 
  update documenentation examples to use new variable names
* Don't load `laeken` and `data.table` in examples
* rewrite `vignettes/TheoryWord.Rmd` into a dynamic markdown file 
  describing the methodology of the package
