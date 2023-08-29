# surveysd 1.3.2

* `rescaled.bootstrap()` has additional parameter `period` which is identical to the one in `draw.bootstrap`. If `period` is not `NULL` the boostraps will be drawn such that in each period and strata/cluster only $\floor{\frac{n}{2}}$ records are drawn. This produces more consisten results and should make calibration afterwards easier.
* Improved numerical weighting with `computeLinearG1` for use with `ipf()`. `computeLinearG1` is now more stable when only numerical variables are used for weighting.
* parameter `numericalWeighting` can be passed to `recalib` and will correctly be passed along to function `ipf`


# surveysd 1.3

* new parameter minMaxTrim in `ipf()` to trim weights
* fix bug for numericalWeighting reported "converged" when it did not converge
* fix bug when checking cluster and strata for bootstrapping


# surveysd 1.2

* Add vignette about `ipf()`
* Resolve bug that appeared if constraints contained only one household or only
  one person ([#17](https://github.com/statistikat/surveysd/issues/17))

# surveysd 1.1.1

* `recalib()` accepts conP and conH in the same way as `ipf`
* `recalib()` has arguments `epsP` and `epsH`, to make convergence limits more transparant
* fixed bug in `ipf()` when supplying `hid="hid"`, see https://github.com/statistikat/surveysd/pull/20. Thanks [@asiripanich](https://github.com/asiripanich)

# surveysd 1.1.0

* uncouple surveysd from simPop
    * use `demo.eusilc()` for all examples and unit tests in this package
    * remove simPop from the "Suggests" field (DESCRIPTION)
* remove hardcoded variable names in default parameters
* allow calls to `recalib()` without specifying `conP.vars` or `conH.vars`

# surveysd 1.0.2

* fix bug with NSE
* update documentation for calc.stError

# surveysd 1.0.1

* fix bug with fpc calculation ([0824834](https://github.com/statistikat/surveysd/commit/0824834))
* fix bug with dummy data ([5adb6b7](https://github.com/statistikat/surveysd/commit/5adb6b7))
* improve readability of code and resolve linters
* resolve issues when datasets have columns with certain names (#7, #10)
* automatize the gh-pages builds
    * css updates
* automated linter-checking
* new badges (code coverage, cran)
* update setup for vignettes
* fix issue when household column is not a factor (#12)

# surveysd 1.0.0

* copy `simPop::ipu2()` and some related functions to surveysd.
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
* Bugfixes for plotting and `...` argument in `calc.stError`
* Documentation updates. Make example plots available in gh-pages

# surveysd 0.2.0

* Automatically pass variable names for household id, weights and
  reference period from `draw.bootstrap` to downstream functions
* Reorganize names of R source files

## Documentation updates
    
* Use markown pre-processing for documentation
* Add github pages, more badges and a logo
* export `draw.bootstrap` and include `prettyNames` argument. 
  update documenentation examples to use new variable names
* Don't load `laeken` and `data.table` in examples
* rewrite `vignettes/TheoryWord.Rmd` into a dynamic markdown file 
  describing the methodology of the package
