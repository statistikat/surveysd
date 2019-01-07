# surveysd 1.0.0

* copy `simPop::ipu2()` and some related functions to surveysd. 
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
