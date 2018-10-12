* Create vignette for error estimation
* Bugfixes: 8c4529a, faed7f9

# surveysd 0.2.0

* Documentation updates
    * Use markown pre-processing for documentation
    * Add github pages, more badges and a logo
    * export `draw.bootstrap` and include `prettyNames` argument. 
      update documenentation examples to use new variable names
    * Don't load `laeken` and `data.table` in examples
    * rewrite `vignettes/TheoryWord.Rmd` into a dynamic markdown file 
      describing the methodology of the package

* Automatically pass variable names for household id, weights and
  reference period from `draw.bootstrap` to downstream functions
* Reorganize names of R source files
