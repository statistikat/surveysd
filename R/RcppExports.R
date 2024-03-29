# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @rdname computeFrac
#' @export
computeLinear <- function(curValue, target, x, w, boundLinear = 10) {
    .Call('_surveysd_computeLinear', PACKAGE = 'surveysd', curValue, target, x, w, boundLinear)
}

#' @rdname computeFrac
#' @export
computeLinearG1_old <- function(curValue, target, x, w, boundLinear = 10) {
    .Call('_surveysd_computeLinearG1_old', PACKAGE = 'surveysd', curValue, target, x, w, boundLinear)
}

#' @rdname computeFrac
#' @export
computeLinearG1 <- function(curValue, target, x, w, boundLinear = 10) {
    .Call('_surveysd_computeLinearG1', PACKAGE = 'surveysd', curValue, target, x, w, boundLinear)
}

#' Calculate mean by factors
#'
#' These functions calculate the arithmetic and geometric mean of the weight for each class. `geometric_mean` and
#' `arithmetic_mean` return a `numeric` vector of the same length as `w` which stores the averaged weight for each
#' observation. `geometric_mean_reference` returns the same value by reference, i.e. the input value `w` gets
#' overwritten by the updated weights. See examples.
#'
#' @md
#' @name cpp_mean
#' @param w An numeric vector. All entries should be positive.
#' @param classes A factor variable. Must have the same length as `w`.
#' @examples
#'
#' \dontrun{
#'
#' ## create random data
#' nobs <- 10
#' classLabels <- letters[1:3]
#' dat = data.frame(
#'   weight = exp(rnorm(nobs)),
#'   household = factor(sample(classLabels, nobs, replace = TRUE))
#' )
#' dat
#'
#' ## calculate weights with geometric_mean
#' geom_weight <- geometric_mean(dat$weight, dat$household)
#' cbind(dat, geom_weight)
#'
#' ## calculate weights with arithmetic_mean
#' arith_weight <- arithmetic_mean(dat$weight, dat$household)
#' cbind(dat, arith_weight)
#'
#' ## calculate weights "by reference"
#' geometric_mean_reference(dat$weight, dat$household)
#' dat
#' }
geometric_mean_reference <- function(w, classes) {
    invisible(.Call('_surveysd_geometric_mean_reference', PACKAGE = 'surveysd', w, classes))
}

geometric_mean <- function(w, classes) {
    .Call('_surveysd_geometric_mean', PACKAGE = 'surveysd', w, classes)
}

arithmetic_mean <- function(w, classes) {
    .Call('_surveysd_arithmetic_mean', PACKAGE = 'surveysd', w, classes)
}

rollMeanC <- function(x, k, type) {
    .Call('_surveysd_rollMeanC', PACKAGE = 'surveysd', x, k, type)
}

rollSumC <- function(x, k, type) {
    .Call('_surveysd_rollSumC', PACKAGE = 'surveysd', x, k, type)
}

#' @name PointEstimates
#' @title Weighted Point Estimates
#'
#' @description Predefined functions for weighted point estimates in package `surveysd`.
#'
#' @param x numeric vector
#' @param w weight vector
#'
#' @details Predefined functions are weighted ratio and weighted sum.
#'
#' @return
#' Each of the functions return a single numeric value
#' @examples
#' x <- 1:10
#' w <- 10:1
#' weightedRatio(x,w)
#' @export
weightedRatio <- function(x, w) {
    .Call('_surveysd_weightedRatio', PACKAGE = 'surveysd', x, w)
}

#' @rdname PointEstimates
#' @examples
#' x <- 1:10
#' w <- 10:1
#' weightedSum(x,w)
#' @export
weightedSum <- function(x, w) {
    .Call('_surveysd_weightedSum', PACKAGE = 'surveysd', x, w)
}

#' Perform one step of iterative proportional updating
#'
#' C++ routines to invoke a single iteration of the Iterative proportional updating (IPU) scheme. Targets and classes
#' are assumed to be one dimensional in the `ipf_step` functions. `combine_factors` aggregates several vectors of
#' type factor into a single one to allow multidimensional ipu-steps. See examples.
#'
#' `ipf_step` returns the adjusted weights. `ipf_step_ref` does the same, but updates `w` by reference rather than
#' returning. `ipf_step_f` returns a multiplicator: adjusted weights divided by unadjusted weights. `combine_factors` is
#' designed to make `ipf_step` work with contingency tables produced by [xtabs].
#'
#' @md
#' @name ipf_step
#' @param w a numeric vector of weights. All entries should be positive.
#' @param classes a factor variable. Must have the same length as `w`.
#' @param targets key figure to target with the ipu scheme. A numeric verctor of the same length as `levels(classes)`.
#' This can also be a `table` produced by `xtabs`. See examples.
#' @examples
#'
#' ############# one-dimensional ipu ##############
#'
#' ## create random data
#' nobs <- 10
#' classLabels <- letters[1:3]
#' dat = data.frame(
#'   weight = exp(rnorm(nobs)),
#'   household = factor(sample(classLabels, nobs, replace = TRUE))
#' )
#' dat
#'
#' ## create targets (same lenght as classLabels!)
#' targets <- 3:5
#'
#' ## calculate weights
#' new_weight <- ipf_step(dat$weight, dat$household, targets)
#' cbind(dat, new_weight)
#'
#' ## check solution
#' xtabs(new_weight ~ dat$household)
#'
#' ## calculate weights "by reference"
#' ipf_step_ref(dat$weight, dat$household, targets)
#' dat
#'
#' ############# multidimensional ipu ##############
#'
#' ## load data
#' factors <- c("time", "sex", "smoker", "day")
#' tips <- data.frame(sex=c("Female","Male","Male"), day=c("Sun","Mon","Tue"),
#' time=c("Dinner","Lunch","Lunch"), smoker=c("No","Yes","No"))
#' tips <- tips[factors]
#'
#' ## combine factors
#' con <- xtabs(~., tips)
#' cf <- combine_factors(tips, con)
#' cbind(tips, cf)[sample(nrow(tips), 10, replace = TRUE),]
#'
#' ## adjust weights
#' weight <- rnorm(nrow(tips)) + 5
#' adjusted_weight <- ipf_step(weight, cf, con)
#'
#' ## check outputs
#' con2 <- xtabs(adjusted_weight ~ ., data = tips)
#' sum((con - con2)^2)
#'
#' @rdname ipf_step
#' @export
ipf_step_ref <- function(w, classes, targets) {
    invisible(.Call('_surveysd_ipf_step_ref', PACKAGE = 'surveysd', w, classes, targets))
}

#' @rdname ipf_step
#' @export
ipf_step <- function(w, classes, targets) {
    .Call('_surveysd_ipf_step', PACKAGE = 'surveysd', w, classes, targets)
}

#' @rdname ipf_step
#' @export
ipf_step_f <- function(w, classes, targets) {
    .Call('_surveysd_ipf_step_f', PACKAGE = 'surveysd', w, classes, targets)
}

