#' Calculate point estimates for weights and bootstrap weights
#'
#' Calculate an estimate for each group and each bootstrap replicate of survey weights. The
#' replicate weihgts are assumed to be in wide format as returned by [recalib] and
#' [draw.bootstrap]
#'
#' @param dat A [data.table] containing bootstrap weights. This will usually be the output of
#'            [recalib] or [draw.bootstrap].
#' @param FUN An aggregation function as returned by [WeightedMean], [WeightedMedian], [WeightedSum],
#'        [WeightedRatio] and [WeightedSize].
#' @param groupVars Column names in `dat` refering to grouping variables. Typically, those
#'        columns are `factor`s or `integer`s
#' @param weights Column names for the survey weights
#' @param b.weights Column names for the boostrap weights
#'
#' @examples
#' eusilc <- demo.eusilc(prettyNames = TRUE)
#' dat_boot <- draw.bootstrap(eusilc, REP = 10, hid = "hid", weights = "pWeight",
#'                            strata = "region", period = "year")
#' dat_boot_calib <- recalib(dat_boot, conP.var = c("gender", "age"), conH.var = c("region", "hsize"))
#'
#' ## poverty rate for each year
#' aggregateReplicates(dat_boot_calib, WeightedRatio("povertyRisk == 1"), c("year"))
#'
#' ## mean income for each year and gender
#' aggregateReplicates(dat_boot_calib, WeightedMean("eqIncome", na.rm = TRUE), c("year", "gender"))
#' @export
aggregateReplicates <- function(dat, FUN, groupVars = NULL,
                                weights = attr(dat, "weights"), b.weights = attr(dat, "b.rep")) {
  b.weights <- c(weights, b.weights)
  res <- dat[, lapply(b.weights, function(w) { FUN(.SD, .SD[[w]]) }), by = groupVars]
  names(res)[length(groupVars) + 1:length(b.weights)] <- b.weights
  attributes(res) <- c(attributes(res), list(call = match.call()))
  res
}

#' Point estimate functions
#'
#' Helper functions to define point estimates for [aggregateReplicates].
#'
#' @param var A variable name in the dataset to be used for estimation
#' @param na.rm logical. Should missing values (including `NaN`) be removed?
#'
#' @return A function taking a subset of the data and a weight vector as input and returning
#'         a single numeric value. For example `function(subset, w){ sum(w) }`
#'
#' @examples
#' df <- data.table::data.table(age = rnorm(100, mean = 30, sd = 10), w = 1)
#'
#' estimator <- WeightedMean("age")
#' estimator(df, df$w)
#' @export
WeightedMean <- function(var, na.rm = FALSE) {
  function(subset, w) {
    mean(subset[[var]] * w, na.rm = na.rm)
  }
}

#' @rdname WeightedMean
#' @export
WeightedMedian <- function(var, na.rm = FALSE) {
  function(subset, w) {
    laeken::weightedMedian(subset[[var]], w, na.rm = na.rm)
  }
}

#' @rdname WeightedMean
#' @export
WeightedSum <- function(var, na.rm = FALSE) {
  function(subset, w) {
    sum(subset[[var]] * w, na.rm = na.rm)
  }
}

#' @rdname WeightedMean
#'
#' @param condition A string specifying a condition that can be evaluated for each observation.
#'                  This condition may contain column names in the dataset.
#'
#' @details
#' `WeightedRatio` returns the relative number of observations in each group
#' satisfying `condition`. For example, a poverty-risk rate for each state.
#'
#' @examples
#'
#' estimator <- WeightedRatio("age>30")
#' estimator(df, df$w)
#'
#' @export
WeightedRatio <- function(condition, na.rm = FALSE) {
  function(subset, w) {
    subset[, sum(w * eval(parse(text = condition)), na.rm = na.rm)/sum(w)]
  }
}

#' @rdname WeightedMean
#' @examples
#'
#' estimator <- WeightedSize("age>30")
#' estimator(df, df$w)
#' @export
WeightedSize <- function(condition, na.rm = FALSE) {
  function(subset, w) {
    subset[, sum(w * eval(parse(text = condition)), na.rm = na.rm)]
  }
}
