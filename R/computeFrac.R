#' Numerical weighting functions
#'
#' Customize weight-updating within factor levels in case of numerical
#' calibration. The functions described here serve as inputs for [ipf].
#'
#' `computeFrac` provides the "standard" IPU updating scheme given as
#'
#' \deqn{f = target/curValue}
#'
#' which means that each weight inside the level will be multtiplied by the same
#' factor when doing the actual update step (`w := f*w`). `computeLinear` on the
#' other hand calculates `f` as
#'
#' \ifelse{html}{
#'   \out{<div style="text-align: center;"> f<sub>i</sub> = a  &middot; x<sub>i</sub> + b </div>}
#' }{\deqn{f_i = ax_i+b}}
#'
#' where `a` and `b` are chosen, so f satisfies the following two equations.
#'
#' \ifelse{html}{
#'   \out{<div style="text-align: center;">&sum; f<sub>i</sub> w<sub>i</sub>
#'        x<sub>i</sub> = target</div>}
#' }{\deqn{\sum f_i * w_i * x_i = target}}
#' \ifelse{html}{
#'   \out{<div style="text-align: center;">&sum; f<sub>i</sub>
#'        w<sub>i</sub> = &sum; w<sub>i</sub></div>}
#' }{\deqn{\sum f_i * w_i = \sum w_i}}
#'
# \eqn{\sum}\out{f<sub>i</sub> w<sub>i</sub> x<sub>i</sub>} = `target`
#' `computeLinearG1` calculates `f` in the same way as `computeLinear`, but if
#'`f_i*w_i<1` `f_i` will be set to `1/w_i`.
#'
#' @md
#' @param curValue Current summed up value. Same as `sum(x*w)`
#' @param target Target value. An element of `conP` in [ipf]
#' @param x Vector of numeric values to be calibrated against
#' @param w Vector of weights
#' @param boundLinear The output `f` will satisfy
#'   `1/boundLinear <= f <= boundLinear`. See `bound` in [ipf]
#'
#' @return A weight multiplier `f`
#'
#' @aliases numericalWeighting
#' @export computeFrac
computeFrac <- function(curValue, target, x, w) {
  target / curValue
}
