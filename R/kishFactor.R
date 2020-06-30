#'
#' Kish Factor
#'
#' Compute the design effect due to unequal weighting.
#'
#' @name kishFactor
#' @param w a numeric vector with weights
#' @return The function will return the the kish factor
#' @author Alexander Kowarik
#' @export kishFactor
#' @details The factor is computed acording to 'Weighting for Unequal P_i', Leslie Kish, Journal of Official Statistics, Vol. 8. No. 2, 1992
#' \deqn{ deff = \sqrt n \sum_j w_j^2 / (\sum_j w_j)^2}

#'
#' @examples
#' kishFactor(rep(1,10))
#' kishFactor(rlnorm(10))
kishFactor <- function(w) {
  if (!is.numeric(w)) {
    stop("The input must be a numeric vector")
  }
  n <- length(w)
  sqrt(n * sum(w ^ 2) / sum(w) ^ 2)
}
