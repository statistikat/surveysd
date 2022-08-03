#'
#' Kish Factor
#'
#' Compute the design effect due to unequal weighting.
#'
#' @name kishFactor
#' @param w a numeric vector with weights
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
#' @return The function will return the the kish factor
#' @author Alexander Kowarik
#' @export kishFactor
#' @details The factor is computed acording to 'Weighting for Unequal P_i', Leslie Kish, Journal of Official Statistics, Vol. 8. No. 2, 1992
#' \deqn{ deff = \sqrt n \sum_j w_j^2 / (\sum_j w_j)^2}

#'
#' @examples
#' kishFactor(rep(1,10))
#' kishFactor(rlnorm(10))
kishFactor <- function(w, na.rm = FALSE) {
  if (!is.numeric(w)) {
    stop("The input must be a numeric vector")
  }
  if(na.rm){
    w <- na.omit(w)
  }
  n <- length(w)
  sqrt(n * sum(w ^ 2) / sum(w) ^ 2)
}
