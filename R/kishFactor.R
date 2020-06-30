#'
#' Kish Factor
#'
#' Compute the kish factor for a specific weight vector
#'
#' @name kishFactor
#' @param w a numeric vector with weights
#' @return The function will return the the kish factor
#' @author Alexander Kowarik
#' @export kishFactor
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
