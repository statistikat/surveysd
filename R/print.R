#' @title Print function for surveysd objects
#'
#' @description
#' Prints the results of a call to [calc.stError]. Shows used variables and
#' function, number of point estiamtes
#' as well as properties of the results.
#'
#' @param x an object of class `'surveysd'`
#' @param ... additonal parameters
#'
#' @export
print.surveysd <- function(x, ...) {

  # get parameter
  col.val <- grepl("^val_*.", colnames(x[["Estimates"]]))
  col.val <- colnames(x[["Estimates"]])[col.val]
  # get number of estimates ~ variables in number of groups using function fun
  #  from package pack
  n.estimates <- nrow(x[["Estimates"]]) * length(col.val)
  # get number of subgroup which fall below size
  n.periods <- unique(x[["Estimates"]][[x[["param"]][["period"]]]])
  n.periods <- length(n.periods[!grepl("-|_", n.periods)])

  if (!is.null(unlist(x[["param"]][["group"]]))) {
    n.groups <- nrow(unique(
      x[["Estimates"]][, .N, by = c(unique(unlist(x[["param"]][["group"]])))]
    ))
  } else {
    n.groups <- NULL
  }

  # get number of missing values in the output due to subgroups with no
  #   observations or subgroups with all observations equal NA
  n.NAs <- sum(is.na(x[["Estimates"]][, mget(col.val)]))

  # print number of estimates calculated as well as functino and variables used
  cat("Calculated point estimates for variable(s)\n\n",
      paste(x[["param"]][["var"]], sep = ","), "\n\n")

  if (!is.null(n.groups)) {
    cat("Results hold", n.estimates, "point estimates for", n.periods,
        "periods in", n.groups, "subgroups\n")
  } else {
    cat("Results hold", n.estimates, "point estimates for", n.periods,
        "periods\n")
  }
  cat("\n")
  if (nrow(x[["smallGroups"]]) != 0) {
    if (nrow(x[["smallGroups"]]) > 10) {
      cat(nrow(x[["smallGroups"]]), "subgroups contained less than",
          x[["param"]][["size.limit"]], "observations\n")
    } else {
      cat("Subgroups with less than", x[["param"]][["size.limit"]],
          "observations\n")
      print(x[["smallGroups"]])
    }
    cat("\n")
  }


  # print number of missing values for point estimates
  if (n.NAs > 0) {
    cat("Point estimates are NAs in", n.NAs, "cases due to either\n no ",
        "observations or only NAs for the variable(s) in the corresponding",
        " subgroups.\n")
    cat("\n")
  }


  # get number of point estimates where sd exceeds cv.limit
  stEtoohigh <- colnames(x[["cvHigh"]])
  stEtoohigh <- stEtoohigh[!stEtoohigh %in% c(
    x[["param"]][["period"]], unique(unlist(x[["param"]][["group"]])))]
  stEtoohigh <- as.matrix(subset(x[["cvHigh"]], select = stEtoohigh))
  if (sum(stEtoohigh, na.rm = TRUE) > 0) {
    cat("Estimated standard error exceeds", x[["param"]][["cv.limit"]],
        "% of the the point estimate in", sum(stEtoohigh, na.rm = TRUE),
        "cases\n")
    cat("\n")
  }
}
