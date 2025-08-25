#' Print method for IPF calibration summary
#'
#' This function provides a concise summary of an IPF (Iterative Proportional Fitting) calibration object. 
#' It automatically extracts the calibration weight from `all_formulas`, computes the Kish factor for the weights, and prints the first 10 rows of any `calib_results_` tables. 
#' Useful for a quick overview of calibration results.
#' Additional details can be explored with `str()` or `names()`.

print.summary.ipf <- function(x, ...) {
  cat("IPF Calibration Summary (Short Version)\n")
  cat("=======================================\n\n")
  
  # --- 1. Automatically extract CalibWeight from all_formulas ---
  if (!is.null(x$all_formulas)) {
    calibWeightName <- trimws(strsplit(x$all_formulas[1], "~")[[1]][1])
  } else {
    calibWeightName <- NULL
  }
  
  # --- 2. Calculate Kish factor ---
  weight_candidates <- names(x[["weighted data"]])
  if (!is.null(calibWeightName) && calibWeightName %in% weight_candidates) {
    w <- x[["weighted data"]][[calibWeightName]]
    kish <- kishFactor(w, na.rm = TRUE)
    cat(sprintf("Kish Factor (based on '%s'): %.4f\n\n", calibWeightName, kish))
  } else {
    cat("Kish Factor: No valid weight column found\n\n")
  }
  
  # --- 3. Display calib_results_ tables (first 10 rows) ---
  calib_names <- grep("^calib_results_", names(x), value = TRUE)
  if (length(calib_names) > 0) {
    cat("Calibration Results (first 10 rows per table):\n")
    for (nm in calib_names) {
      cat("\n---", nm, "---\n")
      print(utils::head(x[[nm]], 10))
    }
  } else {
    cat("No calib_results_ tables found.\n")
  }
  
  cat("\n(See more details with str() or names())\n")
  invisible(x)
}