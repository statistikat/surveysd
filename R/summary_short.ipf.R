#' Generate a concise summary for a ipf calibration
#'
#' @description This function generates a condensed summary of the calibration results, focusing on the most important metrics. It provides a clear overview of how much the original data deviated from the target values and how successfully the calibration process adjusted the weights to meet those targets.
#'
#' @param object An object of class `ipf`.
#' @param ... additional arguments passed to other methods.
#'
#' @return A list of data.tables, containing the following outputs:
#' \itemize{
#'   \item{\strong{`conP_i_summary`}: A data.table for each person-level constraint, showing the original and calibrated relative differences from the target values.}
#'   \item{\strong{`conH_i_summary`}: A data.table for each household-level constraint, showing the original and calibrated relative differences from the target values.}
#'   \item{\strong{`distribution of weights`}: A data.table with key statistics (Coefficient of Variation, min, max, and quotient) for the calibrated weights, grouped by the calibration variables.}
#' }
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # load data
#' eusilc <- demo.eusilc(n = 1, prettyNames = TRUE)
#'
#' # personal constraints
#' conP1 <- xtabs(pWeight ~ age, data = eusilc)
#' conP2 <- xtabs(pWeight ~ gender + region, data = eusilc)
#' conP3 <- xtabs(pWeight*eqIncome ~ gender, data = eusilc)
#'
#' # household constraints
#' conH1 <- xtabs(pWeight ~ hsize + region, data = eusilc)
#'
#' # simple usage ------------------------------------------
#'
#' calibweights1 <- ipf(
#'    eusilc,
#'    conP = list(conP1, conP2, eqIncome = conP3),
#'    bound = NULL,
#'    verbose = TRUE
#' )
#'
#' # Use the new, short summary function
#' output_short <- summary_short.ipf(calibweights1)
#'
#' # The output can still be exported, for example to an Excel file
#' # library(openxlsx)
#' # write.xlsx(output_short, "SummaryShortIPF.xlsx")
#' }

summary_short.ipf <- function(object, ...) {
  terms <- variable <- NULL
  av <- attributes(object)
  w <- av$baseweight
  output <- list()
  
  # Extract variables from the formulas
  vars <- c()
  if (any(names(av) == "conP")) {
    for (i in seq_along(av$formP)) {
      vars <- c(vars, labels(terms(av$formP[[i]])))
    }
  }
  if (any(names(av) == "conH")) {
    for (i in seq_along(av$formH)) {
      vars <- c(vars, labels(terms(av$formH[[i]])))
    }
  }
  vars <- unique(vars)
  
  # Create formulas for original counts
  formPBase <- lapply(av$formP, function(x) {
    x <- as.character(x)
    as.formula(paste(w, "~", x[3]))
  })
  formHBase <- lapply(av$formH, function(x) {
    x <- as.character(x)
    as.formula(paste(w, "~", x[3]))
  })
  
  # Consistency check for person-level control tables
  if (any(names(av) == "conP")) {
    for (i in seq_along(av$conP)) {
      # Convert table objects to data.tables to get column names
      conP_dt <- as.data.table(av$conP[[i]])
      conP_adj_dt <- as.data.table(av$conP_adj[[i]])
      conP_orig_dt <- as.data.table(xtabs(formPBase[[i]], data = object))
      
      # Calculate the relative differences
      rel_diff_original <- as.data.table(round(
        100 * (conP_dt$N - conP_orig_dt$N) / conP_dt$N, 2
      ))
      rel_diff_calib <- as.data.table(round(
        100 * (conP_dt$N - conP_adj_dt$N) / conP_dt$N, 2
      ))
      
      # Restore the original column names and create the final summary table
      original_vars <- labels(terms(av$formP[[i]]))
      summary_dt <- cbind(conP_dt[, .SD, .SDcols = original_vars], rel_diff_original, rel_diff_calib)
      setnames(summary_dt, c(original_vars, "rel_diff_original", "rel_diff_calib"))
      
      output[[paste0("conP_", i, "_summary")]] <- summary_dt
    }
  }
  
  # Consistency check for household-level control tables
  if (any(names(av) == "conH")) {
    for (i in seq_along(av$conH)) {
      # Convert table objects to data.tables to get column names
      conH_dt <- as.data.table(av$conH[[i]])
      conH_adj_dt <- as.data.table(av$conH_adj[[i]])
      conH_orig_dt <- as.data.table(xtabs(formHBase[[i]], data = object))
      
      # Calculate the relative differences
      rel_diff_original <- as.data.table(round(
        100 * (conH_dt$N - conH_orig_dt$N) / conH_dt$N, 2
      ))
      rel_diff_calib <- as.data.table(round(
        100 * (conH_dt$N - conH_adj_dt$N) / conH_dt$N, 2
      ))
      
      # Restore the original column names and create the final summary table
      original_vars <- labels(terms(av$formH[[i]]))
      summary_dt <- cbind(conH_dt[, .SD, .SDcols = original_vars], rel_diff_original, rel_diff_calib)
      setnames(summary_dt, c(original_vars, "rel_diff_original", "rel_diff_calib"))
      
      output[[paste0("conH_", i, "_summary")]] <- summary_dt
    }
  }
  
  # Distribution of the weights (Coefficient of Variation, etc.)
  dist_gew <- list()
  if (!is.null(w)) {
    for (i in unique(vars)) {
      part_i <- object[, list(
        cv = sd(get(w)) / mean(get(w)),
        min = min(get(w)),
        max = max(get(w)),
        quotient = min(get(w)) / max(get(w))
      ), keyby = c(i)]
      # Check if `part_i` is not empty before calling setnames
      if (nrow(part_i) > 0) {
        setnames(part_i, i, paste0("value", 1:length(i)))
        part_i[, variable := paste(i, collapse = " x ")]
        setcolorder(part_i, c("variable", paste0("value", 1:length(i)), "cv", "min", "max", "quotient"))
        dist_gew <- c(dist_gew, list(part_i))
      }
    }
    if (length(dist_gew) > 0) {
      dist_gew <- rbindlist(dist_gew, use.names = TRUE, fill = TRUE)
      output[["distribution of weights"]] <- dist_gew
    }
  }
  
  return(output)
}