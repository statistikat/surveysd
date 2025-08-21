print.summary.ipf <- function(x, ...) {
  cat("IPF Calibration Summary (Short Version)\n")
  cat("=======================================\n\n")
  
  # --- 1. CalibWeight automatisch aus all_formulas extrahieren ---
  if (!is.null(x$all_formulas)) {
    calibWeightName <- trimws(strsplit(x$all_formulas[1], "~")[[1]][1])
  } else {
    calibWeightName <- NULL
  }
  
  # --- 2. Kish-Faktor berechnen ---
  weight_candidates <- names(x[["weighted data"]])
  if (!is.null(calibWeightName) && calibWeightName %in% weight_candidates) {
    w <- x[["weighted data"]][[calibWeightName]]
    kish <- kishFactor(w, na.rm = TRUE)
    cat(sprintf("Kish Factor (basierend auf '%s'): %.4f\n\n", calibWeightName, kish))
  } else {
    cat("Kish Factor: Keine gültige Gewichtsspalte gefunden\n\n")
  }
  
  # --- 3. calib_results_ Tabellen anzeigen (erste 10 Zeilen) ---
  calib_names <- grep("^calib_results_", names(x), value = TRUE)
  if (length(calib_names) > 0) {
    cat("Calibration Results (erste 10 Zeilen pro Tabelle):\n")
    for (nm in calib_names) {
      cat("\n---", nm, "---\n")
      print(utils::head(x[[nm]], 10))
    }
  } else {
    cat("Keine calib_results_-Tabellen gefunden.\n")
  }
  
  cat("\n(Weitere Details mit str() oder names() ansehen)\n")
  invisible(x)
}
