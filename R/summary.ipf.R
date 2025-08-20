#' Generate summary output for a ipf calibration

ipf_summary_calibres <- function(ipf_result) {
  # This function analyzes the results of an IPF (Iterative Proportional Fitting) calibration.It creates summary tables for each constraint used
  
  # 1. INITIAL SETUP AND DATA TYPE DETECTION
  # ----------------------------------------------------------------------
  # Determine the data type based on which calibrated weight column exists.
  # 'gew1' is for quarterly data, 'gew2' is for monthly.
  if ("gew1" %in% names(ipf_result)) {
    cat("Analyzing quarterly data...\n")
    calibWeightName <- "gew1"
    epsP_param_name <- "epsP_lfs"
    epsH_param_name <- "epsH_lfs"
  } else if ("gew2" %in% names(ipf_result)) {
    cat("Analyzing monthly data...\n")
    calibWeightName <- "gew2"
    epsP_param_name <- "epsP_lfs_m"
    epsH_param_name <- "epsH_lfs_m"
  } else {
    stop("Could not determine data type. 'gew1' or 'gew2' column is missing.")
  }
  
  # Extract the attributes of the IPF result object, which contain the constraint data.
  attributes_ipf <- attributes(ipf_result)
  
  # Initialize a list to store the final summary tables for each constraint.
  results_list <- list()
  
  # ----------------------------------------------------------------------
  # 2. PROCESSING THE PERSON CONSTRAINTS (conP)
  # ----------------------------------------------------------------------
  # Retrieve the list of person constraints and their corresponding formulas.
  conP_list <- attributes_ipf$conP
  formP_list <- attributes_ipf$formP
  
  # Get the tolerance ('epsilon') value for person constraints.
  epsP_value <- epsParameters[[epsP_param_name]][[1]]
  
  # Loop through each person constraint.
  for (i in 1:length(conP_list)) {
    # Assign a name to the current constraint for the final table.
    constraint_name <- names(conP_list)[i]
    if (is.null(constraint_name) || is.na(constraint_name) || nchar(constraint_name) == 0) {
      constraint_name <- paste0("conP", i)
    }
    
    # Identify the variables used for grouping from the constraint's dimension names.
    group_vars <- names(attr(conP_list[[i]], "dimnames"))
    
    # Create a copy of the data to perform calculations on.
    tmp_calib_data <- copy(ipf_result)
    # Convert grouping variables to character type to ensure proper merging later.
    tmp_calib_data[, (group_vars) := lapply(.SD, as.character), .SDcols = group_vars]
    
    # Determine the name of the numerical variable from the formula.
    # For a formula like 'eqIncome ~ gender', this extracts 'eqIncome'.
    rhs_vars <- all.vars(formP_list[[i]])
    numericalWeightingVar <- if (length(rhs_vars) > 1) rhs_vars[1] else NULL
    
    # Calculate the calibrated margins for the current constraint.
    calib_results <- tmp_calib_data[
      , .(
        # Count the number of individuals in the sample for each group.
        N = .N,
        CalibMargin = {
          # Check if a numerical weighting variable exists (e.g., eqIncome).
          # If so, calculate the sum of calibrated weights multiplied by that variable.
          # Otherwise, simply sum the calibrated weights.
          if (!is.null(numericalWeightingVar) && numericalWeightingVar %in% names(tmp_calib_data)) {
            sum(get(calibWeightName) * get(numericalWeightingVar))
          } else {
            sum(get(calibWeightName))
          }
        }
      ),
      by = group_vars
    ]
    
    # Process the original population totals ('PopMargin').
    original_xtabs <- conP_list[[i]]
    original_dt <- as.data.table(original_xtabs)
    # Rename columns
    setnames(original_dt, names(original_dt), c(group_vars, "PopMargin"))
    original_dt[, (group_vars) := lapply(.SD, as.character), .SDcols = group_vars]
    
    # Merge the calculated results with the original population totals.
    merged_results <- merge(calib_results, original_dt, by = group_vars, all = TRUE)
    # Fill any missing values with 0.
    merged_results[is.na(N), N := 0]
    merged_results[is.na(CalibMargin), CalibMargin := 0]
    
    # Add the epsilon value and calculate the maximum deviation factor.
    # maxFac indicates how far the calibrated margin is from the target margin.
    merged_results[, epsP := epsP_value]
    merged_results[, maxFac := abs(1 - CalibMargin / PopMargin)]
    
    # Filter out empty groups, order the table, and store it in the results list.
    final_table_p <- merged_results[N > 0][order(N)]
    results_list[[constraint_name]] <- final_table_p
  }
  
  # ----------------------------------------------------------------------
  # 3. PROCESSING THE HOUSEHOLD CONSTRAINTS (conH)
  # ----------------------------------------------------------------------
  # Retrieve the list of household constraints and their corresponding formulas.
  conH_list <- attributes_ipf$conH
  formH_list <- attributes_ipf$formH
  
  # Loop through each household constraint.
  for (i in 1:length(conH_list)) {
    # Assign a name to the current constraint or create a default name.
    constraint_name <- names(conH_list)[i]
    if (is.null(constraint_name) || is.na(constraint_name) || nchar(constraint_name) == 0) {
      constraint_name <- paste0("conH", i)
    }
    
    # Identify the grouping variables.
    group_vars <- names(attr(conH_list[[i]], "dimnames"))
    
    # Create a copy of the data and convert group variables to character.
    tmp_calib_data <- copy(ipf_result)
    tmp_calib_data[, (group_vars) := lapply(.SD, as.character), .SDcols = group_vars]
    
    # Calculate the calibrated margins for the current household constraint.
    calib_results <- tmp_calib_data[
      , .(
        N = .N,
        CalibMargin = {
          # Check if the household size variable ('wg') is one of the grouping variables.
          if ("wg" %in% group_vars) {
            # For household counts, convert person weights to household weights by dividing by the household size ('wg').
            sum(get(calibWeightName) / as.numeric(wg))
          } else {
            # It finds the numerical weighting variable and calculates the weighted sum.
            rhs_vars <- all.vars(formH_list[[i]])
            numericalWeightingVar <- if (length(rhs_vars) > 1) rhs_vars[2] else NULL
            
            if (!is.null(numericalWeightingVar) && numericalWeightingVar %in% names(tmp_calib_data)) {
              sum(get(calibWeightName) * as.numeric(get(numericalWeightingVar)))
            } else {
              sum(get(calibWeightName))
            }
          }
        }
      ),
      by = group_vars
    ]
    
    # Process and format the original population totals ('PopMargin').
    original_xtabs <- conH_list[[i]]
    original_dt <- as.data.table(original_xtabs)
    setnames(original_dt, names(original_dt), c(group_vars, "PopMargin"))
    original_dt[, (group_vars) := lapply(.SD, as.character), .SDcols = group_vars]
    
    # Merge the calculated results with the original population totals.
    merged_results <- merge(calib_results, original_dt, by = group_vars, all = TRUE)
    merged_results[is.na(N), N := 0]
    merged_results[is.na(CalibMargin), CalibMargin := 0]
    
    # Add epsilon values for each group and calculate the max deviation factor.
    epsH_obj <- epsParameters[[epsH_param_name]][[i]]
    
    # Handle the 'epsH' values, which can be stored in different formats (1D, 2D, or higher).
    if (length(group_vars) == 1) {
      epsH_dt_long <- data.table(
        `id_var` = names(epsH_obj),
        epsH = as.numeric(epsH_obj)
      )
      setnames(epsH_dt_long, "id_var", group_vars[1])
      
    } else if (length(group_vars) == 2) {
      id_vars_melt <- group_vars[1]
      epsH_dt <- as.data.table(epsH_obj, keep.rownames = id_vars_melt)
      epsH_dt_long <- melt(
        epsH_dt, 
        id.vars = id_vars_melt, 
        value.name = "epsH"
      )
      setnames(epsH_dt_long, "variable", group_vars[2])
      
    } else {
      epsH_dt_long <- as.data.table(as.table(epsH_obj))
      actual_dim_names <- names(epsH_dt_long)[1:(ncol(epsH_dt_long)-1)]
      setnames(epsH_dt_long, names(epsH_dt_long), c(actual_dim_names, "epsH"))
    }
    
    # This block handles the edge case where the epsilon table is empty.
    if (nrow(epsH_dt_long) > 0) {
      for (col in group_vars) {
        epsH_dt_long[[col]] <- as.character(epsH_dt_long[[col]])
      }
    } else {
      # Creates a new, empty data.table with the correct column names and types.
      empty_dt_list <- vector("list", length(group_vars) + 1)
      names(empty_dt_list) <- c(group_vars, "epsH")
      
      for(j in 1:length(group_vars)) {
        empty_dt_list[[j]] <- character(0)
      }
      empty_dt_list[["epsH"]] <- numeric(0)
      
      epsH_dt_long <- as.data.table(empty_dt_list)
    }
    
    # Merge the results with the epsilon values and calculate the max factor.
    merged_results <- merge(merged_results, epsH_dt_long, by = group_vars, all.x = TRUE)
    merged_results[, maxFac := abs(1 - CalibMargin / PopMargin)]
    
    # Finalize the table, filter, order, and store it.
    final_table_h <- merged_results[N > 0][order(N)]
    results_list[[constraint_name]] <- final_table_h
  }
  
  return(results_list)
}


#' @title Generate Summary Output for IPF Calibration
#'
#' @description Generates a detailed summary of an Iterative Proportional Fitting (IPF) calibration, providing a complete tool for evaluating the calibration's success and the validity of the resulting weights.
#'
#' The output is a list of data.tables for a comprehensive evaluation, including:
#'
#' **Calibration Results:**
#'   - `calib_results_conP_*` and `calib_results_conH_*`: Key diagnostic tables that compare calibrated margins to population targets and assess the goodness of fit via metrics like `maxFac`.
#'
#' **Data and Diagnostics:**
#'   - `weighted data`: An excerpt of the final dataset with the calculated calibration weights.
#'   - `distribution of the weights`: A statistical overview of the weight distribution (min, max, CV).
#'
#' **Detailed Margin Comparisons:**
#'   - `conP_*`, `conH_*`, `*_adjusted`, `*_original`, `*_rel_diff_*`: Tables that compare original sample margins, calibrated margins, and population targets, along with their relative differences.
#'
#' @param object object of class ipf
#' @param ... additional arguments
#'
#' @return a list of the following outputs
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
#'  eusilc,
#'  conP = list(conP1, conP2, eqIncome = conP3),
#'  bound = NULL,
#'  verbose = TRUE
#' )
#' output <- summary(calibweights1)
#' # the output can easily be exported to an Excel file, e.g. with
#' # library(openxlsx)
#' # write.xlsx(output, "SummaryIPF.xlsx")
#' }

summary.ipf <- function(object, ...){
  terms <- variable <- NULL
  dots <- list(...)
  extraExport <- dots$extraExport
  av <- attributes(object)
  w <- av$baseweight
  hid <- av$hid
  
  # --- Start of new table integration ---
  # Calls the generate_ipf_summary function to create the initial tables.
  custom_summary <- ipf_summary_calibres(object)
  
  # Renames the list appropriately to distinguish it from the other tables.
  names(custom_summary) <- paste0("calib_results_", names(custom_summary))
  
  # Initializes the final output with the new tables.
  output <- custom_summary
  # --- End of integration ---
  
  ## ---- FormP/FormH - Variables etc. ---- ##
  vars <- c()
  formulas <- c()
  if(any(names(av) == "conP")){
    calibWeightName <- as.character(av$formP[1][[1]])[2]
    varsP_help <- c()
    for(i in seq_along(av$formP)){
      varsP_help <- c(varsP_help,labels(terms(av$formP[[i]])))
    }
    vars <- c(vars,unique(varsP_help))
    
    #
    formPBase <- lapply(av$formP,function(x){
      x <- as.character(x)
      as.formula(paste(w,"~",x[3]))
    })
    #
    formulas <- c(formulas,as.character(av$formP))
    
  }
  
  if(any(names(av) == "conH")){
    calibWeightName <- as.character(av$formH[1][[1]])[2]
    varsH_help <- c()
    for(i in seq_along(av$formH)){
      varsH_help <- c(varsH_help,labels(terms(av$formH[[i]])))
    }
    vars <- c(vars,unique(varsH_help))
    
    #
    formHBase <- lapply(av$formH,function(x){
      x <- as.character(x)
      as.formula(paste(w,"~",x[3]))
    })
    #
    formulas <- c(formulas,as.character(av$formH))
    
  }
  vars <- as.list(unique(vars))
  
  dist_gew <- list()
  if(!is.null(w)){
    for(i in vars){
      part_i <- object[,list(cv=sd(get(w))/mean(get(w)),
                             min=min(get(w)),
                             max=max(get(w)),
                             quotient=min(get(w))/max(get(w))),keyby=c(i)]
      setnames(part_i,i,paste0("value",1:length(i)))
      part_i[,variable:=paste(i,collapse=" x ")]
      setcolorder(part_i,c("variable",paste0("value",1:length(i)),"cv","min","max","quotient"))
      dist_gew <- c(dist_gew,list(part_i))
    }
    dist_gew <- rbindlist(dist_gew,use.names=TRUE,fill=TRUE)
    setcolorder(dist_gew,c("variable",paste0("value",1:length(i)),"cv","min","max","quotient"))
  }
  
  
  cols <- c(which(colnames(object) %in% c(hid,vars,w,calibWeightName,extraExport)))
  
  # EXTENSION: The existing output is created here.
  original_output <- vector(mode = "list", length = 3)
  names(original_output)[1] <- "weighted data"
  original_output[[1]] <- subset(object,select=cols)
  
  formulaOutput <- data.table(Verteilungen=rep(NA,length(formulas)))
  formulaOutput$Verteilungen[seq_along(formulas)] <- formulas
  
  names(original_output)[2] <- "margin tables"
  original_output[[2]] <- formulaOutput
  
  names(original_output)[3] <- "distribution of the weights"
  original_output[[3]] <- dist_gew
  
  #### ---- ####
  if (any(names(av) == "conP")) {
    for (i in seq_along(av$conP)) {
      output4 <- list(
        "conP_%i" = as.data.table(av$conP[[i]]),
        "conP_%i_adjusted" = as.data.table(av$conP_adj[[i]]),
        "conP_%i_original" = as.data.table(xtabs(formPBase[[i]], data = object)),
        "conP_%i_rel_diff_original" = as.data.table(round(
          100*(av$conP[[i]] - xtabs(formPBase[[i]], data = object)) / av$conP[[i]], 2)),
        "conP_%i_rel_diff_calib" = as.data.table(round(
          100*(av$conP[[i]] - av$conP_adj[[i]]) / av$conP[[i]], 2))
      )
      names(output4) <- sprintf(names(output4), i)
      original_output <- append(original_output, output4)
    }
  }
  
  
  if(any(names(av) == "conH")){
    output5 <- vector(mode = "list", length = length(av$conH)*5)
    for(i in seq_along(av$conH)){
      j <- ((i-1)*5)+1
      names(output5)[(j):(j+4)] <- c(paste0("conH_",i),paste0("conH_",i,"_adjusted"),paste0("conH_",i,"_original"),
                                     paste0("conH_",i,"_rel_diff_original"),paste0("conH_",i,"_rel_diff_calib"))
      output5[[j]] <- as.data.table(av$conH[[i]]) #conH_i
      output5[[j+1]] <- as.data.table(av$conH_adj[[i]]) #conH_i_adjusted
      output5[[j+2]] <- as.data.table(xtabs(formHBase[[i]],data=object)) #conH_i_original
      output5[[j+3]] <- as.data.table(round(100*(av$conH[[i]]-xtabs(formHBase[[i]],data=object))/av$conH[[i]],2)) #conH_i_rel_diff_original
      output5[[j+4]] <- as.data.table(round(100*(av$conH[[i]]-av$conH_adj[[i]])/av$conH[[i]],2)) #conH_i__rel_diff_calib
      
    }
    original_output <- append(original_output,output5)
  }
  
  # Appends the original tables to the custom tables.
  output <- append(output, original_output)
  
  return(output)
}
