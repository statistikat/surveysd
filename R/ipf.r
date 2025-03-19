#' @rdname ipf_step
#' @name ipf_step
#' @md
#'
#' @param dat a `data.frame` containing the factor variables to be combined.
#'
#' @export
#' 

combine_factors <- function(dat, targets) {

  x <- as.data.frame(targets)
  x$ID_ipu <- seq_len(nrow(x))
  x <- merge(dat, x, by = names(dimnames(targets)), sort = FALSE, all.x = TRUE)
  factor(x$ID_ipu, levels = seq_along(targets))
}

getMeanFun <- function(meanHH) {
  if (isTRUE(meanHH))
    meanHH <- "arithmetic"
  if (identical(meanHH, FALSE))
    meanHH <- "none"
  meanfun <- switch(meanHH,
                    arithmetic = arithmetic_mean,
                    geometric = geometric_mean,
                    none = function(x, w) {
                      x
                    }
  )
  if (is.null(meanfun))
    stop("invalid value for meanHH")
  meanfun
}
boundsFak <- function(g1, g0, f, bound = 4, minMaxTrim=NULL) {
  # Berechnet die neuen Gewichte (innerhalb 4, .25 Veraenderungsraten)
  g1 <- g1 * f
  if(!is.null(bound)){
    TF <- which((g1 / g0) > bound)
    TF[is.na(TF)] <- FALSE
    g1[TF] <- bound * g0[TF]
    TF <- which((g1 / g0) < (1 / bound))
    TF[is.na(TF)] <- FALSE
    g1[TF] <- (1 / bound) * g0[TF]
  }
  if(!is.null(minMaxTrim)){
    g1[g1<minMaxTrim[1]] <- minMaxTrim[1]
    g1[g1>minMaxTrim[2]] <- minMaxTrim[2]
  }
  return(g1)
}
boundsFakHH <- function(g1, g0, eps, orig, p, bound = 4, minMaxTrim=NULL) {
  # Berechnet die neuen Gewichte fuer Unter- und Obergrenze (innerhalb 4,
  #   .25 Veraenderungsraten)
  u <- orig * (1 - eps)
  o <- orig * (1 + eps)

  pbo <- which(p > o)
  psu <- which(p < u)
  g1[pbo] <- g1[pbo] * o[pbo] / p[pbo]
  g1[psu] <- g1[psu] * u[psu] / p[psu]
  if(!is.null(bound)){
    TF <- which((g1 / g0) > bound)
    g1[TF] <- bound * g0[TF]
    TF <- which((g1 / g0) < (1 / bound))
    g1[TF] <- (1 / bound) * g0[TF]
  }
  if(!is.null(minMaxTrim)){
    g1[g1<minMaxTrim[1]] <- minMaxTrim[1]
    g1[g1>minMaxTrim[2]] <- minMaxTrim[2]
  }
  return(g1)
}

check_eps <- function(con, eps, type){
  l <- length(con)
  if(is.list(eps)){
    if(length(eps)!=l){
      stop(paste("Provided",type,"eps argument does not fit",type,"constraints."))
    }
    for(i in 1:length(eps)){
      if(is.array(eps[[i]])){
        if(!identical(dim(eps[[i]]),dim(con[[i]]))){
          stop(paste("Provided",type,"eps argument",i,"does not fit in dimension to ",type,"constraints",i,"."))
        }
      }
    }
  }else if(length(eps)>1){
    stop("Individual eps arguments for each constraints must be defined as list.")
  }
}

check_population_totals <- function(con, dat, type = "personal") {
  # do not apply this check for numerical calibration
  if (is.null(names(con))) {
    ind <- seq_along(con)
  } else {
    ind <- which(names(con) == "")
  }

  # do not apply this check for constraints that only cover the population
  #   partially
  ind <- ind[vapply(
    ind,
    function(i) {
      constraint <- con[[i]]
      for (variable in names(dimnames(constraint))) {
        if (!(variable %in% names(dat))){
          stop("variable ", variable, " appears in a constraint but not ",
               "in the dataset")
        }
        vals_dat <- sort(unique(as.character(dat[[variable]])))
        vals_con <- sort(dimnames(constraint)[[variable]])
        if (any(!vals_dat%in%vals_con)) {
          message(type, " constraint ", i, " only covers a subset of the ",
                  "population")
          return(FALSE)
        }
        if(any(!vals_con%in%vals_dat)){
          message(type, " constraint ", i, " covers more combinations of ",
                  paste(names(dimnames(constraint)),collapse=" and "),
                  " then appear in dat")
          return(FALSE)
        }
      }
      return(TRUE)
    },
    TRUE)]

  if (length(ind) == 0)
    return(NULL)

  pop_totals <- vapply(
    ind,
    function(index) {
      sum(con[[index]])
    },
    0)
  rel_errors <- abs(pop_totals - pop_totals[1]) / pop_totals[1]

  # use a 1% tolerance. Maybe it would be better to make the tolerance
  #   dependent on conH?
  if (any(rel_errors > 1e-2))
    stop("population totals for different constraints do not match")
}

calibP <- function(i, dat, error, valueP, pColNames, bound, verbose, calIter,
                   numericalWeighting, numericalWeightingVar, w,
                   cw, minMaxTrim, print_every_n) {
  selectGroupNotConverged <- epsPcur <- maxFac <- OriginalSortingVariable <- V1 <-
    epsvalue <- fVariableForCalibrationIPF <- NULL
  temporary_hvar <- value <-
    wValue <- representativeHouseholdForCalibration <- NULL
  variableKeepingTheBaseWeight <- w
  variableKeepingTheCalibWeight <- cw
  combined_factors <- dat[[paste0("combined_factors_", i)]]
  setnames(dat, valueP[i], "value")
  setnames(dat, paste0("epsP_", i), "epsPcur")
  tmp <- data.table(x = factor(levels(combined_factors)))
  setnames(tmp, "x", paste0("combined_factors_", i))
  con_current <- dat[tmp, on = paste0("combined_factors_", i),
                     mult = "first", value]

  if (!is.null(numericalWeightingVar)) {
    ## numerical variable to be calibrated
    ## use name of conP list element to define numerical variable
    set(dat, j = "fVariableForCalibrationIPF",
        value = ipf_step_f(dat[[variableKeepingTheCalibWeight]] *
                             dat[[numericalWeightingVar]],
                           combined_factors, con_current))
    set(dat, j = "wValue", value = dat[["value"]] /
          dat[["fVariableForCalibrationIPF"]])

    # try to divide the weight between units with larger/smaller values in the
    #   numerical variable linear
    set(dat, j="fVariableForCalibrationIPF", value=1)
    dat[!is.na(value) & value!=0, fVariableForCalibrationIPF := numericalWeighting(
      head(wValue, 1), head(value, 1), get(numericalWeightingVar),
      get(variableKeepingTheCalibWeight)),
      by=c(paste0("combined_factors_", i))]
  
    # adjust weights to hit target 
    # weights might shift due  to applying boundary limits in numerialWeighting()
    
    dat[!is.na(value) & value!=0, fVariableForCalibrationIPF := fVariableForCalibrationIPF*value/sum(fVariableForCalibrationIPF*get(numericalWeightingVar)*get(variableKeepingTheCalibWeight)), by=c(paste0("combined_factors_", i))]
    
    # # do this only if wValue is bigger than target
    # # dat[!is.na(value) & value!=0, value_too_large := value < sum(fVariableForCalibrationIPF*get(numericalWeightingVar)*get(variableKeepingTheCalibWeight)), by=c(paste0("combined_factors_", i))]
    # # dat[!is.na(value) & value!=0 & value_too_large == TRUE, fVariableForCalibrationIPF := fVariableForCalibrationIPF*value/sum(fVariableForCalibrationIPF*get(numericalWeightingVar)*get(variableKeepingTheCalibWeight)), by=c(paste0("combined_factors_", i))]
    
    # result after applying factor
    # dat[!is.na(value) & value!=0, wValue:=sum(fVariableForCalibrationIPF*get(variableKeepingTheCalibWeight)*get(numericalWeightingVar)), by = c(paste0("combined_factors_", i))]
  } else {
    # categorical variable to be calibrated
    set(dat, j = "fVariableForCalibrationIPF", value = ipf_step_f(
      dat[[variableKeepingTheCalibWeight]], combined_factors, con_current))
    set(dat, j = "wValue", value = dat[["value"]] /
          dat[["fVariableForCalibrationIPF"]])
  }
  
  
  dat[, selectGroupNotConverged := (abs(wValue-value)/value)>epsPcur]
  dat[is.na(selectGroupNotConverged),selectGroupNotConverged:=FALSE]
  
  if (dat[,any(selectGroupNotConverged)]) {
    ## sicherheitshalber abs(epsPcur)? Aber es wird schon niemand negative eps
    ##   Werte uebergeben??
    if (verbose && (calIter %% 10 == 0| calIter %% print_every_n == 0)) {
      message(calIter, ":Not yet converged for P-Constraint", i, "\n")
      if (calIter %% print_every_n == 0) {

        tmp <- dat[
          selectGroupNotConverged == TRUE,
          list(
            maxFac = max(abs(1 / fVariableForCalibrationIPF - 1)), .N,
            epsP = head(epsPcur, 1),
            CalibMargin = {
              if (!is.null(numericalWeightingVar)) {
                sum(get(variableKeepingTheCalibWeight) *
                      get(numericalWeightingVar))
              }else{
                sum(get(variableKeepingTheCalibWeight))
              }
            },
            PopMargin = head(value, 1)),
          by = eval(pColNames[[i]])]
       
        print(tmp[order(maxFac, decreasing = TRUE), ])
        message("-----------------------------------------\n")
      }
    }
    
    # variableKeepingTheCalibWeight_helper <- paste0(variableKeepingTheCalibWeight,"_helper")
    # dat[,c(variableKeepingTheCalibWeight_helper):=get(variableKeepingTheCalibWeight)]
    
    if (!is.null(bound) || !is.null(minMaxTrim)) {
      dat[!is.na(fVariableForCalibrationIPF),
          c(variableKeepingTheCalibWeight) :=
            boundsFak(
              get(variableKeepingTheCalibWeight),
              get(variableKeepingTheBaseWeight), fVariableForCalibrationIPF,
              bound = bound, minMaxTrim = minMaxTrim)]
      #,by=eval(pColNames[[i]])]
    } else {
      dat[!is.na(fVariableForCalibrationIPF),
          c(variableKeepingTheCalibWeight) := fVariableForCalibrationIPF *
            get(variableKeepingTheCalibWeight),
          by = eval(paste0("combined_factors_", i))]
    }
    
    # # accept new solution only if it gets closer to target -> computeLinearG1 can diverge in some cases
    # dat[,accept_new:=abs(sum(get(numericalWeightingVar)*get(variableKeepingTheCalibWeight_helper))-value)<
    #       abs(sum(get(numericalWeightingVar)*get(variableKeepingTheCalibWeight))-value),by=eval(paste0("combined_factors_", i))]
    # # dont accept new solution if groups which have already converged and where solution will get worse
    # dat[!(accept_new==FALSE & selectGroupNotConverged==TRUE),c(variableKeepingTheCalibWeight):=get(variableKeepingTheCalibWeight_helper)]
    # dat[,c(variableKeepingTheCalibWeight_helper,"accept_new"):=NULL]
    dat[,c("fVariableForCalibrationIPF"):=NULL]
    error <- TRUE
  }
  dat[,c("selectGroupNotConverged"):=NULL]
  setnames(dat, "value", valueP[i])
  setnames(dat, "epsPcur", paste0("epsP_", i))
  
  # print(dat[combined_factors_1==144,.(kz,calibWeight,fVariableForCalibrationIPF,calibWeight*fVariableForCalibrationIPF)])
  
  # return(list(error=error,dat=dat))
  return(error)
}

calibH <- function(i, dat, error, valueH, hColNames, bound, verbose, calIter,
                   looseH, numericalWeighting, numericalWeightingVar,
                   w, cw, minMaxTrim, print_every_n) {
  variableKeepingTheBaseWeight <- w
  variableKeepingTheCalibWeight <- cw
  selectGroupNotConverged <- epsHcur <- OriginalSortingVariable <- V1 <-
    epsvalue <- fVariableForCalibrationIPF <- NULL
  maxFac <- temporary_hvar <-
    value <- wValue <- representativeHouseholdForCalibration <- NULL

  setnames(dat, valueH[i], "value")
  setnames(dat, paste0("epsH_", i), "epsHcur")

  combined_factors <- dat[[paste0("combined_factors_h_", i)]]
  tmp <- data.table(x = factor(levels(combined_factors)))
  setnames(tmp, "x", paste0("combined_factors_h_", i))
  paste0("combined_factors_h_", i)
  con_current <- dat[tmp, on = paste0("combined_factors_h_", i),
                     mult = "first", value]
  if (!is.null(numericalWeightingVar)) {
    ## numerical variable to be calibrated
    ## use name of conH list element to define numerical variable
    set(dat, j = "fVariableForCalibrationIPF", value = ipf_step_f(
      dat[[variableKeepingTheCalibWeight]] *
        dat[["representativeHouseholdForCalibration"]] *
      dat[[numericalWeightingVar]], combined_factors, con_current))
    set(dat, j = "wValue", value = dat[["value"]] /
          dat[["fVariableForCalibrationIPF"]])

    # try to divide the weight between units with larger/smaller value in the
    #   numerical variable linear
    set(dat, j="fVariableForCalibrationIPF", value=1)
    dat[!is.na(value) & value!=0, fVariableForCalibrationIPF := numericalWeighting(
      head(wValue, 1), head(value, 1), get(numericalWeightingVar),
      get(variableKeepingTheCalibWeight)),
      by = eval(paste0("combined_factors_h_", i))]
  } else {
    # categorical variable to be calibrated
    set(dat, j = "fVariableForCalibrationIPF", value = ipf_step_f(
      dat[[variableKeepingTheCalibWeight]] *
        dat[["representativeHouseholdForCalibration"]],
      combined_factors, con_current))
  }

  set(dat, j = "wValue", value = dat[["value"]] /
        dat[["fVariableForCalibrationIPF"]])

  dat[, selectGroupNotConverged := (abs(wValue-value)/value)>epsHcur]
  dat[is.na(selectGroupNotConverged),selectGroupNotConverged:=FALSE]
  
  if (dat[,any(selectGroupNotConverged)]) {
    if (verbose && (calIter %% 10 == 0| calIter %% print_every_n == 0)) {
      message(calIter, ":Not yet converged for H-Constraint", i, "\n")
      if (calIter %% print_every_n == 0) {
        tmp <- dat[selectGroupNotConverged == TRUE,
          list(maxFac = max(abs(1 / fVariableForCalibrationIPF - 1)), .N,
               epsH = head(epsHcur, 1),
               sumCalibWeight = sum(get(variableKeepingTheCalibWeight) *
                                      representativeHouseholdForCalibration),
               PopMargin = head(value, 1)),
          by = eval(hColNames[[i]])]
        print(tmp[order(maxFac, decreasing = TRUE), ])

        message("-----------------------------------------\n")
      }
    }
    if (!is.null(bound) || !is.null(minMaxTrim)) {
      if (!looseH) {
        set(dat, j = variableKeepingTheCalibWeight, value = boundsFak(
          g1 = dat[[variableKeepingTheCalibWeight]],
          g0 = dat[[variableKeepingTheBaseWeight]],
          f = dat[["fVariableForCalibrationIPF"]],
          bound = bound, minMaxTrim = minMaxTrim))
      }else{
        set(dat, j = variableKeepingTheCalibWeight, value = boundsFakHH(
          g1 = dat[[variableKeepingTheCalibWeight]],
          g0 = dat[[variableKeepingTheBaseWeight]],
          eps = dat[["epsHcur"]], orig = dat[["value"]],
          p = dat[["wValue"]], bound = bound, minMaxTrim = minMaxTrim)
        )
      }
    } else {
      dat[, c(variableKeepingTheCalibWeight) := fVariableForCalibrationIPF *
            get(variableKeepingTheCalibWeight),
          by = eval(paste0("combined_factors_h_", i))]
    }
    dat[,c("fVariableForCalibrationIPF"):=NULL]
    error <- TRUE
  }

  dat[,c("selectGroupNotConverged"):=NULL]
  setnames(dat, "value", valueH[i])
  setnames(dat, "epsHcur", paste0("epsH_", i))
  return(error)
}

## recreate the formula argument to xtabs based on conP, conH
getFormulas <- function(con, w) {
  formOut <- NULL
  for (i in seq_along(con)) {
    lhs <- names(con)[i]
    if (is.null(lhs) || lhs == "") {
      lhs <- w
    } else {
      lhs <- paste(lhs, "*", w)
    }
    rhs <- paste(names(dimnames(con[[i]])), collapse = "+")
    formOut[[i]] <- formula(paste(lhs, "~", rhs), env = .GlobalEnv)
  }
  formOut
}

## enrich dat_original with the calibrated weights and assign attributes

addWeightsAndAttributes <- function(dat, conP, conH, epsP, epsH, dat_original,
                                    maxIter, calIter, returnNA, cw, bw,verbose, looseH, hidVar = NULL) {

  variableKeepingTheCalibWeight <- cw
  representativeHouseholdForCalibration <- OriginalSortingVariable <-
    outTable <- copy(dat_original)

  formP <- getFormulas(conP, w = variableKeepingTheCalibWeight)
  formH <- getFormulas(conH, w = variableKeepingTheCalibWeight)

  # general information
  setattr(outTable, "iterations", min(maxIter, calIter))

  # input constraints
  setattr(outTable, "conP", conP)
  setattr(outTable, "conH", conH)

  # adjusted constraints (conP, conH according to the calibrated weights)
  conP_adj <- lapply(formP, xtabs, dat)
  conH_adj <- lapply(
    formH, xtabs, dat[representativeHouseholdForCalibration == 1])
  setattr(outTable, "conP_adj", conP_adj)
  setattr(outTable, "conH_adj", conH_adj)

  # tolerances
  setattr(outTable, "epsP", epsP)
  setattr(outTable, "epsH", epsH)

  setkey(dat, OriginalSortingVariable)


  # convergence
  conP_converged <- sapply(seq_along(conP), function(i) {
    epsP_current <- switch(is.list(epsP) + 1, epsP, epsP[[i]])
    conP_zero <- conP[[i]]==0
    cond1 <- abs(conP[[i]] - conP_adj[[i]]) <= epsP_current * conP[[i]]
    all(conP_zero==TRUE | cond1==TRUE)
  })
  if(looseH){

    inc <- function(x){
      10^(floor(log10(x))-1) * 9.99
    }

    conH_converged <- sapply(seq_along(conH), function(i) {
      epsH_current <- switch(is.list(epsH) + 1, epsH, epsH[[i]])
      if(verbose){
        message(paste("For looseH=TRUE epsH+",inc(epsH_current),"is allowed as tolerance for the convergence."))
      }
      conH_zero <- conH[[i]] ==0
      cond1 <- abs(conH[[i]] - conH_adj[[i]])/conH[[i]] <= (epsH_current + inc(epsH_current))
      all(conH_zero==TRUE | cond1==TRUE)
    })
  }else{
    conH_converged <- sapply(seq_along(conH), function(i) {
      epsH_current <- switch(is.list(epsH) + 1, epsH, epsH[[i]])
      conH_zero <- conH[[i]] ==0
      cond1 <- abs(conH[[i]] - conH_adj[[i]]) <= epsH_current * conH[[i]]
      all(conH_zero==TRUE | cond1==TRUE)
    })
  }
  converged <- all(conP_converged) && all(conH_converged)
  setattr(outTable, "converged", converged)
  if (verbose) {
    if (converged)
      message("Convergence reached")
    else
      message("No convergence reached")
  }

  # add calibrated weights. Use setkey to make sure the indexes match
  setkey(dat, OriginalSortingVariable)

  if (!converged & returnNA) {
    outTable[, c(variableKeepingTheCalibWeight) := NA]
  } else {
    outTable[, c(variableKeepingTheCalibWeight) :=
               dat[[variableKeepingTheCalibWeight]]]
  }

  # formulas
  setattr(outTable, "formP", formP)
  setattr(outTable, "formH", formH)
  setattr(outTable, "baseweight", bw)
  setattr(outTable, "hid", hidVar)
  # for the summary
  class(outTable) <- c("ipf",class(outTable))

  invisible(outTable)
}


#' Iterative Proportional Fitting
#'
#' Adjust sampling weights to given totals based on household-level and/or
#' individual level constraints.
#'
#' This function implements the weighting procedure described
#' here: \doi{10.17713/ajs.v45i3.120}.
#' Usage examples can be found in the corresponding vignette
#' (`vignette("ipf")`).
#'
#' `conP` and `conH` are contingency tables, which can be created with `xtabs`.
#' The `dimnames` of those tables should match the names and levels of the
#' corresponding columns in `dat`.
#'
#' `maxIter`, `epsP` and `epsH` are the stopping criteria. `epsP` and `epsH`
#' describe relative tolerances in the sense that
#' \deqn{1-epsP < \frac{w_{i+1}}{w_i} < 1+epsP}{1-epsP < w(i+1)/w(i) < 1+epsP}
#' will be used as convergence criterium. Here i is the iteration step and wi is
#' the weight of a specific person at step i.
#'
#' The algorithm
#' performs best if all varables occuring in the constraints (`conP` and `conH`)
#' as well as the household variable are coded as `factor`-columns in `dat`.
#' Otherwise, conversions will be necessary which can be monitored with the
#' `conversion_messages` argument. Setting `check_hh_vars` to `FALSE` can also
#' incease the performance of the scheme.
#'
#' @name ipf
#' @md
#' @aliases ipf
#' @param dat a `data.table` containing household ids (optionally), base
#'   weights (optionally), household and/or personal level variables (numerical
#'   or categorical) that should be fitted.
#' @param hid name of the column containing the household-ids within `dat` or
#'   NULL if such a variable does not exist.
#' @param w name if the column containing the base weights within `dat` or NULL
#'   if such a variable does not exist. In the latter case, every observation
#'   in `dat` is assigned a starting weight of 1.
#' @param conP list or (partly) named list defining the constraints on person
#'   level.  The list elements are contingency tables in array representation
#'   with dimnames corresponding to the names of the relevant calibration
#'   variables in `dat`. If a numerical variable is to be calibrated, the
#'   respective list element has to be named with the name of that numerical
#'   variable. Otherwise the list element shoud NOT be named.
#' @param conH list or (partly) named list defining the constraints on
#'   household level.  The list elements are contingency tables in array
#'   representation with dimnames corresponding to the names of the relevant
#'   calibration variables in `dat`. If a numerical variable is to be
#'   calibrated, the respective list element has to be named with the name of
#'   that numerical variable. Otherwise the list element shoud NOT be named.
#' @param epsP numeric value or list (of numeric values and/or arrays)
#'   specifying the convergence limit(s) for `conP`. The list can contain
#'   numeric values and/or arrays which must appear in the same order as the
#'   corresponding constraints in `conP`. Also, an array must have the same
#'   dimensions and dimnames as the corresponding constraint in `conP`.
#' @param epsH numeric value or list (of numeric values and/or arrays)
#'   specifying the convergence limit(s) for `conH`. The list can contain
#'   numeric values and/or arrays which must appear in the same order as the
#'   corresponding constraints in `conH`. Also, an array must have the same
#'   dimensions and dimnames as the corresponding constraint in `conH`.
#' @param verbose if TRUE, some progress information will be printed.
#' @param bound numeric value specifying the multiplier for determining the
#'   weight trimming boundary if the change of the base weights should be
#'   restricted, i.e. if the weights should stay between 1/`bound`*`w`
#'   and `bound`*\code{w}.
#' @param minMaxTrim numeric vector of length2, first element a minimum value
#'   for weights to be trimmed to, second element a maximum value for weights to
#'   be trimmed to.
#' @param maxIter numeric value specifying the maximum number of iterations
#' that should be performed.
#' @param meanHH if TRUE, every person in a household is assigned the mean of
#'   the person weights corresponding to the household. If `"geometric"`, the
#'   geometric mean is used rather than the arithmetic mean.
#' @param allPthenH if TRUE, all the person level calibration steps are
#'   performed before the houshold level calibration steps (and `meanHH`, if
#'   specified). If FALSE, the houshold level calibration steps (and `meanHH`,
#'   if specified) are performed after everey person level calibration step.
#'   This can lead to better convergence properties in certain cases but also
#'   means that the total number of calibration steps is increased.
#' @param returnNA if TRUE, the calibrated weight will be set to NA in case of
#'   no convergence.
#' @param looseH if FALSE, the actual constraints `conH` are used for
#'   calibrating all the hh weights. If TRUE, only the weights for which the
#'   lower and upper thresholds defined by `conH` and `epsH` are exceeded are
#'   calibrated. They are however not calibrated against the actual constraints
#'   `conH` but against these lower and upper thresholds, i.e.
#'   `conH`-`conH`*`epsH` and `conH`+`conH`*\code{epsH}.
#' @param numericalWeighting See [numericalWeighting]
#' @param check_hh_vars If `TRUE` check for non-unique values inside of a
#'   household for variables in household constraints
#' @param conversion_messages show a message, if inputs need to be reformatted.
#'   This can be useful for speed optimizations if ipf is called several times
#'   with similar inputs (for example bootstrapping)
#' @param nameCalibWeight character defining the name of the variable for the
#'   newly generated calibrated weight.
#' @param print_every_n number of interation steps after which a summary table
#'   is printed. The summary table shows all constraints which are not yet
#'   reached according to `epsP` and `epsH`
#' @return The function will return the input data `dat` with the calibrated
#'   weights `calibWeight` as an additional column as well as attributes. If no
#'   convergence has been reached in `maxIter` steps, and `returnNA` is `TRUE`
#'   (the default), the column `calibWeights` will only consist of `NA`s. The
#'   attributes of the table are attributes derived from the `data.table` class
#'   as well as the following.
#' \tabular{ll}{
#'   `converged` \tab Did the algorithm converge in `maxIter` steps? \cr
#'   `iterations` \tab The number of iterations performed. \cr
#'   `conP`, `conH`, `epsP`, `epsH` \tab See Arguments. \cr
#'   `conP_adj`, `conH_adj` \tab Adjusted versions of `conP` and `conH` \cr
#'   `formP`, `formH` \tab Formulas that were used to calculate `conP_adj` and
#'   `conH_adj` based on the output table.
#' }
#' @export ipf
#' @author Alexander Kowarik, Gregor de Cillia
#' @examples
#' \dontrun{
#'
#' # load data
#' eusilc <- demo.eusilc(n = 1, prettyNames = TRUE)
#'
#' # personal constraints
#' conP1 <- xtabs(pWeight ~ age, data = eusilc)
#' conP2 <- xtabs(pWeight ~ gender + region, data = eusilc)
#' conP3 <- xtabs(pWeight*eqIncome ~ gender, data = eusilc)
#'
#' # household constraints
#' conH1 <- xtabs(pWeight ~ hsize + region, data = eusilc[!duplicated(hid)])
#'
#' # simple usage ------------------------------------------
#'
#' calibweights1 <- ipf(
#'   eusilc,
#'   conP = list(conP1, conP2, eqIncome = conP3),
#'   bound = NULL,
#'   verbose = TRUE
#' )
#'
#' # compare personal weight with the calibweigth
#' calibweights1[, .(hid, pWeight, calibWeight)]
#'
#' # advanced usage ----------------------------------------
#'
#' # use an array of tolerances
#' epsH1 <- conH1
#' epsH1[1:4, ] <- 0.005
#' epsH1[5, ] <- 0.2
#'
#' # create an initial weight for the calibration
#' eusilc[, regSamp := .N, by = region]
#' eusilc[, regPop := sum(pWeight), by = region]
#' eusilc[, baseWeight := regPop/regSamp]
#'
# calibrate
#' calibweights2 <- ipf(
#'   eusilc,
#'   conP = list(conP1, conP2),
#'   conH = list(conH1),
#'   epsP = 1e-6,
#'   epsH = list(epsH1),
#'   bound = 4,
#'   w = "baseWeight",
#'   verbose = TRUE
#' )
#'
#' # show an adjusted version of conP and the original
#' attr(calibweights2, "conP_adj")
#' attr(calibweights2, "conP")
#' }
ipf <- function(
  dat, hid = NULL, conP = NULL, conH = NULL, epsP = 1e-6, epsH = 1e-2,
  verbose = FALSE, w = NULL, bound = 4, maxIter = 200, meanHH = TRUE,
  allPthenH = TRUE, returnNA = TRUE, looseH = FALSE, numericalWeighting =
    computeLinear, check_hh_vars = TRUE, conversion_messages = FALSE,
  nameCalibWeight = "calibWeight", minMaxTrim = NULL, print_every_n = 100) {

  check_population_totals(conP, dat, "personal")
  check_population_totals(conH, dat, "household")
  check_eps(conP, epsP, type = "personal")
  check_eps(conH, epsH, type = "household")
  variableKeepingTheBaseWeight <- w
  variableKeepingTheCalibWeight <- nameCalibWeight
  if ("variableKeepingTheBaseWeight" %in% names(dat))
    stop("The provided dataset must not have a column called",
         " 'variableKeepingTheBaseWeight'")
  if(!is.null(minMaxTrim)){
    if(length(minMaxTrim)!=2)
      stop("minMaxTrim must have exactly 2 elements, a minimum and a maximum.")
    if(!is.numeric(minMaxTrim)){
      stop("minMaxTrim must be a numeric vector of length two.")
    }
    if(minMaxTrim[2]<minMaxTrim[1]){
      stop("minMaxTrim must have a minimum as a first element and a maximum as a second element.
           But in the input the second element is smaller than the first.")
    }
  }
  OriginalSortingVariable <- V1 <- epsvalue <-
    f <- temporary_hvar <-
    value <- wValue <- representativeHouseholdForCalibration <- ..hid <- NULL
  dat_original <- dat
  dat <- copy(dat)
  ## originalsorting is fucked up without this
  dat[, OriginalSortingVariable := .I]
  meanfun <- getMeanFun(meanHH)

  # dat sollte ein data.table sein
  # w ein Name eines Basisgewichts oder NULL
  valueP <- paste0("valueP", seq_along(conP))
  ###fixed target value, should not be changed in iterations
  valueH <- paste0("valueH", seq_along(conH))
  ###Housekeeping of the varNames used
  usedVarNames <- c(valueP, valueH, "value",
                    "representativeHouseholdForCalibration", "wValue")

  if (any(names(dat) %in% usedVarNames)) {
    renameVars <- names(dat)[names(dat) %in% usedVarNames]
    setnames(dat, renameVars, paste0(renameVars, "_safekeeping"))
  }
  ### Treatment of HID, creating 0,1 var for being the first hh member
  #delVars <- c()
  if (is.null(hid)) {
    #delVars <- c("hid")
    hid <- "hid"
    dat[, hid := as.factor(seq_len(nrow(dat)))]
    dat[, representativeHouseholdForCalibration := 1]
  } else {

    if(!hid%in%colnames(dat)){
      stop("dat does not contain column ",hid)
    }

    if(any(is.na(dat[[hid]]))){
      stop("hid contains missing values")
    }

    if (!is.factor(dat[[hid]]))
      data.table::set(dat, NULL, hid, as.factor(dat[[hid]]))
    dat[, representativeHouseholdForCalibration :=
          as.numeric(!duplicated(get(..hid)))]
  }

  ## Names of the calibration variables for Person and household dimension
  pColNames <- lapply(conP, function(x) names(dimnames(x)))
  hColNames <- lapply(conH, function(x) names(dimnames(x)))

  for (i in seq_along(conP)) {
    current_colnames <- pColNames[[i]]

    for (colname in current_colnames) {
      if (!inherits(dat[[colname]], "factor")) {
        if (conversion_messages)
          message("converting column ", colname, " to factor")
        set(
          dat, j = colname,
          value = factor(dat[[colname]],
                         levels = dimnames(conP[[i]])[[colname]])
        )
      }
      else if (!identical(levels(dat[[colname]]),
                          dimnames(conP[[i]])[[colname]])) {
        if (conversion_messages)
          message("correct levels of column ", colname)
        set(
          dat, j = colname, value = factor(
            dat[[colname]], levels = dimnames(conP[[i]])[[colname]])
        )
      }
    }
    combined_factors <- combine_factors(dat, conP[[i]])
    set(dat, j = paste0("combined_factors_", i), value = combined_factors)
    set(dat, j = paste0("valueP", i),
        value = as.vector(conP[[i]][combined_factors]))
  }
  for (i in seq_along(conH)) {
    colnames <- hColNames[[i]]

    ## make sure the columns mentioned in the contingency table are in fact
    ##   factors
    for (colname in colnames) {
      if (!inherits(dat[[colname]], "factor")) {
        if (conversion_messages)
          message("converting column ", colname, " to factor")
        set(
          dat, j = colname, value = factor(
            dat[[colname]], levels = dimnames(conH[[i]])[[colname]])
        )
      }
      else if (!identical(levels(dat[[colname]]),
                          dimnames(conH[[i]])[[colname]])) {
        if (conversion_messages)
          message("correct levels of column ", colname)
        set(
          dat, j = colname, value = factor(
            dat[[colname]], levels = dimnames(conH[[i]])[[colname]])
        )
      }
    }

    combined_factors <- combine_factors(dat, conH[[i]])

    set(dat, j = paste0("combined_factors_h_", i), value = combined_factors)
    set(dat, j = paste0("valueH", i),
        value = as.vector(conH[[i]][combined_factors]))
  }

  if (is.null(variableKeepingTheBaseWeight)) {
    if (!is.null(bound) && is.null(w))
      stop("Bounds are only reasonable if base weights are provided")
    set(dat, j = variableKeepingTheCalibWeight, value = 1)
  } else {
    set(dat, j = variableKeepingTheCalibWeight,
        value = dat[[variableKeepingTheBaseWeight]])
  }

  if (check_hh_vars) {
    ## Check for non-unqiue values inside of a household for variabels used
    ##   in Household constraints
    for (hh in hColNames) {
      for (h in hh) {
        setnames(dat, h, "temporary_hvar")
        if (dat[, length(unique(temporary_hvar)),
                by = c(hid)][, any(V1 != 1)]) {
          stop(paste(h, "has different values inside a household"))
        }
        setnames(dat, "temporary_hvar", h)
      }
    }
  }

  if (is.list(epsP)) {
    for (i in seq_along(epsP)) {
      if (is.array(epsP[[i]])) {
        combined_factors <- dat[[paste0("combined_factors_", i)]]
        set(dat, j = paste0("epsP_", i),
            value = as.vector(epsP[[i]][combined_factors]))
      } else {
        set(dat, j = paste0("epsP_", i), value = epsP[[i]])
      }
    }
  } else {
    for (i in seq_along(conP)) {
      set(dat, j = paste0("epsP_", i), value = epsP)
    }
  }
  if (is.list(epsH)) {
    for (i in seq_along(epsH)) {
      if (is.array(epsH[[i]])) {
        combined_factors <- dat[[paste0("combined_factors_h_", i)]]
        set(dat, j = paste0("epsH_", i),
            value = as.vector(epsH[[i]][combined_factors]))
      } else {
        set(dat, j = paste0("epsH_", i), value = epsH[[i]])
      }
    }
  } else {
    for (i in seq_along(conH)) {
      set(dat, j = paste0("epsH_", i), value = epsH)
    }
  }
  ###Calib
  error <- TRUE
  calIter <- 1
  while (error && calIter <= maxIter) {
    error <- FALSE

    if (allPthenH) {
      ### Person calib
      
      for (i in seq_along(conP)) {
        
        numericalWeightingTmp <- NULL
        if (isTRUE(names(conP)[i] != "")) {
          numericalWeightingTmp <- names(conP)[i]
        }
        error <- calibP(
          i = i, dat = dat, error = error, valueP = valueP,
          pColNames = pColNames, bound = bound, verbose = verbose,
          calIter = calIter, numericalWeighting = numericalWeighting,
          numericalWeightingVar = numericalWeightingTmp,
          w = variableKeepingTheBaseWeight,
          cw = variableKeepingTheCalibWeight, minMaxTrim = minMaxTrim,
          print_every_n = print_every_n)

      }
      
      # ## replace person weight with household average
      set(dat, j = variableKeepingTheCalibWeight,
          value = meanfun(dat[[variableKeepingTheCalibWeight]], dat[[hid]]))
      
      ### Household calib
      for (i in seq_along(conH)) {
        numericalWeightingTmp <- NULL
        if (isTRUE(names(conH)[i] != "")) {
          numericalWeightingTmp <- names(conH)[i]
        }
        error <- calibH(
          i = i, dat = dat, error = error, valueH = valueH,
          hColNames = hColNames, bound = bound, verbose = verbose,
          calIter = calIter, looseH = looseH,
          numericalWeighting = numericalWeighting,
          numericalWeightingVar = numericalWeightingTmp,
          w = variableKeepingTheBaseWeight,
          cw = variableKeepingTheCalibWeight, minMaxTrim = minMaxTrim,
          print_every_n = print_every_n)
        
      }
    } else {
      ### Person calib
      for (i in seq_along(conP)) {
        numericalWeightingTmp <- NULL
        if (isTRUE(names(conP)[i] != "")) {
          numericalWeightingTmp <- names(conP)[i]
        }
        error <- calibP(
          i = i, dat = dat, error = error, valueP = valueP,
          pColNames = pColNames, bound = bound, verbose = verbose,
          calIter = calIter, numericalWeighting = numericalWeighting,
          numericalWeightingVar = numericalWeightingTmp,
          w = variableKeepingTheBaseWeight,
          cw = variableKeepingTheCalibWeight, minMaxTrim = minMaxTrim,
          print_every_n = print_every_n)

        ## replace person weight with household average
        set(dat, j = variableKeepingTheCalibWeight,
            value = meanfun(dat[[variableKeepingTheCalibWeight]], dat[[hid]]))

        ### Household calib
        for (i in seq_along(conH)) {
          numericalWeightingTmp <- NULL
          if (isTRUE(names(conH)[i] != "")) {
            numericalWeightingTmp <- numericalWeighting
          }
          error <- calibH(
            i = i, dat = dat, error = error, valueH = valueH,
            hColNames = hColNames, bound = bound, verbose = verbose,
            calIter = calIter, numericalWeighting = numericalWeighting,
            numericalWeightingVar = numericalWeightingTmp, looseH = looseH,
            w = variableKeepingTheBaseWeight,
            cw = variableKeepingTheCalibWeight, minMaxTrim = minMaxTrim,
            print_every_n = print_every_n)
        }
      }
    }

    if (verbose && !error) {
      message("Iteration stopped after ", calIter, " steps")
    } else if (maxIter == calIter) {
      warning("Not converged in ", maxIter, " steps")
    }
    calIter <- calIter + 1
  }
  # Remove Help Variables
  fVariableForCalibrationIPF <- NULL
  # dat[, fVariableForCalibrationIPF := NULL]
  addWeightsAndAttributes(dat, conP, conH, epsP, epsH, dat_original, maxIter,
                          calIter, returnNA, variableKeepingTheCalibWeight,variableKeepingTheBaseWeight,
                          verbose, looseH, hid)
}
