#' @title Draw bootstrap replicates
#'
#' @description Draw bootstrap replicates from survey data using the rescaled
#'   bootstrap for stratified multistage sampling, presented by Preston, J.
#'   (2009).
#'
#' @param dat either data frame or data table containing the survey sample
#' @param REP integer indicating the number of bootstraps to be drawn
#' @param strata string specifying the column name in `dat` that is used for
#'   stratification. For multistage sampling multiple column names can be
#'   specified by `strata=c("strata1>strata2>strata3")`. See Details for more
#'   information.
#' @param cluster string specifying the column name in `dat` that is used for
#'   clustering. For instance given a household sample the column containing
#'   the household ID should be supplied.
#'   For multistage sampling multiple column names can be specified
#'   by `cluster=c("cluster1>cluster2>cluster3")`.
#' See Details for more information.
#' @param period optional character specifying the name of the column in `dat`
#'   containing he sample periods. This column must be a factor where the factor
#'   levels determine the order of the values in column `period`. See details for
#'   more information
#' @param fpc string specifying the column name in `dat` that contains the
#'   number of PSUs at the first stage. For multistage sampling the number of
#'   PSUs at each stage must be specified by `strata=c("fpc1>fpc2>fpc3")`.
#' @param single.PSU either "merge" or "mean" defining how single PSUs need to
#'   be dealt with. For `single.PSU="merge"` single PSUs at each stage are
#'   merged with the strata or cluster with the next least number of PSUs. If
#'   multiple of those exist one will be select via random draw. For
#'   `single.PSU="mean"` single PSUs will get the mean over all bootstrap
#'   replicates at the stage which did not contain single PSUs.
#' @param return.value either "data" or "replicates" specifying the return value
#'   of the function. For "data" the survey data is returned as class
#'   `data.table`, for "replicates" only the bootstrap replicates are returned
#'   as `data.table`.
#' @param check.input logical, if TRUE the input will be checked before applying
#'   the bootstrap procedure
#'
#' @details For specifying multistage sampling designs the column names in
#' `strata`,`cluster` and `fpc` need to seperated by ">".\cr
#' For multistage sampling the strings are read from left to right meaning that
#' the column name before the first ">" is taken as the column for
#' stratification/clustering/number of PSUs at the first and the column after
#' the last ">" is taken as the column for stratification/clustering/number of
#' PSUs at the last stage.
#' If for some stages the sample was not stratified or clustered one must
#' specify this by "1" or "I", e.g. `strata=c("strata1>I>strata3")` if there was
#' no stratification at the second stage or `cluster=c("cluster1>cluster2>I")`
#' if there were no clusters at the last stage.\cr
#' The number of PSUs at each stage is not calculated internally and must be
#' specified for any sampling design.
#' For single stage sampling using stratification this can usually be done by
#' adding over all sample weights of each PSU by each strata-code.\cr
#' Spaces in each of the strings will be removed, so if column names contain
#' spaces they should be renamed before calling this procedure!\cr 
#' If `period` is supplied the sampling of bootstrap replicates for each period
#' will follow the method proposed by Preston, J. (2009) with the condition
#' that records drawn at the previous period are automatically selected.
#' This can result in more than halve of the records selected
#' for a specific `period`, `strata` and `cluster`. In that case
#' records will be de-selected such that `floor(n/2)` records, with `n` as the 
#' total number of records, are selected for each `period`, `strata` and `cluster`. 
#' If `NULL` (the default), it is assumed that all observations belong#
#' to the same period and this parameter will be ignored.
#'
#' @return returns the complete data set including the bootstrap replicates or
#'   just the bootstrap replicates, depending on `return.value="data"` or
#'   `return.value="replicates"` respectively.
#' @export rescaled.bootstrap
#'
#' @references Preston, J. (2009). Rescaled bootstrap for stratified multistage
#'   sampling. Survey Methodology. 35. 227-234.
#'
#' @author Johannes Gussenbauer, Statistics Austria
#'
#' @examples
#' 
#' library(surveysd)
#' library(data.table)
#' setDTthreads(1)
#' set.seed(1234)
#' eusilc <- demo.eusilc(n = 1,prettyNames = TRUE)
#' 
#' eusilc[,N.households:=uniqueN(hid),by=region]
#' eusilc.bootstrap <- rescaled.bootstrap(eusilc,REP=10,strata="region",
#'                                        cluster="hid",fpc="N.households")
#' 
#' eusilc[,new_strata:=paste(region,hsize,sep="_")]
#' eusilc[,N.housholds:=uniqueN(hid),by=new_strata]
#' eusilc.bootstrap <- rescaled.bootstrap(eusilc,REP=10,strata=c("new_strata"),
#'                                        cluster="hid",fpc="N.households")
#'
#'


rescaled.bootstrap <- function(
  dat, REP = 1000, strata = "DB050>1", cluster = "DB060>DB030", fpc =
    "N.cluster>N.households", single.PSU = c("merge", "mean"), return.value =
    c("data", "replicates"), check.input = TRUE, period = NULL) {

  InitialOrder <- N <- SINGLE_BOOT_FLAG <- SINGLE_BOOT_FLAG_FINAL <- f <-
    n_prev <- n_draw_prev <- sum_prev <- n_draw <- NULL

  dat <- copy(dat)
  
  if (!is.logical(check.input)) {
    stop("check.input can only be logical")
  }
  if (check.input) {
    # check input data
    if (!is.data.table(dat) & !is.data.frame(dat)) {
      stop("dat needs to be a data frame or data table")
    } else {
      dat <- data.table(dat)
    }
  }
  
  # prepare input
  removeCols <- c()
  if (is.null(cluster)) {
    cluster <- generateRandomName(20, colnames(dat))
    removeCols <- c(removeCols, cluster)
    dat[, c(cluster) := 1:.N]
  }
  if (is.null(strata)) {
    strata <- generateRandomName(20, colnames(dat))
    removeCols <- c(removeCols, strata)
    dat[, c(strata) := 1]
  }
  if( is.null(period) ){
    period <- generateRandomName(20, colnames(dat))
    removeCols <- c(removeCols, period)
    dat[,c(period):=1]
  }
  
  input <- c(strata, cluster, fpc)
  input <- gsub("\\s", "", input)
  input <- strsplit(input, ">")

  strata <- input[[1]]
  cluster <- input[[2]]
  fpc <- input[[3]]

  single.PSU <- single.PSU[1]
  return.value <- return.value[1]


  # continue input checks
  if (check.input) {
    
    # check REP
    if (!is.numeric(REP)) {
      stop("REP needs to be numeric")
    } else {
      if (length(REP) > 1) {
        warning("REP has length >1 - First argument will be used!")
        REP <- REP[1]
      }
      if (REP %% 1 != 0) {
        stop("REP cannot have a decimal part")
      }
    }

    # check design variables
    if (length(unique(lapply(input, length))) > 1) {
      stop("strata, cluster, and fpc need to have the same number of arguments",
           " separated with '>'")
    }
    check.values <- unlist(input)
    check.values <- check.values[!check.values %in% c("1", "I")]
    check.values <- check.values[!check.values %in% colnames(dat)]
    if (length(check.values) > 0) {
      stop("dat does not contain the column(s): ", check.values)
    }

    # check if there are missings in fpc
    if (any(is.na(dat[, .SD, .SDcols = c(fpc)]))) {
      stop("missing values not allowed in fpc")
    }

    # check return.value
    if (!return.value %in% c("data", "replicates")) {
      stop("return.value can only take the values 'data' or 'replicates'")
    }

    # check single.PSU
    if (is.null(single.PSU) || !single.PSU %in% c("merge", "mean")) {
      warning("single.PSU was not set to either 'merge' or 'mean'!\n Bootstrap",
              " replicates for single PSUs cases will be missing!")
      single.PSU <- FALSE
    }
    
    # check for each stage that PSUs are not in mutiple strata
    for(i in seq_along(strata)){
      if(!strata[i]%in%c("1","I") & !cluster[i]%in%c("1","I")){
        countMultiple <- dt.eval("dat[,uniqueN(",strata[i],"),by=.(",cluster[i],")][V1>1]")
        if(nrow(countMultiple)>0){
          stop("Some sampling units in ",cluster[i]," occur in multiple strata of ",strata[i])
        }
      }
    }
    
  }

  # check if variable f, N, n are in data.table
  overwrite.names <- c("f", "N", "n", "n_prev", "n_draw", "n_draw_prev")
  overwrite.names <- overwrite.names[overwrite.names %in% colnames(dat)]
  if (length(overwrite.names) > 0) {
    overwrite.names.new <- paste0("ORIGINAL_", overwrite.names)
    setnames(dat, overwrite.names, overwrite.names.new)

    strata[strata %in% overwrite.names] <-
      overwrite.names.new[strata %in% overwrite.names]
    cluster[cluster %in% overwrite.names] <-
      overwrite.names.new[cluster %in% overwrite.names]
    fpc[fpc %in% overwrite.names] <-
      overwrite.names.new[fpc %in% overwrite.names]
  }

  # set index for data to return dat in correct order
  dat[, InitialOrder := .I]

  # calculate bootstrap replicates
  stages <- length(strata)
  

  time_steps <- sort(unique(dat[[period]]))
  n <- nrow(dat)
  # define values for final calculation
  n.calc <- matrix(0, nrow = n, ncol = stages)
  N.calc <- matrix(0, nrow = n, ncol = stages)
  n_draw.calc <- matrix(0, nrow = n, ncol = stages)
  delta.calc <- array(0, dim = c(n, stages, REP))


  for(t in seq_along(time_steps)){
    
    # iterate over time steps
    index_t <-  dt.eval("dat[",period,"==time_steps[t],which=TRUE]")
    dat_t <- dat[index_t,]
    
    if(t==1){
      # save selected deltas for next period
      dati_prev <- list()
    }

    for (i in 1:stages) {
      # iterate over sampling stages
      
      # define by.val
      if (i > 1) {
        by.val <- c(strata[1:(i - 1)], cluster[1:(i - 1)], strata[i])
      } else{
        by.val <- strata[1:i]
      }
      by.val <- by.val[!by.val %in% c("1", "I")]
      
      # define cluster value
      clust.val <- cluster[i]
      if (clust.val %in% c("1", "I")) {
        clust.val <- paste0("ID_help_", i)
        dat_t[, c(clust.val) := .I, by = c(by.val)]
      }
      # check of fpc[i] is well defined
      # fpc[i] can not have different values per each by.val
      check.fpc <- dt.eval("dat_t[,uniqueN(", fpc[i], "),by=c(by.val)][V1>1]")
      if (nrow(check.fpc) > 0) {
        stop("values in ", fpc[i], " do vary in some strata-cluster ",
             "combinations at sampling stage ", i)
      }
      
      singles <- dt.eval("dat_t[", fpc[i], ">1,sum(!duplicated(", clust.val,
                         ")),by=c(by.val)][V1==1]")
      if (nrow(singles) > 0) {
        # if singel.PSU=="merge" change the coding of the current stage of the
        #   single PSU
        # to the coding in the same subgroup, according the next higher stage,
        #   with the smallest number of PSUs
        # if multiple smallest exist choose one per random draw
        higher.stages <- by.val[-length(by.val)]
        by.val.tail <- tail(by.val,1)
        firstStage <- length(higher.stages) == 0
        if (firstStage) {
          higher.stages <- by.val
        }
        singles <- unique(subset(singles, select = higher.stages))
        if (single.PSU == "merge") {
          # save original PSU coding and fpc values to replace changed values
          #   bevore returning the data.table
          if (return.value == "data") {
            dt.eval("dat_t[,", paste0(by.val.tail, "_ORIGINALSINGLES"),
                    ":=", by.val.tail, "]")
            dt.eval("dat_t[,", paste0(fpc[i], "_ORIGINALSINGLES"),
                    ":=", fpc[i], "]")
          }
          
          setkeyv(dat_t, higher.stages)
          
          next.PSU <- dt.eval("dat_t[singles,.(N=sum(!duplicated(", clust.val,
                              "))),by=c(by.val)]")
          
          new.var <- paste0(tail(by.val, 1), "_NEWVAR")
          
          dt.eval("next.PSU[,c(new.var):=next_minimum(N,", by.val.tail,
                  "),by=c(higher.stages)]")
          if(any(is.na(next.PSU[[new.var]]))){
            if(firstStage){
              dt.eval("next.PSU[is.na(",new.var,"),c(new.var):=head(",higher.stages,",1)]")
            }else{
              dt.eval("next.PSU[is.na(",new.var,"),c(new.var):=head(",by.val.tail,",1),by=c(higher.stages)]")
            }
          }
          
          # select singel PSUs and PSUs to join with
          next.PSU <- dt.eval("next.PSU[N == 1 | ",new.var," == ",by.val.tail,"]")
          dat_t <- merge(dat_t, next.PSU[, mget(c(by.val, new.var))],
                       by = c(by.val), all.x = TRUE)
          
          # sum over margins
          fpc.i_ADD <- paste0(fpc[i], "_ADD")
          dt.eval("dat_t[, c(fpc.i_ADD):=", fpc[i], "]")
          dt.eval("dat_t[!is.na(",new.var,"),c(fpc.i_ADD) := 
                sum(",fpc.i_ADD,"[!duplicated(",by.val.tail,")]),
                by=c(new.var)]")
          # assign to new group
          dt.eval("dat_t[!is.na(", new.var, "),c(by.val.tail):=", new.var, "]")
          
          dt.eval("dat_t[,", fpc[i], ":=", fpc[i], "[is.na(", new.var,
                  ")][1],by=c(by.val)]")
          dt.eval("dat_t[!is.na(", new.var, "),", fpc[i], ":=", fpc.i_ADD,"]")
          
          dat[, c(new.var, paste0(fpc[i], "_ADD")) := NULL]
        } else if (single.PSU == "mean") {
          # if single.PSU="mean" flag the observation as well as the all the
          #   observations in the higher group
          singles[, SINGLE_BOOT_FLAG := paste(higher.stages, .GRP, sep = "-"),
                  by = c(higher.stages)]
          
          dat_t <- merge(dat_t, singles, by = c(higher.stages), all.x = TRUE)
          if (!"SINGLE_BOOT_FLAG_FINAL" %in% colnames(dat_t)) {
            dat_t[, SINGLE_BOOT_FLAG_FINAL := SINGLE_BOOT_FLAG]
          } else {
            dat_t[is.na(SINGLE_BOOT_FLAG_FINAL),
                SINGLE_BOOT_FLAG_FINAL := SINGLE_BOOT_FLAG]
          }
          dat_t[, SINGLE_BOOT_FLAG := NULL]
          
        } else {
          message("Single PSUs detected at the following stages:\n")
          print(dt.eval("dat_t[,sum(!duplicated(", clust.val,
                        ")),by=c(by.val)][V1==1,.(",paste(by.val,collapse=","),")]"))
        }
      }
      
      # get Stage
      if (i == 1) {
        dati <- dt.eval(
          "dat_t[,.(N=", fpc[i], "[1],", clust.val, "=unique(",
          clust.val, "),f=1,n_prev=1,n_draw_prev=1,sum_prev=1),by=list(",
          paste(by.val, collapse = ","), ")]")
      } else {
        dati <- dt.eval(
          "dat_t[,.(N=", fpc[i], "[1],", clust.val, "=unique(", clust.val,
          "),f=f[1],n_prev=n_prev[1],n_draw_prev=n_draw_prev[1],",
          "sum_prev=sum_prev[1]),by=list(", paste(by.val, collapse = ","), ")]")
        dat_t[,c("f","n_prev","n_draw_prev","sum_prev") := NULL]

      }
      
      deltai <- paste0("delta_", i, "_", 1:REP)
      dati[, n := .N, by = c(by.val)]
      # determin number of psu to be drawn
      # dati[, n_draw := select.nstar(
      #   n[1], N[1], f[1], n_prev[1], n_draw_prev[1], sum_prev = NULL,
      #   new.method = new.method), by = c(by.val)]
      dati[,n_draw:=floor(n/2),by = c(by.val)]
      
      if (nrow(dati[n_draw == 0]) > 0) {
        stop("Resampling 0 PSUs should not be possible! Please report bug in ",
             "https://github.com/statistikat/surveysd")
      }
      # do bootstrap for i-th stage
      
      if(t == 1){
        # do simple sampling without replacement for first period
        dati[, c(deltai) := as.data.table(
          replicate(REP, draw.without.replacement(n[1], n_draw[1]),
                    simplify = FALSE)),
          by = c(by.val)]
      }else{
        # sample without replacement
        # but consider already selected units from previos period
        dati_selected <- dati_prev[[i]][,.SD,.SDcols=c(by.val,clust.val,deltai)]
        dati[dati_selected,c(deltai):=mget(deltai),on=c(by.val,clust.val)]
        dati[,c(deltai):=lapply(.SD,function(delta,n,n_draw){
          draw.without.replacement(n[1],n_draw[1],delta=delta)
        },n=n,n_draw=n_draw),by=c(by.val),.SDcols=c(deltai)]
        
        dati_check <- dati[,lapply(.SD,function(z,n_draw){
          sum(z)==n_draw[1]
        },n_draw=n_draw),by=c(by.val),.SDcols=c(deltai)]
        
        if(any(!dati_check[,.SD,.SDcols=c(deltai)])){
          stop("Wrong number of units selected! Please report bug in ",
               "https://github.com/statistikat/surveysd")
        }
      }

      
      # merge with data
      dat_t <- merge(dat_t,dati,by=c(by.val,clust.val))
      setorder(dat_t,InitialOrder)
      # extract information from data.table and remove again from data table
      # (less memory intensive)
      # only matrices and arrays needed for final calculation
      n.calc[index_t, i] <- dat_t[, n]
      N.calc[index_t, i] <- dat_t[, N]
      n_draw.calc[index_t, i] <- dat_t[, n_draw]
      delta.calc[index_t, i, ] <- as.matrix(dat_t[, mget(deltai)])
      
      # dat_t[, sum_prev := sum_prev +
      #       (sqrt(n_draw * f * (1 - n / N) / (n - n_draw)) *
      #          sqrt(n_prev / n_draw_prev) * (n / n_draw - 1))]
      dat_t[, f := n / N * f]
      dat_t[, n_prev := n * n_prev]
      dat_t[, n_draw_prev := n_draw * n_draw_prev]
      
      dat_t[, c("n", "N", deltai, "n_draw") := NULL]
      
      # save selected deltas for next period
      if(t==1){
        # at first time step build dati_prev
        dati_prev <- c(dati_prev,list(dati))
      }else{
        # later update dati_prev with previous results
        dati_prev[[i]] <- dati
      }

      rm(dati);gc()
    }
  }
 

  bootRep <- paste0("bootRep", 1:REP)
  dat[, c(bootRep) := as.data.table(calc.replicate(
    n = n.calc, N = N.calc, n_draw = n_draw.calc, delta = delta.calc))]

  if (single.PSU == "mean") {
    dat[!is.na(SINGLE_BOOT_FLAG_FINAL), c(bootRep) := lapply(
      .SD,
      function(z) {
        mean(z, na.rm = TRUE)
      }
      ), by = SINGLE_BOOT_FLAG_FINAL, .SDcols = c(bootRep)]
  }

  setkey(dat, InitialOrder)
  if (length(removeCols) > 0) {
    dat[, c(removeCols) := NULL]
  }

  if (return.value == "data") {
    # get original values for PSUs and fpc - if singles PSUs have been detected
    #   and merged
    if (single.PSU == "merge") {
      c.names <- colnames(dat)
      c.names <- c.names[grepl("_ORIGINALSINGLES", c.names)]
      if (length(c.names) > 0) {
        drop.names <- gsub("_ORIGINALSINGLES", "", c.names)
        dat[, c(drop.names) := NULL]
        setnames(dat, c(c.names), drop.names)
      }
    }
    dat[, c("InitialOrder") := NULL]
    # get original col values
    if (length(overwrite.names) > 0) {
      setnames(dat, overwrite.names.new, overwrite.names)
    }

    return(dat)
  } else if (return.value == "replicates") {
    return(dat[, mget(bootRep)])
  }
}

select.nstar <- function(n, N, f, n_prev, n_draw_prev, lambda_prev,
                         sum_prev = NULL, new.method) {

  if (n == 1) {
    # if only a single unit in strata
    # return missing
    # if single units are
    # not treated missing values are returned
    return(1L)
  }

  if (!is.null(sum_prev)) {
    n_draw <- (sum_prev) ^ 2 / ((1 - (n / N)) * n_prev * f +
                                   (sum_prev) ^ 2) * n
    n_draw <- floor(n_draw)
  } else {
    if (new.method) {
      n_draw <- (n * n_draw_prev) / (n_prev * f * (1 - n / N) + n_draw_prev)
      n_draw <- floor(n_draw)
    } else {
      n_draw <- floor(n / 2)
    }
  }

  return(n_draw)
}

draw.without.replacement <- function(n, n_draw, delta = NULL) {
  
  # if no units have been selected prior
  if(is.null(delta)){
    delta <- rep(c(1.0, 0.0), c(n_draw, n - n_draw))
    if (length(delta) > 1) {
      delta <- sample(delta)
    }
  }else{
    # if units have already been selected
    n_selected <- sum(delta,na.rm=TRUE)
    delta_new_unit <- is.na(delta)
    if(n_draw<n_selected){
      # if more entries have been selected than actually possible
      # due to selection in previous period
      # deselect units
      delta_rest <-  rep(c(1.0, 0.0), c(n_draw, sum(delta_new_unit==FALSE) - n_draw))
      if (length(delta_rest) > 1) {
        delta_rest <- sample(delta_rest)
      }
      
      # reselect n_draw units 
      # from already selected units
      # ~ delta is not missing
      delta[delta_new_unit==FALSE] <- delta_rest
      delta[delta_new_unit==TRUE] <- 0 
    }else{
      n_draw <- n_draw - n_selected
      
      # n_draw can be > sum(delta_new_unit)
      # -> then already assigned 0s need to be 
      # changed to 1s
      change_zero <- max(0,n_draw-sum(delta_new_unit))
      if(change_zero >0){
        n_draw <- min(n_draw,sum(delta_new_unit))
        delta_zero <- which(delta_new_unit==FALSE & delta==0)
        if(length(delta_zero)>1){
          delta_zero <- sample(delta_zero,change_zero)
        }
        delta[delta_zero] <- 1
      }
      
      delta_rest <- rep(c(1.0, 0.0), c(n_draw, sum(delta_new_unit)-n_draw))
      if (length(delta_rest) > 1) {
        delta_rest <- sample(delta_rest)
      }
      delta[delta_new_unit==TRUE] <- delta_rest
    }
  }
  return(delta)
}

calc.replicate <- function(n, N, n_draw, delta) {
  p <- ncol(n)
  # n_draw <- trunc(n/2)
  # n_draw <- floor(n/(2-rowCumprods(n/N))-1)
  # n_draw[n_draw<1] <- 1
  dimdelta <- dim(delta)
  for (i in 1:p) {
    if (i == 1) {
      lambda <- sqrt(n_draw[, 1] *
                       (1 - n[, 1] / N[, 1]) / (n[, 1] - n_draw[, 1]))
      rep_out <- 1 - lambda + lambda * n[, i] / n_draw[, i] * delta[, i, ]
    } else if (i == 2) {
      lambda <- (1 - n[, i] / N[, i]) / (n[, i] - n_draw[, i])
      lambda <- sqrt((n[, i - 1] / N[, i - 1]) * n_draw[, i] * lambda)
      rep_out <- rep_out + lambda *
        (sqrt(n[, i - 1] / n_draw[, i - 1]) * delta[, i - 1, ]) *
        (n[, i] / n_draw[, i] * delta[, i, ] - 1)
    } else {
      lambda <- (1 - n[, i] / N[, i]) / (n[, i] - n_draw[, i])
      lambda <- sqrt(rowProds(n[, 1:(i - 1)] / N[, 1:(i - 1)]) *
                       n_draw[, i] * lambda)
      prod_val <- matrix(0, ncol = dimdelta[3], nrow = dimdelta[1])
      for (r in 1:dimdelta[3]) {
        prod_val[, r] <- rowProds(sqrt(n[, 1:(i - 1)] / n_draw[, 1:(i - 1)]) *
                                    delta[, 1:(i - 1), r])
      }
      # rep_out <- rep_out + lambda*rowProds(sqrt(n[,1:(i-1)]/
      #   n_draw[,1:(i-1)])*delta[,1:(i-1),]) * (n[,i]/n_draw[,i]*delta[,i,]-1)
      rep_out <- rep_out + lambda * prod_val * (n[, i] / n_draw[, i] *
                                                  delta[, i, ] - 1)
    }
  }
  return(rep_out)
}

next_minimum <- function(N, by) {
  N_notOne <- N != 1
  by <- by[N_notOne][which.min(N[N_notOne])]
  if (length(by) > 1) {
    by <- sample(by, 1)
  }
  return(by)
}
