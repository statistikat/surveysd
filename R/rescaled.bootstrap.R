#' @title Draw bootstrap replicates
#'
#' @description Draw bootstrap replicates from survey data using the rescaled
#'   bootstrap for stratified multistage sampling, presented by Preston, J.
#'   (2009).
#'
#' @param dat either data frame or data table containing the survey sample
#' @param method for bootstrap replicates, either preston or rao-wu 
#' @param REP integer indicating the number of bootstraps to be drawn
#' @param strata string specifying the column name in `dat` that is used for
#'   stratification. For multistage sampling multiple column names can be
#'   specified by `strata=c("strata1","strata2","strata3")` or 
#'   `strata=c("strata1>strata2>strata3")`. See Details for more
#'   information.
#' @param cluster string specifying the column name in `dat` that is used for
#'   clustering. For instance given a household sample the column containing
#'   the household ID should be supplied.
#'   For multistage sampling multiple column names can be specified
#'   by `cluster=c("cluster1","cluster2","cluster3")` or 
#'   `cluster=c("cluster1>cluster2>cluster3")`.
#'   See Details for more information.
#' @param fpc string specifying the column name in `dat` that contains the
#'   number of PSUs at the first stage. For multistage sampling the number of
#'   PSUs at each stage must be specified by `strata=c("fpc1","fpc2","fpc3")` 
#'   or `strata=c("fpc1>fpc2>fpc3")`.
#' @param single.PSU either "merge" or "mean" defining how single PSUs need to
#'   be dealt with. For `single.PSU="merge"` single PSUs at each stage are
#'   merged with the strata or cluster with the next least number of PSUs. If
#'   multiple of those exist one will be select via random draw. For
#'   `single.PSU="mean"` single PSUs will get the mean over all bootstrap
#'   replicates at the stage which did not contain single PSUs.
#' @param return.value either "data", "replicates" and/or "selection" 
#'   specifying the return value of the function. For "data" the survey data is 
#'   returned as class `data.table`, for "replicates" only the bootstrap replicates 
#'   are returned as `data.table`. For "selection" list of data.tables with 
#'   length of `length(strata)` is returned containing 1:`REP` 0-1 columns 
#'   indicating if a PSU was selected for each sampling stage.
#' @param run.input.checks logical, if TRUE the input will be checked before applying
#'   the bootstrap procedure
#' @param already.selected list of data.tables or `NULL` where each data.table contains 
#'   columns in `cluster`, `strata` and additionally 1:`REP` columns containing
#'   0-1 values which indicate if a PSU was selected for each bootstrap replicate.
#'   Each of the data.tables corresponds to one of the sampling stages. First entry
#'   in the list corresponds to the first sampling stage and so on.
#' @param seed integer specifying the seed for the random number generator.
#' 
#' @details For specifying multistage sampling designs the column names in
#' `strata`,`cluster` and `fpc` need to be seperated by ">".\cr
#' For multistage sampling the strings are read from left to right meaning that
#' the first vector entry or column name before the first ">" is taken as the column for
#' stratification/clustering/number of PSUs at the first and the last vector entry
#' or column after
#' the last ">" is taken as the column for stratification/clustering/number of
#' PSUs at the last stage.
#' If for some stages the sample was not stratified or clustered one must
#' specify this by "1" or "I", e.g. `strata=c("strata1","I","strata3")` or 
#' `strata=c("strata1>I>strata3")` if there was
#' no stratification at the second stage or 
#' `cluster=c("cluster1","cluster2","I")` respectively 
#' `cluster=c("cluster1>cluster2>I")`
#' if there were no clusters at the last stage.\cr
#' The number of PSUs at each stage is not calculated internally and must be
#' specified for any sampling design.
#' For single stage sampling using stratification this can usually be done by
#' adding over all sample weights of each PSU by each strata-code.\cr
#' Spaces in each of the strings will be removed, so if column names contain
#' spaces they should be renamed before calling this procedure!\cr 
#' If `already.selected` is supplied the sampling of bootstrap replicates 
#' considers if speficif PSUs have already been selected by a previous survey wave.
#' For a specific `strata` and `cluster` this could lead to more than `floor(n/2)` 
#' records selected. In that case records will be de-selected such that `floor(n/2)` records,
#' with `n` as the total number of records, are selected for each 
#' `strata` and `cluster`. This parameter ist mostly used by [draw.bootstrap] in
#' order to consider the rotation of the sampling units over time.
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
    dat, 
    method = c("Preston", "Rao-Wu"), 
    REP = 1000, 
    strata = "DB050>1",   
    cluster = "DB060>DB030",   
    fpc = "N.cluster>N.households",
    single.PSU = c("merge", "mean"), 
    return.value = c("data", "replicates"),
    run.input.checks = TRUE, 
    already.selected = NULL, 
    seed = NULL) {
  
  
  InitialOrder <- N <- SINGLE_BOOT_FLAG <- SINGLE_BOOT_FLAG_FINAL <- f <-
    n_prev <- n_draw_prev <- sum_prev <- n_draw <- NULL  
  
  dat <- copy(dat)
  
  # check run.input.checks
  if (!is.logical(run.input.checks)) {
    stop("run.input.checks can only be logical")
  }
  if (run.input.checks) {
    # check input data
    if (!is.data.table(dat) & !is.data.frame(dat)) {
      stop("dat needs to be a data frame or data table")
    } else {
      dat <- data.table(dat)
    }
    # check seed
    if(!is.null(seed)){
      check.input(seed, input.name = "seed", input.length=1, input.type="numeric")
    }
  }
  set.seed        ################################################################################################# @Johannes: Zahl ? #####################################
  
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
  
  
  # check inputs strata, cluster, fpc
  check.input(strata, input.name = "strata", input.type="character")
  check.input(fpc, input.name = "fpc", input.type="character")
  check.input(cluster, input.name = "cluster", input.type="character")
  
  if(length(strata)==1 && strata %like% ">"){
    strata <- unlist(strsplit(strata, ">"))
    strata <- gsub("\\s", "", strata)
  } 
  
  if(length(cluster)==1 && cluster %like% ">"){
    cluster <- unlist(strsplit(cluster, ">"))
    cluster <- gsub("\\s", "", cluster)
  }
  if(length(fpc)==1 && fpc %like% ">"){
    fpc <- unlist(strsplit(fpc, ">"))
    fpc <- gsub("\\s", "", fpc)
  }
  
  single.PSU <- single.PSU[1]  
  # return.value <- return.value[1]
  
  # continue input checks
  if (run.input.checks) {
    
    # check REP
    check.input(REP, input.name = "REP", input.length=1, input.type="numeric",
                decimal.possible = FALSE)
    
    # check method
    method <- match.arg(method)
    
    # Warning for Rao-Wu method
    if (method == "Rao-Wu") {
      warning("The 'Rao-Wu' method should only be used when the first stage sampling is conducted without replacement. Ensure that your sampling design follows this requirement.")
    }
    
    # check design variables
    # check length of strata, cluster, fpc
    if(length(strata) != length(fpc) | length(strata) != length(cluster)){
      stop("strata, cluster, and fpc need to have the same number of arguments",
           " either separated with '>' or be character vectors of the same length")
    }
    
    # cluster can only be "1" or "I" in the final stage
    if(length(cluster)>1 && any(cluster[1:(length(cluster)-1)] %in% c("I","1"))){
      stop("Specifying 'I' or '1' for parameter 'cluster' 
           can only be set in the final sampling stage!")
    }
    
    check.values <- c(strata, cluster, fpc)
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
    if (any(!return.value %in% c("data", "replicates", "selection"))) {
      stop("return.value can only take the values 'data', 'replicates' or 'selection'")
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
        countMultiple <- dat[,uniqueN(strata_i), by=c(cluster[i]),
                             env = list(strata_i = strata[i])]
        countMultiple <- countMultiple[V1>1]
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
  print("stages: ")
  print(stages)
  
  n <- nrow(dat)
  
  # define values for final calculation
  n.calc <- matrix(0, nrow = n, ncol = stages) # n = sample size
  N.calc <- matrix(0, nrow = n, ncol = stages) #  N = total population
  n_draw.calc <- matrix(0, nrow = n, ncol = stages)  # n_draw = draws in different satges
  delta.calc <- array(0, dim = c(n, stages, REP)) 
  
  delta_selection <- list()
  
  for (i in 1:stages) {  

    # define by.val
    by.val <- strata[i]  # by.val determines the grouping for the current level
    if (i > 1) {
      by.val <- c(strata[1:(i - 1)], cluster[1:(i - 1)], strata[i]) 
    }
    by.val <- by.val[!by.val %in% c("1", "I")] # filter out placeholders
    
    # define cluster value
    clust.val <- cluster[i]
    if (clust.val %in% c("1", "I")) {
      clust.val <- paste0("ID_help_", i)  # creates helper ids if no real cluster is defined, sampling takes place at individual level (because a cluster definition is always required)
      dat[, c(clust.val) := .I, by = c(by.val)]
    }
    
    # check of fpc[i] is well defined
    check.fpc <- dat[,uniqueN(fpc_i), by=c(by.val), 
                     env = list(fpc_i = fpc[i])]  
    check.fpc <- check.fpc[V1>1]
    if (nrow(check.fpc) > 0) {
      stop("values in ", fpc[i], " do vary in some strata-cluster ",
           "combinations at sampling stage ", i)
    }
    
    singles <- dat[fpc_i>1, sum(!duplicated(clust.val)), by=c(by.val),   # finds groups in which only a single psu is present -> this is problematic for bootstrap (cannot be resampled with only one element)
                   env = list(fpc_i = fpc[i],
                              clust.val = clust.val)]
    singles <- singles[V1==1]
    if (nrow(singles) > 0) {
      # if singel.PSU=="merge" change the coding of the current stage of the
      #   single PSU
      # to the coding in the same subgroup, according the next higher stage,
      #   with the smallest number of PSUs
      # if multiple smallest exist choose one per random draw
      higher.stages <- c(strata[1:(i - 1)], cluster[1:(i - 1)])
      by.val.tail <- tail(by.val,1)
      firstStage <- length(higher.stages) == 0
      if (firstStage) {
        higher.stages <- by.val
      }
      singles <- unique(subset(singles, select = higher.stages))
      
      ################ SINGLE BEGIN ################ ##############################################################################################
      
      if (method == "Preston") {
        if (single.PSU == "merge") { # strata with only one PSU are merged with the nearest smaller stratum to avoid issues in the bootstrap procedure.

          if (return.value == "data") {
            by.val.orig <- paste0(by.val.tail, "_ORIGINALSINGLES")
            fpc_orig <- paste0(fpc[i], "_ORIGINALSINGLES")
            
            set(dat, i = by.val.orig, value = dat[[by.val.tail]])
            set(dat, i = fpc_orig, value = dat[[fpc[i]]])
          }
          
          setkeyv(dat, higher.stages)
          next.PSU <- dat[singles, .(N = sum(!duplicated(clust.val))), by = c(by.val),
                          on = c(higher.stages),
                          env = list(clust.val = clust.val)]
          
          if(nrow(next.PSU)==1){
            stop("Only 1 single PSU in first sampling stage! Cannot combine single PSU.\nPlease manually combine sampling STRATA to eliminate single PSU!")
          }
          
          new.var <- paste0(tail(by.val, 1), "_NEWVAR")
          
          next.PSU[,c(new.var) := next_minimum(N, by.val.tail), by = c(higher.stages),
                   env = list(by.val.tail = by.val.tail)]
          if(any(is.na(next.PSU[[new.var]]))){
            if(firstStage){
              next.PUS[is.na(new.var_col), c(new.var) := head(higher.stages, 1),      ################################################################################################# @Johannes: PSU ? #####################################
                       env = list(new.var_col = new_var,
                                  higher.stages = higher.stages)]
            }else{
              next.PUS[is.na(new.var_col), c(new.var) := head(by.val.tail,1),       ################################################################################################# @Johannes: PSU ? #####################################
                       env = list(new.var_col = new_var,
                                  by.val.tail = by.val.tail)]
            }
          }
          
          # select singel PSUs and PSUs to join with
          next.PSU <- next.PSU[N == 1 | new.var == by.val.tail,
                               env = list( new.var = new.var,
                                           by.val.tail = by.val.tail)]
          dat <- merge(dat, next.PSU[, mget(c(by.val, new.var))],
                       by = c(by.val), all.x = TRUE)
          
          # sum over margins
          fpc.i_ADD <- paste0(fpc[i], "_ADD")
          dat[!is.na(new.var), c(fpc.i_ADD) := 
                sum(fpc_i[!duplicated(by.val.tail)]), by = c(new.var),
              env = list(new.var = new.var,
                         fpc_i = fpc[i],
                         by.val.tail = by.val.tail)]
          
          # assign to new group
          dat[!is.na(new.var), c(by.val.tail) := new.var, 
              env = list(new.var = new.var)]
          dat[ c(fpc[i]) := fpc_i[is.na(new.var)][1], by=c(by.val),
               env = list(new.var = new.var)]
          dat[!is.na(new.var), c(fpc[i]) := fpc.i_ADD,
              env = list(new.var = new.var,
                         fpc.i_ADD = fpc.i_ADD)]
          
          dat[, c(new.var, paste0(fpc[i], "_ADD")) := NULL]
        } 
        
        else if (single.PSU == "mean") {  # if single = mean: Marks the affected observations with a SINGLE_BOOT_FLAG in order to simply take the average of other replicates later in the bootstrap process.
          
          singles[, SINGLE_BOOT_FLAG := paste(higher.stages, .GRP, sep = "-"),  
                  by = c(higher.stages)]
          
          dat <- merge(dat, singles, by = c(higher.stages), all.x = TRUE)
          if (!"SINGLE_BOOT_FLAG_FINAL" %in% colnames(dat)) {
            dat[, SINGLE_BOOT_FLAG_FINAL := SINGLE_BOOT_FLAG]
          } else {
            dat[is.na(SINGLE_BOOT_FLAG_FINAL),
                SINGLE_BOOT_FLAG_FINAL := SINGLE_BOOT_FLAG]
          }
          dat[, SINGLE_BOOT_FLAG := NULL]
          
        } 
        
        else {
          message("Single PSUs detected at the following stages:\n")
          dat.print <- dat[,sum(!duplicated(clust.val)), by=c(by.val),
                           env = list(clust.val = clust.val)]
          dat.print <- dat.print[V1==1,.SD,.SDcols=c(by.val)]
        }
      } else if (method == "Rao-Wu") {
        if (single.PSU == "merge") {

          if (return.value == "data") {
            by.val.orig <- paste0(by.val.tail, "_ORIGINALSINGLES")
            fpc_orig <- paste0(fpc[i], "_ORIGINALSINGLES")
            
            set(dat, i = by.val.orig, value = dat[[by.val.tail]])
            set(dat, i = fpc_orig, value = dat[[fpc[i]]])
          }
          
          setkeyv(dat, higher.stages)
          next.PSU <- dat[singles, .(N = sum(!duplicated(clust.val))), by = c(by.val),
                          on = c(higher.stages),
                          env = list(clust.val = clust.val)]
          
          if(nrow(next.PSU)==1){
            warning("Only 1 single PSU in first sampling stage! Results may be unreliable. Consider manually combining strata.")
          }
          
          new.var <- paste0(tail(by.val, 1), "_NEWVAR")
          
          next.PSU[,c(new.var) := next_minimum(N, by.val.tail), by = c(higher.stages),
                   env = list(by.val.tail = by.val.tail)]
          if(any(is.na(next.PSU[[new.var]]))){
            if(firstStage){
              next.PUS[is.na(new.var_col), c(new.var) := head(higher.stages, 1),
                       env = list(new.var_col = new_var,
                                  higher.stages = higher.stages)]
            }else{
              next.PUS[is.na(new.var_col), c(new.var) := head(by.val.tail,1),
                       env = list(new.var_col = new_var,
                                  by.val.tail = by.val.tail)]
            }
          }
          
          # select singel PSUs and PSUs to join with
          next.PSU <- next.PSU[N == 1 | new.var == by.val.tail,
                               env = list( new.var = new.var,
                                           by.val.tail = by.val.tail)]
          dat <- merge(dat, next.PSU[, mget(c(by.val, new.var))],
                       by = c(by.val), all.x = TRUE)
          
          # sum over margins
          fpc.i_ADD <- paste0(fpc[i], "_ADD")  
          dat[!is.na(new.var), c(fpc.i_ADD) := 
                sum(fpc_i[!duplicated(by.val.tail)]), by = c(new.var),
              env = list(new.var = new.var,
                         fpc_i = fpc[i],
                         by.val.tail = by.val.tail)]
          
          # assign to new group
          dat[!is.na(new.var), c(by.val.tail) := new.var, 
              env = list(new.var = new.var)]
          dat[ c(fpc[i]) := fpc_i[is.na(new.var)][1], by=c(by.val),
               env = list(new.var = new.var)]
          dat[!is.na(new.var), c(fpc[i]) := fpc.i_ADD,
              env = list(new.var = new.var,
                         fpc.i_ADD = fpc.i_ADD)]
          
          dat[, c(new.var, paste0(fpc[i], "_ADD")) := NULL]
        } 
        
        else if (single.PSU == "mean") {
          
          dat[!is.na(SINGLE_BOOT_FLAG_FINAL), c(bootRep) := lapply( 
            .SD,
            function(z) {  # Ensure that the single PSU is pulled at least once
              if (all(is.na(z))) {
                avg <- mean(z, na.rm = TRUE)
                
                if (is.na(avg)) avg <- 0  # If all values of a replica column are NA, a default value (e.g. 0) is used -> Rao-Wu method is designed to deal with uncertain data
                return(rep(avg, length(z)))  
              } else {
                return(z)
              }
            }
          ), by = SINGLE_BOOT_FLAG_FINAL, .SDcols = c(bootRep)]
        }
        else {
          message("Single PSUs detected at the following stages:\n")
          dat.print <- dat[,sum(!duplicated(clust.val)), by=c(by.val),
                           env = list(clust.val = clust.val)]
          dat.print <- dat.print[V1==1,.SD,.SDcols=c(by.val)]
          
        }
        
      }
    }
    ################ SINGLE END ################ ##############################################################################################
    
    # get Stage
    if (i == 1) { 
      dati <- dat[,.(N=fpc_i[1], 
                     clust.val = unique(clust.val),
                     f = 1, 
                     n_prev = 1, n_draw_prev = 1, 
                     sum_prev = 1),  
                  by = c(by.val),
                  env = list(fpc_i = fpc[i],
                             clust.val = clust.val)]
      
    } else { # In level 2 or higher: The information from the previous draw is now transferred
      dati <- dat[,.(N=fpc_i[1], 
                     clust.val = unique(clust.val),
                     f = f[1], 
                     n_prev = n_prev[1], 
                     n_draw_prev = n_draw_prev, 
                     sum_prev = sum_prev[1]),
                  by = c(by.val),
                  env = list(fpc_i = fpc[i],
                             clust.val = clust.val)]
    }
    

    deltai <- paste0("delta_", i, "_", 1:REP) 
    dati[, n := .N, by = c(by.val)]  
    
    if (method == "Preston") {
      dati[,n_draw:=floor(n/2),by = c(by.val)]   
    } else {
      dati[, n_draw := n - 1, by = c(by.val)] 
    }
    
    
    if (nrow(dati[n_draw == 0]) > 0) {  
      stop("Resampling 0 PSUs should not be possible! Please report bug in ",
           "https://github.com/statistikat/surveysd")
    }
    
    # do bootstrap for i-th stage
    if(!is.null(already.selected)){ 
      # if already.selected was supplied consider already selected units
      dati_selected <- already.selected[[i]]
      dati[dati_selected,c(deltai):=mget(deltai),on=c(by.val,clust.val)]
      
      dati[,c(deltai):=lapply(.SD,function(delta,n,n_draw){
        
        if (method == "Rao-Wu") {  
          draw.with.replacement(n[1], n_draw[1], delta=delta) 
        } else if (method == "Preston"){
          draw.without.replacement(n[1],n_draw[1], delta=delta)
          print("n[1]:")
          print(n[1])
          
          print("n_draw[1]:")
          print(n_draw)
        }
      },n=n,n_draw=n_draw),by=c(by.val),.SDcols=c(deltai)] 
      
      dati_check <- dati[,lapply(.SD,function(z,n_draw){
        sum(z)==n_draw[1]
      },n_draw=n_draw),by=c(by.val),.SDcols=c(deltai)]
      
      if(any(!dati_check[,.SD,.SDcols=c(deltai)])){
        stop("Wrong number of units selected! Please report bug in ",
             "https://github.com/statistikat/surveysd")
      }
    }else{ 
      # do simple sampling without replacement
      if (method == "Rao-Wu") {
        dati[, c(deltai) := as.data.table(
          replicate(REP, draw.with.replacement(n[1], n_draw[1]),simplify = FALSE)
        ), by = c(by.val)] 
        
      } else if (method == "Preston") {
        print("n[1]:")
        print(n[1])
        
        print("n_draw[1]:")
        print(n_draw)
        
        dati[, c(deltai) := as.data.table(
          replicate(REP, draw.without.replacement(n[1], n_draw[1]), 
                    simplify = FALSE)),
          by = c(by.val)]
      }
    }
    
    
    # merge with data
    if(i>1){ 
      dat[,c("f","n_prev","n_draw_prev","sum_prev"):=NULL]
    }
    
    dat <- merge(dat,dati,by=c(by.val,clust.val)) 
    setorder(dat,InitialOrder) 
    
    # prepare output for return.value %in% "selection"
    delta_selection <- c(delta_selection,
                         list(dat[, .SD, .SDcols = c(cluster, strata, deltai)]))  
    
    # extract information from data.table and remove again from data table
    # (less memory intensive)
    # only matrices and arrays needed for final calculation
    n.calc[, i] <- dat[, n]
    N.calc[, i] <- dat[, N]
    n_draw.calc[, i] <- dat[, n_draw]
    delta.calc[, i, ] <- as.matrix(dat[, mget(deltai)])
    
    dat[, f := n / N * f]    ########################################################################################## Richtig?
    dat[, n_prev := n * n_prev]
    dat[, n_draw_prev := n_draw * n_draw_prev]
    
    dat[, c("n", "N", deltai, "n_draw") := NULL] 
    
    rm(dati);gc() 
  }
  
  print("sampling fraction f: ")
  print (head(dat[,f]))
  # remove names
  dat[,c("f","n_prev","n_draw_prev","sum_prev"):=NULL]
  
  # rename delta_selection for output
  names(delta_selection) <- paste0("SamplingStage",1:stages)
  
  # calculate bootstrap replicate values
  bootRep <- paste0("bootRep", 1:REP) 
  
  dat[, c(bootRep) := as.data.table(calc.replicate( # call function calc.replicate -> calculate bootstrap replicate for each repetition
    method = method,
    n = n.calc, 
    N = N.calc, 
    n_draw = n_draw.calc, 
    delta = delta.calc))] 
  
  if (single.PSU == "mean") {
    dat[!is.na(SINGLE_BOOT_FLAG_FINAL), c(bootRep) := lapply(  
      .SD,
      function(z) {
        mean(z, na.rm = TRUE)
      }
    ), by = SINGLE_BOOT_FLAG_FINAL, .SDcols = c(bootRep)]
  }
  
  print("dat: ")
  print (head(dat))
  
  setkey(dat, InitialOrder)  
  if (length(removeCols) > 0) {
    dat[, c(removeCols) := NULL]
  }
  
  if ("data" %in% return.value) {
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
    
    output_bootstrap <- list(data = dat)
  } else if ("replicates" %in% return.value) {
    output_bootstrap <- list(replicates = dat[, mget(bootRep)]) 
  }
  
  if("selection" %in% return.value){
    if(length(return.value)>1){
      output_bootstrap <- c(output_bootstrap, list(selection = delta_selection))
    }else{
      output_bootstrap <- list(selection = delta_selection)
    }
  }
  
  if(length(output_bootstrap) == 1){
    output_bootstrap <- output_bootstrap[[1]]
  }
  return(output_bootstrap)
}


# Helper Functions 

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

# Sampling WITHOUT replacement (Preston method)
draw.without.replacement <- function(n, n_draw, delta = NULL) {   # delta: 1= selected, 0 = not selected, NA = not yet drawn, NULL = not existent yet

  # if no units have been selected prior
  if(is.null(delta)){
    delta <- rep(c(1.0, 0.0), c(n_draw, n - n_draw)) # if delta is not passed (no units have been created yet), it will be created here. All units are selected with prob of 1.0 or 0.0
    if (length(delta) > 1) {
      delta <- sample(delta)
    }
  }else{
    # if units have already been selected
    n_selected <- sum(delta,na.rm=TRUE)  # if delta exists, i.e. units have already been selected, delta is adjusted: n_selected = number of units selected so far
    delta_new_unit <- is.na(delta) 
    
    # check if only 1 unit needs to be selected
    if(sum(delta_new_unit)==1){
      # deselect one element in delta
      # so that resulting element will have a selection probability
      if(n_selected<n_draw){
        # set one entry in delta==0 to NA
        delta <- set.random2NA(delta,value2NA = 0)  
      }else{
        # set one entry in delta==1 to NA
        delta <- set.random2NA(delta,value2NA = 1) 
      }
      
      n_selected <- sum(delta,na.rm=TRUE)
      delta_new_unit <- is.na(delta)
    }
    
    if(n_draw<=n_selected){
      # if more entries have been selected than actually possible
      # or n_draw units have already been selected
      # due to already selected units
      # set some units with delta==1 to 0
      # such that entries with delta NA have selection probability >0
      add1 <- as.numeric(sum(delta_new_unit)>1) 
      nChanges <- n_selected - n_draw + add1
      
      if(nChanges == sum(delta_new_unit==FALSE)){ 
        delta <- set.random2NA(delta, value2NA=1,
                               nChanges = n_selected - n_draw + add1) 
      }else{
        delta <- change.random.value(delta, changeVal=1,
                                     nChanges = n_selected - n_draw + add1)
      }
      delta_new_unit <- is.na(delta)
      n_selected <- sum(delta,na.rm=TRUE)
    }
    
    if((n_draw-n_selected) >= sum(delta_new_unit)){
      # if new units will all be selected
      # or more units need to be selected than possible
      # set random entries from 0 to 1
      add1 <- as.numeric(sum(delta_new_unit)>1)
      nChanges <- n_draw-n_selected-sum(delta_new_unit) + add1
      delta <- change.random.value(delta, changeVal=0,
                                   nChanges = nChanges)
      n_selected <- sum(delta,na.rm=TRUE)
      delta_new_unit <- is.na(delta)
    }
    
    n_draw <- n_draw - n_selected
    delta_rest <- rep(c(1.0, 0.0), c(n_draw, sum(delta_new_unit)-n_draw))
    if (length(delta_rest) > 1) {
      delta_rest <- sample(delta_rest)
    }
    delta[delta_new_unit==TRUE] <- delta_rest
    
  }
  return(delta)
}


# Sampling WITH replacement (Rao-Wu method)
draw.with.replacement <- function(n, n_draw, delta = NULL) { 

  if (is.null(delta)) {
    delta <- rep(0, n)  # Set all units to 0 (not yet drawn)
  }

  n_selected <- sum(delta, na.rm = TRUE)
  
  # 1. draw additional units if less than n_draw are selected
  if (n_selected < n_draw) {
    n_needed <- n_draw - n_selected
    additional_draws <- sample(1:n, size = n_needed, replace = TRUE) 
    for (i in additional_draws) {
      delta[i] <- delta[i] + 1 # Update delta: Increase the frequency of units drawn
    }
  }
  
  # 2 Reduce excess units if more than n_draw are selected
  if (n_selected > n_draw) {
    over_selected <- n_selected - n_draw 
    for (i in which(delta > 0)[1:over_selected]) {
      delta[i] <- delta[i] - 1
    }
  }
  
  # 3. handle NA values in delta, if available
  if (any(is.na(delta))) {
    delta[is.na(delta)] <- 0     # NA --> 0 
  }
  
  return(delta)
}


change.random.value <- function(delta,changeVal=0,nChanges=1){
  
  changeVal2 <- fifelse(changeVal==0,1,0) 
  set2NA <- which(delta==changeVal & !is.na(delta))
  if(length(set2NA)>1){
    set2NA <- sample(set2NA,min(nChanges,length(set2NA)))
  }
  delta[set2NA] <- changeVal2
  return(delta)
}



set.random2NA <- function(delta,value2NA=0,nChanges=1){
  
  set2NA <- which(delta==value2NA & !is.na(delta))
  if(length(set2NA)>1){
    set2NA <- sample(set2NA,min(nChanges,length(set2NA)))
  }
  delta[set2NA] <- NA
  return(delta)
}


calc.replicate <- function(n, N, n_draw, delta , method = "Preston") {
  
  p <- ncol(n)  
  dimdelta <- dim(delta)
  for (i in 1:p) {
    if (method == "Preston") {
      if (i == 1) { # 
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
          
          print("prod_val:")
          print(prod_val)
        }

        rep_out <- rep_out + lambda * prod_val * (n[, i] / n_draw[, i] *      
                                                    delta[, i, ] - 1)
      }
      
    } else if (method == "Rao-Wu") {

      n_h <- n[, i]               # Number of PSUs drawn 
      m_h <- n_h - 1              # Omit 1 for resampling
      f_h <- n_h / N[, i]         # sample fraction
      print(f_h)
      
      # if (any(f_h > 0.1)) {  
      #   stop("Sampling Fraction is too big for Rao-Wu, choose preston instead") ################################################################################################# @Johannes: Welche Sampling Fraction ? (Recherche sagt maximales f von 0.05 - 0.1) #####################################
      # }
      
      w_hi <- N[, i] / n_h # Designgweights
      
      # Scaling factor λ determines the variance correction
      lambda <- sqrt(m_h * (1 - f_h) / (n_h - 1))
      
      # r_hi_star: Create replicas by drawing (n_h - 1) units without putting them back
      r_hi_star <- delta[, i, ]

      rep_out <- (1 - lambda + lambda * (n_h / m_h) * r_hi_star) * w_hi
      
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

