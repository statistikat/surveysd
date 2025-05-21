#' @title Draw bootstrap replicates
#'
#' @description Draw bootstrap replicates from survey data with rotating panel
#'   design. Survey information, like ID, sample weights, strata and population
#'   totals per strata, should be specified to ensure meaningfull survey
#'   bootstraping.
#'
#' @param dat either data.frame or data.table containing the survey data with
#'   rotating panel design.
#' @param REP integer indicating the number of bootstrap replicates.
#' @param hid character specifying the name of the column in `dat` containing
#'   the household id. If `NULL` (the default), the household structure is not
#'   regarded.
#' @param weights character specifying the name of the column in `dat`
#'   containing the sample weights.
#' @param period character specifying the name of the column in `dat` containing
#'   the sample periods. If `NULL` (the default), it is assumed that all
#'   observations belong to the same period.
#' @param strata character vector specifying the name(s) of the column in `dat`
#'   by which the population was stratified. If `strata` is a vector
#'   stratification will be assumed as the combination of column names contained
#'   in `strata`. Setting in addition `cluster` not NULL stratification will be
#'   assumed on multiple stages, where each additional entry in `strata`
#'   specifies the stratification variable for the next lower stage. see Details
#'   for more information.
#' @param cluster character vector specifying cluster in the data. If not
#'   already specified in `cluster` household ID is taken es the lowest level
#'   cluster.
#' @param totals character specifying the name of the column in `dat` containing
#'   the the totals per strata and/or cluster. Is ONLY optional if `cluster` is
#'   `NULL` or equal `hid` and `strata` contains one columnname! Then the
#'   households per strata will be calcualted using the `weights` argument. If
#'   clusters and strata for multiple stages are specified `totals` needs to be
#'   a vector of `length(strata)` specifying the column on `dat` that contain
#'   the total number of PSUs at each stage. `totals` is interpreted from left
#'   the right, meaning that the first argument corresponds to the number of
#'   PSUs at the first and the last argument to the number of PSUs at the last
#'   stage.
#' @param single.PSU either "merge" or "mean" defining how single PSUs need to
#'   be dealt with. For `single.PSU="merge"` single PSUs at each stage are
#'   merged with the strata or cluster with the next least number of PSUs. If
#'   multiple of those exist one will be select via random draw. For
#'   `single.PSU="mean"` single PSUs will get the mean over all bootstrap
#'   replicates at the stage which did not contain single PSUs.
#' @param boot.names character indicating the leading string of the column names
#'   for each bootstrap replica. If NULL defaults to "w".
#' @param split logical, if TRUE split households are considered using `pid`,
#'   for more information see Details.
#' @param pid column in `dat` specifying the personal identifier. This
#'   identifier needs to be unique for each person throught the whole data set.
#'
#' @return the survey data with the number of REP bootstrap replicates added as
#'   columns.
#'
#' @details `draw.bootstrap` takes `dat` and draws `REP` bootstrap replicates
#' from it.
#' `dat` must be household data where household members correspond to multiple
#' rows with the same household
#' identifier. For most practical applications, the following columns should be
#' available in the dataset
#' and passed via the corresponding parameters:
#'
#' * Column indicating the sample period (parameter `period`).
#' * Column indicating the household ID (parameter `hid`)
#' * Column containing the household sample weights (parameter `weights`);
#' * Columns by which population was stratified during the sampling process
#'   (parameter: `strata`).
#'
#' For single stage sampling design a column the argument `totals` is optional,
#' meaning that a column of the number of PSUs at the first stage does not need
#' to be supplied.
#' For this case the number of PSUs is calculated and added to `dat` using
#' `strata` and `weights`. By setting `cluster` to NULL single stage sampling
#' design is always assumed and
#' if `strata` contains of multiple column names the combination of all those
#' column names will be used for stratification.
#'
#' In the case of multi stage sampling design the argument `totals` needs to be
#' specified and needs to have the same number of arguments as `strata`.
#'
#' If `cluster` is `NULL` or does not contain `hid` at the last stage, `hid`
#' will automatically be used as the final cluster. If, besides `hid`,
#' clustering in additional stages is specified the number of column names in
#' `strata` and `cluster` (including `hid`) must be the same. If for any stage
#' there was no clustering or stratification one can set "1" or "I" for this
#' stage.
#'
#' For example `strata=c("REGION","I"),cluster=c("MUNICIPALITY","HID")` would
#' speficy a 2 stage sampling design where at the first stage the municipalities
#' where drawn stratified by regions
#' and at the 2nd stage housholds are drawn in each municipality without
#' stratification.
#'
#' Bootstrap replicates are drawn for each survey period consecutively (`period`) using the
#' function [rescaled.bootstrap].
#' Bootstrap replicates are drawn consistently in the way that in each `period` and
#' sampling stage always \eqn{\floor{n/2}} clusters are selected in each strata.
#'
#' This ensures that the bootstrap replicates follow the same logic as the
#' sampled households, making the bootstrap replicates more comparable to the
#' actual sample units.
#'
#' If `split` ist set to `TRUE` and `pid` is specified, the bootstrap replicates
#' are carried forward using the personal identifiers instead of the household
#' identifier.
#' This takes into account the issue of a houshold splitting up.
#' Any person in this new split household will get the same bootstrap replicate
#' as the person that has come from an other household in the survey.
#' People who enter already existing households will also get the same bootstrap
#' replicate as the other households members had in the previous periods.
#'
#' @return Returns a data.table containing the original data as well as the
#'   number of `REP` columns containing the bootstrap replicates for each
#'   repetition.\cr
#'   The columns of the bootstrap replicates are by default labeled "w*Number*"
#'   where *Number* goes from 1 to `REP`. If the column names of the bootstrap
#'   replicates should start with a different character or string the parameter
#'   `boot.names` can be used.
#'
#' @seealso [`data.table`][data.table::data.table] for more information on
#'   data.table objects.
#'
#' @author Johannes Gussenbauer, Alexander Kowarik, Statistics Austria
#'
#' @examples
#' 
#' library(surveysd)
#' library(data.table)
#' setDTthreads(1)
#' set.seed(1234)
#' eusilc <- demo.eusilc(n = 3, prettyNames = TRUE)
#'
#' ## draw sample without stratification or clustering
#' dat_boot <- draw.bootstrap(eusilc, REP = 1, weights = "pWeight",
#'                            period = "year")
#'
#' ## use stratification w.r.t. region and clustering w.r.t. households
#' dat_boot <- draw.bootstrap(
#'   eusilc, REP = 1, hid = "hid", weights = "pWeight",
#'   strata = "region", period = "year")
#'
#' ## use multi-level clustering
#' dat_boot <- draw.bootstrap(
#'   eusilc, REP = 1, hid = "hid", weights = "pWeight",
#'   strata = c("region", "hsize"), period = "year")
#'
#'
#' # create spit households
#' eusilc[, pidsplit := pid]
#' year <- eusilc[, unique(year)]
#' year <- year[-1]
#' leaf_out <- c()
#' for(y in year) {
#'   split.person <- eusilc[
#'     year == (y-1) & !duplicated(hid) & !(hid %in% leaf_out),
#'     sample(pid, 20)
#'   ]
#'   overwrite.person <- eusilc[
#'     (year == (y)) & !duplicated(hid) & !(hid %in% leaf_out),
#'     .(pid = sample(pid, 20))
#'   ]
#'   overwrite.person[, c("pidsplit", "year_curr") := .(split.person, y)]
#'
#'   eusilc[overwrite.person, pidsplit := i.pidsplit,
#'          on = .(pid, year >= year_curr)]
#'   leaf_out <- c(leaf_out,
#'                 eusilc[pid %in% c(overwrite.person$pid,
#'                                   overwrite.person$pidsplit),
#'                 unique(hid)])
#' }
#'
#' dat_boot <- draw.bootstrap(
#'   eusilc, REP = 1, hid = "hid", weights = "pWeight",
#'   strata = c("region", "hsize"), period = "year", split = TRUE,
#'   pid = "pidsplit")
#' # split households were considered e.g. household and
#' # split household were both selected or not selected
#' dat_boot[, data.table::uniqueN(w1), by = pidsplit][V1 > 1]
#' 
#'
#' @export draw.bootstrap
#'



draw.bootstrap <- function(
    dat, method = "Preston", REP = 1000, hid = NULL, weights, period = NULL, strata = NULL,
    cluster = NULL, totals = NULL, single.PSU = c("merge", "mean"), boot.names =
      NULL, split = FALSE, pid = NULL, seed = NULL) {
  
  # browser()
  
  occurence_first_period <- NULL
  
  if (method == "Rao-Wu") {
    warning("The 'Rao-Wu' method should only be used when the first stage sampling is conducted without replacement. Ensure that your sampling design follows this requirement.")
  }
  
  # INPUT CHECKING
  if (is.data.frame(dat)) {
    dat <- as.data.table(dat)
  }else if (!is.data.table(dat)) {
    stop("dat must be a data.frame or data.table")
  }
  dat <- copy(dat)

  c.names <- colnames(dat)

  removeCols <- c()

  # check seed
  if(!is.null(seed)){
    check.input(seed, input.name = "seed", input.length=1, input.type="numeric")
    set.seed(seed)
  }

  # check REP
  check.input(REP, input.name = "REP", input.length=1, input.type="numeric",
              decimal.possible = FALSE)

  # check hid
  hidNULL <- is.null(hid)
  if (hidNULL) {
    hid <- generateRandomName(20, colnames(dat))
    dat[, c(hid) := 1:.N]
    removeCols <- c(removeCols, hid)
  }
  check.input(hid, input.name = "hid", input.length=1, input.type="character",
              c.names = c.names)

  # check pid
  pidNULL <- is.null(pid)
  if (pidNULL) {
    pid <- generateRandomName(20, colnames(dat))
    dat[, c(pid) := paste(hid,1:.N, sep="."), by=c(hid)]
    removeCols <- c(removeCols, pid)
  }
  check.input(pid, input.name = "pid", input.length=1, input.type="character",
              c.names = c.names)

  # check weights
  check.input(weights, input.name = "weights", input.length=1, input.type="character",
              c.names = c.names, dat = dat, dat.column.type = "numeric")


  # check period
  periodNULL <- is.null(period)
  if (periodNULL) {
    period <- generateRandomName(20, colnames(dat))
    dat[, c(period) := 1]
    removeCols <- c(removeCols, period)
  }

  check.input(period, input.name = "period", input.length=1, input.type="character",
              c.names = c.names, dat = dat, dat.column.type = "numeric")

  # check design
  strataNULL <- FALSE
  if (is.null(strata)) {
    strata <- "I"
    strataNULL <- TRUE
  }

  clusterNULL <- FALSE
  if (is.null(cluster)) {
    cluster <- hid
    clusterNULL <- TRUE
  } else {
    if (length(cluster) == 1) {
      if (cluster %in% c("1", "I")) {
        cluster <- hid
      }
    }
    if (!hid %in% cluster & !pid %in% cluster) {
      if(pidNULL==FALSE){
        cluster <- c(cluster, pid)
      }else{
        cluster <- c(cluster, hid)
      }
    }
  }

  if (!all(strata[!strata %in% c("1", "I")] %in% c.names)) {
    stop("Not all elements in strata are column names in dat")
  }
  if (any(!cluster[!cluster %in% c("1", "I")] %in% c.names)) {
    stop("Not all names in cluster are column names in dat")
  }
  if (any(!totals %in% c.names)) {
    stop("Not all names in totals are column names in dat")
  }

  # check for missing values
  spec.variables <- c(hid, weights, period, strata, cluster, totals, pid)
  spec.variables <- spec.variables[!spec.variables %in% c("1", "I")]
  dat.na <- dat[, mget(spec.variables)]
  dat.na <- sapply(dat.na, function(z) {
    any(is.na(z))
  })
  if (any(dat.na)) {
    stop("Missing values found in column(s): ",
         paste(names(dat.na[dat.na == TRUE]), collapse = ", "))
  }

  if (length(cluster) > 1) {
    if (length(cluster) != length(strata)) {
      stop("strata and cluster need to have the same number of stages!\n ",
           "Please use either '1' or 'I' if there was no clustering or ",
           "stratification in one of the stages.")
    }
  }else{

    if (length(strata) > 1) {
      if (any(c("1", "I") %in% strata)) {
        stop("When defining multiple strata variables for single stage",
             " sampling design\n none of them can be '1' or 'I'.")
      }
      strata_var_help <- generateRandomName(20, colnames(dat))
      dat[, c(strata_var_help) := do.call(paste, c(.SD, sep = "-")), .SDcols=c(strata)]
      strata <- strata_var_help
      removeCols <- c(removeCols, strata)
    }
  }

  # check single.PSUs
  single.PSU <- single.PSU[1]
  if (is.null(single.PSU) || !single.PSU %in% c("merge", "mean")) {
    warning("single.PSU was not set to either 'merge' or 'mean'!\n Bootstrap",
            " replicates for single PSUs cases will be missing!")
    single.PSU <- FALSE
  }
  
  # check boot.names
  if (!is.null(boot.names)) {
    if (!grepl("^[[:alpha:]]", boot.names)) {
      stop("boot.names must start with an alphabetic character")
    }
  }
  
  # check split and pid
  if (is.null(split)) {
    stop("split needs to be logical")
  }
  if (!is.logical(split)) {
    stop("split needs to be logical")
  }
  if (split) {
    if (!is.character(pid)) {
      stop("when split is TRUE pid needs to be a string")
    } else {
      if (length(pid) > 1) {
        stop("pid can only have length 1")
      } else {
        if (!pid %in% c.names) {
          stop(pid, "is not a column of dat")
        }
      }
    }
    # check if pid is unique in each household and period
    unique.pid <- dat[,uniqueN(pid)==.N, by=c(period, hid),
                      env = list(pid = pid)]
    unique.pid <- unique.pid[V1==FALSE]
    if (nrow(unique.pid) > 0) {
      stop("pid is not unique in each household for each period")
    }
  }
  
  if (is.null(totals)) {
    if (length(cluster) == 1) {
      
      totals <- generateRandomName(20, existingNames = colnames(dat))
      fpc.strata <- strata[!strata %in% c("I", "1")] # nolint
      dat[,c(totals) := sum(weights[!duplicated(hid)]),
          by = c(fpc.strata, period),
          env = list(weights = weights,
                     hid = hid)]
      removeCols <- c(removeCols, totals)
    } else {
      
      stop("For multistage sampling the number of PSUs at each level needs to ",
           "be specified!")
    }
    
  } else {
    
    if (length(totals) != length(strata)) {
      stop("totals must be specified for each stage")
    }
    if (any(!totals %in% c.names)) {
      stop("Not all elements in totals are column names in dat")
    }
    if (!any(unlist(dat[, lapply(.SD, is.numeric), .SDcols = c(totals)]))) {
      stop("Not all elements in totals are numeric columns in dat")
    }
  }
  
  # check for each stage that PSUs are not in mutiple strata
  for(i in seq_along(strata)){
    if(!strata[i]%in%c("1","I") & !cluster[i]%in%c("1","I")){
      countMultiple <- dat[,uniqueN(strata_i), by=c(cluster[i], period),
                           env = list(strata_i = strata[i])]
      countMultiple <- countMultiple[V1>1]
      if(nrow(countMultiple)>0){
        stop("Some sampling units in ",cluster[i]," occur in multiple strata of ",strata[i])
      }
    }
  }
  
  # define sample design
  strata_design <- paste(strata, collapse = ">")
  cluster_design <- paste(cluster, collapse = ">")
  totals_design <- paste(totals, collapse = ">")
  
  if (is.null(boot.names)) {
    w.names <- paste0("w", 1:REP)
  } else{
    w.names <- paste0(boot.names, 1:REP)
  }
  
  # calculate bootstrap replicates
  periods <- sort(unique(dat[[period]]))
  dat_selection_prev <- NULL
  for(p in periods){

      dat_boot <- rescaled.bootstrap(dat[get(period) == p],
                                   method = method, 
                                   REP = REP, strata = strata_design, cluster = cluster_design,
                                   fpc = totals_design, single.PSU = single.PSU, return.value = c("replicates","selection"),
                                   already.selected = dat_selection_prev,
                                   run.input.checks = FALSE)
    dat[period == p, c(w.names) := dat_boot$replicates, env = list(period = period)]
    dat_selection_prev <- dat_boot$selection
    rm(dat_boot)
  }
  
  # respect split households
  if (split) {
    dat <- generate.HHID(dat, period = period, pid = pid, hid = hid)
    
    help_survey_grouping <- function(x){
      x_holes <- which(diff(x)>1)
      x_groups <- diff(c(0,x_holes,length(x)))
      x_groups <- rep(1:length(x_groups),times=x_groups)
      return(x_groups)
    }
    dat[,hid_survey_segment_help := help_survey_grouping(period), by=c(hid),
        env = list(period = period)]
    dat[,occurence_first_period := min(period), by=c(hid,"hid_survey_segment_help"),
        env = list(period = period)]
    # dat[, hid_survey_segment_help := help_survey_grouping(get(period)), by = hid]
    # dat[, occurence_first_period := min(get(period)), 
    #     by = .(hid, hid_survey_segment_help)]
    
    select.first.occurence <- paste0(c(hid,"hid_survey_segment_help", w.names), collapse = ",")
    dat.first.occurence <- unique(dat[period==occurence_first_period, .SD,
                                      .SDcols = c(hid,"hid_survey_segment_help", w.names),
                                      env = list(period = period)])
    # select_cols <- c(hid, "hid_survey_segment_help", w.names)
    # dat.first.occurence <- unique(
    #   dat[get(period) == occurence_first_period, 
    #       ..select_cols] 
    # )
    
    dat[, c(w.names) := NULL]
    dat <- merge(dat, dat.first.occurence, by = c(hid, "hid_survey_segment_help"), all.x = TRUE)
    dat[, c("occurence_first_period","hid_survey_segment_help") := NULL]
    dat[, hid := paste0(hid,"_orig"), env = list(hid = hid)]
    # dat[, (hid) := paste0(get(hid), "_orig")]
    dat[, c(paste0(hid, "_orig")) := NULL]
  }
  
  
  if (length(removeCols) > 0) {
    dat[, c(removeCols) := NULL]
  }
  
  
  if (periodNULL) {
    period <- NULL
  }
  if (hidNULL) {
    hid <- NULL
  }
  if (clusterNULL){
    cluster <- NULL
  }
  if (strataNULL){
    strata <- NULL
  }
  if (pidNULL){
    pid <- NULL
  }
  
  setattr(dat, "weights", weights)
  setattr(dat, "period", period)
  setattr(dat, "b.rep", w.names)
  setattr(dat, "hid", hid)
  setattr(dat, "pid", pid)
  setattr(dat, "cluster", cluster)
  setattr(dat, "strata", strata)
  setattr(dat, "totals", totals)
  
  
  return(dat)
}
