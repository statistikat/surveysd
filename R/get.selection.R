#' @title Get sample selection (~deltas) from drawn bootstrap replicates 
#'
#' @description
#' Reconstruct sample selection, e.g. record was drawn or not drawn (delta = 0/1)
#' in each sampling stage from bootstrap replicates.
#' `get.selection()` needs the `cluster`, `strata` and `hid`/`pid` information (if not `NULL`)
#' to correctly reconstruct if a record was drawn in each sampling stage for each bootstrap replicate.
#' Is only needed if bootstrap replicates are drawn for a survey 
#' with existing bootstrap replicates from a previous period, 
#' see parameter `already.selected` in function [draw.bootstrap()].
#' 
#' @param dat either data.frame or data.table containing the survey data with
#'   rotating panel design. Should contain only survey data from a single time period.
#' @param b.rep character specifying the names of the columns in `dat`
#'   containing bootstrap replicates.
#' @param strata character vector specifying the name(s) of the column in `dat`
#'   by which the population was stratified.
#' @param cluster character vector specifying cluster in the data. 
#' @param hid character specifying the name of the column in `dat` containing
#'   the household id. If `NULL` (the default), the household structure is not
#'   regarded. `hid` and `pid` cannot both be `NULL`.
#' @param pid pid column in `dat` specifying the personal identifier. This
#'   identifier needs to be unique for each person throught the whole data set. 
#'   `hid` and `pid` cannot both be `NULL`.
#'   
#' @returns Returns a list of data.tables. 
#' The length of the list equals the number of sampling stages specified.
#' Each list entry contains a `data.table` with variables for sampling stage and/or
#' `hid`/`pid` as well as `length(attr(dat,"b.rep"))` columns each indicating if 
#' record/cluster was drawn in the respective sampling stage for the i-th boostrap replicate. 
#' @export
#'
#' @examples
#' 
#' library(surveysd)
#' library(data.table)
#' setDTthreads(1)
#' set.seed(1234)
#' eusilc <- demo.eusilc(n = 3, prettyNames = TRUE)
#' 
#' ## draw replicates with stratification
#' dat_boot <- draw.bootstrap(eusilc[year<2012], REP = 3, weights = "pWeight",
#'                            strata = "region", hid = "hid",
#'                            period = "year")
#' 
#' ## get selection matrix for year 2011 
#' dat_selection <- get.selection(dat_boot[year==2011])
#' print(dat_selection)
#' 
#' ## draw bootstrap replicates for year 2012
#' ## respecting already selected units for year 2011 ~ dat_selection
#' ## in order to mimic rotating panel design
#' dat_boot_2012 <- draw.bootstrap(eusilc[year==2012], REP = 3, weights = "pWeight",
#'                                 strata = "region", hid = "hid",
#'                                 period = "year", 
#'                                 already.selected = dat_selection)
#' 
#' 
#' 
get.selection <- function(dat, b.rep = attr(dat,"b.rep"),
                          strata = attr(dat,"strata"),
                          cluster = attr(dat,"cluster"),
                          hid = attr(dat,"hid"),
                          pid = attr(dat,"pid")){
  
  
  . <- ReplicateValues <- selection_col <- selection_prev <- selection_col_prev <-
    problem_col <- ReplicateNames <- c
  
  # -------------------------------
  # check inputs
  
  stages <- length(strata)
  if(is.null(strata)){
    stages <- 1
  }
  
  clusterNULL <- hidNULL <- pidNULL <- strataNULL <- FALSE
  if(is.null(hid)){
    hid <- generateRandomName(existingNames = colnames(dat))
    dat[,c(hid):= 1:.N]
    hidNULL <- TRUE
  }
  if(is.null(pid)){
    pid <- generateRandomName(existingNames = colnames(dat))
    dat[,c(pid):= 1:.N,by=c(hid)]
    pidNULL <- TRUE
  }
  if(is.null(cluster)){
    cluster <- generateRandomName(existingNames = colnames(dat))
    # if cluster is null default to hid
    # if hid is null -> 1:.N
    dat[,c(cluster):= hid, env = list(hid = hid)] 
    clusterNULL <- TRUE
  }
  if(is.null(strata)){
    strata <- generateRandomName(existingNames = colnames(dat))
    dat[,c(strata):= 1]
    strataNULL <- TRUE
  }

  if(all(hidNULL,pidNULL)){
    stop("hid and pid should not both be NULL on order to 
    link the selection variables to the original data!
         \nPlease specify either hid or pid")
  }
  
  
  # --------------------------------
  # melt data for getting deltas
  dat_melt <- melt(dat, id.vars=unique(c(hid,pid,strata,cluster)), measure.vars = b.rep,
                   variable.name = "ReplicateNames", value.name = "ReplicateValues")
  help_selection <- function(v, selection_prev=NULL){
    
    v_out <- rep(0,length(v))
    if(!is.null(selection_prev)){
      v_cutoff <- v[selection_prev==1]
      v_out[v>=v_cutoff] <- 1
      
    }else{
      v_out <- max(v)>1
    }
    
    return(v_out)
  }
  dat_selection <- list()
  selection_cols <- paste0("SelectionStage_",1:stages,"_",generateRandomName(existingNames = colnames(dat_melt), nchar = 10))
  problem_cols <- paste0("ProblemStage_",1:stages,"_",generateRandomName(existingNames = colnames(dat_melt), nchar = 10))
  
  for(i in 1:stages){
    
    # per stage define selected and not selected PSUs
    # get mean of unqiue values 
    # greater than 1 -> selected
    # lower than 1 -> not selected
    by_col <- c(strata[1:i],cluster[1:i],"ReplicateNames")
    selection_i <- selection_cols[i]
    problem_i <- problem_cols[i]
    if(i==1){
      dat_melt_i <- dat_melt[,.(ReplicateValues=unique(ReplicateValues)), by=c(by_col)]
      dat_melt_i[, c(selection_i) := ReplicateValues>1]
      
      dat_melt_i[,c(problem_i) := floor(.N/2)!=sum(selection_col), by=c(strata[i],"ReplicateNames"), 
                 env=list(selection_col = selection_i)]
    }else{
      dat_melt_i <- dat_melt[,.(ReplicateValues=unique(ReplicateValues)), by=c(by_col,selection_cols[1:(i-1)])]
      dat_melt_i[, c(selection_i) := 0]
      dat_melt_i[selection_prev>0 , c(selection_i) := ReplicateValues>mean(ReplicateValues), 
                 by=c(strata[1:(i-1)],cluster[1:(i-1)]), 
                 env = list(selection_prev = selection_cols[i-1])]
      
      dat_melt_i[selection_col_prev == 1,c(problem_i) := floor(.N/2)!=sum(selection_col), 
                 by=c(strata[i],"ReplicateNames"), 
                 env=list(selection_col = selection_i,
                          selection_col_prev = selection_cols[i-1])]
    }
    
    dat_problem <- dat_melt_i[problem_col == TRUE, env=list(problem_col = problem_i)]
    if(nrow(dat_problem)>0){
      # should not happen - faulty input
      # print error
      print_cols <- c(cluster[i][clusterNULL==FALSE],
                      strata[i][strataNULL==FALSE])
      rep_number <- dat_problem[,unique(ReplicateNames)]
      if(is.null(print_cols)){
        example_rep <- paste(rep_number[1:min(5,length(rep_number))],collapse=", ")
        if(length(rep_number)>5){
          example_rep <- paste0(example_rep,", ...")
        }
        stop("Uneven amount of unique replicate values for\n",
             uniqueN(rep_number)," replicate weights (",example_rep,")\n",
             "This function only accepts replicate weights drawn in the same fashion as 
             surveysd::draw.bootstrap(..., method = 'Preston', ...) or surveysd::rescaled.bootstrap(..., method = 'Preston', ...) !")
      }else{
        dat_problem <- unique(dat_problem[,.SD,.SDcols=c(print_cols,"ReplicateNames")])
        stop("Uneven amount of unique replicate values found in \n",
             paste(capture.output(print(dat_problem, row.names = FALSE)),collapse="\n"),
             "\nThis function only accepts replicate weights drawn in the same fashion as 
             surveysd::draw.bootstrap(..., method = 'Preston', ...) or surveysd::rescaled.bootstrap(..., method = 'Preston', ...) !")
      }
    }
    
    # merge dat_melt_i to dat_melt
    dat_melt[dat_melt_i,c(selection_i):=selection_col, on=c(by_col),
             env = list(selection_col = selection_i)]
    
    id_cols <- c()
    if(strataNULL==FALSE){
      id_cols <- c(id_cols, strata)
    }
    if(clusterNULL==FALSE){
      id_cols <- c(id_cols, cluster)
    }
    if(hidNULL==FALSE){
      id_cols <- c(id_cols, hid)
    }
    if(pidNULL==FALSE){
      id_cols <- c(id_cols, pid)
    }
    id_cols <- unique(id_cols)

    form_cast <- as.formula(paste0(paste(id_cols,collapse="+"),"~ReplicateNames"))
    dat_melt <- unique(dat_melt, by=c(id_cols, "ReplicateNames"))
    dat_selection_i <- dcast(dat_melt, formula = form_cast, value.var=selection_i)
    setnames(dat_selection_i, b.rep, paste0("delta_",i,"_",1:length(b.rep)))
    
    dat_selection <- c(dat_selection, list(dat_selection_i))
  }
  
  names(dat_selection) <- paste0("SamplingStage",1:stages)
  return(dat_selection)
}



