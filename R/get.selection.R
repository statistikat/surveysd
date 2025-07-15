#' Title
#'
#' @param dat
#' @param bootRep
#' @param cluster
#' @param strata
#' @param pid
#' @param hid
#'
#' @returns
#' @export
#'
#' @examples
#'
get.selection <- function(dat, bootRep = attr(dat,"b.rep"),
                          cluster = attr(dat,"cluster"),
                          strata = attr(dat,"strata"),
                          pid = attr(dat,"pid"),
                          hid = attr(dat,"hid") ){


  # -------------------------------
  # check inputs

  stages <- length(strata)
  clusterNULL <- hidNULL <- pidNULL <- strataNULL <- FALSE
  if(is.null(cluster)){
    cluster <- generateRandomName(existingNames = colnames(dat))
    dat[,c(cluster):= 1]
    clusterNULL <- TRUE
  }
  if(is.null(strata)){
    strata <- generateRandomName(existingNames = colnames(dat))
    dat[,c(strata):= 1]
    strataNULL <- TRUE
  }
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


  # --------------------------------
  # melt data for getting deltas
  dat_melt <- melt(dat, id.vars=c(hid,pid,strata,cluster), measure.vars = bootRep,
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
      dat_melt_i <- dat_melt[,.(ReplicateValues=max(ReplicateValues)), by=c(by_col)]
      dat_melt_i[, c(selection_i) := ReplicateValues>1]

      dat_melt_i[,c(problem_i) := floor(.N/2)!=sum(selection_col), by=c(strata[i],"ReplicateNames"),
                 env=list(selection_col = selection_i)]
    }else{
      dat_melt_i <- dat_melt[,.(ReplicateValues=max(ReplicateValues)), by=c(by_col,selection_cols[1:(i-1)])]
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
             "This function only accepts replicate weights drawn in the same fashion as surveysd::draw.bootstrap() or surveysd::rescaled.bootstrap() !")
      }else{
        dat_problem <- unique(dat_problem[,.SD,.SDcols=c(print_cols,"ReplicateNames")])
        stop("Uneven amount of unique replicate values found in \n",
             paste(capture.output(print(dat_problem, row.names = FALSE)),collapse="\n"),
             "\nThis function only accepts replicate weights drawn in the same fashion as surveysd::draw.bootstrap() or surveysd::rescaled.bootstrap() !")
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
    dat_selection_i <- dcast(dat_melt, formula = form_cast, value.var=selection_i)
    setnames(dat_selection_i, bootRep, paste0("delta_",i,"_",1:length(bootRep)))

    dat_selection <- c(dat_selection, list(dat_selection_i))
  }

  names(dat_selection) <- paste0("SamplingStage",1:stages)
  return(dat_selection)
}



