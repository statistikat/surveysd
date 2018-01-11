#' @title Draw bootstrap replicates
#'
#' @description Draw bootstrap replicates from survey data with rotating panel design.
#' Survey information, like ID, sample weights, strata and population totals per strata, should be specified to ensure meaningfull survey bootstraping.
#'
#' @usage draw.bootstrap(dat,REP=1000,hid="hid",weights="hgew",strata="bundesld",
#'                      year="jahr",totals=NULL,boot.names=NULL)
#'
#' @param dat either data.frame or data.table containing the survey data with rotating panel design.
#' @param REP integer indicating the number of bootstrap replicates.
#' @param hid character specifying the name of the column in \code{dat} containing the household ID.
#' @param weights character specifying the name of the column in \code{dat} containing the sample weights.
#' @param strata character vector specifying the name of the column in \code{dat} by which the population was stratified.
#' @param year character specifying the name of the column in \code{dat} containing the sample years.
#' @param country character specifying the name of the column in \code{dat} containing the country name. Is only used if \code{dat} contains data from multiple countries.
#' In this case the bootstep procedure will be applied on each country seperately. If \code{country=NULL} the household identifier must be unique for each household.
#' @param cluster character vector specifying cluster in the data. If \code{NULL} household ID is taken es the lowest level cluster.
#' @param totals (optional) character specifying the name of the column in \code{dat} containing the the totals per strata and/or cluster. If totals and cluster is \code{NULL}, the households per strata will be calcualted using the \code{weights} argument and named 'fpc'.
#' If clusters are specified then totals need to be supplied by the user, otherwise they will be set to \code{NULL}. When multiple cluster and or strata are specified totals needs to contain multiple argument each corresponding to a column name in \code{dat}.
#' Each column needs to contains the total number if units in the population regarding the subsequent level. The vector is interpreted from left to right meaning that the most left value of \code{totals}
#' specifies the column names with the number of units in the population at the highest level and the most right value specifies the column names with the number of units in the population at the lowest level.
#' This argument will be passed onto the function \code{svydesign()} from package \code{survey} through the argument \code{fpc}.
#' @param boot.names character indicating the leading string of the column names for each bootstrap replica. If NULL defaults to "w".
#' @return the survey data with the number of REP bootstrap replicates added as columns.
#'
#' @details \code{draw.bootstrap} takes \code{dat} and draws \code{REP} bootstrap replicates from it.
#' \code{dat} must be household data where household members correspond to multiple rows with the same household identifier. The data should at least containt the following columns:
#' \itemize{
#'   \item Column indicating the sample year;
#'   \item Column indicating the household ID;
#'   \item Column containing the household sample weights;
#'   \item Columns by which population was stratified during the sampling process.
#' }
#' A column for the totals in each strat can be included, but is only optional. If it is not included, e.g \code{totals=NULL}, this column will be calculated and added to \code{dat} using \code{strata} and \code{weights}.\cr
#' The bootstrap replicates are drawn for each survey year (\code{year}) using the function \code{\link[survey]{as.svrepdesign}} from the package \code{survey}.
#' Afterwards the bootstrap replicates for each household are carried forward from the first year the household enters the survey to all the censecutive years it stays in the survey.\cr
#' This ensures that the bootstrap replicates follow the same logic as the sampled households, making the bootstrap replicates more comparable to the actual sample units.
#'
#' @return Returns a data.table containing the original data as well as the number of \code{REP} columns containing the bootstrap replicates for each repetition.\cr
#' The columns of the bootstrap replicates are by default labeled "w\emph{Number}" where \emph{Number} goes from 1 to \code{REP}.
#' If the column names of the bootstrap replicates should start with a different character or string the parameter \code{boot.names} can be used.
#'
#' @seealso \code{\link[data.table]{data.table}} for more information on data.table objects.\cr
#' \code{\link[survey]{svydesign}} for more information on how to create survey-objects.\cr
#' \code{\link[survey]{as.svrepdesign}} for more information on how bootstrap replicates are drawn from survey-objects.
#'
#' @author Johannes Gussenbauer, Alexander Kowarik, Statistics Austria
#'
#' @examples
#' # read in data (must be changed..)
#' dat <- data.table(read_sas("PATH"))
#'
#' # create 20 bootstrap replicates using the column "bundesld" as strata
#' dat_boot <- draw.bootstrap(dat=copy(dat),REP=20,hid="hid",weights="hgew",
#'                           strata="bundesld",year="jahr")
#'
#' # do the same with more strata
#' dat_boot <- draw.bootstrap(dat=copy(dat),REP=20,hid="hid",weights="hgew",
#'                           strata=c("bundesld","sex","hsize"),year="jahr")
#'
#' # change column names for bootstrap replicates
#' dat_boot <- draw.bootstrap(dat=copy(dat),REP=20,hid="hid",weights="hgew",
#'                           strata=c("bundesld"),year="jahr",boot.names="replicate")
#'
#' # save bootstrap replicates as .RData
#' save(dat_boot,file="dat_replicates.RData")
#' # or .csv-file
#' write.csv2(dat_boot,file="dat_replicates.csv",row.names=FALSE)
#'
#' @export draw.bootstrap
#' @import survey data.table


draw.bootstrap <- function(dat,REP=1000,hid,weights,strata=NULL,year,country=NULL,cluster=NULL,totals=NULL,boot.names=NULL){

  ##########################################################
  # INPUT CHECKING
  if(class(dat)[1]=="data.frame"){
    dat <- as.data.table(dat)
  }else if(class(dat)[1]!="data.table"){
    stop("dat must be a data.frame or data.table")
  }


  c.names <- colnames(dat)

  # check REP
  if(length(REP)!=1){
    stop("REP must have length 1")
  }
  if(!is.numeric(REP)){
    stop("REP must contain one numeric value")
  }
  if(REP%%1!=0){
    stop("REP cannot have a decimal part")
  }

  # check hid
  if(length(hid)!=1){
    stop("hid must have length 1")
  }
  if(!hid%in%c.names){
    stop(paste0(hid," is not a column in dat"))
  }

  # check weights
  if(length(weights)!=1){
    stop("weights must have length 1")
  }
  if(!weights%in%c.names){
    stop(paste0(weights," is not a column in dat"))
  }
  if(!is.numeric(dt.eval("dat[,",weights,"]"))){
    stop(paste0(weights," must be a numeric column"))
  }

  # check year
  if(length(year)!=1){
    stop("year must have length 1")
  }
  if(!year%in%c.names){
    stop(paste0(year," is not a column in dat"))
  }

  # check strata
  if(!is.null(strata)){
    if(!all(strata%in%c.names)){
      stop("Not all elements in strata are column names in dat")
    }
  }

  # check design
  if(!is.null(cluster)){
    if(any(!cluster%in%c.names)){
      stop("Not all names in cluster are column names in dat")
    }
  }

  # check country
  if(!is.null(country)){
    if(length(country)!=1){
      stop("country must have length 1")
    }
    if(!country%in%c.names){
      stop(paste0(country," is not a column in dat"))
    }
  }

  # check boot.names
  if(!is.null(boot.names)){
    if(!grepl("^[[:alpha:]]",boot.names)){
      stop("boot.names must start with an alphabetic character")
    }
  }


  # check totals
  # if clusters are specified the finite population correction factors must be user specified (at least for now)
  # check input for totals
  # if no totals are specified then leave them NULL
  if(is.null(totals)){
    if(is.null(cluster)){
      # if no clusters are specified the use calculate number of households in each strata
      totals <- "fpc"
      dt.eval("dat[,fpc:=sum(",weights,"[!duplicated(",hid,")]),by=list(",paste(c(strata,country),collapse=","),")]")
      add.totals <- TRUE
      optwarn <- FALSE
    }else{
      # else leave totals NULL
      add.totals <- FALSE
      cat("Number of Clusters at each level are not specified\nDesign is sampled with replacement\n")
      options(warn=-1)
      optwarn <- TRUE
      # totals <- "fpc"
      # dt.eval("dat[,fpc:=sum(",weights,"[!duplicated(",hid,")]),by=list(",paste(c(strata,cluster,country),collapse=","),")]")
      # add.totals <- TRUE
      # optwarn <- FALSE
    }
  }else{

    if(length(totals)!=length(c(strata,cluster))){
      stop("totals must specified for each stage")
    }
    if(any(!totals%in%c.names)){
      stop("Not all elements in totals are column names in dat")
    }
    if(!any(unlist(dat[,lapply(.SD,is.numeric),.SDcols=c(totals)]))){
      stop("Not all elements in totals are numeric columns in dat")
    }

    add.totals <- FALSE
    optwarn <- FALSE
  }
  ##########################################################

  # make arguments usable for survey package
  strata <- as.formula(paste0("~",paste(strata,collapse=":")))
  weights <- as.formula(paste0("~",weights))
  if(!is.null(totals)){
    totals <- as.formula(paste0("~",paste(totals,collapse="+")))
  }
  # define sample design
  if(is.null(cluster)){
    cluster <- as.formula(paste0("~",hid))
  }else{
    cluster <- as.formula(paste0("~",paste(c(cluster,hid),collapse = "+")))
  }

  if(is.null(boot.names)){
    w.names <- paste0("w",1:REP)
  }else{
    w.names <- paste0(boot.names,1:REP)
  }


  # calculate bootstrap replicates
  dat[,c(w.names):=gen.boot(.SD,REP=REP,cluster=cluster,weights=weights,strata=strata,totals=totals),by=c(year,country)]
  if(optwarn){
    options(warn=0)
  }

  # keep bootstrap replicates of first year for each household
  w.names.c <- paste0("'",paste(w.names,collapse="','"),"'")
  by.c <- paste(c(hid,country),collapse=",")
  dt.eval("dat[,c(",w.names.c,"):=.SD[",year,"==min(",year,"),.(",paste(w.names,collapse=","),")][1],by=list(",by.c,")]")

  # remove columns
  if(add.totals){
    dt.eval("dat[,",gsub('~','',totals)[2],":=NULL]")
  }
  return(dat)
}


gen.boot <- function(dat,REP=1000,cluster="~hid",weights="hgew",strata="bundesld",totals="fpc"){

  dat.svy <- svydesign(ids=cluster,weights=weights,fpc=totals,strata=strata,data=dat)
  dat.svyboot <- as.svrepdesign(dat.svy,type="mrbbootstrap",replicates=REP,compress=FALSE)$repweights

  return(data.table(dat.svyboot))
}

