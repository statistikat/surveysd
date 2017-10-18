#' @title Draw bootstrap replicates
#'
#' @description Draw bootstrap replicates from survey data with rotating panel design.
#' Survey information, like ID, sample weights, strata and population totals per strata, should be specified to ensure meaningfull survey bootstraping.
#'
#' @usage bootstrap.rep(dat,REP=1000,hid="hid",weights="hgew",strata="bundesld",
#'                      year="jahr",totals=NULL,boot.names=NULL)
#'
#' @param dat either data.frame or data.table containing the survey data with rotating panel design.
#' @param REP integer indicating the number of bootstrap replicates.
#' @param hid string specifying the name of the column in \code{dat} containing the household ID.
#' @param weights string specifying the name of the column in \code{dat} containing the sample weights.
#' @param strata string vector specifying the name of the column in \code{dat} by which the population was stratified.
#' @param year string specifying the name of the column in \code{dat} containing the sample years.
#' @param totals (optional) string specifying the name of the column in \code{dat} containing the the totals per strata. If totals is \code{NULL}, the sum of weights per strata will be calcualted and named 'fpc'.
#' @param boot.names character indicating the leading string of the column names for each bootstrap replica. If NULL defaults to "w".
#' @return the survey data with the number of REP bootstrap replicates added as columns.
#'
#' @details \code{bootstrap.rep} takes \code{dat} and draws \code{REP} bootstrap replicates from it.
#' \code{dat} contains household survey data, where each row corresponds to one household. In addition the following columns should be included in \code{dat}:
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
#' library(data.table)
#' dat <- data.table(read_sas("O:/B/3-AP/Analyse/sonstiges/
#'                             bundesländerschätzungen 2008-2018/daten/bldaten0816.sas7bdat"))
#'
#' # create 20 bootstrap replicates using the column "bundesld" as strata
#' dat_boot <- bootstrap.rep(dat=copy(dat),REP=20,hid="hid",weights="hgew",
#'                           strata="bundesld",year="jahr")
#'
#' # do the same with more strata
#' dat_boot <- bootstrap.rep(dat=copy(dat),REP=20,hid="hid",weights="hgew",
#'                           strata=c("bundesld","sex","hsize"),year="jahr")
#'
#' # change column names for bootstrap replicates
#' dat_boot <- bootstrap.rep(dat=copy(dat),REP=20,hid="hid",weights="hgew",
#'                           strata=c("bundesld"),year="jahr",boot.names="replicate")
#'
#' @export bootstrap.rep
#' @import survey data.table


bootstrap.rep <- function(dat,REP=1000,hid="hid",weights="hgew",strata="bundesld",year="jahr",totals=NULL,boot.names=NULL){

  if(class(dat)[1]=="data.frame"){dat <- as.data.table(dat)}

 	if(is.null(totals)){
		totals <- "fpc"
   	dt.eval("dat[,fpc:=sum(",weights,"),by=list(",paste(strata,collapse=","),")]")
	}

	# make arguments usable for survey package
	hid <- as.formula(paste0("~",hid))
	strata <- as.formula(paste0("~",paste(strata,collapse=":")))
	weights <- as.formula(paste0("~",weights))
	totals <- as.formula(paste0("~",totals))

	if(is.null(boot.names)){
		w.names <- paste0("w",1:REP)
	}else{
		w.names <- paste0(boot.names,1:REP)
	}

	# calculate bootstrap replicates
	dat[,c(w.names):=gen.boot(.SD,REP=REP,hid=hid,weights=weights,strata=strata,totals=totals),by=c(year)]

	# keep bootstrap replicates of first year for each household
	dat[,c(w.names):=.SD[year==min(year),mget(w.names)],by=c(gsub("~","",hid)[2])]

	return(dat)
}


gen.boot <- function(dat,REP=1000,hid="hid",weights="hgew",strata="bundesld",totals="fpc"){

	dat.svy <- svydesign(ids=hid,weights=weights,fpc=totals,strata=strata,data=dat)
	dat.svyboot <- as.svrepdesign(dat.svy,type="mrbbootstrap",replicates=REP,compress=FALSE)$repweights

	return(data.table(dat.svyboot))
}

