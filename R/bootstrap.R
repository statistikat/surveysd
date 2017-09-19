#' bootstrap.rep
#'
#' Draw bootstrap replicates from survey data taken in various years.
#' Survey information, like ID, sample weights, strata and population totals per strata, should be specified to ensure meaningfull survey bootstraping.
#'
#'
#' @param dat either data.frame or data.table containing the sample survey for various years.
#' @param REP integer indicating the number of bootstrap replicates.
#' @param hid string specifying the name of the column in \code{dat} containing the household ID.
#' @param weights string specifying the name of the column in \code{dat} containing the sample weights.
#' @param strata string vector specifying the name of the column in \code{dat} by which the population was stratified.
#' @param year string specifying the name of the column in \code{dat} containing the sample years.
#' @param totals string specifying the name of the column in \code{dat} containing the the totals per strata. If totals is NULL, the sum of weights per strata will be calcualted and named 'fpc'.
#' @param boot.names character indicating the leading string of the column names for each bootstrap replica. If NULL defaults to "w".
#' @return the survey data with the number of REP bootstrap replicates added as columns.
#' @author Johannes Gussenbauer, Statistik Austria
#'
#' @examples
#'
#' @export bootstrap.rep
#'


bootstrap.rep <- function(dat,REP=1000,hid="hid",weights="hgew",strata="bundesld",year="jahr",totals=NULL,boot.names=NULL){

  if(class(dat)[1]=="data.frame"){dat <- as.data.table(dat)}

 	if(is.null(totals)){
		totals <- "fpc"
    #dat[,fpc:=sum(hgew),by=c(strata)]
		dt.eval("dat[,fpc:=sum(",weights,"),by=list(",paste(strata,collapse=","),")]")
		#eval(parse(text=paste0("dat[,fpc:=sum(",weights,"),by=list(",paste(strata,collapse=","),")]")))
		#test <- dt.eval("dat[,sum(",weights,"),by=list(",paste(strata,collapse=","),")]")

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
	# optional:
	#	dt.eval("dat[,c(w.names):=gen.boot(.SD,REP=",REP,",hid=",hid,",weights=",weights,",strata=",strata,",totals=",totals,"),by=c(year)]")

	# keep bootstrap replicates of first year for each household
	dat[,c(w.names):=.SD[year==min(year),mget(w.names)],by=c(gsub("~","",hid)[2])]

	return(dat)
}


gen.boot <- function(dat,REP=1000,hid="hid",weights="hgew",strata="bundesld",totals="fpc"){

	dat.svy <- svydesign(ids=hid,weights=weights,fpc=totals,strata=strata,data=dat)
	dat.svyboot <- as.svrepdesign(dat.svy,type="mrbbootstrap",replicates=REP,compress=FALSE)$repweights

	return(data.table(dat.svyboot))
}

