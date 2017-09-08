#' recalib
#' 
#' Recalibrate bootstrap weights by using iterative proportional updating to match population totals on various household and personal levels.
#'   
#' 
#' @param dat either data.frame or data.table containing the sample survey for various years.
#' @param hid string specifying the name of the column in \code{dat} containing the household ID. 
#' @param weights string specifying the name of the column in \code{dat} containing the sample weights. 
#' @param b.weights string specifying the names of the columns in \code{dat} containing bootstrap weights which should be recalibratet 
#' @param year string specifying the name of the column in \code{dat} containing the sample years. 
#' @param conP.var character vector containig person-specific variables for which contingency tables for the population tables are calculatet per year. 
#' @param conH.var character vector containig household-specific variables for which contingency tables for households are calculatet per year.
#' @param ... additional arguments passed on to function \code{\link{ipu2}} from the simPop package.
#' @return the survey data with the number of REP bootstrap replicates added as columns. 
#' @author Johannes Gussenbauer, Statistik Austria
#' 
#' @examples
#'
#' @export recalib
#'

recalib <- function(dat,hid="hid",weights="hgew",b.weights=paste0("w",1:1000),year="jahr",conP.var=c("ksex","kausl","al","erw","pension"),
										conH.var=c("bundesld","hsize","recht"),...){
	
	# define default values for 
	getElement("verbose",TRUE,...)
	getElement("epsP",1e-2,...)
	getElement("epsH",5e-2,...)
	getElement("bound",4,...)
	getElement("maxIter",50,...)
	getElement("meanHH",TRUE,...)
	
	# calculate contingency tables
	if(!is.null(conP.var)){
		conP <- lapply(conP.var,function(z){
			form.z <- paste0("V1~",paste(year,z,sep="+"))
			dt.eval("xtabs(",form.z,",data=dat[,sum(",weights,"),by=list(",year,",",z,")])")
		})
	}else{
		conP <- NULL
	}
	if(!is.null(conH.var)){
		
		dat[,ind_p:=c(1L,rep(0,.N-1)),by=c(hid,year)]
		
		conH <- lapply(conH.var,function(z){
			form.z <- paste0("V1~",paste(year,z,sep="+"))
			dt.eval("xtabs(",form.z,",data=dat[,sum(ind_p*",weights,"),by=list(",year,",",z,")])")
		})
	}else{
		conH <- NULL
	}
	
	
	# define new Index
	new_id <- paste(hid,year,sep=",")
	dt.eval("dat[,hidf:=paste0(",new_id,")]")
	
	# calibrate weights to conP and conH
	select.var <- c("hidf",weights,year,conP.var,conH.var)
	for(g in b.weights){
		set(dat,j=g,value=dt.eval("dat[,",g,"*",weights,"]"))
		set(dat,j=g,value=ipu2(dat=copy(dat[,mget(c(g,select.var))]),conP=conP,
										 conH=conH,verbose=verbose,epsP=epsP,epsH=epsH,
										 w=g,bound=bound,maxIter=maxIter,meanHH=,meanHH,hid="hidf")[,calibWeight])
	}
	
	return(dat)
}

