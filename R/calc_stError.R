#' calc.stError
#'
#' Calcualte estimates of variables per year as well as their 3 year rolling mean.
#' In addition the standard deviation per year and for each three consecutive years are calculated using the bootstrap weights.
#'
#' @param dat either data.frame or data.table containing the sample survey for various years.
#' @param weights character specifying the name of the column in \code{dat} containing the sample weights. Used to calculate yearly and 3-year estimates (defined by \code{var} and \code{fun}).
#' @param b.weights character specifying the names of the columns in \code{dat} containing bootstrap weights. Used to calculate standard error of yearly and 3-year estimates (defined by \code{var} and \code{fun}).
#' @param year character specifying the name of the column in \code{dat} containing the sample years.
#' @param var character vector containing variable names in \code{dat} on which \code{fun} shall be applied for each \code{year}. Must have the same length as \code{fun}.
#' @param fun character vector containig function-names which are to be applied on \code{var} for each \code{year}. Must have the same length as \code{var}.
#' Possible arguments are \code{weightedRatio,weightedSum,sampSize,popSize} as well as any other function which returns a double or integer and has \code{weight} as its second argument.
#' @param cross_var list of character vectors containig variables in \code{dat}. For each list entry \code{dat} will be split in subgroups according to the containing variables as well as \code{year}
#' and \code{var} and \code{fun} are calculated for each of those groups.
#' @param year.diff character vectors, defining the years for which the differenzes as well as the standard deviation of the differenzes should be computed. Each entry must have the form of \code{"year1 - year2"}. Can be NULL
#' @param year.mean integer, defining the range of years over which mean of estimates will be defined. Can be NULL.
#' @param add.arg list containing strings for additional function arguments. Must be the same length as \code{var} oder \code{fun}. Can be a list of \code{NULL}s.
#' @return data.table containing yearly and 3-year estimates of \code{fun} applied to \code{var} as well as the corresponding standard errors, which were calculated using the bootstrap weights.
#' @author Johannes Gussenbauer, Statistik Austria
#'
#' @examples
#' var <- c("povmd60","povmd60","aroseCopy","aroseCopy")
#' fun <- c("weightedRatio","weightedSum","weightedRatio","weightedSum")
#'
#' cross_var <- list(NULL,c("bundesld"),c("sex"),c("sex","bundesld"))
#'
#' tab(datbw,weights="hgew,b.weights=paste0("w",1:1000),year="jahr",var=var,fun=fun,cross_var=cross_var)
#'
#' @export calc.stError
#'


# wrapper-function to apply fun to var using weights (~weights, b.weights)
# and calculating standard devation (using the bootstrap replicates) per year and für 3-year rolling means
calc.stError <- function(dat,weights="hgew",b.weights=paste0("w",1:1000),year="jahr",var="povmd60",fun="weightedRatio",cross_var=list(NULL),year.diff=NULL,year.mean=3,add.arg=vector("list",length(var))){


	# define columns in which NAs are present (will be discarded for the evaluation)
	col_cross <- unique(unlist(cross_var))
	no.na <- unlist(dat[,lapply(.SD,function(z){any(is.na(z))}),.SDcols=col_cross])
	no.na <- colnames(dat)[colnames(dat)%in%col_cross][!no.na]

	if(!is.null(year.diff)){
	  year.diff <- strsplit(year.diff,"-")
	}


	outx <- list()
	for(i in 1:length(var)){
		out_i <- help.stError(dat=dat,year=year,var=var[i],weights=weights,b.weights=b.weights,fun=fun[i],cross_var=cross_var,year.diff=year.diff,year.mean=year.mean,no.na=no.na,add.arg=add.arg[[i]])
		out_i[,variable:=paste(var[i],fun[i],sep="_")]
		outx <- c(outx,list(out_i))
	}

	outx <- rbindlist(outx,use.names=TRUE,fill=TRUE)
	outx.names <- colnames(outx)
	form <- as.formula(paste(paste(outx.names[!outx.names%in%c("val","st.Error","variable")],collapse="+"),"variable",sep="~"))

	outx <- dcast(outx,form,value.var=c("val","st.Error"),fill=NA)
	return(outx)

}


# function to apply fun to var using weights (~weights, b.weights)
# and calculating standard devation (using the bootstrap replicates) per year and für 3-year rolling means
help.stError <- function(dat,year,var,weights,b.weights=paste0("w",1:1000),fun,cross_var,year.diff=NULL,year.mean=NULL,no.na,add.arg=NULL){

	# define names for estimates for each weight (normal weights and boostrap weights)
	# makes it easier to (sort of) verctorize expressions
	res.names=paste0("r",1:length(c(weights,b.weights)))
	# formulate function and arguments
	if(fun=="sampSize"){
		eval.fun <- c(res.names,"=.N)")
	}else if(fun=="popSize"){
		eval.fun <- paste0(res.names,"=sum(",c(weights,b.weights),")")
	}else{
	  if(!is.null(add.arg)){
	    eval.fun <- paste0(res.names,"=",fun,"(",paste(var,c(weights,b.weights),add.arg,sep=","),")")
	  }else{
	    eval.fun <- paste0(res.names,"=",fun,"(",paste(var,c(weights,b.weights),sep=","),")")
	  }

	  eval.fun <- paste0(".(",paste(eval.fun,collapse=","),")")
	}

	years <- dt.eval("dat[,unique(",year,")]")
	# formulate 3 consecutive years
	yearsList <- unlist(lapply(years[1:c(length(years)-2)],function(z){
		paste(z:c(z+2),collapse="_")
	}))

	year.diff.b <- !is.null(year.diff)

	# apply function to all elemnts of cross_var
	# apply also mean and standard deviation for estimates
	out <- lapply(cross_var,function(z){

		na.check <- z[!z%in%no.na]
		if(length(na.check)>0){
			na.eval <- paste(paste0("(!is.na(",na.check,"))"),collapse="&")
		}else{
			na.eval <- NULL
		}
		by.eval <- paste(c(year,z),collapse=",")

		# calcualte estimate
		var.est <- dt.eval("dat[",na.eval,",",eval.fun,",by=list(",by.eval,")]")
		var.est <- melt(unique(var.est,by=c(year,z)),id.vars=c(year,z),measure.vars=res.names,value.name = "V1")
		var.est[,ID:=substring(variable,2)]
		var.est[,variable:=NULL]

		# calculate mean of estimate over 'year.mean'- years
		if(!is.null(year.mean)){
		  roll.est <- var.est[,list(V1=rollMeanC(V1,k=year.mean,type="c"),V2=yearsList),by=c("ID",z)]
		  setnames(roll.est,"V2",year)
		}

		if(year.diff.b){
		  by.diff <- paste(c("ID",z),collapse=",")
		  diff.est <- lapply(year.diff,function(y){
		    y_cond <- paste(year,paste0("c(",paste(y,collapse=","),")"),sep="%in%")
		    diff.y <- dt.eval("var.est[",y_cond,",V1[",year,"==",y[1],"]-V1[",year,"==",y[2],"],by=list(",by.diff,")]")
		    diff.y[,c(year):=paste(y,collapse="-")]
		    return(diff.y)
		  })
		  diff.est <- rbindlist(diff.est)
		}

		if(!is.null(year.mean)){
		  var.est <- rbind(var.est,roll.est)
		}
		if(year.diff.b){
		  var.est <- rbind(var.est,diff.est)
		}

		sd.est <- var.est[ID!=1,.(st.Error=sd(V1)),by=c(year,z)]

		#out.z <- dt.eval("var.est[ID==1,.(",by.eval,",V1)]")
		#dt.eval("out.z[,c(year):=as.character(",year,")]")
		out.z <- merge(subset(var.est,ID==1,select=c(year,z,"V1")),sd.est,by=c(year,z))
		return(out.z)

	})

	out <- rbindlist(out,fill=TRUE,use.names=TRUE)
	setnames(out,"V1","val")
	return(out)
}


# additional help function
weightedRatio <- function(x,weightvar,x2=NULL){
	if(is.null(x2)){
		return(sum(weightvar[x==1],na.rm=TRUE)/sum(weightvar,na.rm=TRUE))
	}else{
		if(sum(weightvar[x==1],na.rm=TRUE)==0&sum(weightvar[x2==1],na.rm=TRUE)==0){
			return(0)
		}
		return(sum(weightvar[x==1],na.rm=TRUE)/sum(weightvar[x2==1],na.rm=TRUE))
	}
}

weightedSum <- function(x,w){
	sum(as.numeric(x)*w,na.rm=TRUE)
}

povmd <- function(epinc,weightvar){
	md <- weightedMedian(epinc,weightvar)*0.6
	pmd60 <- epinc<md
	sum(weightvar[pmd60])/sum(weightvar)
}
kish_fact <- function(w){
	n <- length(w)
	sqrt(n*sum(w^2)/sum(w)^2)
}
