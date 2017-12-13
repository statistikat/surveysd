#' @title Calcualte point estimates and their standard errors using bootstrap weights.
#'
#' @description
#' Calculate point estimates as well as standard errors of variables in surveys. Standard errors are estimated using bootstrap weights (see \code{\link{draw.bootstrap}} and \code{\link{recalib}}).
#' In addition the standard error of an estimate can be calcualted using the survey data for 3 or more consecutive years, which results in a reduction of the standard error.
#'
#' @usage calc.stError(dat,weights="hgew",b.weights=paste0("w",1:1000),year="jahr",var="povmd60",
#'                     fun="weightedRatio",cross_var=NULL,year.diff=NULL,
#'                     year.mean=3,bias=FALSE,add.arg=NULL,size.limit=20,cv.limit=10)
#'
#'
#' @param dat either data.frame or data.table containing the survey data. Surveys can be a panel survey or rotating panel survey, but does not need to be. For rotating panel survey bootstrap weights can be created using \code{\link{draw.bootstrap}} and \code{\link{recalib}}.
#' @param weights character specifying the name of the column in \code{dat} containing the original sample weights. Used to calculate point estimates.
#' @param b.weights character vector specifying the names of the columns in \code{dat} containing bootstrap weights. Used to calculate standard errors.
#' @param year character specifying the name of the column in \code{dat} containing the sample years.
#' @param var character vector containing variable names in \code{dat} on which \code{fun} shall be applied for each sample year.
#' @param fun character specifying the function which will be applied on \code{var} for each sample year.
#' Possible arguments are \code{weightedRatio,weightedSum,sampSize,popSize} as well as any other function which returns a double or integer and uses weights as its second argument.
#' @param cross_var character vectors or list of character vectors containig variables in \code{dat}. For each list entry \code{dat} will be split in subgroups according to the containing variables as well as \code{year}.
#' The pointestimates are then estimated for each subgroup seperately. If \code{cross_var=NULL} the data will split into sample years by default.
#' @param year.diff character vectors, defining years for which the differences in the point estimate as well it's standard error is calculated. Each entry must have the form of \code{"year1 - year2"}. Can be NULL
#' @param year.mean integer, defining the range of years over which the sample mean of point estimates is additionally calcualted.
#' @param bias boolean, if TRUE the sample mean over the point estimates of the bootstrap weights is returned.
#' @param add.arg character specifying additional arguments for \code{fun}. Can be \code{NULL}.
#' @param size.limit integer defining a lower bound on the number of observations on \code{dat} in each group defined by \code{year} and the entries in \code{cross_var}.
#' Warnings are returned if the number of observations in a subgroup falls below \code{size.limit}. In addition the concerned groups are available in the function output.
#' @param cv.limit non-negativ value defining a upper bound for the standard error in relation to the point estimate. If this relation exceed \code{cv.limit}, for a point estimate, they are flagged and available in the function output.
#'
#' @details \code{calc.stError} takes survey data (\code{dat}) and returns point estimates as well as their standard Errors defined by \code{fun} and \code{var} for each sample year in \code{dat}.
#' \code{dat} must be household data where household members correspond to multiple rows with the same household identifier. The data should at least containt the following columns:
#' \itemize{
#'   \item Column indicating the sample year;
#'   \item Column indicating the household ID;
#'   \item Column containing the household sample weights;
#'   \item Columns which contain the bootstrap weights (see output of \code{\link{recalib}});
#'   \item Columns listed in \code{var} as well as in \code{cross_var}
#' }
#' For each variable in \code{var} as well as sample year the function \code{fun} is applied using the original as well as the bootstrap sample weights.\cr
#' The point estimate is then selected as the result of \code{fun} when using the original sample weights and it's standard error is estimated with the result of \code{fun} using the bootstrap sample weights. \cr
#' \cr
#' \code{fun} can be any function which returns a double or integer and uses sample weights as it's second argument. The predifined options are \code{weightedRatio,weightedSum,sampSize} and \code{popSize}, for wich \code{sampSize} and \code{popSize} indicate the sample and population size respectively.\cr
#' \cr
#' For the option \code{weightedRatio} a weighted ratio (in \%) of \code{var} is calculated for \code{var} equal to 1, e.g \code{sum(weight[var==1])/sum(weight[!is.na(var)])*100}.\cr
#' \cr
#' If \code{cross_var} is not \code{NULL} but a vector of variables from \code{dat} then \code{fun} is applied on each subset of \code{dat} defined by all combinations of values in \code{cross_var}.\cr
#' For instance if \code{cross_var = "sex"} with "sex" having the values "Male" and "Female" in \code{dat} the point estimate and standard error is calculated on the subsets of \code{dat} with only "Male" or "Female" value for "sex". This is done for each value of \code{year}.
#' For variables in \code{cross_var} which have \code{NA}s in \code{dat} the rows containing the missings will be discarded. \cr
#' When \code{cross_var} is a list of character vectors, subsets of \code{dat} and the following estimation of the point estimate, including the estimate for the standard error, are calculated for each list entry.\cr
#' \cr
#' When defining \code{year.diff} the difference of point estimates between years as well their standard errors are calculated.\cr
#' The entries in \code{year.diff} must have the form of \code{"year1 - year2"} which means that the results of the point estimates for \code{year2} will be substracted from the results of the point estimates for \code{year1}.\cr
#' \cr
#' Specifying \code{year.mean} leads to an improvement in standard error by averaging the results for the point estimates, using the bootstrap weights, over \code{year.mean} years.
#' Setting, for instance, \code{year.mean = 3} the results in averaging these results over each consecutive set of 3 years.\cr
#' Estimating the standard error over these averages gives an improved estimate of the standard error for the central year, which was used for averaging.\cr
#' The averaging of the results is also applied in differences of point estimates. For instance defining \code{year.diff = "2015-2009"} and \code{year.mean = 3}
#' the differences in point estimates of 2014 and 2008, 2015 and 2009 as well as 2016 and 2010 are calcualated and finally the average over these 3 differences is calculated.
#' \cr
#' Setting \code{bias} to \code{TRUE} returns the calculation of a mean over the results from the bootstrap replicates. In  the output the corresponding columns is labeled \emph{_mean} at the end.\cr
#' \cr
#' If \code{fun} needs more arguments they can be set in add.arg.\cr
#' \cr
#' The parameter \code{size.limit} indicates a lower bound of the sample size for subsets in \code{dat} created by \code{cross_var}. If the sample size of a subset falls below \code{size.limit} a warning will be displayed.\cr
#' In addition all subsets for which this is the case can be selected from the output of \code{calc.stError} with \code{$smallGroups}.\cr
#' With the parameter \code{cv.limit} one can set an upper bound on the coefficient of variantion. Estimates which exceed this bound are flagged with \code{TRUE} and are available in the function output with \code{$cvHigh}.
#' \code{cv.limit} must be a positive integer and is treated internally as \%, e.g. for \code{cv.limit=1} the estimate will be flagged if the coefficient of variantion exceeds 1\%.\cr
#' \cr
#' When specifying \code{year.mean}, the decrease in standard error for choosing this method is internally calcualted and a rough estimate for an implied increase in sample size is available in the output with \code{$stEDecrease}.
#' The rough estimate for the increase in sample size uses the fact that for a sample of size \eqn{n} the sample estimate for the standard error of most point estimates converges with a factor \eqn{1/\sqrt{n}} against the true standard error \eqn{\sigma}.
#'
#' @return Returns a list containing:
#' \itemize{
#'   \item \code{Estimates}: data.table containing yearly, differences and/or k year averages for estimates of \code{fun} applied to \code{var} as well as the corresponding standard errors, which are calculated using the bootstrap weights.
#'   \item \code{smallGroups}: data.table containing groups for which the number of observation falls below \code{size.limit}.
#'   \item \code{cvHigh}: data.table containing a boolean variable which indicates for each estimate if the estimated standard error exceeds \code{cv.limit}.
#'   \item \code{stEDecrease}: data.table indicating for each estimate the theoretical increase in sample size which is gained when averaging over k years. Only returned if \code{year.mean} is not \code{NULL}.
#' }
#'
#' @seealso \code{\link{draw.bootstrap}} \cr
#' \code{\link{recalib}}
#' @keywords survey manip
#' @author Johannes Gussenbauer, Alexander Kowarik, Statistics Austria
#'
#' @examples
#' # read in and prepare data
#' library(data.table)
#' dat <- data.table(read_sas("O:/B/3-AP/Analyse/sonstiges/
#'                             bundesländerschätzungen 2008-2018/daten/bldaten0816.sas7bdat"))
#'
#' dat <- draw.bootstrap(dat,REP=20,hid="hid",weights="hgew",strata="bundesld",
#'                      year="jahr",totals=NULL,boot.names=NULL)
#' dat <- recalib(dat,hid="hid",weights="hgew",b.rep=paste0("w",1:20),
#'                year="jahr",conP.var=c("ksex","kausl","al","erw","pension"),
#'                conH.var=c("bundesld","hsize","recht"))
#'
#' # or load file with calibrated bootstrap weights
#' # load("dat_calibweight.RData")
#'
#' # estimate weightedRatio for povmd60 per year
#' err.est <- calc.stError(dat,weights="hgew",b.weights=paste0("w",1:20),year="jahr",var="povmd60",
#'                        fun="weightedRatio",cross_var=NULL,year.diff=NULL,year.mean=NULL)
#'
#' # estimate weightedRatio for povmd60 per year and sex
#' cross_var <- "sex"
#' err.est <- calc.stError(dat,weights="hgew",b.weights=paste0("w",1:20),
#'                         year="jahr",var="povmd60",fun="weightedRatio",
#'                         cross_var=cross_var,year.diff=NULL,year.mean=NULL)
#'
#' # use average over 3 years for standard error estimation
#' err.est <- calc.stError(dat,weights="hgew",b.weights=paste0("w",1:20),year="jahr",var="povmd60",
#'                         fun="weightedRatio",cross_var=cross_var,year.diff=NULL,year.mean=3)
#'
#' # get estimate for difference of year 2016 and 2013
#' year.diff <- c("2016-2013")
#' err.est <- calc.stError(dat,weights="hgew",b.weights=paste0("w",1:20),year="jahr",var="povmd60",
#'                        fun="weightedRatio",cross_var=cross_var,year.diff=year.diff,year.mean=3)
#'
#' # apply function to multiple variables and define different subsets
#' var <- c("povmd60","arose")
#' cross_var <- list("sex","bundesld",c("sex","bundesld"))
#' err.est <- calc.stError(dat,weights="hgew",b.weights=paste0("w",1:20),year="jahr",var=var,
#'                        fun="weightedRatio",cross_var=cross_var,year.diff=year.diff,year.mean=3)
#'
#' # use a function from an other package that has sampling weights as its second argument
#' # for example ging() from laeken
#' library(laeken)
#'
#' # set up help function that returns only the gini index
#' help_gini <- function(x,w){
#'  return(gini(x,w)$value)
#' }
#'
#' # exporting data
#' # get point estimates
#' results <- err.est$Estimates
#' write2.csv(results,file="My_Results.csv",row.names=FALSE)
#'
#' err.est <- calc.stError(dat,weights="hgew",b.weights=paste0("w",1:20),year="jahr",var="epinc_real",
#'                        fun="help_gini",cross_var=cross_var,year.diff=year.diff,year.mean=3)
#'
#' @export calc.stError
#' @export print.surveysd
#' @useDynLib surveysd
#' @importFrom Rcpp sourceCpp
#' @import simPop data.table Rcpp


# wrapper-function to apply fun to var using weights (~weights, b.weights)
# and calculating standard devation (using the bootstrap replicates) per year and für 3-year rolling means
calc.stError <- function(dat,weights,b.weights=paste0("w",1:1000),year,var,fun="weightedRatio",cross_var=NULL,year.diff=NULL,year.mean=3,bias=FALSE,add.arg=NULL,size.limit=20,cv.limit=10,p=NULL){

  ##########################################################
  # INPUT CHECKING
  if(class(dat)[1]=="data.frame"){
    dat <- as.data.table(dat)
  }else if(class(dat)[1]!="data.table"){
    stop("dat must be a data.frame or data.table")
  }

  c.names <- colnames(dat)

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

  # check b.weights
  if(!all(b.weights%in%c.names)){
    stop("Not all elements in b.rep are column names in dat")
  }
  if(any(!grepl("^[[:alpha:]]",b.weights))){
    stop("Column names of bootstrap replicates must start with alphabetic character")
  }
  if(any(!unlist(lapply(dat[,mget(b.weights)],is.numeric)))){
    stop("Columns containing bootstrap replicates must be numeric")
  }

  # check year
  if(length(year)!=1){
    stop("year must have length 1")
  }
  if(!year%in%c.names){
    stop(paste0(year," is not a column in dat"))
  }

  # check var
  if(any(!var%in%c.names)){
    stop("Not all elements in var are column names in dat")
  }
  if(any(!unlist(lapply(dat[,mget(b.weights)],is.numeric)))){
    stop(paste0("Columns containing ",paste(var,collapse=",")," must all be numeric"))
  }

  # check fun
  if(length(fun)!=1){
    stop("fun must have length 1")
  }
  if(!fun%in%c("weightedRatio","weightedSum","sampSize","popSize"))
  if(length(find(fun))==0){
    stop(paste0("Function ",fun," is undefined"))
  }

  # check cross_var
  if(is.null(cross_var)){
    cross_var <- list(NULL)
  }
  if(class(cross_var)!="list"){
    cross_var <- as.list(cross_var)
  }
  if(!any(unlist(lapply(cross_var,is.null)))){
    cross_var <- c(list(NULL),cross_var)
  }
  if(any(!unlist(cross_var)%in%c.names)){
    stop("Not all elements on cross_var are column names in dat")
  }

  # check year.mean
  if(length(year.mean)!=1){
    stop("year.mean must have length 1")
  }
  if(!is.numeric(year.mean)){
    stop("year.mean must contain one numeric value")
  }
  if(year.mean%%1!=0){
    stop("year.mean cannot have a decimal part")
  }

  # check bias
  if(length(bias)!=1){
    stop("year.mean must have length 1")
  }
  if(!is.logical(bias)){
    stop("bias can only be TRUE of FALSE")
  }

  # check size.limit
  if(length(size.limit)!=1){
    stop("size.limit must have length 1")
  }
  if(!is.numeric(size.limit)){
    stop("size.limit must contain one numeric value")
  }

  # check cv.limit
  if(length(cv.limit)!=1){
    stop("cv.limit must have length 1")
  }
  if(!is.numeric(cv.limit)){
    stop("cv.limit must contain one numeric value")
  }

  # check p
  if(!is.null(p)){
    if(!is.numeric(p)){
      stop("p must be a numeric vector")
    }
    if(any(!p%between%c(0,1))){
      stop("Values in p must be between 0 and 1")
    }
  }

  # check year.diff
  if(!is.null(year.diff)){
    years.dat <- dt.eval("dat[,unique(",year,")]")
    year.diff <- strsplit(year.diff,"-")

    rm.index <- rep(0,length(year.diff))
    for(i in 1:length(year.diff)){
      if(any(!year.diff[[i]]%in%years.dat)){
        cat("Removing",paste(year.diff[[i]],collapse="-"),"from year.diff - year(s) not present in column",year,"\n")
        rm.index[i] <- 1
      }
    }
    if(all(rm.index==1)){
      cat("No differences will be calculated\n")
      year.diff <- NULL
    }else{
      year.diff <- year.diff[rm.index==0]
    }
  }


  ##########################################################

  ##########################################################
  # setup parameters
	# define columns in which NAs are present (will be discarded for the evaluation)

	col_cross <- unique(unlist(cross_var))
	no.na <- unlist(dat[,lapply(.SD,function(z){all(!is.na(z))}),.SDcols=col_cross])
	no.na <- colnames(dat)[colnames(dat)%in%col_cross][!no.na]

	if(!is.null(p)){
    p.names <- paste0("p",p)
	}else{
	  p.names <- NULL
	}

	# calculate point estimates
  outx <- help.stError(dat=dat,year=year,var=var,weights=weights,b.weights=b.weights,fun=fun,cross_var=cross_var,year.diff=year.diff,year.mean=year.mean,bias=bias,no.na=no.na,add.arg=add.arg,size.limit=size.limit,p=p)

  outx.names <- colnames(outx)
  outx.names <- outx.names[!outx.names%in%c("val","N","est_type","stE","mean","size",p.names)]
  # get meta data like stE_high - size - increase in effektive sample size
  # flag stE if values are especially high
  outx[,stE_high:=((stE/val)*100)>cv.limit]

  # create bool matrix for stE_high
  sd_bool <- subset(outx,select=c("stE_high",outx.names))
  form <- as.formula(paste(paste(outx.names[outx.names!="est"],collapse="+"),"est",sep="~"))
  sd_bool <- dcast(sd_bool,form,value.var="stE_high")

  # create matrix for increase of sample size
  if(!is.null(year.mean)){
    # estimate (roughly) the effektive sample size per
    samp_eff <- outx[est_type=="roll"]
    setnames(samp_eff,year,paste0(year,"_roll"))
    dt.eval("samp_eff[,",year,":=tstrsplit(",year,"_roll,'_',keep=2)]")
    if(bias){
      samp_eff[,"mean":=NULL]
    }
    samp_eff[,c("size","val","est_type","N"):=NULL]
    setnames(samp_eff,"stE","stE_roll")
    same_names <- intersect(colnames(samp_eff),colnames(outx))

    samp_eff <- merge(samp_eff,outx[,mget(c("stE","N",same_names))],by=same_names)
    samp_eff[,N_inc:=((stE/stE_roll)^2-1)*N]
    samp_eff[,c(paste0(year,"_roll"),"stE_roll","stE","N"):=NULL]
  }else{
    samp_eff <- NULL
  }

  # create Matrix for groups which have small sizes
  size_group <- unique(subset(outx[size==TRUE],select=c(outx.names[outx.names!="est"],"N")))

  val.var <- c("val","stE",p.names)
  if(bias){
    val.var <- c(val.var,"mean")
  }

  outx <- dcast(outx,form,value.var=val.var,fill=NA)
  # reorder output
  outx.names <- colnames(outx)
  outx.names.sub <- gsub(paste0(paste0(val.var,"_"),collapse="|"),"",outx.names)
  ord.names <- c()
  for(i in var){
    ord.names <- c(ord.names,which(outx.names.sub==i))
  }

  outx.names <- outx.names[c(1:(min(ord.names)-1),ord.names)]
  outx <- outx[,mget(outx.names)]

  # specify parameters for output
  param <- list(number.bweights=length(b.weights),year=year,var=var,fun=fun,package=find(fun),cross_var=cross_var,year.diff=year.diff,year.mean=year.mean,
                bias=bias,add.arg=add.arg,size.limit=size.limit,cv.limit=cv.limit)

  output <- list(Estimates=outx,smallGroups=size_group,cvHigh=sd_bool,stEDecrease=samp_eff,param=param)

  class(output) <- c("surveysd", class(output))

	return(output)

}


# function to apply fun to var using weights (~weights, b.weights)
# and calculating standard devation (using the bootstrap replicates) per year and für 3-year rolling means
help.stError <- function(dat,year,var,weights,b.weights=paste0("w",1:1000),fun,cross_var,year.diff=NULL,year.mean=NULL,bias=FALSE,no.na,add.arg=NULL,size.limit=20,p=NULL){

  # use c++ implementation for weightedRatio
  if(fun=="weightedRatio"){
    fun <- "weightedRatioC"
  }

	# define names for estimates for each weight (normal weights and boostrap weights)
	# makes it easier to (sort of) verctorize expressions
  res.names <- c(t(outer(var, 1:length(c(weights,b.weights)), paste_)))

  #res.names <- paste0("r",1:length(c(weights,b.weights)))
	# formulate function and arguments
	if(fun=="sampSize"){
		eval.fun <- c(res.names,"=.N)")
	}else if(fun=="popSize"){
		eval.fun <- paste0(res.names,"=sum(",c(weights,b.weights),")")
	}else{
	  if(!is.null(add.arg)){
	    add.arg <- paste(add.arg,collapse=",")
	    eval.fun <- paste0(res.names,"=",fun,"(",paste(c(t(outer(var,c(weights,b.weights), paste_c))),add.arg,sep=","),")")
	  }else{
	    eval.fun <- paste0(res.names,"=",fun,"(",c(t(outer(var,c(weights,b.weights), paste_c))),")")
	  }

	  eval.fun <- paste0(".(.N,",paste(eval.fun,collapse=","),")")
	}

  # define additional parameters for mean over consecutive years
	years <- dt.eval("dat[,unique(",year,")]")
	# formulate k consecutive years
	if(length(years)>=year.mean){
	  yearsList <- unlist(lapply(years[1:c(length(years)-year.mean+1)],function(z){
	    if(!((z+year.mean-1)>max(years))){
	      paste(z:c(z+year.mean-1),collapse="_")
	    }
	  }))
	}else{
	  yearsList <- NULL
	  cat(paste0("Not enough years present in data to calculate mean over ",year.mean," years.\n"))
	}

	# define parameter for quantile calculation
	if(!is.null(p)){
	  p.names <- paste0("p",p)
	  np <- length(p.names)
	}

	# get years for k year mean over yearly differences
  # get for each difference a list of vectors that correspond to the differences needed to use for the mean over differences
	if(!is.null(year.diff)){
    year.diff.b <- TRUE
    year.diff.mean <- lapply(year.diff,function(z){
      z <- as.numeric(z)
      if(z[1]+(year.mean-1)<=max(years)){
        lapply(1:year.mean,function(s){
          z+s
        })
      }else{
        cat(paste0("Can not calcualte differences between years ",z[1]," and ",z[2]," over ",year.mean," years.\n"))
        NULL
      }
    })
  }else{
    year.diff.b <- FALSE
    year.diff.mean <- NULL
  }

	# apply function to all elemnts of cross_var
	# apply also mean and standard deviation for estimates
	out <- lapply(cross_var,function(z){

	  # use only unique values for grouping (duplicates are discarded)
	  # if year in z also discard -> always group by year
	  z <- unique(z[!z==year])
		na.check <- z[!z%in%no.na]
		if(length(na.check)>0){
			na.eval <- paste(paste0("(!is.na(",na.check,"))"),collapse="&")
		}else{
			na.eval <- NULL
		}
		by.eval <- paste(c(year,z),collapse=",")

		# calcualte estimate
		var.est <- dt.eval("dat[",na.eval,",",eval.fun,",by=list(",by.eval,")]")

		if(nrow(var.est[N<size.limit])>0){
		  small.group <- var.est[N<size.limit,mget(c(year,z))]
		  cat(paste0("For grouping by ",paste(c(year,z),collapse="~"),": \n"))
		  if(nrow(var.est[N<size.limit])>10){
		    cat(paste0("Sample size lower ",size.limit," for ",nrow(var.est[N<size.limit]) ," groups \n"))
		  }else{
		    cat(paste0("Sample size lower ",size.limit," for groups \n"))
		    print(small.group)
		  }
		}

		var.est <- melt(var.est,id.vars=c(year,z,"N"),measure.vars=res.names,value.name = "V1")
		var.est[,c("est","ID"):=tstrsplit(variable,"\\.")]
		var.est[,est_type:="norm"]
		var.est[,variable:=NULL]

		# add groups which are not in var.est
		# because they were not present in sample
		if(!is.null(year.mean)|year.diff.b){
		  # if group did not exist in year
		  # add group with NA
		  roll.miss <- unique(dt.eval("dat[",na.eval,",list(",by.eval,")]"))
		  roll.miss <- lapply(roll.miss,function(l){unique(na.omit(l))})
		  roll.miss <- c(roll.miss,list(est=var.est[,unique(est)]))
		  roll.est <- data.table(expand.grid(c(roll.miss,list(ID=unique(var.est$ID)))))

		  if(nrow(roll.est)>nrow(unique(var.est,by=c(year,z)))){
		    var.est <- merge(roll.est,var.est,by=c(year,z,"ID","est"),all.x=TRUE)
		    setkeyv(var.est,year)
		    var.est[is.na(V1),N:=0]
		  }
		}

		# calculate mean of estimate over 'year.mean'- years
		if(!is.null(year.mean)){
		  roll.est <- var.est[,list(V1=rollMeanC(x=V1,k=year.mean,type="c"),V2=yearsList),by=c("ID",z,"est")]
		  setnames(roll.est,"V2",year)
		  roll.est[,est_type:="roll"]
		}

		if(year.diff.b){
		  # calcualte differences between years
		  by.diff <- paste(c("ID","est",z),collapse=",")
		  diff.est <- lapply(year.diff,function(y){
		    y_cond <- paste(year,paste0("c(",paste(y,collapse=","),")"),sep="%in%")
		    diff.y <- dt.eval("var.est[",y_cond,",V1[",year,"==",y[1],"]-V1[",year,"==",y[2],"],by=list(",by.diff,")]")
		    diff.y[,c(year):=paste(y,collapse="-")]
		    return(diff.y)
		  })
		  diff.est <- rbindlist(diff.est)
		  diff.est[,est_type:="diff"]

		  # calcualte differences between years and mean over consecutive differences
		  if(!is.null(unlist(year.diff.mean))){
		    i.diff.mean <- (year.mean+1)/2
		    diff.mean.est <- lapply(year.diff.mean,function(d){
		      # calculate differences for all pairwise years in d
		      d.m.est <- lapply(d,function(y){
		        y_cond <- paste(year,paste0("c(",paste(y,collapse=","),")"),sep="%in%")
		        diff.y <- dt.eval("var.est[",y_cond,",V1[",year,"==",y[1],"]-V1[",year,"==",y[2],"],by=list(",by.diff,")]")
		        diff.y[,c(year):=paste(y,collapse="-")]
		        return(diff.y)
		      })
		      # calcualte mean over consecutive differences
		      d.m.est <- rbindlist(d.m.est)
		      d.m.est <- dt.eval("d.m.est[,mean(V1),by=list(",by.diff,")]")
		      d.m.est[,c(year):=paste0(paste(d[[i.diff.mean]],collapse="-"),"-mean")]
		    })
		    diff.mean.est <- rbindlist(diff.mean.est)
		  }
		}

		# add results to one data.table
		if(!is.null(year.mean)){
		  var.est <- rbind(var.est,roll.est,fill=TRUE)
		}
		if(year.diff.b){
		  var.est <- rbind(var.est,diff.est,fill=TRUE)

		  if(!is.null(unlist(year.diff.mean))){
		    var.est <- rbind(var.est,diff.mean.est,fill=TRUE)
		  }
		}

		if(!is.null(p)){
		  sd.est <- var.est[ID!=1,as.list(c(stE=sd(V1),quantileNA(V1,probs=p,p.names=p.names,np=np))),by=c(year,z,'est')]
		}else{
		  sd.est <- var.est[ID!=1,.(stE=sd(V1)),by=c(year,z,'est')]
		}


		out.z <- merge(subset(var.est,ID==1,select=c(year,z,"V1","N","est","est_type")),sd.est,by=c(year,z,"est"))

		if(bias){
		  bias.est <- var.est[ID!=1,.(mean=mean(V1)),by=c(year,z)]
		  out.z <- merge(out.z,bias.est,by=c(year,z))
		}

		# define size groups - groups with zero size do not fall under size-output
		# for groups with zero size the resutling estimatese will be NAs
		out.z[(!is.na(N))&N>0,size:=N<size.limit]

		return(out.z)

	})

	out <- rbindlist(out,fill=TRUE,use.names=TRUE)
	setnames(out,"V1","val")
	return(out)
}


# additional help function
weightedRatio <- function(x,weightvar){
  sum(weightvar[x==1],na.rm=TRUE)/sum(weightvar[!is.na(x)],na.rm=TRUE)*100
}


# weightedRatio <- function(x,weightvar,x2=NULL){
# 	if(is.null(x2)){
# 		return(sum(weightvar[x==1],na.rm=TRUE)/sum(weightvar[!is.na(x)],na.rm=TRUE)*100)
# 	}else{
# 		if(sum(weightvar[x==1],na.rm=TRUE)==0&sum(weightvar[x2==1],na.rm=TRUE)==0){
# 			return(0)
# 		}
# 		return(sum(weightvar[x==1],na.rm=TRUE)/sum(weightvar[x2==1],na.rm=TRUE))
# 	}
# }

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

# helpfunction for quantile calcultion with missings
quantileNA <- function(x,probs,p.names,np=length(probs)){

  if(any(is.na(x))){
    out <- rep(NA_real_,np)
  }else{
    out <- quantile(x,probs=probs)
  }
  names(out) <- p.names
  return(out)
}


# print function for surveysd objects
print.surveysd <- function(sd.result){

  # get parameter
  col.val <- grepl("^val_*.",colnames(sd.result[["Estimates"]]))
  col.val <- colnames(sd.result[["Estimates"]])[col.val]
  # get number of estimates ~ variables in number of groups using function fun from package pack
  n.estimates <- nrow(sd.result[["Estimates"]])*length(col.val)
  # get number of subgroup which fall below size
  n.years <- unique(sd.result[["Estimates"]][[sd.result[["param"]][["year"]]]])
  n.years <- length(n.years[!grepl("-|_",n.years)])
  n.groups <- nrow(unique(sd.result[["Estimates"]][,.N,by=c(unique(unlist(sd.result[["param"]][["cross_var"]])))]))
  # get number of missing values in the output due to subgroups with no observations or subgroups with all observations equal NA
  n.NAs <- sum(is.na(sd.result[["Estimates"]][,mget(col.val)]))

  # print number of estimates calculated as well as functino and variables used
  if(!is.na(sd.result[["param"]][["package"]][1])){
    cat("Calculated point estimates for variable(s)\n\n",paste(sd.result[["param"]][["var"]],sep=","),"\n\nusing function",sd.result[["param"]][["fun"]],"from",sd.result[["param"]][["package"]][1],"\n\n")
  }else{
    cat("Calculated point estimates for variable(s)\n\n",paste(sd.result[["param"]][["var"]],sep=","),"\n\nusing function",sd.result[["param"]][["fun"]],"\n\n")
  }

  cat("Results hold",n.estimates,"point estimates for",n.years,"years in",n.groups,"subgroups\n")
  cat("\n")
  if(nrow(sd.result[["smallGroups"]])>10|nrow(sd.result[["smallGroups"]])==0){
    cat(nrow(sd.result[["smallGroups"]]),"subgroups contained less than",sd.result[["param"]][["size.limit"]],"observations\n")
  }else{
    cat("Subgroups with less than",sd.result[["param"]][["size.limit"]],"observations\n")
    print(sd.result[["smallGroups"]])
  }
  cat("\n")
  # print number of missing values for point estimates
  cat("Point estimates are NAs in",n.NAs,"cases due to either\n no observations or only NAs for the variable(s) in the corresponding subgroups.\n")

  cat("\n")
  # get number of point estimates where sd exceeds cv.limit
  stEtoohigh <- colnames(sd.result[["cvHigh"]])
  stEtoohigh <- stEtoohigh[!stEtoohigh%in%c(sd.result[["param"]][["year"]],unique(unlist(sd.result[["param"]][["cross_var"]])))]
  stEtoohigh <- as.matrix(subset(sd.result[["cvHigh"]],select=stEtoohigh))
  cat("Estimted standard error exceeds",sd.result[["param"]][["cv.limit"]],"% of the the point estimate in",sum(stEtoohigh,na.rm=TRUE),"cases\n")
  cat("\n")

}
