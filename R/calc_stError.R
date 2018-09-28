#' @title Calcualte point estimates and their standard errors using bootstrap weights.
#'
#' @description
#' Calculate point estimates as well as standard errors of variables in surveys. Standard errors are estimated using bootstrap weights (see \code{\link{draw.bootstrap}} and \code{\link{recalib}}).
#' In addition the standard error of an estimate can be calcualted using the survey data for 3 or more consecutive periods, which results in a reduction of the standard error.
#'
#'
#' @param dat either data.frame or data.table containing the survey data. Surveys can be a panel survey or rotating panel survey, but does not need to be. For rotating panel survey bootstrap weights can be created using \code{\link{draw.bootstrap}} and \code{\link{recalib}}.
#' @param weights character specifying the name of the column in \code{dat} containing the original sample weights. Used to calculate point estimates.
#' @param b.weights character vector specifying the names of the columns in \code{dat} containing bootstrap weights. Used to calculate standard errors.
#' @param period character specifying the name of the column in \code{dat} containing the sample periods.
#' @param var character vector containing variable names in \code{dat} on which \code{fun} shall be applied for each sample period.
#' @param fun function which will be applied on \code{var} for each sample period.
#' Predefined functions are \code{weightedRatio,weightedSum}, but can also take any other function which returns a double or integer and uses weights as its second argument.
#' @param national boolean, if TRUE point estimates resulting from fun will be divided by the point estimate at the national level.
#' @param group character vectors or list of character vectors containig variables in \code{dat}. For each list entry \code{dat} will be split in subgroups according to the containing variables as well as \code{period}.
#' The pointestimates are then estimated for each subgroup seperately. If \code{group=NULL} the data will split into sample periods by default.
#' @param fun.adjust.var can be either \code{NULL} or a function. This argument can be used to apply a function for each \code{period} and bootstrap weight to the data. The resulting estimates will be passed down to \code{fun}. See details for more explanations.
#' @param adjust.var can be either \code{NULL} or a character specifying the first argument in \code{fun.adjust.var}.
#' @param period.diff character vectors, defining periods for which the differences in the point estimate as well it's standard error is calculated. Each entry must have the form of \code{"period1 - period2"}. Can be NULL
#' @param period.mean odd integer, defining the range of periods over which the sample mean of point estimates is additionally calcualted.
#' @param bias boolean, if TRUE the sample mean over the point estimates of the bootstrap weights is returned.
#' @param size.limit integer defining a lower bound on the number of observations on \code{dat} in each group defined by \code{period} and the entries in \code{group}.
#' Warnings are returned if the number of observations in a subgroup falls below \code{size.limit}. In addition the concerned groups are available in the function output.
#' @param cv.limit non-negativ value defining a upper bound for the standard error in relation to the point estimate. If this relation exceed \code{cv.limit}, for a point estimate, they are flagged and available in the function output.
#' @param p numeric vector containing values between 0 and 1. Defines which quantiles for the distribution of \code{var} are additionally estimated.
#' @param ... additional arguments which will be passed to fun.
#'
#' @details \code{calc.stError} takes survey data (\code{dat}) and returns point estimates as well as their standard Errors defined by \code{fun} and \code{var} for each sample period in \code{dat}.
#' \code{dat} must be household data where household members correspond to multiple rows with the same household identifier. The data should at least containt the following columns:
#' \itemize{
#'   \item Column indicating the sample period;
#'   \item Column indicating the household ID;
#'   \item Column containing the household sample weights;
#'   \item Columns which contain the bootstrap weights (see output of \code{\link{recalib}});
#'   \item Columns listed in \code{var} as well as in \code{group}
#' }
#' For each variable in \code{var} as well as sample period the function \code{fun} is applied using the original as well as the bootstrap sample weights.\cr
#' The point estimate is then selected as the result of \code{fun} when using the original sample weights and it's standard error is estimated with the result of \code{fun} using the bootstrap sample weights. \cr
#' \cr
#' \code{fun} can be any function which returns a double or integer and uses sample weights as it's second argument. The predifined options are \code{weightedRatio} and \code{weightedSum}.\cr
#' \cr
#' For the option \code{weightedRatio} a weighted ratio (in \%) of \code{var} is calculated for \code{var} equal to 1, e.g \code{sum(weight[var==1])/sum(weight[!is.na(var)])*100}.\cr
#' Additionally using the option \code{national=TRUE} the weighted ratio (in \%) is divided by the weighted ratio at the national level for each \code{period}.
#' \cr
#' If \code{group} is not \code{NULL} but a vector of variables from \code{dat} then \code{fun} is applied on each subset of \code{dat} defined by all combinations of values in \code{group}.\cr
#' For instance if \code{group = "sex"} with "sex" having the values "Male" and "Female" in \code{dat} the point estimate and standard error is calculated on the subsets of \code{dat} with only "Male" or "Female" value for "sex". This is done for each value of \code{period}.
#' For variables in \code{group} which have \code{NA}s in \code{dat} the rows containing the missings will be discarded. \cr
#' When \code{group} is a list of character vectors, subsets of \code{dat} and the following estimation of the point estimate, including the estimate for the standard error, are calculated for each list entry.\cr
#' \cr
#' The optional parameters \code{fun.adjust.var} and \code{adjust.var} can be used if the values in \code{var} are dependent on the \code{weights}. As is for instance the case for the poverty thershhold calculated from EU-SILC.
#' In such a case an additional function can be supplied using \code{fun.adjust.var} as well as its first argument \code{adjust.var}, which needs to be part of the data set \code{dat}. Then, before applying \code{fun} on variable \code{var}
#' for all \code{period} and groups, the function \code{fun.adjust.var} is applied to \code{adjust.var} using each of the bootstrap weights seperately (NOTE: weight is used as the second argument of \code{fun.adjust.var}).
#' Thus creating i=1,...,\code{length(b.weights)} additional variables.
#' For applying \code{fun} on \code{var} the estimates for the bootstrap replicate will now use each of the corresponding new additional variables. So instead of
#' \deqn{fun(var,weights,...),fun(var,b.weights[1],...),fun(var,b.weights[2],...),...}
#' the function \code{fun} will be applied in the way
#' \deqn{fun(var,weights,...),fun(var.1,b.weights[1],...),fun(var.2,b.weights[2],...),...}
#' where \code{var.1}, \code{var.2},... correspond to the estimates resulting from \code{fun.adjust.var} and \code{adjust.var}.
#' NOTE: This procedure is especially usefull if the \code{var} is dependent on \code{weights} and \code{fun} is applied on subgroups of the data set. Then it is not possible to capture this
#' procedure with \code{fun} and \code{var}, see examples for a more hands on explanation.
#' \cr
#' When defining \code{period.diff} the difference of point estimates between periods as well their standard errors are calculated.\cr
#' The entries in \code{period.diff} must have the form of \code{"period1 - period2"} which means that the results of the point estimates for \code{period2} will be substracted from the results of the point estimates for \code{period1}.\cr
#' \cr
#' Specifying \code{period.mean} leads to an improvement in standard error by averaging the results for the point estimates, using the bootstrap weights, over \code{period.mean} periods.
#' Setting, for instance, \code{period.mean = 3} the results in averaging these results over each consecutive set of 3 periods.\cr
#' Estimating the standard error over these averages gives an improved estimate of the standard error for the central period, which was used for averaging.\cr
#' The averaging of the results is also applied in differences of point estimates. For instance defining \code{period.diff = "2015-2009"} and \code{period.mean = 3}
#' the differences in point estimates of 2015 and 2009, 2016 and 2010 as well as 2017 and 2011 are calcualated and finally the average over these 3 differences is calculated.
#' The periods set in \code{period.diff} are always used as starting periods from which \code{period.mean-1} consecutive periods are used to build the average.
#' \cr
#' Setting \code{bias} to \code{TRUE} returns the calculation of a mean over the results from the bootstrap replicates. In  the output the corresponding columns is labeled \emph{_mean} at the end.\cr
#' \cr
#' If \code{fun} needs more arguments they can be supplied in \code{...}\cr
#' \cr
#' The parameter \code{size.limit} indicates a lower bound of the sample size for subsets in \code{dat} created by \code{group}. If the sample size of a subset falls below \code{size.limit} a warning will be displayed.\cr
#' In addition all subsets for which this is the case can be selected from the output of \code{calc.stError} with \code{$smallGroups}.\cr
#' With the parameter \code{cv.limit} one can set an upper bound on the coefficient of variantion. Estimates which exceed this bound are flagged with \code{TRUE} and are available in the function output with \code{$cvHigh}.
#' \code{cv.limit} must be a positive integer and is treated internally as \%, e.g. for \code{cv.limit=1} the estimate will be flagged if the coefficient of variantion exceeds 1\%.\cr
#' \cr
#' When specifying \code{period.mean}, the decrease in standard error for choosing this method is internally calcualted and a rough estimate for an implied increase in sample size is available in the output with \code{$stEDecrease}.
#' The rough estimate for the increase in sample size uses the fact that for a sample of size \eqn{n} the sample estimate for the standard error of most point estimates converges with a factor \eqn{1/\sqrt{n}} against the true standard error \eqn{\sigma}.
#'
#' @return Returns a list containing:
#' \itemize{
#'   \item \code{Estimates}: data.table containing period differences and/or k period averages for estimates of \code{fun} applied to \code{var} as well as the corresponding standard errors, which are calculated using the bootstrap weights.
#'   In addition the sample size, \code{n}, and poplutaion size for each group is added to the output.
#'   \item \code{smallGroups}: data.table containing groups for which the number of observation falls below \code{size.limit}.
#'   \item \code{cvHigh}: data.table containing a boolean variable which indicates for each estimate if the estimated standard error exceeds \code{cv.limit}.
#'   \item \code{stEDecrease}: data.table indicating for each estimate the theoretical increase in sample size which is gained when averaging over k periods. Only returned if \code{period.mean} is not \code{NULL}.
#' }
#'
#' @seealso \code{\link{draw.bootstrap}} \cr
#' \code{\link{recalib}}
#' @keywords survey manip
#' @author Johannes Gussenbauer, Alexander Kowarik, Statistics Austria
#'
#' @examples
#' \dontrun{
#' library(surveysd)
#' library(laeken)
#' library(data.table)
#'
#' eusilc <- surveysd:::demo.eusilc()
#'
#' dat_boot <- draw.bootstrap(eusilc,REP=10,hid="db030",weights="rb050",strata=c("db040"),
#'                            period="year")
#'
#' # calibrate weight for bootstrap replicates
#' dat_boot_calib <- recalib(copy(dat_boot),hid="db030",weights="rb050",
#'                           period="year",b.rep=paste0("w",1:10),conP.var=c("rb090"),
#'                           conH.var = c("db040"))
#'
#' # estimate weightedRatio for povmd60 per period
#' err.est <- calc.stError(dat_boot_calib,weights="rb050",b.weights=paste0("w",1:10),
#'                         period="year",var="povmd60",fun=weightedRatio,
#'                         group=NULL,period.diff=NULL,period.mean=NULL)
#'
#' # estimate weightedRatio for povmd60 per period and rb090
#' group <- "rb090"
#' err.est <- calc.stError(dat_boot_calib,weights="rb050",b.weights=paste0("w",1:10),
#'                         period="year",var="povmd60",fun=weightedRatio,
#'                         group=group,period.diff=NULL,period.mean=NULL)
#'
#'
#' # estimate weightedRatio for povmd60 per period and rb090, db040 and combination of both
#' group <- list("rb090","db040",c("rb090","db040"))
#' err.est <- calc.stError(dat_boot_calib,weights="rb050",b.weights=paste0("w",1:10),
#'                         period="year",var="povmd60",fun=weightedRatio,
#'                         group=group,period.diff=NULL,period.mean=NULL)
#'
#' err.est <- calc.stError(dat_boot_calib,weights="rb050",b.weights=paste0("w",1:10),
#'                         period="year",var="povmd60",fun=weightedRatio,
#'                         group=group,period.diff=NULL,period.mean=NULL)
#'
#' # use average over 3 periods for standard error estimation
#' err.est <- calc.stError(dat_boot_calib,weights="rb050",b.weights=paste0("w",1:10),
#'                         period="year",var="povmd60",fun=weightedRatio,
#'                         group=group,period.diff=NULL,period.mean=3)
#'
#' # get estimate for difference of period 2016 and 2013
#' period.diff <- c("2015-2011")
#' err.est <- calc.stError(dat_boot_calib,weights="rb050",b.weights=paste0("w",1:10),
#'                         period="year",var="povmd60",fun=weightedRatio,
#'                         group=group,period.diff=period.diff,period.mean=3)
#'
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
#' err.est <- calc.stError(dat_boot_calib,weights="rb050",b.weights=paste0("w",1:10),
#'                         period="year",var="povmd60",fun=help_gini,group=group,
#'                         period.diff=period.diff,period.mean=3)
#'
#'
#' # using fun.adjust.var and adjust.var to estimate povmd60 indicator
#' # for each period and bootstrap weight before applying the weightedRatio point estimate
#'
#' # this function estimates the povmd60 indicator with x as income vector
#' # and w as weight vector
#' povmd <- function(x,w){
#'  md <- laeken::weightedMedian(x,w)*0.6
#'  pmd60 <- x<md
#'  return(as.integer(pmd60))
#' }
#'
#' # set adjust.var="eqIncome" so the income vector ist used to estimate
#' # the povmd60 indicator for each bootstrap weight
#' # and the resultung indicators are passed to function weightedRatio
#'
#' err.est <- calc.stError(dat_boot_calib,weights="rb050",b.weights=paste0("w",1:10),
#'                         period="year",var="povmd60",fun=weightedRatio,
#'                         group=group,fun.adjust.var=povmd,adjust.var="eqIncome")
#'
#' # why fun.adjust.var and adjust.var are needed (!!!):
#' # one could also use the following function
#' # and set fun.adjust.var=NULL,adjust.var=NULL
#' # and set fun = povmd, var = "eqIncome"
#'
#' povmd2 <- function(x,w){
#'  md <- laeken::weightedMedian(x,w)*0.6
#'  pmd60 <- x<md
#'  # weighted ratio is directly estimated inside my function
#'  return(sum(w[pmd60])/sum(w)*100)
#' }
#'
#' # but this results in different results in subgroups
#' # compared to using fun.adjust.var and adjust.var
#'
#' err.est.different <- calc.stError(dat_boot_calib,weights="rb050",b.weights=paste0("w",1:10),
#'                         period="year",var="eqIncome",fun=povmd2,
#'                         group=group,fun.adjust.var=NULL,adjust.var=NULL)
#'
#' # results are equal for yearly estimates
#' all.equal(err.est.different$Estimates[is.na(rb090)&is.na(db040)],
#' err.est$Estimates[is.na(rb090)&is.na(db040)],
#' check.attributes = FALSE)
#'
#' # but for subgroups (rb090, db040) results vary
#' all.equal(err.est.different$Estimates[!(is.na(rb090)&is.na(db040))],
#' err.est$Estimates[!(is.na(rb090)&is.na(db040))],
#' check.attributes = FALSE)
#' }
#'
#' @export calc.stError
#'


# wrapper-function to apply fun to var using weights (~weights, b.weights)
# and calculating standard devation (using the bootstrap replicates) per period and for k-period rolling means
calc.stError <- function(dat,weights,b.weights=paste0("w",1:1000),period,var,
                         fun=weightedRatio,national=FALSE,group=NULL,fun.adjust.var=NULL,adjust.var=NULL,
                         period.diff=NULL,period.mean=3,bias=FALSE,size.limit=20,cv.limit=10,p=NULL,...){

  stE_high <- stE <- val <- as.formula <- est_type <- n_inc <- stE_roll <- n <- size <- NULL

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
    stop(weights," is not a column in dat")
  }
  if(!is.numeric(dt.eval("dat[,",weights,"]"))){
    stop(weights," must be a numeric column")
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

  # check period
  if(length(period)!=1){
    stop("period must have length 1")
  }
  if(!period%in%c.names){
    stop(paste0(period," is not a column in dat"))
  }

  # check var
  if(any(!var%in%c(c.names))){
    stop("Not all elements in var are column names in dat")
  }
  if(any(!unlist(lapply(dat[,mget(b.weights)],is.numeric)))){
    stop("Columns containing ",paste(var,collapse=",")," must all be numeric")
  }

  # check national
  if(class(national)!="logical"){
    stop("national can only be logical")
  }

  # check fun and ...
  if(class(fun)!="function"){
    stop("fun can only be a function")
  }
  if(!"..."%in%formalArgs(fun)){
    ellipsNames <- names(list(...))
    if(any(!ellipsNames%in%formalArgs(fun))){
      stop(ellipsNames[!ellipsNames%in%formalArgs(fun)]," not argument(s) of supplied function.")
    }
  }

  test.val <- dt.eval("dat[,fun(",var,",",weights,",...)]")
  if(!is.numeric(test.val)&!is.integer(test.val)){
    stop("Function in fun does not return integer or numeric value")
  }
  if(length(test.val)>1){
    stop("Function in fun does return more than one value. Only functions which return a single value are allowed.")
  }

  # check fun.adjust.var
  if(!is.null(fun.adjust.var)){
    if(class(fun.adjust.var)!="function"){
      stop("fun.adjust.var can only be a function or NULL")
    }

    test.val <- dt.eval("dat[,fun.adjust.var(",var,",",weights,",...)]")
    if(!is.numeric(test.val)&!is.integer(test.val)){
      stop("Function in fun.adjust.var does not return integer or numeric value")
    }
  }
  # check adjust.var
  if(!is.null(adjust.var)){
    if(class(adjust.var)!="character"){
      stop("adjust.var needs to be a character")
    }
    if(length(adjust.var)>1){
      stop("adjust.var can only be a single variable name")
    }
    if(!adjust.var%in%c.names){
      stop("adjust.var must be a column name in dat")
    }
  }

  # check group
  if(is.null(group)){
    group <- list(NULL)
  }
  if(class(group)!="list"){
    group <- as.list(group)
  }
  if(!any(unlist(lapply(group,is.null)))){
    group <- c(list(NULL),group)
  }
  if(any(!unlist(group)%in%c.names)){
    stop("Not all elements on group are column names in dat")
  }

  # check period.mean
  if(!is.null(period.mean)){
    if(length(period.mean)!=1){
      stop("period.mean must have length 1")
    }
    if(!is.numeric(period.mean)){
      stop("period.mean must contain one numeric value")
    }
    if(period.mean%%1!=0){
      stop("period.mean cannot have a decimal part")
    }
    if(period.mean%%2==0){
      warning("period.mean must be odd - mean over periods will not be calculated")
      period.mean <- NULL
    }else{
      if(period.mean==1){
        period.mean <- NULL
      }
    }
  }


  # check bias
  if(length(bias)!=1){
    stop("period.mean must have length 1")
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

  # check period.diff
  if(!is.null(period.diff)){
    periods.dat <- dt.eval("dat[,unique(",period,")]")
    period.diff <- strsplit(period.diff,"-")

    rm.index <- rep(0,length(period.diff))
    for(i in 1:length(period.diff)){
      if(any(!period.diff[[i]]%in%periods.dat)){
        warning("Removing ",paste(period.diff[[i]],collapse="-")," from period.diff - period(s) not present in column ",period,"\n")
        rm.index[i] <- 1
      }
    }
    if(all(rm.index==1)){
      warning("No differences will be calculated\n")
      period.diff <- NULL
    }else{
      period.diff <- period.diff[rm.index==0]
    }
  }


  ##########################################################

  ##########################################################
  # setup parameters
	# define columns in which NAs are present (will be discarded for the evaluation)


  col_cross <- unique(unlist(group))
  if(!is.null(col_cross)){
    no.na <- unlist(dat[,lapply(.SD,function(z){all(!is.na(z))}),.SDcols=col_cross])
    no.na <- names(no.na)[!no.na]

    if(length(no.na)>0){
      warning("Missing values found in column name(s)",no.na,"\n Cells with missing values are discarded for the calculation!\n")
    }

  }else{
    no.na <- NULL
  }

	if(!is.null(p)){
    p.names <- paste0("p",p)
	}else{
	  p.names <- NULL
	}

	# calculate point estimates
  outx <- help.stError(dat=dat,period=period,var=var,weights=weights,b.weights=b.weights,fun=fun,national=national,group=group,
                       fun.adjust.var=fun.adjust.var,adjust.var=adjust.var,period.diff=period.diff,period.mean=period.mean,bias=bias,no.na=no.na,size.limit=size.limit,p=p,...)

  outx.names <- colnames(outx)
  outx.names <- outx.names[!outx.names%in%c("val","est_type","stE","mean","size",p.names)]
  # get meta data like stE_high - size - increase in effektive sample size
  # flag stE if values are especially high
  outx[,stE_high:=((stE/val)*100)>cv.limit]

  # create bool matrix for stE_high
  sd_bool <- subset(outx,select=c("stE_high",outx.names[!outx.names%in%c("N","n")]))
  form <- as.formula(paste(paste(outx.names[!outx.names%in%c("N","n","est")],collapse="+"),"est",sep="~"))
  sd_bool <- dcast(sd_bool,form,value.var="stE_high")

  # create matrix for increase of sample size
  if(nrow(outx[est_type=="roll"])>0){
    # estimate (roughly) the effektive sample size per
    samp_eff <- outx[est_type=="roll"]
    setnames(samp_eff,period,paste0(period,"_roll"))
    dt.eval("samp_eff[,",period,":=tstrsplit(",period,"_roll,'_',keep=2)]")
    if(bias){
      samp_eff[,"mean":=NULL]
    }
    samp_eff[,c("size","val","est_type","N","n"):=NULL]
    setnames(samp_eff,"stE","stE_roll")
    same_names <- intersect(colnames(samp_eff),colnames(outx))

    samp_eff <- merge(samp_eff,outx[,mget(c("stE","n",same_names))],by=same_names)
    samp_eff[,n_inc:=((stE/stE_roll)^2-1)*n]
    samp_eff[,c(paste0(period,"_roll"),"stE_roll","stE","n"):=NULL]
  }else{
    samp_eff <- NULL
  }

  # create Matrix for groups which have small sizes
  size_group <- unique(subset(outx[size==TRUE],select=c(outx.names[!outx.names%in%c("est","N")])))

  val.var <- c("val","stE",p.names)
  if(bias){
    val.var <- c(val.var,"mean")
  }

  form <- as.formula(paste(paste(outx.names[outx.names!="est"],collapse="+"),"est",sep="~"))
  outx <- dcast(outx,form,value.var=val.var,fill=NA)
  # reorder output
  col.order <- as.vector(outer(paste0(val.var,"_"),var,FUN="paste0"))
  outx.names <- colnames(outx)
  col.order <- c(outx.names[!outx.names%in%col.order],col.order)
  setcolorder(outx,col.order)

  param <- list(number.bweights=length(b.weights),period=period,var=var,fun=fun,fun.adjust.var=fun.adjust.var,adjust.var=adjust.var,group=group,period.diff=period.diff,period.mean=period.mean,
                bias=bias,size.limit=size.limit,cv.limit=cv.limit,...)

  output <- list(Estimates=outx,smallGroups=size_group,cvHigh=sd_bool,stEDecrease=samp_eff,param=param)

  class(output) <- c("surveysd", class(output))

	return(output)

}


# function to apply fun to var using weights (~weights, b.weights)
# and calculating standard devation (using the bootstrap replicates) per period and for k-period rolling means
help.stError <- function(dat,period,var,weights,b.weights=paste0("w",1:1000),fun,national,group,
                         fun.adjust.var,adjust.var,period.diff=NULL,period.mean=NULL,bias=FALSE,no.na,size.limit=20,p=NULL,...){

  N <- variable <- est_type <- est <- V1 <- n <- ID <- sd <- . <- size <-

  # create point estimate for subnational result in relation to national level
  if(national){
    add.arg <- c("Nat[1]")
    dt.eval("dat[,Nat:=fun(",var,",",weights,"),by=list(",period,")]")

    # create new functions which divides by national level
    fun_original <- fun
    fun <- function(x,w,additionalArgument,...){
      fun_original(x,w,...)/additionalArgument*100
    }
  }

  # define names for estimates for each weight (normal weights and boostrap weights)
  # makes it easier to (sort of) verctorize expressions
  if(!is.null(fun.adjust.var)){
    varnew <- paste0("'",var,".",2:(length(b.weights)+1),"'",collapse=",")
    varnew <- paste0("c(",varnew,")")
    eval.fun.adjust <- paste0("fun.adjust.var(",adjust.var,",",c(b.weights),")",collapse=",")

    dt.eval("dat[,",varnew,":=.(",eval.fun.adjust,"),by=list(",period,")]")

    res.names <- c(t(outer(var, 1:length(c(weights,b.weights)), paste_)))

    varnew <- c(var,paste0(var,".",2:(length(b.weights)+1)))

    if(national){
      eval.fun <- paste0(res.names,"=fun(",varnew,",",c(weights,b.weights),",",add.arg,",...)")
    }else{
      eval.fun <- paste0(res.names,"=fun(",varnew,",",c(weights,b.weights),",...)")
    }
  }else{

    res.names <- c(t(outer(var, 1:length(c(weights,b.weights)), paste_)))
    if(national){
      eval.fun <- paste0(res.names,"=fun(",paste(c(t(outer(var,c(weights,b.weights), paste_c))),add.arg,sep=","),",...)")
    }else{
      eval.fun <- paste0(res.names,"=fun(",c(t(outer(var,c(weights,b.weights), paste_c))),",...)")
    }

  }



  eval.fun <- paste0(".(n=.N,N=sum(",weights,"),",paste(eval.fun,collapse=","),")")

  # define parameter for quantile calculation
  if(!is.null(p)){
    p.names <- paste0("p",p)
    np <- length(p.names)
  }

  # define additional parameters for mean over consecutive periods
  periods <- sort(dt.eval("dat[,unique(",period,")]"))

  if(!is.null(period.mean)){
    # formulate k consecutive periods
    if(length(periods)>=period.mean&period.mean%%2==1){
      periodsList <- unlist(lapply(periods[1:c(length(periods)-period.mean+1)],function(z){
        if(!((z+period.mean-1)>max(periods))){
          paste(z:c(z+period.mean-1),collapse="_")
        }
      }))
    }else{
      periodsList <- NULL
      warning("Not enough periods present in data to calculate mean over ",period.mean," periods.\n")
      period.mean <- NULL
    }

    # get periods for k period mean over period-differences
    # get for each difference a list of vectors that correspond to the differences needed to use for the mean over differences
    if(!is.null(period.diff)){
      period.diff.b <- TRUE
      period.diff.mean <- lapply(period.diff,function(z){
        z <- as.numeric(z)
        steps <- (period.mean-1)/2
        z_upper <- (z[1]+steps):(z[1]-steps)
        z_lower <- (z[2]+steps):(z[2]-steps)
        diff.feasable <- all(c(z_lower,z_upper)%in%periods)
        if(diff.feasable){
          lapply(1:period.mean,function(s){
            c(z_upper[s],z_lower[s])
          })
        }else{
          warning("Cannot calculate differences between periods ",z[1]," and ",z[2]," over ",period.mean," periods.\n")
          NULL
        }
      })
      period.diff.mean[unlist(lapply(period.diff.mean,is.null))] <- NULL
    }else{
      period.diff.b <- FALSE
      period.diff.mean <- NULL
    }
  }else{
    period.diff.b <- FALSE
    period.diff.mean <- NULL
    periodsList <- NULL
  }

	# apply function to all elemnts of group
	# apply also mean and standard deviation for estimates
	out <- lapply(group,function(z){

	  # use only unique values for grouping (duplicates are discarded)
	  # if period in z also discard -> always group by period
	  z <- unique(z[!z==period])
		na.check <- z[z%in%no.na]
		if(length(na.check)>0){
			na.eval <- paste(paste0("(!is.na(",na.check,"))"),collapse="&")
		}else{
			na.eval <- NULL
		}
		by.eval <- paste(c(period,z),collapse=",")

		# calcualte estimate
		var.est <- dt.eval("dat[",na.eval,",",eval.fun,",by=list(",by.eval,")]")

		if(nrow(var.est[N<size.limit])>0){
		  small.group <- var.est[N<size.limit,mget(c(period,z))]
		  cat(paste0("For grouping by ",paste(c(period,z),collapse="~"),": \n"))
		  if(nrow(var.est[N<size.limit])>10){
		    cat(paste0("Sample size lower ",size.limit," for ",nrow(var.est[N<size.limit]) ," groups \n"))
		  }else{
		    cat(paste0("Sample size lower ",size.limit," for groups \n"))
		    print(small.group)
		  }
		}

		var.est <- melt(var.est,id.vars=c(period,z,"n","N"),measure.vars=res.names,value.name = "V1")
		var.est[,c("est","ID"):=tstrsplit(variable,"\\.")]
		var.est[,est_type:="norm"]
		var.est[,variable:=NULL]

		# add groups which are not in var.est
		# because they were not present in sample
		if(!is.null(period.mean)|period.diff.b){
		  # if group did not exist in period
		  # add group with NA
		  roll.miss <- unique(dt.eval("dat[",na.eval,",list(",by.eval,")]"))
		  roll.miss <- lapply(roll.miss,function(l){unique(na.omit(l))})
		  roll.miss <- c(roll.miss,list(est=var.est[,unique(est)]))
		  roll.est <- data.table(expand.grid(c(roll.miss,list(ID=unique(var.est$ID)))))

		  if(nrow(roll.est)>nrow(var.est)){
		    var.est <- merge(roll.est,var.est,by=c(period,z,"ID","est"),all.x=TRUE)
		    setkeyv(var.est,period)
		    var.est[is.na(V1),N:=0]
		    var.est[is.na(V1),n:=0]
		  }
		}

		# calculate mean of estimate over 'period.mean'- periods
		if(!is.null(periodsList)){
		  roll.est <- var.est[,list(V1=rollMeanC(x=V1,k=period.mean,type="c"),V2=periodsList),by=c("ID",z,"est")]
		  setnames(roll.est,"V2",period)
		  roll.est[,est_type:="roll"]
		  if(!is.null(z)){
		    roll.Nn <- var.est[ID==1&est==var[1],list(N=rollMeanC(x=N,k=period.mean,type="c"),
		                                              n=rollMeanC(x=n,k=period.mean,type="c"),V2=periodsList),by=c(z)]
		  }else{
		    roll.Nn <- var.est[ID==1&est==var[1],list(N=rollMeanC(x=N,k=period.mean,type="c"),
		                                              n=rollMeanC(x=n,k=period.mean,type="c"),V2=periodsList)]
		  }

		  setnames(roll.Nn,"V2",period)

		  # merge results
		  roll.est <- merge(roll.est,roll.Nn,by=c(z,period),all.x=TRUE)

		}

		if(period.diff.b){
		  # calcualte differences between periods
		  by.diff <- paste(c("ID","est",z),collapse=",")
		  diff.est <- lapply(period.diff,function(y){
		    y_cond <- paste(period,paste0("c(",paste(y,collapse=","),")"),sep="%in%")
		    diff.y <- dt.eval("var.est[",y_cond,",V1[",period,"==",y[1],"]-V1[",period,"==",y[2],"],by=list(",by.diff,")]")
		    diff.y[,c(period):=paste(y,collapse="-")]
		    return(diff.y)
		  })
		  diff.est <- rbindlist(diff.est)
		  diff.est[,est_type:="diff"]

		  # calcualte N and n for groups and diff
		  diff.Nn <- lapply(period.diff,function(y){
		    y_cond <- paste(period,paste0("c(",paste(y,collapse=","),")"),sep="%in%")
		    diff.y <- dt.eval("var.est[ID==1&est==var[1]&",y_cond,",.(n=mean(n),N=mean(N)),by=list(",by.diff,")]")
		    diff.y[,c(period):=paste(y,collapse="-")]
		    diff.y[,c("est","ID"):=NULL]
		    return(diff.y)
		  })
		  diff.Nn <- rbindlist(diff.Nn)

		  # merge results
		  diff.est <- merge(diff.est,diff.Nn,by=c(z,period),all.x=TRUE)

		  # calcualte differences between periods and mean over consecutive differences
		  if(!is.null(unlist(period.diff.mean))){
		    # i.diff.mean <- (period.mean+1)/2
		    i.diff.mean <- 1
		    diff.mean.est <- lapply(period.diff.mean,function(d){
		      # calculate differences for all pairwise periods in d
		      d.m.est <- lapply(d,function(y){
		        y_cond <- paste(period,paste0("c(",paste(y,collapse=","),")"),sep="%in%")
		        diff.y <- dt.eval("var.est[",y_cond,",V1[",period,"==",y[1],"]-V1[",period,"==",y[2],"],by=list(",by.diff,")]")
		        diff.y[,c(period):=paste(y,collapse="-")]
		        return(diff.y)
		      })
		      # calcualte mean over consecutive differences
		      d.m.est <- rbindlist(d.m.est)
		      d.m.est <- dt.eval("d.m.est[,mean(V1),by=list(",by.diff,")]")
		      d.m.est[,c(period):=paste0(paste(d[[i.diff.mean]],collapse="-"),"-mean")]
		    })
		    diff.mean.est <- rbindlist(diff.mean.est)
		    diff.mean.est[,est_type:="diff_mean"]

		    # calcualte N and n for groups and diff
		    diff.roll.Nn <- lapply(period.diff.mean,function(y){
		      y_cond <- paste(period,paste0("c(",paste(unlist(y),collapse=","),")"),sep="%in%")
		      diff.y <- dt.eval("var.est[ID==1&est==var[1]&",y_cond,",.(n=mean(n),N=mean(N)),by=list(",by.diff,")]")
		      diff.y[,c(period):=paste0(paste(y[[i.diff.mean]],collapse="-"),"-mean")]
		      diff.y[,c("est","ID"):=NULL]
		      return(diff.y)
		    })
		    diff.roll.Nn <- rbindlist(diff.roll.Nn)

		    # merge results
		    diff.mean.est <- merge(diff.mean.est,diff.roll.Nn,by=c(z,period),all.x=TRUE)
		  }
		}

		# add results to one data.table
		if(!is.null(period.mean)){
		  var.est <- rbind(var.est,roll.est,fill=TRUE)
		}
		if(period.diff.b){
		  var.est <- rbind(var.est,diff.est,fill=TRUE)

		  if(!is.null(unlist(period.diff.mean))){
		    var.est <- rbind(var.est,diff.mean.est,fill=TRUE)
		  }
		}

		if(!is.null(p)){
		  sd.est <- var.est[ID!=1,as.list(c(stE=sd(V1),quantileNA(V1,probs=p,p.names=p.names,np=np))),by=c(period,z,'est')]
		}else{
		  sd.est <- var.est[ID!=1,.(stE=sd(V1)),by=c(period,z,'est')]
		}


		out.z <- merge(var.est[ID==1,-"ID"],sd.est,by=c(period,z,"est"))

		if(bias){
		  bias.est <- var.est[ID!=1,.(mean=mean(V1)),by=c(period,z)]
		  out.z <- merge(out.z,bias.est,by=c(period,z))
		}

		# define size groups - groups with zero size do not fall under size-output
		# for groups with zero size the resutling estimatese will be NAs
		out.z[(!is.na(n))&n>0,size:=n<size.limit]

		return(out.z)

	})

	out <- rbindlist(out,fill=TRUE,use.names=TRUE)
	setnames(out,"V1","val")
	return(out)
}
