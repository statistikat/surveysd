#' Draw bootstrap replicates
#'
#' Draw bootstrap replicates from survey data using the rescaled bootstrap for stratified multistage sampling, presented by Preston, J. (2009).
#'
#' @usage rescaled.boot(dat,REP=1000,strata="DB050>1",cluster="DB060>DB030",fpc="N.cluster>N.households",
#'                       check.input=TRUE,single.PSU=c("merge"), return.value=c("data"))
#'
#' @param dat either data frame or data table containing the survey sample
#' @param REP integer indicating the number of bootstraps to be drawn
#' @param strata string specifying the column name in \code{dat} that is used for stratification. For multistage sampling multiple column names can be specified by \code{strata=c("strata1>strata2>strata3")}.
#' See Details for more information.
#' @param cluster string specifying the column name in \code{dat} that is used for clustering. For multistage sampling multiple column names can be specified by \code{cluster=c("cluster1>cluster2>cluster3")}.
#' See Details for more information.
#' @param fpc string specifying the column name in \code{dat} that contains the number of PSUs at the first stage. For multistage sampling the number of PSUs at each stage must be specified by \code{strata=c("fpc1>fpc2>fpc3")}.
#' @param single.PSU either "merge" or "mean" defining how single PSUs need to be dealt with.
#' For \code{single.PSU="merge"} single PSUs at each stage are merged with the strata or cluster with the next least number of PSUs. If multiple of those exist one will be select via random draw
#' For \code{single.PSU="mean"} single PSUs will get the mean over all bootstrap replicates at the stage which did not contain single PSUs.
#' @param return.value either "data" or "replicates" specifying the return value of the function. For "data" the survey data is returned as class \code{data.table}, for "replicates" only the bootstrap replicates are returned as \code{data.table}.
#' @param check.input logical, if TRUE the input will be checked before applying the bootstrap procedure
#'
#' @details For specifying multistage sampling designs the column names in \code{strata},\code{cluster} and \code{fpc} need to seperated by ">".\cr
#' For multistage sampling the strings are read from left to right meaning that the column name before the first ">" is taken as the column for stratification/clustering/number of PSUs at the first and the column after the last ">" is taken as the column for stratification/clustering/number of PSUs at the last stage.
#' If for some stages the sample was not stratified or clustered one must specify this by "1" or "I", e.g. \code{strata=c("strata1>I>strata3")} if there was no stratification at the second stage or \code{cluster=c("cluster1>cluster2>I")} if there were no clusters at the last stage.\cr
#' The number of PSUs at each stage is not calculated internally and must be specified for any sampling design.
#' For single stage sampling using stratification this can usually be done by adding over all sample weights of each PSU by each strata-code.\cr
#' Spaces in each of the strings will be removed, so if column names contain spaces they should be renamed before calling this procedure!
#'
#' @return returns the complete data set including the bootstrap replicates or just the bootstrap replicates, depending on \code{return.value="data"} or \code{return.value="replicates"} respectively.
#' @export bootstrap
#'
#' @references Preston, J. (2009). Rescaled bootstrap for stratified multistage sampling. Survey Methodology. 35. 227-234.
#'
#' @author Johannes Gussenbauer, Statistics Austria
#'
#' @examples
#' @import matrixStats


bootstrap <- function(dat,REP=1000,strata="DB050>1",cluster=" DB060>DB030",fpc=" N.cluster>N.households",
                      single.PSU=c("merge","mean"), return.value=c("data","replicates"),check.input=TRUE){

  # prepare input
  input <- c(strata,cluster,fpc)
  input <- gsub("\\s","",input)
  input <- strsplit(input,">")

  strata <- input[[1]]
  cluster <- input[[2]]
  fpc <- input[[3]]

  single.PSU <- single.PSU[1]
  return.value <- return.value[1]

  if(!is.logical(check.input)){
    stop("check.input can only be logical")
  }
  # check input
  if(check.input){
    # check input data
    if(!any(class(dat)%in%c("data.table","data.frame"))){
      stop("dat needs to be a data frame or data table")
    }else{
      dat <- data.table(dat)
    }

    # check REP
    if(!is.numeric(REP)){
      stop("REP needs to be numeric")
    }else{
      if(length(REP)>1){
        warning("REP has length >1 - First argument will be used!")
        REP <- REP[1]
      }
      if(REP%%1!=0){
        stop("REP cannot have a decimal part")
      }
    }

    # check design variables
    if(length(unique(lapply(input,length)))>1){
      stop("strata, cluster, and fpc need to have the same number of arguments separated with '>'")
    }
    check.values <- unlist(input)
    check.values <- check.values[!by.val%in%c("1","I")]
    check.values <- check.values[!check.values%in%colnames(dat)]
    if(length(check.values)>0){
      stop("dat does not contain the column(s)",check.values)
    }

    # check return.value
    if(!return.value%in%c("data","replicates")){
      stop("return.value can only take the values 'data' or 'replicates'")
    }

    # check single.PSU
    if(is.null(single.PSU)){
      warning("single.PSU was not set to either 'merge' or 'mean'!\n Bootstrap replicates for single PSUs cases will be missing!")
    }else{
      if(!single.PSU%in%c("merge","mean")){
        warning("single.PSU was not set to either 'merge' or 'mean'!\n Bootstrap replicates for single PSUs cases will be missing!")
      }
    }
  }

  # set index for data to return dat in correct order
  dat[,InitialOrder:=.I]

  # calculate bootstrap replicates
  stages <- length(strata)
  n <- nrow(dat)
  # define values for final calculation
  n.calc <- matrix(0,nrow=n,ncol=stages)
  N.calc <- matrix(0,nrow=n,ncol=stages)
  delta.calc <- array(0,dim=c(n,stages,REP))

  for(i in 1:stages){

    # define by.val
    if(i>1){
      by.val <- c(strata[1:(i-1)],cluster[1:(i-1)],strata[i])
    }else{
      by.val <- strata[1:i]
    }
    by.val <- by.val[!by.val%in%c("1","I")]

    # define cluster value
    clust.val <- cluster[i]
    if(clust.val%in%c("1","I")){
      clust.val <- paste0("ID_help_",i)
      dat[,c(clust.val):=.I,by=c(by.val)]
    }

    singles <- dt.eval("dat[,sum(!duplicated(",clust.val,")),by=c(by.val)][V1==1]")
    if(nrow(singles)>0){
      # if singel.PSU=="merge" change the coding of the current stage of the single PSU
      # to the coding in the same subgroup, according the next higher stage, with the smallest number of PSUs
      # if multiple smallest exist choose one per random draw
      higher.stages <- by.val[-length(by.val)]
      if(length(higher.stages)==0){
        higher.stages <- by.val
      }
      singles <- unique(subset(singles,select=higher.stages))
      if(single.PSU=="merge"){

        setkeyv(dat,higher.stages)
        next.PSU <- dt.eval("dat[singles,.(N=sum(!duplicated(",clust.val,"))),by=c(by.val)]")

        new.var <- paste0(tail(by.val,1),"_NEWVAR")
        dt.eval("next.PSU[,c(new.var):=next_minimum(N,",tail(by.val,1),"),by=c(higher.stages)]")

        next.PSU <- next.PSU[N==1]
        dat <- merge(dat,next.PSU[,mget(c(by.val,new.var))],by=c(by.val),all.x=TRUE)
        dt.eval("dat[,c(tail(by.val,1)):=",new.var,"]")
        dt.eval[,c(new.var):=NULL]
      }else if(single.PSU=="mean"){
        # if single.PSU="mean" flag the observation as well as the all the observations in the higher group
        singles[,SINGLE_BOOT_FLAG:=paste(higher.stages,.GRP,sep="-"),by=c(higher.stages)]

        dat <- merge(dat,singles,by=c(higher.stages),all.x=TRUE)
        if(!"SINGLE_BOOT_FLAG_FINAL"%in%colnames(dat)){
          dat[,SINGLE_BOOT_FLAG_FINAL:=SINGLE_BOOT_FLAG]
        }else{
          dat[is.na(SINGLE_BOOT_FLAG_FINAL),SINGLE_BOOT_FLAG_FINAL:=SINGLE_BOOT_FLAG]
        }
        dat[,SINGLE_BOOT_FLAG:=NULL]

      }else{
        warning("Single PSUs detected at the following stages:")
        print(dt.eval("dat[,sum(!duplicated(",clust.val,")),by=c(by.val)][V1==1]"))
      }
    }

    # get Stage
    dati <- dt.eval("dat[,.(N=",fpc[i],"[1],",clust.val,"=unique(",clust.val,")),by=list(",paste(by.val,collapse=","),")]")

    deltai <- paste0("delta_",i,"_",1:REP)
    dati[,n:=.N,by=c(by.val)]
    # do bootstrap for i-th stage
    dati[,c(deltai):=as.data.table(replicate(REP,draw.without.replacement(n[1]),simplify = FALSE)),by=c(by.val)]

    # merge with data
    dat <- merge(dat,dati,by=c(by.val,clust.val))

    # extract information from data.table and remove again from data table (less memory intensive)
    # only matrices and arrays needed for final calculation
    n.calc[,i] <- dat[,n]
    N.calc[,i] <- dat[,N]
    delta.calc[,i,] <- as.matrix(dat[,mget(deltai)])

    dat[,c("n","N",deltai):=NULL]

  }

  bootRep <- paste0("bootRep",1:REP)
  dat[,c(bootRep):=as.data.table(calc.replicate(n=n.calc,N=N.calc,delta=delta.calc))]

  if(single.PSU=="mean"){
    dat[!is.na(SINGLE_BOOT_FLAG_FINAL),c(bootRep):=lapply(.SD,function(z){mean(z,na.rm=TRUE)}),by=SINGLE_BOOT_FLAG_FINAL,.SDcols=c(bootRep)]
  }

  setkey(dat,InitialOrder)
  if(return.value=="data"){
    return(dat)
  }else if(return.value=="replicates"){
    return(dat[,mget(bootRep)])
  }
}

draw.without.replacement <- function(n){
  n_draw <- trunc(n/2)
  delta <- rep(c(1,0),c(n_draw,n-n_draw))
  delta <- sample(delta)
  return(delta)
}

calc.replicate <- function(n,N,delta){
  p <- ncol(n)
  ndraw <- trunc(n/2)
  dimdelta <- dim(delta)
  for(i in 1:p){
    if(i==1){
      lambda <- sqrt(ndraw[,1]*(1-n[,1]/N[,1])/(n[,1]-ndraw[,1]))
      rep_out <- 1-lambda+lambda*n[,i]/ndraw[,i]*delta[,i,]
    }else if(i==2){
      lambda <- (1-n[,i]/N[,i])/(n[,i]-ndraw[,i])
      lambda <- sqrt((n[,i-1]/N[,i-1])*ndraw[,i]*lambda)
      rep_out <- rep_out + lambda*(sqrt(n[,i-1]/ndraw[,i-1])*delta[,i-1,]) * (n[,i]/ndraw[,i]*delta[,i,]-1)
    }else{
      lambda <- (1-n[,i]/N[,i])/(n[,i]-ndraw[,i])
      lambda <- sqrt(rowProds(n[,1:(i-1)]/N[,1:(i-1)])*ndraw[,i]*lambda)
      prod_val <- matrix(0,ncol=dimdelta[3],nrow=dimdelta[1])
      for(r in 1:dimdelta[3]){
        prod_val[,r] <- rowProds(sqrt(n[,1:(i-1)]/ndraw[,1:(i-1)])*delta[,1:(i-1),r])
      }
      # rep_out <- rep_out + lambda*rowProds(sqrt(n[,1:(i-1)]/ndraw[,1:(i-1)])*delta[,1:(i-1),]) * (n[,i]/ndraw[,i]*delta[,i,]-1)
      rep_out <- rep_out + lambda*prod_val * (n[,i]/ndraw[,i]*delta[,i,]-1)
    }
  }
  return(rep_out)
}

next_minimum <- function(N,by){
  N_notOne <- N!=1
  by <- by[n_notOne][which.min(N[N_notOne])]
  if(length(by)>1){
    by <- sample(by,1)
  }
  return(by)
}
