#' @title Draw bootstrap replicates
#'
#' @description Draw bootstrap replicates from survey data using the rescaled bootstrap for stratified multistage sampling, presented by Preston, J. (2009).
#'
#'
#' @param dat either data frame or data table containing the survey sample
#' @param REP integer indicating the number of bootstraps to be drawn
#' @param strata string specifying the column name in `dat` that is used for stratification. For multistage sampling multiple column names can be specified by `strata=c("strata1>strata2>strata3")`.
#' See Details for more information.
#' @param cluster string specifying the column name in `dat` that is used for clustering. For multistage sampling multiple column names can be specified by `cluster=c("cluster1>cluster2>cluster3")`.
#' See Details for more information.
#' @param fpc string specifying the column name in `dat` that contains the number of PSUs at the first stage. For multistage sampling the number of PSUs at each stage must be specified by `strata=c("fpc1>fpc2>fpc3")`.
#' @param single.PSU either "merge" or "mean" defining how single PSUs need to be dealt with.
#' For `single.PSU="merge"` single PSUs at each stage are merged with the strata or cluster with the next least number of PSUs. If multiple of those exist one will be select via random draw
#' For `single.PSU="mean"` single PSUs will get the mean over all bootstrap replicates at the stage which did not contain single PSUs.
#' @param return.value either "data" or "replicates" specifying the return value of the function. For "data" the survey data is returned as class `data.table`, for "replicates" only the bootstrap replicates are returned as `data.table`.
#' @param check.input logical, if TRUE the input will be checked before applying the bootstrap procedure
#' @param new.method logical, if TRUE bootstrap replicates will never be negative even if in some strata the whole population is in the sample. WARNING: This is still experimental and resulting standard errors might be underestimated! Use this if for some strata the whole population is in the sample!
#'
#' @details For specifying multistage sampling designs the column names in `strata`,`cluster` and `fpc` need to seperated by ">".\cr
#' For multistage sampling the strings are read from left to right meaning that the column name before the first ">" is taken as the column for stratification/clustering/number of PSUs at the first and the column after the last ">" is taken as the column for stratification/clustering/number of PSUs at the last stage.
#' If for some stages the sample was not stratified or clustered one must specify this by "1" or "I", e.g. `strata=c("strata1>I>strata3")` if there was no stratification at the second stage or `cluster=c("cluster1>cluster2>I")` if there were no clusters at the last stage.\cr
#' The number of PSUs at each stage is not calculated internally and must be specified for any sampling design.
#' For single stage sampling using stratification this can usually be done by adding over all sample weights of each PSU by each strata-code.\cr
#' Spaces in each of the strings will be removed, so if column names contain spaces they should be renamed before calling this procedure!
#'
#' @return returns the complete data set including the bootstrap replicates or just the bootstrap replicates, depending on `return.value="data"` or `return.value="replicates"` respectively.
#' @export rescaled.bootstrap
#'
#' @references Preston, J. (2009). Rescaled bootstrap for stratified multistage sampling. Survey Methodology. 35. 227-234.
#'
#' @author Johannes Gussenbauer, Statistics Austria
#'
#' @examples
#' data(eusilc, package = "laeken")
#' data.table::setDT(eusilc)
#'
#' eusilc[,N.households:=sum(db090[!duplicated(db030)]),by=db040]
#' eusilc.bootstrap <- rescaled.bootstrap(eusilc,REP=100,strata="db040",
#'                                        cluster="db030",fpc="N.households")
#'
#' eusilc[,new_strata:=paste(db040,rb090,sep="_")]
#' eusilc[,N.housholds:=sum(db090[!duplicated(db030)]),by=new_strata]
#' eusilc.bootstrap <- rescaled.bootstrap(eusilc,REP=100,strata=c("new_strata"),
#'                                        cluster="db030",fpc="N.households")
#'
#' eusilc.bootstrap <- rescaled.bootstrap(eusilc,REP=100,strata=c("new_strata"),
#'                                        cluster="db030",fpc="N.households")
#'
#'


rescaled.bootstrap <- function(dat,REP=1000,strata="DB050>1",cluster="DB060>DB030",fpc="N.cluster>N.households",
                               single.PSU=c("merge","mean"), return.value=c("data","replicates"),check.input=TRUE,new.method=FALSE){
  
  InitialOrder <- N <- SINGLE_BOOT_FLAG <- SINGLE_BOOT_FLAG_FINAL <- f <- n_prev <- n_draw_prev <- sum_prev <- n_draw <- NULL
  
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
    if (!is.data.table(dat) & !is.data.frame(dat)) {
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
    check.values <- check.values[!check.values%in%c("1","I")]
    check.values <- check.values[!check.values%in%colnames(dat)]
    if(length(check.values)>0){
      stop("dat does not contain the column(s): ",check.values)
    }
    
    # check if there are missings in fpc
    if(any(is.na(dat[,.SD,.SDcols=c(fpc)]))){
      stop("missing values not allowed in fpc")
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
  dat <- copy(dat)
  # check if variable f, N, n are in data.table
  overwrite.names <- c("f","N","n","n_prev","n_draw","n_draw_prev")
  overwrite.names <- overwrite.names[overwrite.names%in%colnames(dat)]
  if(length(overwrite.names)>0){
    overwrite.names.new <- paste0("ORIGINAL_",overwrite.names)
    setnames(dat,overwrite.names,overwrite.names.new)
    
    strata[strata%in%overwrite.names] <- overwrite.names.new[strata%in%overwrite.names]
    cluster[cluster%in%overwrite.names] <- overwrite.names.new[cluster%in%overwrite.names]
    fpc[fpc%in%overwrite.names] <- overwrite.names.new[fpc%in%overwrite.names]
  }
  
  # set index for data to return dat in correct order
  dat[,InitialOrder:=.I]
  
  # calculate bootstrap replicates
  stages <- length(strata)
  n <- nrow(dat)
  # define values for final calculation
  n.calc <- matrix(0,nrow=n,ncol=stages)
  N.calc <- matrix(0,nrow=n,ncol=stages)
  n_draw.calc <- matrix(0,nrow=n,ncol=stages)
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
    # check of fpc[i] is well defined
    # fpc[i] can not have different values per each by.val
    check.fpc <- dt.eval("dat[,uniqueN(",fpc[i],"),by=c(by.val)][V1>1]")
    if(nrow(check.fpc)>0){
      stop("values in ", fpc[i]," do vary in some strata-cluster combinations at stage ",i)
    }
    
    singles <- dt.eval("dat[",fpc[i],">1,sum(!duplicated(",clust.val,")),by=c(by.val)][V1==1]")
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
        # save original PSU coding and fpc values to replace changed values bevore returning the data.table
        if(return.value=="data"){
          dt.eval("dat[,",paste0(tail(by.val,1),"_ORIGINALSINGLES"),":=",tail(by.val,1),"]")
          dt.eval("dat[,",paste0(fpc[i],"_ORIGINALSINGLES"),":=",fpc[i],"]")
        }
        
        setkeyv(dat,higher.stages)
        next.PSU <- dt.eval("dat[singles,.(N=sum(!duplicated(",clust.val,"))),by=c(by.val)]")
        
        new.var <- paste0(tail(by.val,1),"_NEWVAR")
        dt.eval("next.PSU[,c(new.var):=next_minimum(N,",tail(by.val,1),"),by=c(higher.stages)]")
        
        next.PSU <- next.PSU[N==1]
        dat <- merge(dat,next.PSU[,mget(c(by.val,new.var))],by=c(by.val),all.x=TRUE)
        # sum over margins
        dt.eval("dat[,",paste0(fpc[i],"_ADD"),":=",fpc[i],"]")
        # assign to new group
        dt.eval("dat[!is.na(",new.var,"),c(tail(by.val,1)):=",new.var,"]")
        dt.eval("dat[,",fpc[i],":=",fpc[i],"[is.na(",new.var,")][1],by=c(by.val)]")
        dt.eval("dat[!is.na(",new.var,"),",fpc[i],":=",fpc[i],"+",paste0(fpc[i],"_ADD"),"]")
        dt.eval("dat[,",fpc[i],":=max(",fpc[i],"),by=c(by.val)]")
        
        dat[,c(new.var,paste0(fpc[i],"_ADD")):=NULL]
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
    if(i==1){
      dati <- dt.eval("dat[,.(N=",fpc[i],"[1],",clust.val,"=unique(",clust.val,"),f=1,n_prev=1,n_draw_prev=1,sum_prev=1),by=list(",paste(by.val,collapse=","),")]")
    }else{
      dati <- dt.eval("dat[,.(N=",fpc[i],"[1],",clust.val,"=unique(",clust.val,"),f=f[1],n_prev=n_prev[1],n_draw_prev=n_draw_prev[1],sum_prev=sum_prev[1]),by=list(",paste(by.val,collapse=","),")]")
      dat[,f:=NULL]
      dat[,n_prev:=NULL]
      dat[,n_draw_prev:=NULL]
      dat[,sum_prev:=NULL]
    }
    
    deltai <- paste0("delta_",i,"_",1:REP)
    dati[,n:=.N,by=c(by.val)]
    # determin number of psu to be drawn
    dati[,n_draw:=select.nstar(n[1],N[1],f[1],n_prev[1],n_draw_prev[1],sum_prev=NULL,new.method=new.method),by=c(by.val)]
    # dati[,n_draw_old:=select.nstar(n[1],N[1],f[1],n_prev[1],n_draw_prev[1],sum_prev=NULL,new.method=FALSE),by=c(by.val)]
    # dati[,n_draw_new:=select.nstar(n[1],N[1],f[1],n_prev[1],n_draw_prev[1],sum_prev=sum_prev[1],new.method=FALSE),by=c(by.val)]
    if(nrow(dati[n_draw==0])>0){
      stop("Resampling 0 PSUs should not be possible! Please report bug in https://github.com/statistikat/surveysd")
    }
    # do bootstrap for i-th stage
    dati[,c(deltai):=as.data.table(replicate(REP,draw.without.replacement(n[1],n_draw[1]),simplify = FALSE)),by=c(by.val)]
    
    # merge with data
    dat <- merge(dat,dati,by=c(by.val,clust.val))
    
    # extract information from data.table and remove again from data table (less memory intensive)
    # only matrices and arrays needed for final calculation
    n.calc[,i] <- dat[,n]
    N.calc[,i] <- dat[,N]
    n_draw.calc[,i] <- dat[,n_draw]
    delta.calc[,i,] <- as.matrix(dat[,mget(deltai)])
    
    dat[,sum_prev:=sum_prev + (sqrt(n_draw*f*(1-n/N)/(n-n_draw)) * sqrt(n_prev/n_draw_prev) * (n/n_draw-1))]
    dat[,f:=n/N*f]
    dat[,n_prev:=n*n_prev]
    dat[,n_draw_prev:=n_draw*n_draw_prev]
    
    dat[,c("n","N",deltai,"n_draw"):=NULL]
    
  }
  
  bootRep <- paste0("bootRep",1:REP)
  dat[,c(bootRep):=as.data.table(calc.replicate(n=n.calc,N=N.calc,n_draw=n_draw.calc,delta=delta.calc))]
  
  if(single.PSU=="mean"){
    dat[!is.na(SINGLE_BOOT_FLAG_FINAL),c(bootRep):=lapply(.SD,function(z){mean(z,na.rm=TRUE)}),by=SINGLE_BOOT_FLAG_FINAL,.SDcols=c(bootRep)]
  }
  
  setkey(dat,InitialOrder)
  if(return.value=="data"){
    # get original col values
    if(length(overwrite.names)>0){
      setnames(dat,overwrite.names.new,overwrite.names)
    }
    # get original values for PSUs and fpc - if singles PSUs have been detected and merged
    if(single.PSU=="merge"){
      c.names <- colnames(dat)
      c.names <- c.names[grepl("_ORIGINALSINGLES",c.names)]
      if(length(c.names)>0){
        drop.names <- gsub("_ORIGINALSINGLES","",c.names)
        dat[,c(drop.names):=NULL]
        setnames(dat,c(c.names),drop.names)
      }
    }
    dat[,c("f","n_prev","n_draw_prev","InitialOrder","sum_prev"):=NULL]
    return(dat)
  }else if(return.value=="replicates"){
    return(dat[,mget(bootRep)])
  }
}

select.nstar <- function(n,N,f,n_prev,n_draw_prev,lambda_prev,sum_prev=NULL,new.method){
  
  if(N==1){
    # if only a single unit in strata
    # then always select this unit
    return(1.0)
  }
  
  if(!is.null(sum_prev)){
    n_draw <- (sum_prev)^2/((1-(n/N))*n_prev*f+(sum_prev)^2)*n
    n_draw <- floor(n_draw)
  }else{
    if(new.method){
      n_draw <- (n*n_draw_prev)/(n_prev*f*(1-n/N)+n_draw_prev)
      n_draw <- floor(n_draw)
    }else{
      n_draw <- floor(n/2)
    }
  }
  
  return(n_draw)
}

draw.without.replacement <- function(n,n_draw){
  delta <- rep(c(1.0,0.0),c(n_draw,n-n_draw))
  if(length(delta)>1){
    delta <- sample(delta)
  }
  return(delta)
}

calc.replicate <- function(n,N,n_draw,delta){
  p <- ncol(n)
  # n_draw <- trunc(n/2)
  # n_draw <- floor(n/(2-rowCumprods(n/N))-1)
  # n_draw[n_draw<1] <- 1
  dimdelta <- dim(delta)
  for(i in 1:p){
    if(i==1){
      lambda <- sqrt(n_draw[,1]*(1-n[,1]/N[,1])/(n[,1]-n_draw[,1]+1))
      rep_out <- 1-lambda+lambda*n[,i]/n_draw[,i]*delta[,i,]
    }else if(i==2){
      lambda <- (1-n[,i]/N[,i])/(n[,i]-n_draw[,i]+1)
      lambda <- sqrt((n[,i-1]/N[,i-1])*n_draw[,i]*lambda)
      rep_out <- rep_out + lambda*(sqrt(n[,i-1]/n_draw[,i-1])*delta[,i-1,]) * (n[,i]/n_draw[,i]*delta[,i,]-1)
    }else{
      lambda <- (1-n[,i]/N[,i])/(n[,i]-n_draw[,i]+1)
      lambda <- sqrt(rowProds(n[,1:(i-1)]/N[,1:(i-1)])*n_draw[,i]*lambda)
      prod_val <- matrix(0,ncol=dimdelta[3],nrow=dimdelta[1])
      for(r in 1:dimdelta[3]){
        prod_val[,r] <- rowProds(sqrt(n[,1:(i-1)]/n_draw[,1:(i-1)])*delta[,1:(i-1),r])
      }
      # rep_out <- rep_out + lambda*rowProds(sqrt(n[,1:(i-1)]/n_draw[,1:(i-1)])*delta[,1:(i-1),]) * (n[,i]/n_draw[,i]*delta[,i,]-1)
      rep_out <- rep_out + lambda*prod_val * (n[,i]/n_draw[,i]*delta[,i,]-1)
    }
  }
  return(rep_out)
}

next_minimum <- function(N,by){
  N_notOne <- N!=1
  by <- by[N_notOne][which.min(N[N_notOne])]
  if(length(by)>1){
    by <- sample(by,1)
  }
  return(by)
}
