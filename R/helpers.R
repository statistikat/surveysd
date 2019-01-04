##############################################################
# Helper-Functions
#
#' @import ggplot2
#' @import data.table Rcpp
#' @importFrom "graphics" "plot"
#' @importFrom "stats" "as.formula" "na.omit" "quantile" "sd" "xtabs" "formula"
#' @importFrom "utils" "data" "find" "tail" "head"
#' @importFrom "matrixStats" "rowProds"
#' @importFrom "laeken" "weightedMedian"
#' @importFrom "methods" "formalArgs"
#' @useDynLib surveysd




dt.eval <- function(...,env=parent.frame()){

	expressions <- paste0(...)
	if(length(expressions)>1){
		return(lapply(expressions,function(z){eval(parse(text=z),envir=env)}))
	}else{
		return(eval(parse(text=paste0(...)),envir=env))
	}
}

dt.eval2 <- function(...){
	return(eval(parse(text=paste0(...)),envir=parent.frame()))
}

getEllipsis <- function(element, default, ell){
  ifelse(is.null(ell[[element]]), default, ell[[element]])
}

paste_ <- function(a,b){
  paste(a,b,sep=".")
}
paste_c <- function(a,b){
  paste(a,b,sep=",")
}
paste_addarg <- function(a,b){

  a <- tstrsplit(a,",")

  return(paste(a[[1]],b,paste(a[2:length(a)],collapse=","),sep=","))
}

# helpfunctions for point estimates
weightedRatioNat <- function(x,w,N){
  weightedRatio(x,w)/N*100
}
weightedRatioR <- function(x,w){
  sum(w[x==1],na.rm=TRUE)/sum(w[!is.na(x)],na.rm=TRUE)*100
}
weightedSumR <- function(x,w){
  sum(as.numeric(x)*w,na.rm=TRUE)
}

povmd <- function(x,w){
  md <- laeken::weightedMedian(x,w)*0.6
  pmd60 <- x<md
  return(as.integer(pmd60))
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

randomInsert <- function(x, y, n = 20){
  if(length(x) < 20 | length(y) < 20){
    stop("n must be smaller than length(x) and length(y)")
  }

  x.indices <- sample(length(x), n)
  y.values <- sample(y, n)
  x[x.indices] <- y.values
  return(x)
}
