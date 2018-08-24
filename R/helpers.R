##############################################################
# Helper-Functions
#
#' @import ggplot2 matrixStats laeken
#' @import simPop data.table Rcpp
#' @importFrom "graphics" "plot"
#' @importFrom "stats" "as.formula" "na.omit" "quantile" "sd"
#' @importFrom "utils" "data" "find" "tail"
#' @useDynLib surveysd




dt.eval <- function(...,env=parent.frame()){

	expressions <- paste0(...)
	if(length(expressions)>1){
		return(lapply(expressions,function(z){eval(parse(text=z))}))
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


# helpfunction to generate multiple years of eusilc data
demo.eusilc <- function(y=7){
  
  db030 <- rb030 <- povmd60 <- eqincome <- db090 <- age <- hsize <- . <- NULL
  
  data("eusilc")
  setDT(eusilc)
  # generate yearly data for y years
  # 25% drop out from 1 year to the other
  eusilc[,year:=2010]
  eusilc.i <- copy(eusilc)
  nsamp <- round(eusilc[,uniqueN(db030)]*.25)
  hhincome <- eusilc[!duplicated(db030)][["eqIncome"]]
  nextIDs <- (1:nsamp)+eusilc[,max(db030)]
  for(i in 1:y){
    eusilc.i[db030%in%sample(unique(eusilc.i$db030),nsamp),
             c("db030","eqIncome"):=.(nextIDs[.GRP],sample(hhincome,.N)),
             by=db030]
    eusilc.i[,year:=year+1]
    eusilc <- rbind(eusilc,eusilc.i)
    nextIDs <- (1:nsamp)+eusilc[,max(db030)]
  }
  
  eusilc[,rb030:=as.integer(paste0(db030,"0",1:.N)),by=list(year,db030)]
  eusilc[,povmd60:=as.numeric(eqIncome<.6*laeken::weightedMedian(eqIncome[!duplicated(db030)],w=db090[!duplicated(db030)])),by=year]
  eusilc[,age:=cut(age,c(-Inf,16,25,45,65,Inf))]
  eusilc[,hsize:=cut(hsize,c(0:5,Inf))]
  
  return(eusilc)
}

randomInsert <- function(x,y,n=20){
  if(length(x)<20|length(y)<20){
    stop("n must be smaller than length(x) and length(y)")
  }

  x.indices <- sample(length(x),n)
  y.values <- sample(y,n)
  x[x.indices] <- y.values
  return(x)
}

