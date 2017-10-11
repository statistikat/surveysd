##############################################################
# Helper-Functions
#

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
