##############################################################
# Helper-Functions
#

dt.eval <- function(...,env=parent.frame(2)){
	
	expressions <- paste0(...)
	if(length(expressions)>1){
		return(lapply(expressions,function(z){eval(parse(text=z))}))
	}else{
		return(eval(parse(text=paste0(...)),envir=env))
	}
}

dt.eval2 <- function(...){
	return(eval(parse(text=paste0(...)),envir=parent.frame(2)))
}

getElement <- function(element, default, list) {
	ifelse(is.null(list[[element]]), default, list[[element]])
}

