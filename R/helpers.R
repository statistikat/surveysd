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

# print function for surveysd objects
print.surveysd <- function(sd.result){

  # print number of estimates ~ variables in number of groups using function fun from package pack
  n.estimates <- nrow(sd.result[["Estimates"]])*sum(grepl("^val_*.",colnames(sd.result[["Estimates"]])))
  # print number of subgroup which fall below size
  n.years <- unique(sd.result[["Estimates"]][[sd.result[["param"]][["year"]]]])
  n.years <- length(n.years[!grepl("-|_",n.years)])
  n.groups <- nrow(unique(sd.result[["Estimates"]][,.N,by=c(unique(unlist(sd.result[["param"]][["cross_var"]])))]))

  cat("Calculated point estimates for variable(s)\n\n",paste(sd.result[["param"]][["var"]],sep=","),"\n\nusing function",sd.result[["param"]][["fun"]],"from",sd.result[["param"]][["package"]][1],"\n\n")

  cat("Results hold",n.estimates,"point estimates for",n.years,"years in",n.groups,"subgroups\n")
  cat("\n")
  if(nrow(sd.result[["smallGroups"]])>10|nrow(sd.result[["smallGroups"]])==0){
    cat(nrow(sd.result[["smallGroups"]]),"subgroups contained less than",sd.result[["param"]][["size.limit"]],"observations\n")
  }else{
    cat("Subgroups with less than",sd.result[["param"]][["size.limit"]],"observations\n")
    print(sd.result[["smallGroups"]])
  }
  cat("\n")
  stEtoohigh <- colnames(sd.result[["cvHigh"]])
  stEtoohigh <- stEtoohigh[!stEtoohigh%in%c(sd.result[["param"]][["year"]],unique(unlist(sd.result[["param"]][["cross_var"]])))]
  stEtoohigh <- as.matrix(subset(sd.result[["stEHigh"]],select=stEtoohigh))
  cat("Estimted standard error exceeds",sd.result[["param"]][["stE.limit"]],"% of the the point estimate in",sum(stEtoohigh),"cases\n")
  cat("\n")

}



