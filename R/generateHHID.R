#' Generate new houshold ID for survey data with rotating panel design taking into account split households
#'
#' Generating a new houshold ID for survey data using a houshold ID and a personal ID.
#' For surveys with rotating panel design containing housholds, houshold members can move from an existing household to a new one, that was not originally in the sample.
#' This leads to the creation of so called split households. Using a peronal ID (that stays fixed over the whole survey) an indicator for different time steps and a houshold ID
#' a new houshold ID is assigned to the original and the split houshold.
#'
#' @param dat data table of data frame containing the survey data
#' @param time.step column name of \code{dat} containing an indicator for the rotations, e.g years, quarters, months, ect...
#' @param pid column name of \code{dat} containing the personal identifier. This needs to be fixed for an indiviual throught the whole survey
#' @param hid column name of \code{dat} containing the household id. This needs to for a household throught the whole survey
#'
#' @return the survey data \code{dat} as data.table object containing a new and an old houhsehold ID. The new houshold ID which considers the split households is now named \code{hid} and the original household ID
#' has a trailing "_orig".
#' @export generate.HHID
#'
#' @examples
#'

generate.HHID <- function(dat,time.step="RB010",pid="RB030",hid="DB030"){

  # check input
  if(!any(class(dat)%in%c("data.frame","data.table"))){
    stop("dat needs to be a data frame or data table")
  }
  if(class(dat)[1]!="data.table"){
    dat <- data.table(dat)
  }
  c.names <- colnames(dat)

  #
  if(!is.character(time.step)){
    stop("time.stop must be a string")
  }else{
    if(length(time.step)>1){
      stop("time.step must have length 1")
    }
  }

  if(!is.character(pid)){
    stop("pid must be a string")
  }else{
    if(length(pid)>1){
      stop("pid must have length 1")
    }
  }
  if(!is.character(hid)){
    stop("hid must be a string")
  }else{
    if(length(hid)>1){
      stop("hid must have length 1")
    }
  }

  if(!time.step%in%c.names){
    stop(time.step," is not a column of dat")
  }
  if(!pid%in%c.names){
    stop(pid," is not a column of dat")
  }
  if(!hid%in%c.names){
    stop(hid," is not a column of dat")
  }


  # create lookup table starting from first time.step
  ID_lookup <- dt.eval("dat[",time.step,"==min(",time.step,"),.(",pid,",ID_orig=",hid,")]")
  ID_lookup[,ID_new:=.GRP,by=ID_orig]
  ID_lookup[,ID_orig:=NULL]

  time.steps <- sort(dt.eval("dat[,unique(",time.step,")]"))

  for(i in time.steps[-1]){

    ID_lookup <- merge(ID_lookup,dt.eval("dat[",time.step,"==",i,",.(",pid,",ID_orig=",hid,")]"),by=pid,all=TRUE)
    ID_lookup[!is.na(ID_orig),ALL_NEW:=all(is.na(ID_new)),by=ID_orig]
    ID_lookup[ALL_NEW==FALSE,ID_new:=na.omit(ID_new)[1],by=ID_orig]
    ID_next <- ID_lookup[,max(ID_new,na.rm=TRUE)]
    ID_lookup[ALL_NEW==TRUE,ID_new:=.GRP+ID_next,by=ID_orig]
    ID_lookup[,c("ID_orig","ALL_NEW"):=NULL]
  }
  dat <- merge(dat,ID_lookup,by=pid)
  # if ID not unique by hid and year
  # leave original grouping for this year
  # this happens if household splits up and people move to already existing households
  group_broke <- dat[,length(unique(ID_new)),by=c(time.step,hid)][V1>1,mget(c(time.step,hid))]
  if(nrow(group_broke)>0){
    setkeyv(dat,c(time.step,hid))
    dt.eval("dat[group_broke,ID_new_help:=paste0(head(",hid,",1),'_1'),by=list(ID_new)]")
    dat[is.na(ID_new_help),ID_new_help:=as.character(ID_new)]
    dat[,ID_new:=.GRP,by=ID_new_help]
    dat[,ID_new_help:=NULL]
  }

  setnames(dat,hid,paste0(hid,"_orig"))
  setnames(dat,"ID_new",hid)
  return(dat)
}
