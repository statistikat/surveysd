#' Generate new houshold ID for survey data with rotating panel design taking into account split households
#'
#' Generating a new houshold ID for survey data using a houshold ID and a personal ID.
#' For surveys with rotating panel design containing housholds, houshold members can move from an existing household to a new one, that was not originally in the sample.
#' This leads to the creation of so called split households. Using a peronal ID (that stays fixed over the whole survey), an indicator for different time steps and a houshold ID,
#' a new houshold ID is assigned to the original and the split household.
#'
#' @param dat data table of data frame containing the survey data
#' @param period column name of \code{dat} containing an indicator for the rotations, e.g years, quarters, months, ect...
#' @param pid column name of \code{dat} containing the personal identifier. This needs to be fixed for an indiviual throught the whole survey
#' @param hid column name of \code{dat} containing the household id. This needs to for a household throught the whole survey
#'
#' @return the survey data \code{dat} as data.table object containing a new and an old household ID. The new household ID which considers the split households is now named \code{hid} and the original household ID
#' has a trailing "_orig".
#' @export generate.HHID
#'
#' @examples
#' 
#' library(surveysd)
#' library(laeken)
#' library(data.table)
#'
#' eusilc <- surveysd:::demo.eusilc() 
#' eusilc[,rb030split:=rb030]
#' year <- eusilc[,unique(year)]
#' year <- year[-1]
#' leaf_out <- c()
#' for(y in year){
#'   split.person <- eusilc[year==(y-1)&!duplicated(db030)&!db030%in%leaf_out,
#'                          sample(rb030,20)]
#'   overwrite.person <- eusilc[year==(y)&!duplicated(db030)&!db030%in%leaf_out,
#'                              .(rb030=sample(rb030,20))]
#'   overwrite.person[,c("rb030split","year_curr"):=.(split.person,y)]
#'   
#'   eusilc[overwrite.person,rb030split:=i.rb030split,on=.(rb030,year>=year_curr)]
#'   leaf_out <- c(leaf_out,
#'                 eusilc[rb030%in%c(overwrite.person$rb030,overwrite.person$rb030split),
#'                 unique(db030)])
#' }
#' 
#' # create spit households
#' eusilc[year>min(year)&!duplicated(db030),
#'        rb030split:=surveysd:::randomInsert(rb030split,eusilc[year==(unlist(.BY)-1)]$rb030,20),
#'        by=year]
#' 
#' # pid which are in split households
#' eusilc[,.(uniqueN(db030)),by=list(rb030split)][V1>1]
#' 
#' eusilc.new <- generate.HHID(eusilc,period="year",pid="rb030split",hid="db030")
#' 
#' # no longer any split households in the data
#' eusilc.new[,.(uniqueN(db030)),by=list(rb030split)][V1>1]
#' 
#'


generate.HHID <- function(dat,period="RB010",pid="RB030",hid="DB030"){
  
  ID_new <- ID_orig <- ALL_NEW <- ID_new_help <- NULL
  
  
  ID_new <- ID_orig <- ALL_NEW <- na.omit
  # check input
  if(!any(class(dat)%in%c("data.frame","data.table"))){
    stop("dat must be a data frame or data table")
  }
  if(class(dat)[1]!="data.table"){
    dat <- data.table(dat)
  }
  dat <- copy(dat)
  c.names <- colnames(dat)
  
  #
  if(!is.character(period)){
    stop("period must be a string")
  }else{
    if(length(period)>1){
      stop("period must have length 1")
    }
    if(!period%in%c.names){
      stop(period," is not a column of dat")
    }
    if(!class(dat[[period]])%in%c("numeric","integer")){
      stop(period," must be an integer or numeric vector")
    }
  }

  if(!is.character(pid)){
    stop("pid must be a string")
  }else{
    if(length(pid)>1){
      stop("pid must have length 1")
    }
    if(!pid%in%c.names){
      stop(pid," is not a column of dat")
    }
  }
  if(!is.character(hid)){
    stop("hid must be a string")
  }else{
    if(length(hid)>1){
      stop("hid must have length 1")
    }
    if(!hid%in%c.names){
      stop(hid," is not a column of dat")
    }
  }

  # create lookup table starting from first period
  ID_lookup <- dt.eval("dat[",period,"==min(",period,"),.(",pid,",ID_orig=",hid,")]")
  ID_lookup[,ID_new:=.GRP,by=ID_orig]
  ID_lookup[,ID_orig:=NULL]

  periods <- sort(dt.eval("dat[,unique(",period,")]"))

  for(i in periods[-1]){

    ID_lookup <- merge(ID_lookup,dt.eval("dat[",period,"==",i,",.(",pid,",ID_orig=",hid,")]"),by=pid,all=TRUE)
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
  group_broke <- dt.eval("dat[,length(unique(ID_new)),by=list(",period,",",hid,")][V1>1,.(",period,",",hid,")]")
  if(nrow(group_broke)>0){
    setkeyv(dat,c(period,hid))
    dt.eval("dat[group_broke,ID_new_help:=paste0(head(",hid,",1),'_1'),by=list(ID_new)]")
    dat[is.na(ID_new_help),ID_new_help:=as.character(ID_new)]
    dat[,ID_new:=.GRP,by=ID_new_help]
    dat[,ID_new_help:=NULL]
  }

  setnames(dat,hid,paste0(hid,"_orig"))
  setnames(dat,"ID_new",hid)
  return(dat)
}
