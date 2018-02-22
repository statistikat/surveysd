#' Title
#'
#' @param dat
#' @param year
#' @param pid
#' @param hid
#' @param ID_name
#'
#' @return
#' @export generate.HHID
#'
#' @examples
#'

generate.HHID <- function(dat,year="RB010",pid="RB030",hid="db030",ID_name="HHID_new"){


  # create lookup table starting from first year
  ID_lookup <- dt.eval("dat[",year,"==min(",year,"),.(",pid,",ID_orig=",hid,")]")
  ID_lookup[,ID_new:=.GRP,by=ID_orig]
  ID_lookup[,ID_orig:=NULL]

  years <- sort(dt.eval("dat[,unique(",year,")]"))

  for(i in years[-1]){

    ID_lookup <- merge(ID_lookup,dt.eval("dat[",year,"==",i,",.(",pid,",ID_orig=",hid,")]"),by=pid,all=TRUE)
    ID_lookup[!is.na(ID_orig),ALL_NEW:=all(is.na(ID_new)),by=ID_orig]
    ID_lookup[ALL_NEW==FALSE,ID_new:=na.omit(ID_new)[1],by=ID_orig]
    ID_next <- ID_lookup[,max(ID_new,na.rm=TRUE)]
    ID_lookup[ALL_NEW==TRUE,ID_new:=.GRP+ID_next,by=ID_orig]
    ID_lookup[,c("ID_orig","ALL_NEW"):=NULL]
  }
  dat <- merge(dat,ID_lookup,by=pid)
  setnames(dat,hid,paste0(hid,"_orig"))
  setnames(dat,"ID_new",hid)
  return(dat)
}
