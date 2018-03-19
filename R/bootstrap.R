#' @title Draw bootstrap replicates
#'
#' @description Draw bootstrap replicates from survey data with rotating panel design.
#' Survey information, like ID, sample weights, strata and population totals per strata, should be specified to ensure meaningfull survey bootstraping.
#'
#' @usage draw.bootstrap(dat,REP=1000,hid="DB030",weights="RB050",strata="DB040",
#'                      year="RB010",totals=NULL,boot.names=NULL)
#'
#' @param dat either data.frame or data.table containing the survey data with rotating panel design.
#' @param REP integer indicating the number of bootstrap replicates.
#' @param hid character specifying the name of the column in \code{dat} containing the household ID.
#' @param weights character specifying the name of the column in \code{dat} containing the sample weights.
#' @param year character specifying the name of the column in \code{dat} containing the sample years.
#' @param strata character vector specifying the name of the column in \code{dat} by which the population was stratified.
#' If \code{strata} is a vector stratification will be assumed as the combination of column names contained in \code{strata}.
#' Setting in addition \code{cluter} not NULL stratification will be assumed on multiple stages, where each additional entry in \code{strata} specifies the stratification variable for the next lower stage. see Details for more information.
#' @param cluster character vector specifying cluster in the data. If \code{NULL} household ID is taken es the lowest level cluster.
#' @param totals character specifying the name of the column in \code{dat} containing the the totals per strata and/or cluster. Is ONLY optional if \code{cluster} is \code{NULL} or equal \code{hid} and \code{strata} contains one columnname!
#' Then the households per strata will be calcualted using the \code{weights} argument. If clusters and strata for multiple stages are specified \code{totals} needs to be a vector of \code{length(strata)} specifying the column on \code{dat}
#' that contain the total number of PSUs at each stage. \code{totals} is interpreted from left the right, meaning that the first argument corresponds to the number of PSUs at the first and the last argument to the number of PSUs at the last stage.
#' @param boot.names character indicating the leading string of the column names for each bootstrap replica. If NULL defaults to "w".
#' @param country character specifying the name of the column in \code{dat} containing the country name. Is only used if \code{dat} contains data from multiple countries.
#' In this case the bootstep procedure will be applied on each country seperately. If \code{country=NULL} the household identifier must be unique for each household.
#' @param split logical, if TRUE split households are considered using \code{pid}, for more information see Details.
#' @param pid column in \code{dat} specifying the personal identifier. This identifier needs to be unique for each person throught the whole data set.
#' @return the survey data with the number of REP bootstrap replicates added as columns.
#'
#' @details \code{draw.bootstrap} takes \code{dat} and draws \code{REP} bootstrap replicates from it.
#' \code{dat} must be household data where household members correspond to multiple rows with the same household identifier. The data should at least containt the following columns:
#' \itemize{
#'   \item Column indicating the sample year;
#'   \item Column indicating the household ID;
#'   \item Column containing the household sample weights;
#'   \item Columns by which population was stratified during the sampling process.
#' }
#' For single stage sampling design a column the argument \code{totals} is optional, meaning that a column of the number of PSUs at the first stage does not need to be supplied.
#' For this case the number of PSUs is calculated and added to \code{dat} using \code{strata} and \code{weights}. By setting \code{cluster} to NULL single stage sampling design is always assumed and
#' if \code{strata} contains of multiple column names the combination of all those column names will be used for stratification.\cr
#' In the case of multi stage sampling design the argument \code{totals} needs to be specified and needs to have the same number of arguments as \code{strata}.\cr
#'
#' If \code{cluster} is \code{NULL} or does not contain the \code{hid} at the last stage \code{hid} it will automatically be used as the final cluster. If, besides \code{hid}, clustering in additional stages is specified the number of column names in
#' \code{strata} and \code{cluster} (including \code{hid}) must be the same. If for any stage there was no clustering or stratification one can set "1" or "I" for this stage.\cr
#' For example \code{strata=c("REGION","I"),cluster=c("MUNICIPALITY","HID")} would speficy a 2 stage sampling design where at the first stage the municipalities where drawn stratified by regions
#' and at the 2nd stage housholds are drawn in each municipality without stratification.\cr
#'
#' The bootstrap replicates are drawn for each survey year (\code{year}) using the function \code{\link{bootstrap}}.
#' Afterwards the bootstrap replicates for each household are carried forward from the first year the household enters the survey to all the censecutive years it stays in the survey.\cr
#' This ensures that the bootstrap replicates follow the same logic as the sampled households, making the bootstrap replicates more comparable to the actual sample units.\cr
#' If \code{split} ist set to \code{TRUE} and \code{pid} is specified, the bootstrap replicates are carried forward using the personal identifiers instead of the houshold identifier.
#' This takes into account the issue of a houshold splitting up.
#' Any person in this new split household will get the same bootstrap replicate as the person that has come from an other household in the survey.
#' People who enter already existing households will also get the same bootstrap replicate as the other households members had in the previous years.
#'
#' @return Returns a data.table containing the original data as well as the number of \code{REP} columns containing the bootstrap replicates for each repetition.\cr
#' The columns of the bootstrap replicates are by default labeled "w\emph{Number}" where \emph{Number} goes from 1 to \code{REP}.
#' If the column names of the bootstrap replicates should start with a different character or string the parameter \code{boot.names} can be used.
#'
#' @seealso \code{\link[data.table]{data.table}} for more information on data.table objects.
#'
#' @author Johannes Gussenbauer, Alexander Kowarik, Statistics Austria
#'
#' @examples
#' library(data.table)
#'
#' # run on UDB SILC-data
#'
#' # example for SILC-data for Spain
#' # dat_es <- fread("path//to//spanish//data.csv")
#'
#' # approximate Number of clusters if not known
#' strata <- dat_es[,.(STRATA_sum=sum(RB050[!duplicated(DB030)])),by=list(DB050,RB010)]
#' strata[,STRATA_ratio:=STRATA_sum/sum(STRATA_sum),by=RB010]
#' strata[,N.cluster:=random_round(STRATA_ratio*35917),by=RB010]
#' strata[,N.households:=STRATA_sum/N.cluster]#'
#' dat_es <- merge(dat_es,strata[,.(DB050,RB010,N.cluster,N.households)],by=c("DB050","RB010"))
#'
#' dat_boot <- draw.bootstrap(dat=dat_es,REP=250,hid="DB030",weights="RB050",strata=c("DB050","I"),cluster="DB060",
#'                            year="RB010",totals=c("N.cluster","N.households"),split=TRUE,pid="RB030")
#'
#' # example for SILC-data for Austria
#' # dat_at <- fread("path//to//austrian//data.csv")
#' dat_boot <- draw.bootstrap(dat=dat_at,REP=250,hid="DB030",weights="RB050",strata="DB040",
#'                            year="RB010",split=TRUE,pid="RB030")
#'
#'
#' # save bootstrap replicates as .RData
#' save(dat_boot,file="dat_replicates.RData")
#' # or .csv-file
#' write.csv2(dat_boot,file="dat_replicates.csv",row.names=FALSE)
#'
#' @export draw.bootstrap
#' @import data.table


draw.bootstrap <- function(dat,REP=1000,hid,weights,year,strata=NULL,cluster=NULL,totals=NULL,
                           single.PSU=c("merge","mean"),boot.names=NULL,country=NULL,split=FALSE,pid=NULL){

  ##########################################################
  # INPUT CHECKING
  if(class(dat)[1]=="data.frame"){
    dat <- as.data.table(dat)
  }else if(class(dat)[1]!="data.table"){
    stop("dat must be a data.frame or data.table")
  }


  c.names <- colnames(dat)

  # check REP
  if(length(REP)!=1){
    stop("REP must have length 1")
  }
  if(!is.numeric(REP)){
    stop("REP must contain one numeric value")
  }
  if(REP%%1!=0){
    stop("REP cannot have a decimal part")
  }

  # check hid
  if(length(hid)!=1){
    stop("hid must have length 1")
  }
  if(!hid%in%c.names){
    stop(paste0(hid," is not a column in dat"))
  }

  # check weights
  if(length(weights)!=1){
    stop("weights must have length 1")
  }
  if(!weights%in%c.names){
    stop(paste0(weights," is not a column in dat"))
  }
  if(!is.numeric(dt.eval("dat[,",weights,"]"))){
    stop(paste0(weights," must be a numeric column"))
  }

  # check year
  if(length(year)!=1){
    stop("year must have length 1")
  }
  if(!year%in%c.names){
    stop(paste0(year," is not a column in dat"))
  }

  # check design
  if(is.null(strata)){
    strata <- "I"
  }
  if(is.null(cluster)){
    cluster <- hid
  }else{
    if(!hid%in%cluster){
      cluster <- c(cluster,hid)
    }
  }
  if(!all(strata[!strata%in%c("1","I")]%in%c.names)){
    stop("Not all elements in strata are column names in dat")
  }
  if(any(!cluster[!cluster%in%c("1","I")]%in%c.names)){
    stop("Not all names in cluster are column names in dat")
  }

  if(length(cluster)>1){
    if(length(cluster)!=length(strata)){
      stop("strata and cluster need to have the same number of stages!\n Please use either '1' or 'I' if there was no clustering or stratification in one of the stages.")
    }
  }else{

    if(length(strata)>1){
      if(any(c("1","I")%in%strata)){
        stop("When defining multiple strata variables for single stage sampling design\n none of them can be '1' or 'I'.")
      }

      dt.eval("dat[,STRATA_VAR_HELP:=paste(",paste0(strata,collapse=","),",sep='-')]")
      strata <- "STRATA_VAR_HELP"
    }
  }


  # check single.PSUs
  single.PSU <- single.PSU[1]
  if(is.null(single.PSU)){
    warning("single.PSU was not set to either 'merge' or 'mean'!\n Bootstrap replicates for single PSUs cases will be missing!")
  }else{
    if(!single.PSU%in%c("merge","mean")){
      warning("single.PSU was not set to either 'merge' or 'mean'!\n Bootstrap replicates for single PSUs cases will be missing!")
    }
  }


  # check country
  if(!is.null(country)){
    if(length(country)!=1){
      stop("country must have length 1")
    }
    if(!country%in%c.names){
      stop(paste0(country," is not a column in dat"))
    }
  }

  # check boot.names
  if(!is.null(boot.names)){
    if(!grepl("^[[:alpha:]]",boot.names)){
      stop("boot.names must start with an alphabetic character")
    }
  }

  # check split and pid
  if(is.null(split)){
    stop("split needs to be logical")
  }
  if(!is.logical(split)){
    stop("split needs to be logical")
  }
  if(split){
    if(!is.character(pid)){
      stop("when split is TRUE pid needs to be a string")
    }else{
      if(length(pid)>1){
        stop("pid can only have length 1")
      }else{
        if(!pid%in%c.names){
          stop(pid,"is not a column of dat")
        }
      }
    }
  }

  # check totals
  # if clusters are specified the finite population correction factors must be user specified (at least for now)
  # check input for totals
  # if no totals are specified then leave them NULL
  if(is.null(totals)){
    if(length(cluster)==1){
      # if no clusters are specified calculate number of households in each strata
      totals <- "fpc"
      dt.eval("dat[,fpc:=sum(",weights,"[!duplicated(",hid,")]),by=list(",paste(c(strata,country),collapse=","),")]")
    }else{
      # else leave totals NULL
      # if(length(cluster)>1){
      #   stop("If sample ist clusterd at multiple stages the number of Clusters at each stage must be specified!\n")
      # }
      #
      # warning("Number of Clusters is not specified and will therefor be roughly estimated.
      #         \n Resulting bootstrap replicates might be biased. To avoid this define number of clusters in each strata through parameter 'totals'")
      stop("For multistage sampling the number of PSUs at each level needs to be specified!")
    }
    add.totals <- TRUE
  }else{

    if(length(totals)!=length(strata)){
      stop("totals must specified for each stage")
    }
    if(any(!totals%in%c.names)){
      stop("Not all elements in totals are column names in dat")
    }
    if(!any(unlist(dat[,lapply(.SD,is.numeric),.SDcols=c(totals)]))){
      stop("Not all elements in totals are numeric columns in dat")
    }

    add.totals <- FALSE

  }
  ##########################################################

  # define sample design
  strata <- paste(strata,collapse=">")
  cluster <- paste(cluster,collapse=">")
  totals <- paste(totals,collapse=">")

  if(is.null(boot.names)){
    w.names <- paste0("w",1:REP)
  }else{
    w.names <- paste0(boot.names,1:REP)
  }


  # calculate bootstrap replicates
  dat[,c(w.names):=rescaled.bootstrap(dat=copy(.SD),REP=REP,strata=strata,cluster=cluster,fpc=totals,single.PSU = single.PSU,return.value="replicates",check.input=FALSE),by=c(year,country)]

  # keep bootstrap replicates of first year for each household
  if(split){
    dat <- generate.HHID(dat,time.step=year,pid=pid,hid=hid)
  }
  by.c <- paste(c(hid,country),collapse=",")
  dt.eval("dat[,occurence_first_year :=min(",year,"),by=list(",by.c,")]")
  dat.first.occurence <- unique(subset(dt.eval("dat[",year,"==occurence_first_year]"),select=c(hid,w.names)),by=hid)
  dat[,c(w.names):=NULL]
  dat <- merge(dat,dat.first.occurence,by=hid,all.x=TRUE)
  dat[,occurence_first_year:=NULL]

  # remove columns
  if(split){
    dt.eval("dat[,",hid,":=",paste0(hid,"_orig"),"]")
    dat[,c(paste0(hid,"_orig")):=NULL]
  }
  if(add.totals){
    dt.eval("dat[,",totals,":=NULL]")
  }
  if("STRATA_VAR_HELP"%in%colnames(dat)){
    dat[,STRATA_VAR_HELP:=NULL]
  }

  return(dat)
}


