################################################################
# TEST CODE

library(data.table)
library(haven)
library(surveysd)
library(simPop)

dat <- fread("/mnt/meth/Gussenbauer/surveysd/udb_short.csv")
dat <- dat[RB020!="CZ"]
colnames(dat) <- tolower(colnames(dat))
dat[,rb050:=as.numeric(gsub(",",".",rb050))]
dat[,db030_neu:=paste(rb020,db030,sep="_")]

#dat_boot <- draw.bootstrap(dat,REP=10,hid="db030_neu",weights="rb050",strata="db040",
#                          year="rb010",totals=NULL,boot.names=NULL)

dat_boot <- draw.bootstrap(dat,REP=3,hid="db030",weights="rb050",strata=c("db040"),
                           year="rb010",country="rb020",totals=NULL,boot.names=NULL)


dat_boot <- dat_boot[!is.na(hx080)]
dat_boot[,hx080:=factor(hx080)]
dat_boot_test <- copy(dat_boot)

t <- Sys.time()
dat_boot_calib <- recalib(dat=copy(dat_boot),hid="db030",weights="rb050",
                          year="rb010",country="rb020",b.rep=paste0("w",1),conP.var=c("rb090"),conH.var = c("db040","hx080"),conversion_messages=TRUE)
Sys.time()-t

str(dat_boot_calib)

microbenchmark(HID_FACTOR=recalib(dat=copy(dat_boot),hid="db030",weights="rb050",
                       year="rb010",country="rb020",b.rep=paste0("w",1),conP.var=c("rb090"),conH.var = c("db040","hx080"),conversion_messages=TRUE),
               HID_NOFACTOR=recalib(dat=copy(dat_boot),hid="db030",weights="rb050",
                       year="rb010",country="rb020",b.rep=paste0("w",1),conP.var=c("rb090"),conH.var = c("db040","hx080"),conversion_messages=TRUE,hidf_factor = FALSE),times=20)



dat_boot_test[,rb090:=factor(rb090)]
dat_boot_test[,db040:=factor(db040)]
dat_boot_test[,rb010:=factor(rb010)]
t <- Sys.time()
dat_boot_calib <- recalib(dat=copy(dat_boot_test),hid="db030",weights="rb050",
                          year="rb010",country="rb020",b.rep=paste0("w",1),conP.var=c("rb090"),conH.var = c("db040","hx080"),hidf_factor = FALSE,conversion_messages=TRUE)
Sys.time()-t

microbenchmark(HID_FACTOR=recalib(dat=copy(dat_boot_test),hid="db030",weights="rb050",
                                  year="rb010",country="rb020",b.rep=paste0("w",1),conP.var=c("rb090"),conH.var = c("db040","hx080")),
               HID_NOFACTOR=recalib(dat=copy(dat_boot_test),hid="db030",weights="rb050",
                                    year="rb010",country="rb020",b.rep=paste0("w",1),conP.var=c("rb090"),conH.var = c("db040","hx080"),hidf_factor = FALSE))



dat_boot_calib

erg <- calc.stError(dat=copy(dat_boot_calib),weights="rb050",year="rb010",b.weights=paste0("w",1:110),
                    var="hx080",cross_var=list(c("db040","db100")),year.diff=c("2016-2008"),
                    p=c(.01,.05,.1,.9,.95,.99))

erg

erg$smallGroups


year.diff=c("84984-51515","36235-2525","134525-2525")


a <- data.table(1:10,LETTERS[1:10])
b <- data.table(20:19,LETTERS[1:10])
b[,V2:=factor(V2)]
a[,V1:=factor(V1)]
merge(b,a,by="V2")
b[,V2:=as.numeric(factor(V2))]

b[,V2:=levels(V2)[V2]]

a <- data.table(1:10,LETTERS[1:10])
b <- data.table(0:9,LETTERS[11:20])
b[,V1:=factor(V1)]
b[,V1:=as.numeric(as.character(V1))]

merge(b,a,by="V1")

# bsp1:
a <- data.table(1:10,LETTERS[1:10])
b <- data.table(20:19,LETTERS[1:10])
b[,V2:=factor(V2)]
setkey(a, V2); setkey(b, V2)
b[a]

# bsp2
a <- data.table(1:10,LETTERS[1:10])
b <- data.table(20:19,LETTERS[1:10])
b[,V1:=factor(V1)]
setkey(a, V1); setkey(b, V1)
b[a]

#####################################################
# TEST INPUT CHECKING
a <- list(1:1000)
recalib(dat=copy(dat_boot),hid="db030",weights=c("rb050"),
        year="rb010",country="rb020",b.rep=paste0("w",1:10),conP.var=c("rb090"),conH.var = c("db040","hx080"))



# TEST C++ code
library(Rcpp)
sourceCpp("src/helpersC.cpp")

comp <- rep(FALSE,10000)
for(i in 1:10000){
  x <- sample(c(1,0,NA_real_),200,prob=c(.45,.45,.1),replace=TRUE)
  w <- sample(30:400,200,replace=TRUE)

  r_res <- weightedRatio(x,w)
  c_res <- weightedRatioC(x,w)

  comp[i] <- abs(r_res-c_res)<1e-10
}

comp <- rep(FALSE,10000)
for(i in 1:10000){
  x <- sample(c(1,0,NA_real_),200,prob=c(.45,.45,.1),replace=TRUE)
  w <- sample(30:400,200,replace=TRUE)

  r_res <- weightedSum(x,w)
  c_res <- weightedSumC(x,w)

  comp[i] <- abs(r_res-c_res)<1e-10
}



#########################################
# teste clustering
#

library(data.table)
library(haven)
library(surveysd)
library(simPop)
library(survey)

dat <- fread("/mnt/obdatenaustausch/NETSILC3/udb_short_new.csv")
dat[,RB050:=gsub(",","\\.",RB050)]
dat[,RB050:=as.numeric(RB050)]
# DB060 CLUSTER
# DB050 STRATA
# CLUSTER WERDEN IN STRATA GEZOGEN

# check was mit lonely PSU passiert

dat[,.N,by=RB020]

dat_es <- dat[RB020=="ES"]
dat_es[is.na(db050)]
dat_es[is.na(DB060)]

dat_es[,fpc1:=length(unique(DB060))/.05,by=list(db050,RB010)]
dat_es[,fpc2:=sum(RB050[!duplicated(db030)]),by=list(DB060,db050,RB010)]

set.seed(1234)
dat_boot_1 <- draw.bootstrap(dat=copy(dat_es),REP=100,hid="db030",weights="RB050",strata="db050",cluster=c("DB060"),
                           year="RB010",totals=c("fpc1","fpc2"),boot.names=NULL)
set.seed(1234)
dat_boot_2 <- draw.bootstrap(dat=copy(dat_es),REP=100,hid="db030",weights="RB050",strata="db050",cluster=c("DB060"),
                           year="RB010",boot.names=NULL)


dat_boot <- dat_boot[!is.na(HX080)]
dat_boot[,hx080:=factor(HX080)]

dat_boot_calib <- recalib(dat=copy(dat_boot),hid="db030",weights="RB050",
                          year="RB010",b.rep=paste0("w",1:1000),conP.var=c("RB090"),conH.var = c("DB040"))

erg <- calc.stError(dat=copy(dat_boot_calib),fun="weightedRatioNat",weights="RB050",year="RB010",b.weights=paste0("w",1:250),
                    var="HX080",cross_var=list(c("DB040","DB100")),year.diff=c("2016-2008"),
                    p=c(.01,.05,.1,.9,.95,.99))


set.seed(1234)
dat_boot_1 <- draw.bootstrap(dat=copy(dat_es),REP=100,hid="db030",weights="RB050",strata="db050",cluster=c("DB060"),
                             year="RB010",totals=c("fpc1","fpc2"),boot.names=NULL)
set.seed(1234)
dat_boot_2 <- draw.bootstrap(dat=copy(dat_es),REP=100,hid="db030",weights="RB050",strata="db050",cluster=c("DB060"),
                             year="RB010",boot.names=NULL)




# TESTE ipu2
library(data.table)
library(surveysd)
dat <- fread("/mnt/obdatenaustausch/NETSILC3/udb_short_new.csv")
dat[,RB050:=gsub(",","\\.",RB050)]
dat[,RB050:=as.numeric(RB050)]

dat_es <- dat[RB020=="ES"]
dat_es[,.(RB010,RB030,DB040,arose,hsize,HX040,db050)]

dat_es[,totals1:=round(sum(RB050[!duplicated(db030)])/400),by=list(RB010,db050)]
dat_es[,totals2:=400]

# define stratified 1-Stage cluster sample
set.seed(1234)
dat_boot <- draw.bootstrap(dat=copy(dat_es),REP=150,hid="db030",weights="RB050",strata="DB040",cluster="DB060",
                           year="RB010",totals=c("totals1","totals2"))

dat_boot_calib <- recalib(dat=copy(dat_boot),hid="db030",weights="RB050",
                          year="RB010",b.rep=paste0("w",1:150),conP.var=c("RB090"),conH.var = c("DB040","HX080"))
save(dat_boot_calib, file="/mnt/meth/Gussenbauer/surveysd/TEST.RData")
load("/mnt/meth/Gussenbauer/surveysd/TEST.RData")
erg <- calc.stError(dat=copy(dat_boot_calib),weights="RB050",year="RB010",b.weights=paste0("w",1:150),
                    var="HX080",cross_var=list("DB040","RB090","hsize",c("DB040","hsize"),c("DB040","RB090")),year.diff=c("2014-2008"),
                    p=c(.01,.05,.1,.9,.95,.99))

# cluster(DB060) in strata(db050) bestehen aus ca 400 Einheiten
# Anzahl cluster in strate = HH_STRATA / 400
# test plot device
class(erg)

plot(erg)

plot(erg,type="grouping",groups=c("RB090","hsize"))

#####################################################
#
library(data.table)
library(surveysd)
dat <- fread("/mnt/obdatenaustausch/NETSILC3/udb_short_new.csv")

dat <- generate.HHID(dat)


#############################
# cluster schätzen
library(data.table)
library(surveysd)
dat <- fread("/mnt/obdatenaustausch/NETSILC3/udb_short_new.csv")
dat[,RB050:=gsub(",","\\.",RB050)]
dat[,RB050:=as.numeric(RB050)]

dat_es <- dat[RB020=="ES"][!duplicated(paste(RB030,RB010,sep="_"))]

dat_es[,db050_old:=db050]

# definiere strate für Spanien
for(i in 2013:2009){
  strat_lookup <- unique(na.omit(dat_es[RB010>=i,.(DB060,db050_neu=db050)]))

  dat_es_i <- dat_es[RB010==(i-1)]
  dat_es_i <- merge(dat_es_i,strat_lookup,by=c("DB060"),all.x=TRUE)
  dat_es_i[,db050_neu:=na.omit(db050_neu)[1],by=list(DB060)]

  na.group <- unique(dat_es_i[is.na(db050_neu),.(DB040,db050)])
  if(nrow(na.group)>0){
    setkeyv(dat_es_i,c("DB040","db050"))
    choose.group <- dat_es_i[na.group,length(unique(na.omit(db050_neu))),by=list(DB040,db050)][V1==1,.(DB040,db050)]
    dat_es_i[choose.group,db050_neu:=na.omit(db050_neu)[1],by=list(DB040,db050)]
  }

  dat_es <- merge(dat_es,dat_es_i[,.(RB030,RB010,db050_neu)],by=c("RB030","RB010"),all.x=TRUE)
  dat_es[RB010==(i-1),db050:=db050_neu]
  dat_es[,db050_neu:=NULL]
}
dat_es[is.na(db050)]
dcast(dat_es[,sum(RB050),by=list(RB010,db050)],db050~RB010)

# Cluster aufteilen
# Anteil census section 35917

random_round <- function(x){
  set.seed(1234)
  x_off <- sum(x-floor(x))
  up_down <- rep(FALSE,length(x))
  if(x_off>0){
    up_down[1:x_off] <- TRUE
    up_down <- sample(up_down)
    x[up_down] <- ceiling(x[up_down])
    x[!up_down] <- floor(x[!up_down])
  }
return(x)
}

strata <- dat_es[,.(STRATA_sum=sum(RB050[!duplicated(db030)])),by=list(db050,RB010)]
strata[,STRATA_ratio:=STRATA_sum/sum(STRATA_sum),by=RB010]
strata[,N.cluster:=random_round(STRATA_ratio*35917),by=RB010]
strata[,N.households:=STRATA_sum/N.cluster]

dat_es <- merge(dat_es,strata[,.(db050,RB010,N.cluster,N.households)],by=c("db050","RB010"))

# define stratified 1-Stage cluster sample
set.seed(1234)
dat_boot <- draw.bootstrap(dat=copy(dat_es),REP=250,hid="db030",weights="RB050",strata="db050",cluster="DB060",
                           year="RB010",totals=c("N.cluster","N.households"))

dat_boot_calib <- recalib(dat=copy(dat_boot),hid="db030",weights="RB050",
                          year="RB010",b.rep=paste0("w",1:250),conP.var=c("RB090"),conH.var = c("db050","DB100"))

write.csv2(dat_boot_calib,file="/mnt/obdatenaustausch/NETSILC3/udb_ES_calib.csv")



# compare survey and own calibration bootstrap sample
#


#####################################################
#
library(mountSTAT)
library(data.table)
library(surveysd)
library(haven)

require(Hmisc)

pfad_meth <- mountWinShare("DatenREG","REG_METHODIK","meth")[1]

# W16 <- data.table(read_spss(file=paste0(pfad_meth,"/Gussenbauer/surveysd/indicators_var.sav")))
W16 <- data.table(spss.get(paste0(pfad_meth,"/Gussenbauer/surveysd/indicators_var.sav"), use.value.labels = FALSE))

set.seed(123)
W16_boot1 <- draw.bootstrap(dat=copy(W16),REP=2000,hid="Hid",weights="hgew",strata="bundesld",
                           year="JAHR",split=TRUE,pid="pid")

set.seed(123)
W16_boot2 <- draw.bootstrap(dat=copy(W16),REP=250,hid="Hid",weights="hgew",strata="bundesld",
                           year="JAHR",split=FALSE,pid="pid")

set.seed(123)
W16_boot3 <- draw.bootstrap(dat=copy(W16),REP=250,hid="Hid",weights="hgew",strata="bundesld",cluster="Hid",
                            year="JAHR",split=FALSE,pid="pid")

setkey(W16_boot1,Hid)
setkey(W16_boot2,Hid)
setkey(W16_boot3,Hid)

all.equal(W16_boot1[,c(paste0("w",1:250))],W16_boot2[,c(paste0("w",1:250))])
all.equal(W16_boot2[,c(paste0("w",1:250))],W16_boot3[,c(paste0("w",1:250))])
#head(W16_boot)

dim(W16_boot)

W16_boot_calib1 <- recalib(dat=copy(W16_boot1),hid="Hid",weights="hgew",
                          year="JAHR",b.rep=paste0("w",1:2000),conP.var=c("sex"),conH.var = c("bundesld"))

W16_boot_calib2 <- recalib(dat=copy(W16_boot2),hid="Hid",weights="hgew",
                          year="JAHR",b.rep=paste0("w",1:250),conP.var=c("sex"),conH.var = c("bundesld"))

W16_boot_calib3 <- recalib(dat=copy(W16_boot3),hid="Hid",weights="hgew",
                           year="JAHR",b.rep=paste0("w",1:250),conP.var=c("sex"),conH.var = c("bundesld"))


help_wRatio <- function(x,w){
  sum(x*w)/sum(w)*100
}

W16_err.est1 <- calc.stError(dat=copy(W16_boot_calib1), weights="hgew",b.weights=paste0("w",1:250),year="JAHR",var="povmd60",
                            fun="help_wRatio",cross_var=NULL,year.diff=NULL,year.mean=NULL)
W16_err.est1$Estimates

W16_err.est2 <- calc.stError(dat=copy(W16_boot_calib2), weights="hgew",b.weights=paste0("w",1:250),year="JAHR",var="povmd60",
                            fun="help_wRatio",cross_var=NULL,year.diff=NULL,year.mean=NULL)
W16_err.est2$Estimates

W16_err.est3 <- calc.stError(dat=copy(W16_boot_calib3), weights="hgew",b.weights=paste0("w",1:250),year="JAHR",var="povmd60",
                             fun="help_wRatio",cross_var=NULL,year.diff=NULL,year.mean=NULL)
W16_err.est3$Estimates

library(laeken)
help_AROP <- function(x,w){
  x.med <- 0.6*laeken::weightedMedian(x, weights=w)
  AROP <- rep(0,length(x))
  AROP[x < x.med] <- 1
  sum(AROP*w)/sum(w)*100
}

W16_err.est <- calc.stError(dat=copy(W16_boot_calib1), weights="hgew",b.weights=paste0("w",1:2000),year="JAHR",var="epinc",
                            fun="help_AROP",cross_var=NULL,year.diff=NULL,year.mean=NULL)
W16_err.est$Estimates



# set up help function that returns only the gini index
help_gini <- function(x,w){
  return(gini(x,w)$value)
}

help_wRatio <- function(x,w){
  sum(x*w)/sum(w)*100
}
#help_wRatio(W16_boot_calib$povmd60,W16_boot_calib$hgew)

W16_err.est <- calc.stError(dat=copy(W16_boot_calib), weights="hgew",b.weights=paste0("w",1:250),year="JAHR",var="epinc",
                            fun="help_gini",cross_var=NULL,year.diff=NULL,year.mean=NULL)
W16_err.est$Estimates

W16_err.est <- calc.stError(dat=copy(W16_boot_calib), weights="hgew",b.weights=paste0("w",1:250),year="JAHR",var="povmd60",
                            fun="help_wRatio",cross_var=NULL,year.diff=NULL,year.mean=NULL)
W16_err.est$Estimates

W16_err.est <- calc.stError(dat=copy(W16_boot_calib), weights="hgew",b.weights=paste0("w",1:250),year="JAHR",var="povmd60",
                            fun="weightedRatio",cross_var=NULL,year.diff=NULL,year.mean=NULL)
W16_err.est$Estimates


W16 <- fread(paste0(pfad_meth,"/Kowarik/Projekte/SILC/Projekt2016/bldaten0816.csv"))

W16_boot <- draw.bootstrap(dat=W16,REP=250,hid="hid",weights="hgew",strata="bundesld",
                           year="jahr",split=TRUE,pid="pid")
#head(W16_boot)

dim(W16_boot)

W16_boot_calib <- recalib(dat=copy(W16_boot),hid="Hid",weights="hgew",
                          year="JAHR",b.rep=paste0("w",1:250),conP.var=c("sex"),conH.var = c("bundesld"))
#Error in vars.class[[vars[i]]] : subscript out of bounds  (?!)





test <- function(fun){
  if(exists(fun,mode="function",inherits=FALSE)){
    stop(paste0("Function ",fun," is undefined"))
  }else{
    exists(fun,mode="function")
  }
}



############################
# teste veste
#
library(data.table)
library(mountSTAT)
library(surveysd)
pfad_meth <- mountWinShare("DatenREG","REG_METHODIK","meth")[1]

load(paste0(pfad_meth,"/Gussenbauer/veste14.RData"))

dat <- copy(stppers12)
REP=10
strata="sn>1"
cluster="kz_z>key_e"
fpc=" gnUnt>anzbtot"
single.PSU=c("merge","mean")
return.value=c("data","replicates")
check.input=TRUE

# dat <- dat[sn%in%stppers12[!duplicated(kz_z),.N!=gnUnt,by=sn][V1==TRUE]$sn]

dat_boot <- rescaled.bootstrap(dat,REP=200,strata=strata,cluster=cluster,
                               fpc=fpc,check.input=TRUE,single.PSU=c("merge"),
                               return.value=c("data"))

any(dat_boot[,mget(paste0("bootRep",1:200))]<0)


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
  if(!any(class(dat)%in%c("data.table","data.frame"))){
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
    stop("dat does not contain the column(s)",check.values)
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
  
  singles <- dt.eval("dat[,sum(!duplicated(",clust.val,")),by=c(by.val)][V1==1]")
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
    dati <- dt.eval("dat[,.(N=",fpc[i],"[1],",clust.val,"=unique(",clust.val,"),f=1,n_prev=1,n_draw_prev=1),by=list(",paste(by.val,collapse=","),")]")
  }else{
    dati <- dt.eval("dat[,.(N=",fpc[i],"[1],",clust.val,"=unique(",clust.val,"),f=f,n_prev=n_prev,n_draw_prev=n_draw_prev),by=list(",paste(by.val,collapse=","),")]")
    dat[,f:=NULL]
    dat[,n_prev:=NULL]
    dat[,n_draw_prev:=NULL]
  }
  
  deltai <- paste0("delta_",i,"_",1:REP)
  dati[,n:=.N,by=c(by.val)]
  # determin number of psu to be drawn
  dati[,n_draw:=select.nstar(n[1],N[1],f[1],n_prev[1],n_draw_prev[1]),by=c(by.val)]
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
  
  dat[,f:=n/N*f]
  dat[,n_prev:=n*n_prev]
  dat[,n_draw_prev:=n_draw*n_draw_prev]
  
  dat[,c("n","N",deltai,"n_draw"):=NULL]
  
}


bootRep <- paste0("bootRep",1:REP)
dat[,c(bootRep):=as.data.table(calc.replicate(n=n.calc,N=N.calc,n_draw=n_draw.calc,delta=delta.calc))]

dat[bootRep10<0]

n=n.calc
N=N.calc
delta=delta.calc

index <- dat[kz_z=="Z018O730K",which=TRUE]
n <- n.calc[index,]
N <- N.calc[index,]
n_draw <- n_draw.calc[index,]
delta <- delta.calc[index,,]


N <- dati[n_draw==0]$N
n <- dati[n_draw==0]$n
f <- dati[n_draw==0]$f
n_prev <- dati[n_draw==0]$n_prev
n_draw_prev <- dati[n_draw==0]$n_draw_prev


neg.index <- c(592,   595,   596,   597,   603,   604)
neg.index

n1 <- n[neg.index,1]
N1 <- N[neg.index,1]
n2 <- n[neg.index,2]
N2 <- N[neg.index,2]

n.star2 <- n_draw[neg.index,2]
n.star1 <- n_draw[neg.index,1]
n.star1 <- floor(n1/(2-(n1/N1))-1)
n.star1 <- rep(1,length(n.star1))

n.star2 <- floor(n2/(2-(n2/N2)*(n1/N1))-1)
n.star2 <- ceiling((n.star1*n2)/(n1/N1*n1*(1-n2/N2)+n.star1)-1)
lambda <- (1-n2/N2)/(n2-n.star2)
lambda <- sqrt((n1/N1)*n.star2*lambda)

(sqrt(n1/n.star1))*lambda

floor(n2/(2-(n2/N2)*n1/N1)-1)


select.nstar <- function(n,N,f,n_prev,n_draw_prev){
  n_draw <- (n*n_draw_prev)/(n_prev*f*(1-n/N)+n_draw_prev)
  n_draw <- max(1,ceiling(n_draw-1))
  if(n_draw<0){
    n_draw <- 1
  }
  return(n_draw)
}

draw.without.replacement <- function(n,n_draw){
  delta <- rep(c(1,0),c(n_draw,n-n_draw))
  delta <- sample(delta)
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




#
# testen
#
n <- n.calc
N <- N.calc
delta <- delta.calc

neg.index <- c(376,   378,   380,   381,   384)
fac.1 <- lambda*(sqrt(n[,i-1]/n_draw[,i-1])*delta[,i-1,])
fac.1[neg.index,]
fac.2 <- (n[,i]/n_draw[,i]*delta[,i,]-1)
fac.2[neg.index,]
rep_out[neg.index,]
n[neg.index,]
N[neg.index,]
n_draw[neg.index,i]

n1 <- n[neg.index,1]
N1 <- N[neg.index,1]
n2 <- n[neg.index,2]
N2 <- N[neg.index,2]
n.star <- n_draw[neg.index,1]

n1 <- 2
N1 <- 200
n2 <- 5
N2 <- 500
n.star <- ceiling(n1/(2-n1/N1)-1)
n.star2 <- floor(n2/(2-n2/N2)-1)


n.star2 <- ceiling((n.star*n2)/(n1*n1/N1*(1-n2/N2)+n.star)-1)
sqrt(n.star2*n1/N1*((1-n2/N2)/(n2-n.star2))) * sqrt(n1/n.star)
lambda[neg.index]



######
# test funktion
x <- 1:1000
y <- sqrt(x/(1000-x+1))
plot(x,y)
