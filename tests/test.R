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

W16 <- data.table(read_spss(file=paste0(pfad_meth,"/Gussenbauer/surveysd/indicators_var.sav")))
W16 <- data.table(spss.get(paste0(pfad_meth,"/Gussenbauer/surveysd/indicators_var.sav"), use.value.labels = FALSE))
W16_boot <- draw.bootstrap(dat=W16,REP=250,hid="Hid",weights="hgew",strata="bundesld",
                           year="JAHR",split=TRUE,pid="pid")
#head(W16_boot)

dim(W16_boot)

W16_boot_calib <- recalib(dat=copy(W16_boot),hid="Hid",weights="hgew",
                          year="JAHR",b.rep=paste0("w",1:250),conP.var=c("sex"),conH.var = c("bundesld"))
#Error in vars.class[[vars[i]]] : subscript out of bounds  (?!)



W16 <- fread(paste0(pfad_meth,"/Kowarik/Projekte/SILC/Projekt2016/bldaten0816.csv"))

W16_boot <- draw.bootstrap(dat=W16,REP=250,hid="hid",weights="hgew",strata="bundesld",
                           year="jahr",split=TRUE,pid="pid")
#head(W16_boot)

dim(W16_boot)

W16_boot_calib <- recalib(dat=copy(W16_boot),hid="Hid",weights="hgew",
                          year="JAHR",b.rep=paste0("w",1:250),conP.var=c("sex"),conH.var = c("bundesld"))
#Error in vars.class[[vars[i]]] : subscript out of bounds  (?!)




