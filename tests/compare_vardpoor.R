###################################################
# teste vardpoor package
#
library(data.table)
library(vardpoor)
library(surveysd)

dat <- fread("/mnt/obdatenaustausch/NETSILC3/udb_short_new.csv")
dat[,RB050:=gsub(",","\\.",RB050)]
dat[,RB050:=as.numeric(RB050)]

dat_at <- dat[rb020=="AT"][!duplicated(paste(RB030,RB010,sep="_"))]



dat_at[,IDd:=paste0("V",.I)]
erg2 <- varpoord(Y="HX080",H="DB040",PSU="DB030",ID_level1 = "DB030",N_h=NULL,
                 ID_level2 = "IDd",w_final="RB050",period="RB010",dataset=dat_at,type="linarpr")


dat_boot_calib <- fread("/mnt/obdatenaustausch/NETSILC3/udb_AT_calib.csv")
dat_boot_calib[,c(paste0("w",1:250)):=lapply(.SD,function(z){as.numeric(gsub(",","\\.",z))}),.SDcols=c(paste0("w",1:250))]
dat_boot_calib[,RB050:=as.numeric(gsub(",","\\.",RB050))]

dat_boot_calib[,HX080_neu:=ifelse(HX080==1,0,1)]
erg_sd <- calc.stError(dat=copy(dat_boot_calib),weights="RB050",year="RB010",b.weights=paste0("w",1:250),
                    var="HX080",cross_var=list(c("DB040")),
                    p=c(.025,.975))
erg_sd_neu <- calc.stError(dat=copy(dat_boot_calib),weights="RB050",year="RB010",b.weights=paste0("w",1:250),
                       var="HX080_neu",cross_var=list(c("DB040")),
                       p=c(.025,.975))
erg_sd_neu$Estimates[RB010==2008]

erg2$all_result[RB010==2008,.(RB010,DB040,value,se)]

erg_comp <- merge(erg_sd_neu$Estimates[RB010==2008,.(RB010, DB040, val_HX080_surveysd=val_HX080_neu, stE_HX080_surveysd=stE_HX080_neu)],
                  erg2$all_result[RB010==2008,.(RB010,DB040,val_HX080_vardpoor=value,stE_HX080_vardpoor=se)],by=c("RB010","DB040"),all=TRUE)

###################################################
# DATA ES
library(data.table)
library(surveysd)
library(vardpoor)
library(mountSTAT)

dat <- fread(paste0(mountO(),"/B/Datenaustausch/NETSILC3/udb_short_new.csv"))
dat[,RB050:=gsub(",","\\.",RB050)]
dat[,RB050:=as.numeric(RB050)]

dat_es <- dat[rb020=="ES"][!duplicated(paste(RB030,RB010,sep="_"))]

dat_es[,DB050_old:=DB050]
dat_es <- dat_es[!is.na(DB030)]

# definiere strate für Spanien
for(i in 2013:2009){
  strat_lookup <- unique(na.omit(dat_es[RB010>=i,.(DB060,DB050_neu=DB050)]))

  dat_es_i <- dat_es[RB010==(i-1)]
  dat_es_i <- merge(dat_es_i,strat_lookup,by=c("DB060"),all.x=TRUE)
  dat_es_i[,DB050_neu:=na.omit(DB050_neu)[1],by=list(DB060)]

  na.group <- unique(dat_es_i[is.na(DB050_neu),.(DB040,DB050)])
  if(nrow(na.group)>0){
    setkeyv(dat_es_i,c("DB040","DB050"))
    choose.group <- dat_es_i[na.group,length(unique(na.omit(DB050_neu))),by=list(DB040,DB050)][V1==1,.(DB040,DB050)]
    dat_es_i[choose.group,DB050_neu:=na.omit(DB050_neu)[1],by=list(DB040,DB050)]
  }

  dat_es <- merge(dat_es,dat_es_i[,.(RB030,RB010,DB050_neu)],by=c("RB030","RB010"),all.x=TRUE)
  dat_es[RB010==(i-1),DB050:=DB050_neu]
  dat_es[,DB050_neu:=NULL]
}
#dat_es[is.na(db050)]
#dcast(dat_es[,sum(RB050),by=list(RB010,db050)],db050~RB010)

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

strata <- dat_es[,.(STRATA_sum=sum(RB050[!duplicated(DB030)])),by=list(DB050,RB010)]
strata[,STRATA_ratio:=STRATA_sum/sum(STRATA_sum),by=RB010]
strata[,N.cluster:=random_round(STRATA_ratio*35917),by=RB010]
strata[,N.households:=STRATA_sum/N.cluster]

dat_es <- merge(dat_es,strata[,.(DB050,RB010,N.cluster,N.households)],by=c("DB050","RB010"))
dat_es[,IDd:=paste0("V",.I)]

N_h <- dat_es[,sum(RB050[!duplicated(DB030)]),by=list(RB010,DB050)]
erg2 <- varpoord(Y="HX080",H="DB050",PSU="DB030",Dom="DB040",ID_level1 = "DB030",N_h=N_h,
                 ID_level2 = "IDd",w_final="RB050",period="RB010",dataset=dat_es,type="linarpr")

dat_boot_calib <- fread(paste0(mountO(),"/B/Datenaustausch/NETSILC3/udb_ES_calib_250_nocluster.csv"))
select.column <- copy(colnames(dat_boot_calib))
# dat_boot_calib[,c(paste0("w",1:500)):=lapply(.SD,function(z){as.numeric(gsub(",","\\.",z))}),.SDcols=c(paste0("w",1:500))]
# dat_boot_calib[,RB050:=as.numeric(gsub(",","\\.",RB050))]

dat_boot_calib[,HX080_neu:=ifelse(HX080==1,0,1)]
erg_sd_neu <- calc.stError(dat=copy(dat_boot_calib),weights="RB050",year="RB010",b.weights=paste0("w",1:250),
                           var="HX080_neu",cross_var=list(c("DB040")),
                           p=c(.025,.975))

dat_boot_calib <- fread(paste0(mountO(),"/B/Datenaustausch/NETSILC3/udb_ES_calib.csv"))
# dat_boot_calib <- dat_boot_calib[,mget(select.column)]
# dat_boot_calib[,c(paste0("w",1:500)):=lapply(.SD,function(z){as.numeric(gsub(",","\\.",z))}),.SDcols=c(paste0("w",1:500))]
# dat_boot_calib[,RB050:=as.numeric(gsub(",","\\.",RB050))]

dat_boot_calib[,HX080_neu:=ifelse(HX080==1,0,1)]
erg_sd_cluster <- calc.stError(dat=copy(dat_boot_calib),weights="RB050",year="RB010",b.weights=paste0("w",1:1000),
                           var="HX080_neu",cross_var=list(c("DB040")),
                           p=c(.025,.975))

erg_comp <- merge(erg2$all_result[RB010==2008,.(RB010,DB040,stE_HX080_vardpoor=se)],
                  erg_sd_neu$Estimates[RB010==2008,.(RB010, DB040, val_HX080_neu,stE_HX080_surveysdNOCLUSTER=stE_HX080_neu)],
                  by=c("RB010","DB040"),all=TRUE)

erg_comp <- merge(erg_comp,erg_sd_cluster$Estimates[RB010==2008,.(RB010, DB040, stE_HX080_surveysdCLUSTER=stE_HX080_neu)],
                  by=c("RB010","DB040"),all=TRUE)
erg_comp[,c(1,2,4,3,5,6),with=FALSE]


dat_boot_calib <- fread("/mnt/obdatenaustausch/NETSILC3/udb_ES_calib_1000.csv")
dat_boot_calib[,c(paste0("w",1:249)):=lapply(.SD,function(z){as.numeric(gsub(",","\\.",z))}),.SDcols=c(paste0("w",1:249))]
dat_boot_calib[,RB050:=as.numeric(gsub(",","\\.",RB050))]

dat_boot_calib[,HX080_neu:=ifelse(HX080==1,0,1)]
erg_sd <- calc.stError(dat=copy(dat_boot_calib),weights="RB050",year="RB010",b.weights=paste0("w",1:249),
                           var="HX080_neu",cross_var=list(c("DB040")),
                           p=c(.025,.975))

dat_boot_calib <- fread("/mnt/obdatenaustausch/NETSILC3/udb_ES_calib_NOCLUSTER.csv")
dat_boot_calib[,c(paste0("w",1:249)):=lapply(.SD,function(z){as.numeric(gsub(",","\\.",z))}),.SDcols=c(paste0("w",1:249))]
dat_boot_calib[,RB050:=as.numeric(gsub(",","\\.",RB050))]

dat_boot_calib[,HX080_neu:=ifelse(HX080==1,0,1)]
erg_sd_noclust <- calc.stError(dat=copy(dat_boot_calib),weights="RB050",year="RB010",b.weights=paste0("w",1:249),
                       var="HX080_neu",cross_var=list(c("DB040")),
                       p=c(.025,.975))

erg_comp <- merge(erg_sd$Estimates[,.(RB010, DB040, val_HX080_surveysd=val_HX080_neu, stE_HX080_surveysd=stE_HX080_neu)],
                  erg_sd_noclust$Estimates[,.(RB010, DB040, val_HX080_noclust=val_HX080_neu, stE_HX080_noclust=stE_HX080_neu)],by=c("RB010","DB040"),all=TRUE)
erg_comp[,summary(stE_HX080_surveysd-stE_HX080_noclust)]



