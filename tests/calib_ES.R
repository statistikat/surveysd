#############################
# cluster schätzen
library(data.table)
library(surveysd)

dat <- fread("/mnt/obdatenaustausch/NETSILC3/udb_short_new.csv")
dat[,RB050:=gsub(",","\\.",RB050)]
dat[,RB050:=as.numeric(RB050)]

dat_es <- dat[RB020=="ES"][!duplicated(paste(RB030,RB010,sep="_"))]

dat_es <- generate.HHID(dat_es)

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

strata <- dat_es[,.(STRATA_sum=sum(RB050[!duplicated(db030)])),by=list(db050,RB010)]
strata[,STRATA_ratio:=STRATA_sum/sum(STRATA_sum),by=RB010]
strata[,N.cluster:=random_round(STRATA_ratio*35917),by=RB010]
strata[,N.households:=STRATA_sum/N.cluster]

dat_es <- merge(dat_es,strata[,.(db050,RB010,N.cluster,N.households)],by=c("db050","RB010"))

# define stratified 1-Stage cluster sample
set.seed(1234)
dat_boot <- draw.bootstrap(dat=copy(dat_es),REP=250,hid="db030",weights="RB050",strata="db050",cluster="DB060",
                           year="RB010",totals=c("N.cluster","N.households"))

write.csv2(dat_boot,file="/mnt/obdatenaustausch/NETSILC3/udb_ES_bootweight.csv")

dat_boot_calib[,db030:=db030_orig]
dat_boot_calib[,db030_orig:=NULL]
dat_boot_calib <- recalib(dat=copy(dat_boot),hid="db030",weights="RB050",
                          year="RB010",b.rep=paste0("w",1:250),conP.var=c("RB090"),conH.var = c("db050"))

write.csv2(dat_boot_calib,file="/mnt/obdatenaustausch/NETSILC3/udb_ES_calib.csv")
