#############################
# cluster schätzen
library(data.table)
library(surveysd)

dat <- fread("/mnt/obdatenaustausch/NETSILC3/udb_short_new.csv")
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

# define stratified 1-Stage cluster sample
set.seed(12345)
dat_boot <- draw.bootstrap(dat=copy(dat_es),REP=1000,hid="DB030",weights="RB050",strata="DB050",cluster="DB060",
                           year="RB010",totals=c("N.cluster","N.households"),split=TRUE,pid="RB030")

write.csv2(dat_boot,file="/mnt/obdatenaustausch/NETSILC3/udb_ES_bootweight_1000.csv")

# dat_boot <- fread("/mnt/obdatenaustausch/NETSILC3/udb_ES_bootweight.csv")
# dat_boot[,c(paste0("w",1:1000)):=lapply(.SD,function(z){as.numeric(gsub(",","\\.",z))}),.SDcols=c(paste0("w",1:1000))]
#
# dat_boot <- merge(dat_es,dat_boot[,c("RB010","RB030",paste0("w",1:1000))],by=c("RB010","RB030"),all.x=TRUE)

# dat_boot[,ind1:=paste(hsize,DB040,sep="-")]
# dat_boot[,ind2:=paste(agex,RB090,DB040,sep="-")]
#
# dat_boot_calib <- recalib(dat=copy(dat_boot),hid="DB030",weights="RB050",
#                           year="RB010",b.rep=paste0("w",1:2),conP.var=c("ind2"),conH.var = c("ind1"),maxIter=200,bound=6,epsP=0.025)


dat_boot[,ind3:=paste(agex,RB090,sep="-")]
dat_boot_calib <- recalib(dat=copy(dat_boot),hid="DB030",weights="RB050",
                          year="RB010",b.rep=paste0("w",1:1000),conP.var=c("ind3"),conH.var = c("hsize","DB040"),maxIter=200)

write.csv2(dat_boot_calib,file="/mnt/obdatenaustausch/NETSILC3/udb_ES_calib_1000.csv")
