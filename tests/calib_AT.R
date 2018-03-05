#############################
# cluster schätzen
library(data.table)
library(surveysd)
library(mountSTAT)

dat <- fread(paste0(mountO(),"/B/Datenaustausch/NETSILC3/udb_short_new.csv"))
dat[,RB050:=gsub(",","\\.",RB050)]
dat[,RB050:=as.numeric(RB050)]

dat_at <- dat[rb020=="AT"][!duplicated(paste(RB030,RB010,sep="_"))]

# define stratified 1-Stage cluster sample
set.seed(1234)

dat_boot <- draw.bootstrap(dat=copy(dat_at),REP=1000,hid="DB030",weights="RB050",strata="DB040",
                           year="RB010",split=TRUE,pid="RB030")

# write.csv2(dat_boot,file="/mnt/obdatenaustausch/NETSILC3/udb_AT_bootweight_250.csv")

dat_boot_calib <- recalib(dat=copy(dat_boot),hid="DB030",weights="RB050",
                          year="RB010",b.rep=paste0("w",1:1000),conP.var=c("RB090"),conH.var = c("DB040"))

fwrite(dat_boot_calib,file=paste0(mountO(),"/B/Datenaustausch/NETSILC3/udb_AT_calib.csv"))

