############################################################
# TESTE MATCHING VON PID FÜR VERGANGENE JAHRE
#

library(data.table)
library(VIM)

load("/mnt/meth/Kowarik/Projekte/SILC/StpZ/silcFrameMerged.RData")
#silc <- unique(silc,by="hid")
# mittels VIM funktion
setnames(silcM,"HHID_ENCODE","hid")

setkey(silcM,hid,ID)
silcM[,pid:=1:.N,by=list(hid,jahr)]
silcM[,hhsize:=.N,by=list(hid,jahr)]

dat_pid <- silcM[jahr==2016]
dat_pid_miss <- silcM[jahr%in%c(2013:2015)]

# gehe jahr für jahr durch, verwende Haushalt specifische Merkmale
dat_pid_miss[,hid_alt:=.GRP,by=list(hid,jahr)]
setnames(dat_pid_miss,"hid","hid_true")

setkeyv(dat_pid,c("BDL","GESCHL","FAMST","EDU_HAB_NAT","EC_URTYP"))

dat_pid_miss_1 <- dat_pid_miss[jahr==2015&HH_SIZE==1,.(hid_alt,hid_true,pid,BDL,GESCHL,FAMST,EDU_HAB_NAT,EC_URTYP)]

dat_pid_comp <- dat_pid[HH_SIZE==1,.(hid,pid,BDL,GESCHL,FAMST,EDU_HAB_NAT,EC_URTYP)]

dat_merge <- merge(dat_pid_comp,dat_pid_miss_1,by=c("BDL","GESCHL","FAMST","EDU_HAB_NAT","EC_URTYP"),
                   all=TRUE,allow.cartesian = TRUE)

dat <- silcM[jahr==2016]
dat_miss <- silcM[jahr%in%c(2013:2015)]
dat_miss <- silcM[jahr==2015]

dat_matched <- match_PID(dat=dat,dat_miss=dat_miss,hvars = c("BDL","EC_URTYP"), pvars = c("GESCHL","FAMST","EDU_HAB_NAT"), max.size=5)



####################################################
# lese daten ein
#
library(foreign)
library(data.table)
library(haeven)

silc <- read.spss("/mnt/obdatenaustausch/NETSILC3/udbworkfile07_16.sav",to.data.frame=TRUE)
silc <- data.table(silc)

colnames(silc)

#silc <- silc[,.(RB010,RB020,RB030,RB050,DB040)]
#silc[,unique(DB040)]

# PB140 -> Geburtsjahr
# RB010  -> jahr
# RB020 -> Land
# db030 -> HID
# rb050 -> gewicht
# RB090

# bestimme neue HID die Haushalte über Jahre verknüpft

# definiere variablen
# year <- rb010
# country <- RB020
# region
# year of birth yob
# quarter of birth qob
# sex
# basic activity status bas
# currenc economic status <- econ
# highest ISCED level <- isced
# ID

y2016 <-  c("RB010","RB020","DB040","RB080","RB070","RB090","RB210","PL031","PE040","RB030")
y2015 <-  y2016
y2014 <-  y2016

silc_part <- silc[RB010>2013&RB020=="ES",mget(y2016)]
setnames(silc_part,colnames(silc_part),c("jahr","land","region","gebjahr","gebquartal","geschl","bas","econ","isced","pid"))
silc_part[,hid:=as.numeric(substr(pid,1,nchar(pid)-2))]
setkeyv(silc_part,c("jahr","gebjahr","gebquartal","geschl"))
silc_part[,pid:=NULL]
silc_part[,pid:=1:.N,by=list(hid,jahr)]
silc_part[,hhsize:=.N,by=list(hid,jahr)]

dat <- silc_part[jahr==2015]
dat_miss <- silc_part[jahr==2014]


dat_matched <- match_PID(dat=dat,dat_miss=dat_miss,hvars = "land", pvars = c("gebquartal","geschl","bas","econ"), max.size=5,hsize="hhsize")

# function to match IDs given household and personal variables
match_PID <- function(dat,dat_miss, hvars = c("BDL","EC_URTYP"), pvars = c("GESCHL","FAMST","EDU_HAB_NAT"), max.size=5,hsize="hhsize"){

  # gehe jahr für jahr durch, verwende haushalt und personen Merkmale
  setnames(dat_miss,"hid","hid_orig")

  hsize <- dat[hhsize<=max.size,unique(hhsize)]

  dat_matched <- list()

  for(k in 1:length(hsize)){

    if(hsize[k]>=max.size){
      dat_pid_miss_k <- dat_miss[hhsize>=max.size&pid<=max.size,mget(c("hid_orig","pid",hvars,pvars))]
      dat_pid_comp <- dat[hhsize>=max.size&pid<=max.size,mget(c("hid","pid",hvars,pvars))]
    }else{
      dat_pid_miss_k <- dat_miss[hhsize==hsize[k],mget(c("hid_orig","pid",hvars,pvars))]
      dat_pid_comp <- dat[hhsize==hsize[k],mget(c("hid","pid",hvars,pvars))]
    }

    dcast.form <- as.formula(paste(paste(c(hvars,"hid_orig"),collapse="+"),"pid",sep="~"))
    dat_pid_miss_k <- dcast(dat_pid_miss_k,dcast.form,value.var=pvars)

    dcast.form <- as.formula(paste(paste(c(hvars,"hid"),collapse="+"),"pid",sep="~"))
    dat_pid_comp <- dcast(dat_pid_comp,dcast.form,value.var=pvars)

    merge.by <- colnames(dat_pid_miss_k)
    merge.by <- merge.by[!merge.by%in%c("hid","hid_orig")]

    dat_merge <- merge(dat_pid_comp,dat_pid_miss_k,by=merge.by,
                       all=TRUE,allow.cartesian = TRUE)
    n_matched <- dat_merge[!is.na(hid_orig)&!is.na(hid),unique(hid_orig)]
    n_total <- nrow(dat_pid_miss_k)
    if(n_total*.69<length(n_matched)){
      sf <- sample(seq(.69,.71,by=.01),1)
      n_matched <- sample(n_matched,length(n_matched)*sf)
      dat_merge <- dat_merge[hid_orig%in%n_matched]
    }

    dat_merge[,GRUPPE:=.GRP,by=merge.by]
    dat_merge_neu <- dat_merge[,match_groups(hid,hid_orig),by=c("GRUPPE")]

    #hid <-  dat_merge[GRUPPE%in%unique(GRUPPE)]$hid
    #hid_alt <- dat_merge[GRUPPE%in%unique(GRUPPE)]$hid_alt
    #dat_merge_neu[,.N,by=list(hid_alt,hid)][N>1]

    #dat_compare <- merge(dat_merge_neu,unique(dat_merge[,.(hid_orig)]),by="hid_orig")

    #dat_compare[,hhsize:=hsize[k]]
    dat_matched <- c(dat_matched,list(dat_merge_neu))
  }

  dat_matched <- rbindlist(dat_matched)

  return(na.omit(dat_matched))
}


silc[,hid:=as.numeric(substr(RB030,1,nchar(RB030)-2))]

unique(silc,c("hid","RB010"))[,.N,by=list(hid,RB010)][,table(N)]
library(ggplot2)

ggplot(unique(silc_part,by=c("hid","jahr"))[,.N,by=list(hid)],aes(N))+
  geom_bar(position="dodge")


years <- silc_part[,unique(jahr)]

row_overlapp <- rep(0,length(years)-1)
for(i in length(years):2){
  dat <- silc_part[jahr==years[i]]
  dat_miss <- silc_part[jahr==years[i-1]]

  row_overlapp[i] <- nrow(dat_miss[hid%in%dat$hid])
}

dat <- silc_part[jahr==2015]
dat_miss <- silc_part[jahr==2014]

# betrachte overlapp 2016 und 2015 - wo liegen unterschiede

dat <- silc_part[jahr==2016]
dat_miss <- silc_part[jahr==2015]

nrow(unique(dat_miss[!hid%in%dat$hid],by="hid"))/nrow(unique(dat_miss,by="hid"))


dat_hh <- unique(dat[,.(hid,hhsize)])
dat_miss_hh <- unique(dat_miss[,.(hid,hhsize)])
compare_hh <- merge(dat_hh,dat_miss_hh,by="hid")
compare_hh[hhsize.x!=hhsize.y]

hsize <- dat_miss[,unique(hhsize)]

for(i in 1:length(hsize)){

  dat_h <- dcast(dat[hhsize==hsize[i]],land+hid~pid,value.var=c("gebjahr","gebquartal","geschl","bas","econ","isced"))
  dat_miss_h <- dcast(dat_miss[hhsize==hsize[i]],land+hid~pid,value.var=c("gebjahr","gebquartal","geschl","bas","econ","isced"))

  hh_same <- dat_h[hid%in%dat_miss_h$hid,hid]
  for(j in hh_same){
    dat_h[hid==j]
  }

}





silc_part <- silc[RB010>2013&RB020=="AT",mget(y2016)]
setnames(silc_part,colnames(silc_part),c("jahr","land","region","gebjahr","gebquartal","geschl","bas","econ","isced","pid"))
silc_part[,hid:=as.numeric(substr(pid,1,nchar(pid)-2))]
setkeyv(silc_part,c("jahr","gebjahr","gebquartal","geschl"))
silc_part[,pid:=NULL]
silc_part[,pid:=1:.N,by=list(hid,jahr)]
silc_part[,hhsize:=.N,by=list(hid,jahr)]
unique(silc_part,by="hid")[,summary(hid),by=jahr]

year <- 2015
dat <- silc_part[jahr==year]
dat_miss <- silc_part[jahr==year-1]

dat[hid%in%dat_miss$hid]

