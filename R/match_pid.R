############################################################
# MATCH ID FOR SPAIN AND AUSTRIA
#
#
# lese daten ein
#
library(foreign)
library(data.table)
source("R/match_pid_help.R")

silc <- read.spss("/mnt/obdatenaustausch/NETSILC3/udbworkfile07_16.sav",to.data.frame=TRUE)
silc <- data.table(silc)

country <- "ES"

# definiere Spalten für "jahr","land","region","gebjahr","gebquartal","geschl","bas","econ","isced","pid","hid
if(country=="ES"){
  start.year <- 2015 # jahr ab dem IDs noch übereinstimmem
}else if(country=="AT"){
  start.year <- 2014 # jahr ab dem IDs noch übereinstimmem
}

pvars <- c("gebjahr","gebquartal","geschl") # personen Variablen an denen gematched wird
hvars <- c("region") # household variablen an denen gematched wird

y2016 <-  c("RB010","RB020","DB040","RB080","RB070","RB090","RB210","PL031","PE040","db030","pid")

#silc_part1 <- silc[RB010>2007&RB010<2010&RB020==country,mget(y2008)]
#silc_part2 <- silc[RB010>2009&RB020==country,mget(y2016)]
#silc_part <- rbind(silc_part1,silc_part2,use.names=FALSE)
silc_part <- silc[RB010>2007&RB020==country,mget(y2016)]
setnames(silc_part,colnames(silc_part),c("jahr","land","region","gebjahr","gebquartal","geschl","bas","econ","isced","hid","pid"))
# definiere hid über PID
#silc_part[,hid_neu:=as.numeric(substr(ID,1,nchar(ID)-2))]
#silc_part[hid_neu!=hid]
#silc_part[,pid:=as.numeric(substr(pid,nchar(pid)-1,nchar(pid)))]
silc_part[,pid:=NULL]
silc_part[,hhsize:=.N,by=list(hid,jahr)]
key_ID <- silc_part[,as.numeric(.GRP),by=list(hid,jahr)]
key_ID[jahr>start.year-1,V1:=hid]
silc_part[jahr<start.year,hid:=as.numeric(.GRP),by=list(hid,jahr)]

silc_part[,summary(hid),by=jahr]

go_back <- start.year:2009
leafe_hid <- c()
for(y in 1:length(go_back)){

  if(length(leafe_hid)>0){
    dat <- silc_part[!hid%in%leafe_hid&jahr==go_back[y]]
  }else{
    dat <- silc_part[jahr==go_back[y]]
  }

  dat_miss <- silc_part[jahr==go_back[y]-1]
  dat_miss[gebjahr==min(gebjahr),gebjahr:=dat[,min(gebjahr)]]
  setkeyv(dat_miss,c("jahr","gebjahr","gebquartal","geschl"))
  setkeyv(dat,c("jahr","gebjahr","gebquartal","geschl"))
  dat_miss[,pid:=1:.N,by=list(hid,jahr)]
  dat[,pid:=1:.N,by=list(hid,jahr)]

  # matche IDs von jahr go_back[y] zu go_back[y]-1
  dat_matched <- match_PID(dat=dat,dat_miss=dat_miss,hvars = "region", pvars = c("gebjahr","gebquartal","geschl"), max.size=5,hsize="hhsize")

  dat_matched[,jahr:=go_back[y]-1]
  setnames(dat_matched,c("hid"),c("hid_new"))
  silc_part <- merge(silc_part,dat_matched[,.(jahr,hid_orig,hid_new)],by.x=c("jahr","hid"),by.y=c("jahr","hid_orig"),all.x=TRUE)
  silc_part[!is.na(hid_new),hid:=hid_new]
  silc_part[,hid_new:=NULL]

  silc_rotn <- unique(silc_part[,.(jahr,hid)],by=c("jahr","hid"))
  setkeyv(silc_rotn,c("jahr","hid"))
  silc_rotn[jahr>=go_back[y]-1,rotn:=1:.N,by=hid]
  leafe_hid <- c(leafe_hid,silc_rotn[rotn==4,hid])

}

# speichere Ergebnis ab
erg <- na.omit(merge(unique(key_ID[,.(hid,V1)]),unique(silc_part[,.(jahr,hid)],by=c("jahr","hid")),by.x="V1",by.y="hid",all.x=TRUE))
setnames(erg,"V1","hid_new")
write.csv2(erg,file=paste0("/mnt/obdatenaustausch/NETSILC3/",country,"_matchedID.csv"))


