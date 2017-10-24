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


# resample function definieren -> ?sample
resample <- function(x, ...){
  x[sample.int(length(x), ...)]
}

# assign hid to hid_alt via random sampling
# done in each group defined by household and personal variables
match_groups <- function(hid,hid_alt){
  hid <- unique(na.omit(hid))
  hid_alt <- unique(na.omit(hid_alt))

  nhid <- length(hid)
  nhid_alt <- length(hid_alt)

  if(nhid>0&nhid_alt>0){
    if(nhid>nhid_alt){
      id_merge <- cbind(hid_alt,resample(hid)[1:nhid_alt])
    }else if(nhid<nhid_alt){
      id_merge <- cbind(resample(hid_alt)[1:nhid],hid)
    }else{
      id_merge <- cbind(resample(hid_alt),hid)
    }
  }else{
    if(nhid==0){
      id_merge <- cbind(hid_alt,NA_real_)
    }else{
      id_merge <- cbind(NA_real_,hid)
    }
  }

  colnames(id_merge) <- c("hid_alt","hid")
  return(data.table(id_merge))
}

# function to match IDs given household and personal variables
match_PID <- function(dat,dat_miss, hvars = c("BDL","EC_URTYP"), pvars = c("GESCHL","FAMST","EDU_HAB_NAT"), max.size=5,hsize="hhsize"){

  # gehe jahr für jahr durch, verwende haushalt und personen Merkmale
  dat_miss[,hid_alt:=.GRP,by=list(hid,jahr)]
  setnames(dat_miss,"hid","hid_true")

  hsize <- silcM[hhsize<=max.size,unique(hhsize)]

  dat_matched <- list()

  for(k in 1:length(hsize)){

    if(hsize[k]>=max.size){
      dat_pid_miss_k <- dat_pid_miss[jahr==2015&hhsize>=max.size&pid<=max.size,mget(c("hid_alt","hid_true","pid",hvars,pvars))]
      dat_pid_comp <- dat_pid[hhsize>=max.size&pid<=max.size,mget(c("hid","pid",hvars,pvars))]
    }else{
      dat_pid_miss_k <- dat_pid_miss[jahr==2015&hhsize==hsize[k],mget(c("hid_alt","hid_true","pid",hvars,pvars))]
      dat_pid_comp <- dat_pid[hhsize==hsize[k],mget(c("hid","pid",hvars,pvars))]
    }

    dcast.form <- as.formula(paste(paste(c(hvars,"hid_alt","hid_true"),collapse="+"),"pid",sep="~"))
    dat_pid_miss_k <- dcast(dat_pid_miss_k,dcast.form,value.var=pvars)

    dcast.form <- as.formula(paste(paste(c(hvars,"hid"),collapse="+"),"pid",sep="~"))
    dat_pid_comp <- dcast(dat_pid_comp,dcast.form,value.var=pvars)

    merge.by <- colnames(dat_pid_miss_k)
    merge.by <- merge.by[!merge.by%in%c("hid","hid_alt","hid_true")]

    dat_merge <- merge(dat_pid_comp,dat_pid_miss_k,by=merge.by,
                       all=TRUE,allow.cartesian = TRUE)

    dat_merge[,GRUPPE:=.GRP,by=merge.by]
    dat_merge_neu <- dat_merge[,match_groups(hid,hid_alt),by=c("GRUPPE")]

    hid <-  dat_merge[GRUPPE%in%unique(GRUPPE)[3]]$hid
    hid_alt <- dat_merge[GRUPPE%in%unique(GRUPPE)[3]]$hid_alt
    dat_merge_neu[,.N,by=list(hid_alt,hid)][N>1]

    dat_compare <- merge(dat_merge_neu,unique(dat_merge[,.(hid_alt,hid_true)]),by="hid_alt")

    dat_compare[,hhsize:=hsize[k]]
    dat_matched <- c(dat_matched,list(dat_compare))
  }

  dat_matched <- rbindlist(dat_matched)

  return(na.omit(dat_matched))
}


dat <- silcM[jahr==2016]
dat_miss <- silcM[jahr%in%c(2013:2015)]
dat_miss <- silcM[jahr==2015]

dat_matched <- match_PID(dat=dat,dat_miss=dat_miss,hvars = c("BDL","EC_URTYP"), pvars = c("GESCHL","FAMST","EDU_HAB_NAT"), max.size=5)
