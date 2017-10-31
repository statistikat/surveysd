
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

  colnames(id_merge) <- c("hid_orig","hid")
  return(data.table(id_merge))
}

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
    dat_merge <- dat_merge[!is.na(hid_orig)&!is.na(hid)]
    n_matched <- dat_merge[,unique(hid_orig)]
    n_total <- nrow(dat_pid_miss_k)
    if(n_total*.79<length(n_matched)){
      sf <- sample(seq(.72,.75,by=.01),1)
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
