############################################################
# TESTE MATCHING VON PID FÜR VERGANGENE JAHRE
#

library(data.table)
library(VIM)

silc <- fread("/mnt/meth/Kowarik/Projekte/SILC/Projekt2016/bldaten0816.csv")
silc <- unique(silc,by="hid")
# mittels VIM funktion
silc[is.na,bildung:=0]

dat_pid <- silc[jahr==2016]
dat_pid_miss <- silc[jahr%in%c(2013:2015)]

# gehe jahr für jahr durch, verwende Haushalt specifische Merkmale
dat_pid_miss[,pid:=NULL]
dat_pid_miss[,hid_alt:=.GRP,by=hid]
setnames(dat_pid_miss,"hid","hid_true")

setkeyv(dat_pid,c("bundesld","sex","bildung","region","alter","hsize"))

dat_pid_miss_1 <- dat_pid_miss[jahr==2015&hsize==1,.(hid_alt,hid_true,bundesld,sex,bildung,region,alter)]

dat_pid_comp <- dat_pid[hsize==1,.(hid,bundesld,sex,bildung,region)]

dat_merge <- merge(dat_pid_comp,dat_pid_miss_1,by=c("bundesld","sex","bildung","region"),all=TRUE,allow.cartesian = TRUE)

dat_merge

