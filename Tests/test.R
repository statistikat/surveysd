################################################################
# TEST CODE

library(data.table)
library(haven)
library(surveysd)
library(simPop)

dat <- fread("/mnt/meth/Gussenbauer/surveysd/udb_short.csv")
colnames(dat) <- tolower(colnames(dat))
dat[,rb050:=as.numeric(gsub(",",".",rb050))]
dat[,db030_neu:=paste(rb020,db030,sep="_")]

dat_boot <- draw.bootstrap(dat,REP=10,hid="db030_neu",weights="rb050",strata="db040",
                          year="rb010",totals=NULL,boot.names=NULL)


dat_boot <- dat_boot[!is.na(hx080)]
dat_boot[,hx080:=factor(hx080)]
dat_boot[,rb090:=factor(rb090)]

dat_boot_calib <- recalib(dat=copy(dat_boot),hid="db030_neu",weights="rb050",
                          year="rb010",b.rep=paste0("w",1:10),conP.var=c("rb090"),conH.var = c("db040","hx080"))

dat_boot_calib

erg <- calc.stError(dat=copy(dat_boot_calib[rb020!="CZ"]),weights="rb050",year="jahr",b.weights=paste0("w",1:10),var="hx080",cross_var=list(c("db040","db100")),year.diff=c("2016-2008"))

erg

erg$smallGroups


