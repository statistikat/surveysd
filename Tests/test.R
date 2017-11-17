################################################################
# TEST CODE

library(data.table)
library(haven)
library(surveysd)
library(simPop)

dat <- fread("/mnt/meth/Gussenbauer/surveysd/udb_short.csv")
dat <- dat[RB020!="CZ"]
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

erg <- calc.stError(dat=copy(dat_boot_calib[rb020=="ES"]),weights="rb050",year="rb010",b.weights=paste0("w",1:10),var="hx080",cross_var=list(c("db040","db100")),year.diff=c("2015-2009"))

erg

erg$smallGroups


a <- data.table(1:10,LETTERS[1:10])
b <- data.table(20:19,LETTERS[1:10])
b[,V2:=factor(V2)]
a[,V1:=factor(V1)]
merge(b,a,by="V2")
b[,V2:=as.numeric(factor(V2))]

b[,V2:=levels(V2)[V2]]

a <- data.table(1:10,LETTERS[1:10])
b <- data.table(0:9,LETTERS[11:20])
b[,V1:=factor(V1)]
b[,V1:=as.numeric(as.character(V1))]

merge(b,a,by="V1")

# bsp1:
a <- data.table(1:10,LETTERS[1:10])
b <- data.table(20:19,LETTERS[1:10])
b[,V2:=factor(V2)]
setkey(a, V2); setkey(b, V2)
b[a]

# bsp2
a <- data.table(1:10,LETTERS[1:10])
b <- data.table(20:19,LETTERS[1:10])
b[,V1:=factor(V1)]
setkey(a, V1); setkey(b, V1)
b[a]




