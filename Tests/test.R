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

#dat_boot <- draw.bootstrap(dat,REP=10,hid="db030_neu",weights="rb050",strata="db040",
#                          year="rb010",totals=NULL,boot.names=NULL)

dat_boot <- draw.bootstrap(dat,REP=10,hid="db030",weights="rb050",strata=c("db040"),
                           year="rb010",country="rb020",totals=NULL,boot.names=NULL)


dat_boot <- dat_boot[!is.na(hx080)]
dat_boot_test <- copy(dat_boot)

t <- Sys.time()
dat_boot_calib <- recalib(dat=copy(dat_boot),hid="db030",weights="rb050",
                          year="rb010",country="rb020",b.rep=paste0("w",1),conP.var=c("rb090"),conH.var = c("db040","hx080"))
Sys.time()-t

dat_boot_test[,hx080:=factor(hx080)]
dat_boot_test[,rb090:=factor(rb090)]
dat_boot_test[,db040:=factor(db040)]
dat_boot_test[,rb010:=factor(rb010)]
t <- Sys.time()
dat_boot_calib <- recalib(dat=copy(dat_boot_test),hid="db030",weights="rb050",
                          year="rb010",country="rb020",b.rep=paste0("w",1),conP.var=c("rb090"),conH.var = c("db040","hx080"))
Sys.time()-t

dat_boot_calib

erg <- calc.stError(dat=copy(dat_boot_calib),weights="rb050",year="rb010",b.weights=paste0("w",1:110),
                    var="hx080",cross_var=list(c("db040","db100")),year.diff=c("2016-2008"),
                    p=c(.01,.05,.1,.9,.95,.99))

erg

erg$smallGroups


year.diff=c("84984-51515","36235-2525","134525-2525")


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

#####################################################
# TEST INPUT CHECKING
a <- list(1:1000)
recalib(dat=copy(dat_boot),hid="db030",weights=c("rb050"),
        year="rb010",country="rb020",b.rep=paste0("w",1:10),conP.var=c("rb090"),conH.var = c("db040","hx080"))



# TEST C++ code
library(Rcpp)
sourceCpp("src/helpersC.cpp")

comp <- rep(FALSE,10000)
for(i in 1:10000){
  x <- sample(c(1,0,NA_real_),200,prob=c(.45,.45,.1),replace=TRUE)
  w <- sample(30:400,200,replace=TRUE)

  r_res <- weightedRatio(x,w)
  c_res <- weightedRatioC(x,w)

  comp[i] <- abs(r_res-c_res)<1e-10
}

comp <- rep(FALSE,10000)
for(i in 1:10000){
  x <- sample(c(1,0,NA_real_),200,prob=c(.45,.45,.1),replace=TRUE)
  w <- sample(30:400,200,replace=TRUE)

  r_res <- weightedSum(x,w)
  c_res <- weightedSumC(x,w)

  comp[i] <- abs(r_res-c_res)<1e-10
}







