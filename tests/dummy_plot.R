########################################################
# create dummy data to plot output of calc.stError
#
#

library(data.table)
library(treemp)
library(d3treeR)

dat <- data.table(year=2008:2016)
# add states
dat.state <- dat[,.(state=1:5),by=year]
# add regions
dat.state.region <- dat.state[,.(region=sample(1:100,sample(3:20,1))),by=list(year,state)]

dat <- rbindlist(list(dat,dat.state,dat.state.region),fill=TRUE,use.names = TRUE)

# generate values
dat[,val:=sin(rnorm(.N,mean=0,sd=.7))^2]
dat[,sd:=sin(rnorm(.N,mean=0,sd=.2))^2]


# add ethnicity
dat <- dat[,.(ethnicity=1:5),by=list(year,state,region)]
# add gender
dat <- dat[,.(gender=1:2),by=list(year,state,region,ethnicity)]





library(treemap)
library(d3treeR)

# example 1 from ?treemap
data(GNI2014)
d3tree2(
  treemap(
    GNI2014
    ,index=c("continent", "iso3")
    ,vSize="population"
    ,vColor="GNI"
    ,type="value"
  )
  , rootname = "World"
)


library(surveysd)
library(laeken)
library(data.table)

data("eusilc")
eusilc <- data.table(eusilc)

years <- 2013:2018
n <- eusilc[!duplicated(db030),.N]
eusilc[,year:=2013]
id.max <- eusilc[,max(db030)]
for(i in years){
  
  eusilc.old <- eusilc[year==i&db030%in%eusilc[!duplicated(db030),sample(db030,ceiling(2/3*n))]]
  
  eusilc.new <- eusilc[year==i&db030%in%eusilc[!duplicated(db030),sample(db030,floor(1/3*n))]]
  eusilc.new[!duplicated(db030),db030:=1:.N+id.max]
  
  id.max <- eusilc.new[,max(db030)]
  
  eusilc.year <- rbind(eusilc.old,eusilc.new)
  eusilc.year[,year:=i+1]
  
  eusilc <- rbind(eusilc,eusilc.year)
}




eusilc.boot <- draw.bootstrap(dat=eusilc,hid="db030",year="year",REP=1000,weights="rb050",strata="db040")


# calibrate weight for bootstrap replicates
eusilc.boot.calib <- recalib(dat=copy(eusilc.boot),hid="db030",weights="rb050",
                          year="year",b.rep=paste0("w",1:1000),conP.var=NULL,conH.var = c("db040"))
