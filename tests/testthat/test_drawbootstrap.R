#################################
# test rescaled.bootstrap()
#

context("rescaled.bootstrap()")
library(surveysd)
library(laeken)
library(data.table)

data("eusilc")
setDT(eusilc)
eusilc[,N.households:=sum(db090[!duplicated(db030)]),by=db040]
eusilc[!duplicated(db030),N.households.error:=sum(db090),by=db040]
eusilc[,N.households.all:=sum(db090[!duplicated(db030)])]

# generate yearly data for y years
# 25% drop out from 1 year to the other
y <- 7
eusilc[,year:=2010]
eusilc.i <- copy(eusilc)
nsamp <- round(eusilc[,uniqueN(db030)]*.25)
nextIDs <- (1:nsamp)+eusilc[,max(db030)]
for(i in 1:7){
  eusilc.i[db030%in%sample(unique(eusilc.i$db030),nsamp),db030:=nextIDs[.GRP],by=db030]
  eusilc.i[,year:=year+1]
  eusilc <- rbind(eusilc,eusilc.i)
  nextIDs <- (1:nsamp)+eusilc[,max(db030)]
}
eusilc[,rb030:=as.integer(paste0(db030,"0",.N)),by=list(year,db030)]

# test input parameter
test_that("test para - data",{
  expect_error(draw.bootstrap(as.matrix(eusilc),REP=10,hid="db030",weights="db090",period="year",strata="db040"))
  expect_success(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040"))
})

test_that("test para - REP",{
  expect_error(draw.bootstrap(eusilc,REP="a",hid="db030",weights="db090",period="year",strata="db040"))
  expect_error(draw.bootstrap(eusilc,REP=1:10,hid="db030",weights="db090",period="year",strata="db040"))
})

test_that("test para - hid, weights and period"){
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030s",weights="db090",period="year",strata="db040"))
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090s",period="year",strata="db040"))
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="years",strata="db040"))
  
  eusilc[,year.char :=as.character(year)]
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year.char",strata="db040"))
}

test_that("test para - strata, cluster and totals",{
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db0400"))
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040",cluster="a"))
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040","rb090"),cluster=c("db040","db030"),totals=c("N.households")))
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040"),cluster=c("db040","db030")))
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040"),cluster=c("db030"),totals=c("N.something")))
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040"),cluster=c("db030"),totals=c("N.households.error")))
  
  expect_success(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040"),cluster=c("1")))
  expect_success(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040"),cluster=c("I")))
  expect_success(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=NULL))
  expect_success(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="I"))
  expect_success(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="1"))
  
})


test_that("test para - check.input",{
  expect_success(rescaled.bootstrap(eusilc,REP=10,strata="db040",cluster="db030",fpc="N.households",check.input=FALSE))
  expect_error(rescaled.bootstrap(eusilc,REP=10,strata="db040",cluster="db030",fpc="N.households",check.input="FALSE"))
})

 
test_that("test para - single.PSU",{
  expect_warning(rescaled.bootstrap(eusilc,REP=10,strata="db040",cluster="db030",fpc="N.households",single.PSU="something"))
})

test_that("test return",{
  dat.boot <- rescaled.bootstrap(eusilc,REP=10,strata="db040",cluster="db030",fpc="N.households")
  expect_true(ncol(dat.boot)==(10+ncol(eusilc)))
  dat.unique <- dat.boot[,lapply(.SD,uniqueN),by="db030",.SDcols=c(paste0("bootRep",1:10))]
  dat.unique[,db030:=NULL]
  expect_true(all(dat.unique==1))
  expect_false(any(unlist(dat.boot[,
                               lapply(.SD,function(z){any(is.infinite(z))}),
                               .SDcols=c(paste0("bootRep",1:10))])))
  expect_false(any(is.na(dat.boot[,.SD,.SDcols=c(paste0("bootRep",1:10))])))
  expect_true(all(dat.boot[,.SD,.SDcols=c(paste0("bootRep",1:10))]>0))
  
  dat.boot <- rescaled.bootstrap(eusilc,REP=10,strata="db040",cluster="db030",fpc="N.households",return.value = "replicates")
  expect_true(ncol(dat.boot)==10)
})











