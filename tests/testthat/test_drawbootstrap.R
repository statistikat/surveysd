#################################
# test draw.bootstrap()
#

context("draw.bootstrap()")
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
eusilc[,rb030:=as.integer(paste0(db030,"0",1:.N)),by=list(year,db030)]

# test input parameter
test_that("test para - data",{
  expect_error(draw.bootstrap(as.matrix(eusilc),REP=10,hid="db030",weights="db090",period="year",strata="db040"),
               "dat must be a data.frame or data.table")
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040"),NA)
})

test_that("test para - REP",{
  expect_error(draw.bootstrap(eusilc,REP="a",hid="db030",weights="db090",period="year",strata="db040"),
               "REP must contain one numeric value")
  expect_error(draw.bootstrap(eusilc,REP=1:10,hid="db030",weights="db090",period="year",strata="db040"),
               "REP must have length 1")
})

test_that("test para - hid, weights and period"{
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030s",weights="db090",period="year",strata="db040"),
               "db030s is not a column in dat")
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090s",period="year",strata="db040"),
               "db090s is not a column in dat")
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="years",strata="db040"),
               "years is not a column in dat")
  
  eusilc[,year.char :=as.character(year)]
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year.char",strata="db040"),
               "year.char is not an integer or numeric column")
})

test_that("test para - strata, cluster and totals",{
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db0400"),
               "Not all elements in strata are column names in dat")
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040",cluster="a"),
               "Not all names in cluster are column names in dat")
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040","rb090"),cluster=c("db040","db030"),totals=c("N.households")),
               "totals must be specified for each stage")
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040"),cluster=c("db040","db030")),
               "strata and cluster need to have the same number of stages!\n Please use either '1' or 'I' if there was no clustering or stratification in one of the stages.",
               fixed=TRUE)
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040"),cluster=c("db030"),totals=c("N.something")),
               "Not all names in totals are column names in dat")
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040"),cluster=c("db030"),totals=c("N.households.error")),
               "Missing values found in column(s): N.households.error",fixed=TRUE)
  
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040"),cluster=c("1")),NA)
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=c("db040"),cluster=c("I")),NA)
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata=NULL),NA)
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="I"),NA)
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="1"),NA)
  
})

test_that("test para - bootnames, split and pid"){
  
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040",split="FALSE"),
               "split needs to be logical")
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040",split=TRUE),
               "when split is TRUE pid needs to be a string")
  
  eusilc[,rb030error:=rb030[1],by=list(year,db030)]
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040",split=TRUE,pid="rb030error"),
               "pid is not unique in each household for each period")
  
  eusilc[,rb030new:=rb030]
  eusilc[year>min(year)][!duplicated(db030),rb030new:=sample(eusilc[year==(unlist(.BY)-1)]$rb030,20),by=year]
  
  eusilc.split <- draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040",split=TRUE,pid="rb030new")
  eusilc.split <- eusilc.split[,lapply(.SD,uniqueN),by=rb030new,.SDcols=paste0("w",1:10)]
  expect_true(all(eusilc.split[,.SD,.SDcols=paste0("w",1:10)]==1))
  
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040",boot.names="1"),
               "boot.names must start with an alphabetic character")
  expect_error(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040",boot.names="weight"),NA)
}

test_that("test para - single.PSU",{
  expect_warning(draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040",single.PSU="something"),
                 "single.PSU was not set to either 'merge' or 'mean'!\n Bootstrap replicates for single PSUs cases will be missing!",
                 fixed=TRUE)
})

test_that("test return",{
  dat.boot <- draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040")
  expect_true(ncol(dat.boot)==(10+ncol(eusilc)))
  dat.unique <- dat.boot[,lapply(.SD,uniqueN),by="db030",.SDcols=c(paste0("w",1:10))]
  dat.unique[,db030:=NULL]
  expect_true(all(dat.unique==1))
  expect_false(any(unlist(dat.boot[,
                               lapply(.SD,function(z){any(is.infinite(z))}),
                               .SDcols=c(paste0("w",1:10))])))
  expect_false(any(is.na(dat.boot[,.SD,.SDcols=c(paste0("w",1:10))])))
  expect_true(all(dat.boot[,.SD,.SDcols=c(paste0("w",1:10))]>0))
  
})




