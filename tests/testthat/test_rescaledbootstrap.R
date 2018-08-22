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

# test input parameter
test_that("test para - data",{
  expect_error(rescaled.bootstrap(as.matrix(eusilc),REP="a",strata="db040",cluster="db030",fpc="N.households"))
  expect_success(rescaled.bootstrap(eusilc,REP=100,strata="db040",cluster="db030",fpc="N.households"))
})

test_that("test para - REP",{
  expect_error(rescaled.bootstrap(eusilc,REP="a",strata="db040",cluster="db030",fpc="N.households"))
  expect_warning(rescaled.bootstrap(eusilc,REP=c(1:10),strata="db040",cluster="db030",fpc="N.households"))
})

test_that("test para - strata, cluster and fpc",{
  expect_error(rescaled.bootstrap(eusilc,REP=100,strata="db0400",cluster="db030",fpc="N.households"))
  expect_error(rescaled.bootstrap(eusilc,REP=100,strata="db040>db030",cluster="db030",fpc="N.households"))
  expect_error(rescaled.bootstrap(eusilc,REP=100,strata="db040",cluster="10",fpc="N.households"))
  expect_error(rescaled.bootstrap(eusilc,REP=100,strata="db040",cluster="db040>db030",fpc="N.households"))
  expect_error(rescaled.bootstrap(eusilc,REP=100,strata="db040>db030",cluster="db040>db030",fpc="N.households>N.something"))
  expect_error(rescaled.bootstrap(eusilc,REP=100,strata="I",cluster="db030",fpc="N.households"))
  expect_error(rescaled.bootstrap(eusilc,REP=100,strata="1",cluster="db030",fpc="N.households"))
  expect_error(rescaled.bootstrap(eusilc,REP=100,strata="db040",cluster="db030",fpc="N.households.error"))
  
  expect_success(rescaled.bootstrap(eusilc,REP=100,strata="db040",cluster="I",fpc="N.households"))
  expect_success(rescaled.bootstrap(eusilc,REP=100,strata="db040",cluster="1",fpc="N.households"))
  expect_success(rescaled.bootstrap(eusilc,REP=100,strata="1",cluster="db030",fpc="N.households.all"))
  expect_success(rescaled.bootstrap(eusilc,REP=100,strata="I",cluster="db030",fpc="N.households.all"))
})


test_that("test para - check.input",{
  expect_success(rescaled.bootstrap(eusilc,REP=100,strata="db040",cluster="db030",fpc="N.households",check.input=FALSE))
  expect_error(rescaled.bootstrap(eusilc,REP=100,strata="db040",cluster="db030",fpc="N.households",check.input="FALSE"))
})

 
test_that("test para - single.PSU",{
  expect_warning(rescaled.bootstrap(eusilc,REP=100,strata="db040",cluster="db030",fpc="N.households",single.PSU="something"))
})

test_that("test return",{
  dat.boot <- rescaled.bootstrap(eusilc,REP=100,strata="db040",cluster="db030",fpc="N.households")
  expect_true(ncol(dat.boot)==(100+ncol(eusilc)))
  dat.unique <- dat.boot[,lapply(.SD,uniqueN),by="db030",.SDcols=c(paste0("bootRep",1:100))]
  dat.unique[,db030:=NULL]
  expect_true(all(dat.unique==1))
  expect_false(any(unlist(dat.boot[,
                               lapply(.SD,function(z){any(is.infinite(z))}),
                               .SDcols=c(paste0("bootRep",1:100))])))
  expect_false(any(is.na(dat.boot[,.SD,.SDcols=c(paste0("bootRep",1:100))])))
  expect_true(all(dat.boot[,.SD,.SDcols=c(paste0("bootRep",1:100))]>0))
  
  dat.boot <- rescaled.bootstrap(eusilc,REP=100,strata="db040",cluster="db030",fpc="N.households",return.value = "replicates")
  expect_true(ncol(dat.boot)==100)
})











