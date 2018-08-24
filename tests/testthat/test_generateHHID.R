#################################
# test generate.HHID()
#

context("generate.HHID()")
library(surveysd)
library(laeken)
library(data.table)

eusilc <- surveysd:::demo.eusilc()
eusilc[,rb030split:=rb030]
# create spit households

eusilc[year>min(year)&!duplicated(db030),
       rb030split:=surveysd:::randomInsert(rb030split,eusilc[year==(unlist(.BY)-1)]$rb030,20),
       by=year]

# test input parameter
test_that("test para - data",{
  expect_error(generate.HHID(as.matrix(eusilc),period="year",pid="rb030",hid="db030"),
               "dat must be a data.frame or data.table")
  expect_error(generate.HHID(eusilc,period="year",pid="rb030",hid="db030"),NA)
})


test_that("test para - hid, pid and period",{
  expect_error(generate.HHID(eusilc,period="years",pid="rb030",hid="db030"),
               "years is not a column of dat")
  expect_error(generate.HHID(eusilc,period="year",pid="rb030s",hid="db030"),
               "rb030s is not a column of dat")
  expect_error(generate.HHID(eusilc,period="year",pid="rb030",hid="db030s"),
               "db030s is not a column of dat")
  
  eusilc[,year.char:=as.character(year)]
  expect_error(generate.HHID(eusilc,period="year.char",pid="rb030",hid="db030"),
               "year.char must be an integer or numeric vector")
})

test_that("test return",{
  dat.HHID <- generate.HHID(eusilc,period="year",pid="rb030split",hid="db030")
  dat.HHID <- dat.HHID[,uniqueN(db030),by=rb030split][V1>1]
  expect_true(nrow(dat.HHID)==0)
})




