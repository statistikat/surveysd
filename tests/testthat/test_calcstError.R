#################################
# test calc.stError()
#

context("calc.stError()")
library(surveysd)
library(laeken)
library(data.table)

data("eusilc")
setDT(eusilc)
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
eusilc[,povmd60:=as.numeric(eqIncome<.6*laeken::weightedMedian(eqIncome[!duplicated(db030)],w=db090[!duplicated(db030)])),by=year]
eusilc[,age:=cut(age,c(-Inf,16,25,45,65,Inf))]
eusilc[,hsize:=cut(hsize,c(0:5,Inf))]

eusilc <- draw.bootstrap(eusilc,REP=10,hid="db030",weights="db090",period="year",strata="db040")
eusilc <- recalib(eusilc,hid="db030",weights="db090",b.rep=paste0("w",1:10),period="year",
                  conP.var="rb090",conH.var="db040")

# test input parameter
test_that("test para - data",{
  expect_error(calc.stError(as.matrix(eusilc),weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040")),"dat must be a data.frame or data.table")
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040")),NA)
})

test_that("test para - weights, b.weights, year and group",{
  expect_error(calc.stError(eusilc,weights="db090",b.weights="a",period="year",var="povmd60",
                            group=c("rb090","db040")),"Not all elements in b.rep are column names in dat")
  expect_error(calc.stError(eusilc,weights="db090",b.weights=1:10,period="year",var="povmd60",
                            group=c("rb090","db040")),"Not all elements in b.rep are column names in dat")
  
  expect_error(calc.stError(eusilc,weights="db090s",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040")),"db090s is not a column in dat")
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="years",var="povmd60",
                            group=c("rb090","db040")),"years is not a column in dat")
  
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090s","db040")),"Not all elements on group are column names in dat")
  
  eusilc.est <- calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                             group=list("rb090","db040",c("rb090","db040"),c("hsize","age")))
  ngroups <- eusilc[,uniqueN(rb090)+uniqueN(db040)+uniqueN(paste(rb090,db040))+uniqueN(paste(hsize,age))]
  expect_true(nrow(unique(eusilc.est$Estimates,by=c("rb090","db040","hsize","age")))==ngroups+1)
  
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=list("rb090","db040",c("rb090","db040"),c("hsize","age"))),NA)
})


myfun <- function(x,w){
  return(sum(w*x))
}
myfun.char <- function(x,w){
  return(as.character(sum(w*x)))
}
myfun.mulval <- function(x,w){
  return(w*x)
}
test_that("test para -  var and fun",{
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60s",
                            group=c("rb090","db040")),"Not all elements in var are column names in dat")
  
  eusilc[sample(.N,100),povmd60NA:=NA]
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60NA",
                            group=c("rb090","db040")),NA)
  
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            fun="myfun.undefined",group=c("rb090","db040")),"Function myfun.undefined is undefined")
  
  
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60NA",
                            fun="myfun",group=c("rb090","db040")),NA)
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            fun="myfun",group=c("rb090","db040")),NA)

  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            fun="myfun.char",group=c("rb090","db040")),"Function myfun.char does not return integer or numeric value")

  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            fun="myfun.mulval",group=c("rb090","db040")),"Function myfun.mulval does return more than one value. Only functions which return a single value are allowed.")

})

test_that("test para - period.diff, period.mean",{
 
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),period.mean=NULL),NA)
  expect_warning(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),period.mean=4),"period.mean must be odd - mean over periods will not be calculated")
  
  diff.warning <- capture_warnings(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),period.diff = "2015-2008"))
  expect_true(diff.warning[1]=="Removing 2015-2008 from period.diff - period(s) not present in column year\n")
  expect_true(diff.warning[2]=="No differences will be calculated\n")
  
  expect_warning(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                              group=c("rb090","db040"),period.diff = c("2015-2008","2016-2011")),"Removing 2015-2008 from period.diff - period(s) not present in column year",fixed=TRUE)
  
  expect_warning(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                              group=c("rb090","db040"),period.diff = c("2015-2010","2016-2011")),"Cannot calculate differences between periods 2015 and 2010 over 3 periods.")
  
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),period.diff = c("2015-2011","2016-2012"),period.mean=3),NA)
})

test_that("test para - bias, size.limit, cv.limit, p",{
  
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),bias="FALSE"),"bias can only be TRUE of FALSE")
  eusilc.bias <- calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                               group=c("rb090","db040"),bias=TRUE)
  expect_true("mean_povmd60"%in%colnames(eusilc.bias$Estimates))
  
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),size.limit="10"),"size.limit must contain one numeric value")
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),size.limit=1:10),"size.limit must have length 1")
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),size.limit=50),NA)
  
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),cv.limit=1:10),"cv.limit must have length 1")
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),cv.limit="1"),"cv.limit must contain one numeric value")
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),cv.limit=20),NA)
  
  
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),p=".5"),"p must be a numeric vector")
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),p=c(.1,.7,1.2)),"Values in p must be between 0 and 1")
  expect_error(calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                            group=c("rb090","db040"),p=c(.1,.7,.9)),NA)
})

test_that("test return",{
  eusilc.est <- calc.stError(eusilc,weights="db090",b.weights=paste0("w",1:10),period="year",var="povmd60",
                             group=c("rb090","db040"))
  eusilc.comp <- rbindlist(list(eusilc[,.(V1=surveysd:::weightedRatioC(povmd60,db090),N_true=sum(db090)),by=year],
                                eusilc[,.(V1=surveysd:::weightedRatioC(povmd60,db090),N_true=sum(db090)),by=list(year,rb090)],
                                eusilc[,.(V1=surveysd:::weightedRatioC(povmd60,db090),N_true=sum(db090)),by=list(year,db040)]),use.names=TRUE,fill=TRUE)
  eusilc.comp[,year:=as.character(year)]
  eusilc.comp <- merge(eusilc.comp, eusilc.est$Estimates[,.(year,rb090,db040,N,val_povmd60)])
  expect_true(nrow(eusilc.comp[V1!=val_povmd60])==0)
  expect_true(nrow(eusilc.comp[N_true!=N])==0)
})



