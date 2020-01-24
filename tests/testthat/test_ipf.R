context("ipf")
library(surveysd)
library(data.table)
eusilc <- demo.eusilc(n = 1, prettyNames = TRUE)

# treat households as a factor variable
eusilc[, hid := as.factor(hid)]

## example for base weights assuming a simple random sample of households
## stratified per region
eusilc[, regSamp := .N, by = region]
eusilc[, regPop := sum(pWeight), by = region]
eusilc[, baseWeight := regPop / regSamp]
eusilc[, pIncome := eqIncome / .N, by = hid]

## constraints on person level
# age
conP1 <- xtabs(pWeight ~ age, data = eusilc)
# gender by region
conP2 <- xtabs(pWeight ~ gender + region, data = eusilc)
# personal net income by gender
conP3 <- xtabs(pWeight * pIncome ~ gender, data = eusilc)

## constraints on household level
conH1 <- xtabs(pWeight ~ hsize + region, data = eusilc,
               subset = !duplicated(hid))
## constraints on household level netIncome
conH2 <- xtabs(pWeight * eqIncome ~ region, data = eusilc,
               subset = !duplicated(hid))
# array of convergence limits for conH1
epsH1 <- conH1
epsH1[1:4, ] <- 0.005
epsH1[5, ] <- 0.2


test_that("ipf with a numerical variable works as expected - computeLinear", {
  # without array epsP1
  calibweights1 <- ipf(
    eusilc, hid = "hid",
    conP = list(conP1, conP2, pIncome = conP3),
    conH = list(conH1),
    epsP = list(1e-06, 1e-06, 1e-03),
    epsH = 0.01,
    bound = NULL, verbose = FALSE,  maxIter = 200,
    numericalWeighting = computeLinear)
  conP3_adj <- xtabs(calibWeight * pIncome ~ gender, data = calibweights1)
  expect_true(abs(sum(conP3_adj) - sum(conP3)) / sum(conP3) < .01)
  expect_true(all(abs(conP3_adj - conP3) / conP3 < .01))

  err <- max(c(
    max(abs(xtabs(calibWeight ~ age, data = calibweights1) -
              conP1) / conP1),
    max(abs(xtabs(calibWeight ~ gender + region, data = calibweights1) -
              conP2) / conP2),
    max(abs(xtabs(calibWeight ~ hsize + region, data = calibweights1,
                  subset = !duplicated(hid)) - conH1) / conH1)))
  expect_true(err < .01)
})

test_that("ipf with a numerical variable works as expected - computeLinearG1", {
  # without array epsP1
  calibweights1 <- ipf(
    eusilc, hid = "hid",
    conP = list(conP1, conP2, pIncome = conP3),
    conH = list(conH1),
    epsP = list(1e-06, 1e-06, 1e-03),
    epsH = 0.01,
    bound = NULL, verbose = FALSE, maxIter = 200,
    numericalWeighting = computeLinearG1)
  conP3_adj <- xtabs(calibWeight * pIncome ~ gender, data = calibweights1)
  expect_true(abs(sum(conP3_adj) - sum(conP3)) / sum(conP3) < .01)
  expect_true(all(abs(conP3_adj - conP3) / conP3 < .01))

  err <- max(c(
    max(abs(xtabs(calibWeight ~ age, data = calibweights1) -
              conP1) / conP1),
    max(abs(xtabs(calibWeight ~ gender + region, data = calibweights1) -
              conP2) / conP2),
    max(abs(xtabs(calibWeight ~ hsize + region, data = calibweights1,
                  subset = !duplicated(hid)) - conH1) / conH1)))
  expect_true(err < .01)
})


test_that("ipf with a numerical variable in households  as expected", {
  # without array epsP1
  calibweights1 <- ipf(
    eusilc, hid = "hid",
    conP = list(conP1, conP2),
    conH = list(conH1, eqIncome = conH2),
    epsP = list(1e-06, 1e-06),
    epsH = list(0.01, 0.01),
    bound = NULL, verbose = FALSE,  maxIter = 50,
    numericalWeighting = computeFrac)
  conP3_adj <- xtabs(calibWeight * pIncome ~ gender, data = calibweights1)
  expect_true(abs(sum(conP3_adj) - sum(conP3)) / sum(conP3) < .01)
  expect_true(all(abs(conP3_adj - conP3) / conP3 < .01))

  err <- max(c(
    max(abs(xtabs(calibWeight ~ age, data = calibweights1) -
              conP1) / conP1),
    max(abs(xtabs(calibWeight ~ gender + region, data = calibweights1) -
              conP2) / conP2),
    max(abs(xtabs(calibWeight ~ hsize + region, data = calibweights1,
                  subset = !duplicated(hid)) - conH1) / conH1)))
  expect_true(err < .01)
})

test_that("ipf works as expected", {
  # with array epsP1, base weights and bound
  calibweights2 <- ipf(
    eusilc, hid = "hid", conP = list(conP1, conP2), conH = list(conH1),
    epsP = 1e-06, epsH = list(epsH1), w = "baseWeight", bound = 4,
    verbose = FALSE, maxIter = 200)
  err <- max(c(
    max(abs(xtabs(calibWeight ~ age, data = calibweights2) - conP1) /
          conP1),
    max(
      abs(xtabs(calibWeight ~ gender + region, data = calibweights2) - conP2) /
        conP2),
    max(abs(xtabs(calibWeight ~ hsize + region, data = calibweights2,
                  subset = !duplicated(hid)) - conH1) / conH1)))
  expect_true(err < .01)
})

test_that("ipf works as expected calibWeight renamed", {
  # with array epsP1, base weights and bound
  calibweights2 <- ipf(
    eusilc, hid = "hid", conP = list(conP1, conP2), conH = list(conH1),
    epsP = 1e-06, epsH = list(epsH1), w = "baseWeight", bound = 4,
    verbose = FALSE, maxIter = 200, nameCalibWeight = "calibWeightNew")
  err <- max(c(
    max(abs(xtabs(calibWeightNew ~ age, data = calibweights2) - conP1) /
          conP1),
    max(abs(
      xtabs(calibWeightNew ~ gender + region, data = calibweights2) - conP2) /
          conP2),
    max(abs(xtabs(calibWeightNew ~ hsize + region, data = calibweights2,
                  subset = !duplicated(hid)) - conH1) / conH1)))
  expect_true(err < .01)
})



test_that("ipf errors work as expected", {
  
  expect_error(ipf(
    eusilc, hid = "hid", conP = list(conP1, conP2),
    w = "baseWeight9"),
    "Base weight baseWeight9 is not a column name in dat")
  
  expect_error(ipf(
    eusilc, hid = "hid", conP = list(conP1, eqincome = conP3),
    w = "baseWeight"),
    "Numerical constraints must be named by variables in dat")
  
  setNA <- eusilc[,sample(hid,1),by=.(hsize,region)]
  eusilc[.(hid=setNA$V1),eqIncome:=NA,on=.(hid)]
  
  expect_error(ipf(
    eusilc, hid = "hid", conP = list(conP1, eqIncome = conP3),
    w = "baseWeight"),
    "Numeric variable eqIncome contains missing values")
  
  setNA <- eusilc[,sample(hid,1),by=.(hsize,region)]
  eusilc[.(hid=setNA$V1),baseWeight:=NA,on=.(hid)]
  
  expect_error(ipf(
    eusilc, hid = "hid", conP = list(conP1, eqIncome = conP3),
    w = "baseWeight"),
    "Base weight baseWeight contains missing values")
})

