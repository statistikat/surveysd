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

test_that("ipf stops  as expected when dimensionality of eps does not fit", {
  # with array epsP1, base weights and bound
  expect_error(ipf(
    eusilc, hid = "hid", conP = list(conP1, conP2), conH = list(conH1,conH1,conH1),
    epsP = c(1e-06), epsH = list(epsH1,1e-06), w = "baseWeight", bound = 4,
    verbose = FALSE, maxIter = 200),
    regexp = "Provided household eps argument does not fit household constraints.")

  expect_error(ipf(
    eusilc, hid = "hid", conP = list(conP1, conP2), conH = list(conH1),
    epsP = c(1e-06), epsH = list(epsH1[,1:8]), w = "baseWeight", bound = 4,
    verbose = FALSE, maxIter = 200),
    regexp = "Provided household eps argument 1 does not fit in dimension to  household constraints 1 .")

    expect_error(ipf(
    eusilc, hid = "hid", conP = list(conP1, conP2), conH = list(conH1,conH1,conH1),
    epsP = c(1e-06,22), epsH = list(epsH1,1e-06,1e-6), w = "baseWeight", bound = 4,
    verbose = FALSE, maxIter = 200),
    regexp = "Individual eps arguments for each constraints must be defined as list.")

})


test_that("ipf works as expected with minMaxTrim only P", {
  # with array epsP1, base weights and bound
  calibweights2 <- ipf(
    eusilc, hid = "hid", conP = list(conP1, conP2),
    epsP = 1e-06, w = "baseWeight", bound = 4,
    verbose = FALSE, maxIter = 200, minMaxTrim = c(450,700))
  err <- max(c(
    max(abs(xtabs(calibWeight ~ age, data = calibweights2) - conP1) /
          conP1),
    max(
      abs(xtabs(calibWeight ~ gender + region, data = calibweights2) - conP2) /
        conP2)))
  expect_true(err < .01)
  expect_true(all(calibweights2$calibWeight<=700)&&all(calibweights2$calibWeight>=450))
})

test_that("ipf works as expected with minMaxTrim P and HH", {
  # with array epsP1, base weights and bound
  calibweights2 <- ipf(
    eusilc, hid = "hid", conP = list(conP1, conP2), conH = list(conH1),
    epsP = 1e-06, epsH = list(epsH1), w = "baseWeight", bound = 4,
    verbose = FALSE, maxIter = 200, minMaxTrim = c(340,870))
  err <- max(c(
    max(abs(xtabs(calibWeight ~ age, data = calibweights2) - conP1) /
          conP1),
    max(
      abs(xtabs(calibWeight ~ gender + region, data = calibweights2) - conP2) /
        conP2),
    max(abs(xtabs(calibWeight ~ hsize + region, data = calibweights2,
                  subset = !duplicated(hid)) - conH1) / conH1)))
  expect_true(err < .01)
  expect_true(all(calibweights2$calibWeight<=870)&&all(calibweights2$calibWeight>=340))
})

test_that("summary.ipf", {
  # with array epsP1, base weights and bound
  calibweights2 <- ipf(
    eusilc, hid = "hid", conP = list(conP1, conP2), conH = list(conH1),
    epsP = 1e-06, epsH = list(epsH1), w = "baseWeight", bound = 4,
    verbose = FALSE, maxIter = 200, minMaxTrim = c(340,870))
  out <- summary(calibweights2)
  expect_length(out, 18)
  expect_length(out$`distribution of the weights`,6)
})


test_that("summary.ipf without base weight", {
  # with array epsP1, base weights and bound
  calibweights2 <- ipf(
    eusilc, hid = "hid", conP = list(conP1, conP2), conH = list(conH1),
    epsP = 1e-06, epsH = list(epsH1), bound = NULL)
  out <- summary(calibweights2)
  expect_length(out, 18)
  expect_length(out$`distribution of the weights`,0)
})
