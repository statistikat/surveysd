#################################
# test calc.stError()
#
#Sys.unsetenv("R_TESTS")

context("calc.stError()")
library(surveysd)
library(laeken)
library(data.table)

eusilc <- surveysd:::demo.eusilc()
eusilc <- eusilc[!db040 %in% c("Vienna", "Lower Austria", "Upper Austria")]

eusilc <- draw.bootstrap(eusilc, REP = 2, hid = "db030", weights = "db090",
                         period = "year", strata = "db040")
eusilc <- recalib(eusilc, conP.var = "rb090", conH.var = "db040")

# test input parameter
test_that("test para - data", {
  expect_error(calc.stError(as.matrix(eusilc),
                            var = "povmd60", group = c("rb090", "db040")),
               "dat must be a data.frame or data.table")
  expect_error(
    calc.stError(eusilc, var = "povmd60", group = c("rb090", "db040")),
    NA
  )
})

test_that("test para - weights, b.weights, year and group", {
  expect_error(
    calc.stError(
      eusilc, b.weights = "a", var = "povmd60", group = c("rb090", "db040")),
    "Not all elements in b.rep are column names in dat"
  )
  expect_error(
    calc.stError(eusilc, b.weights = 1:2, var = "povmd60",
                 group = c("rb090", "db040")),
    "Not all elements in b.rep are column names in dat"
  )

  expect_error(
    calc.stError(
      eusilc, weights = "db090s", var = "povmd60", group = c("rb090", "db040")),
    "db090s is not a column in dat")
  expect_error(
    calc.stError(
      eusilc, period = "years", var = "povmd60", group = c("rb090", "db040")),
    "years is not a column in dat")

  expect_error(
    calc.stError(
      eusilc, var = "povmd60", group = c("rb090s", "db040")),
    "Not all elements on group are column names in dat")

  eusilc.est <- calc.stError(
    eusilc, var = "povmd60", group = list("rb090", "db040", c("rb090", "db040"),
                                          c("hsize", "age")))
  ngroups <- eusilc[, uniqueN(rb090) + uniqueN(db040) +
                      uniqueN(paste(rb090, db040)) +
                      uniqueN(paste(hsize, age))]
  expect_true(
    nrow(
      unique(eusilc.est$Estimates, by = c("rb090", "db040", "hsize", "age"))
    ) == ngroups + 1)

  expect_error(calc.stError(
    eusilc, weights = "db090", b.weights = paste0("w", 1:2), period = "year",
    var = "povmd60", group = list("rb090", "db040", c("rb090", "db040"),
                                  c("hsize", "age"))), NA)
})


test_that("test para - period.diff, period.mean", {

  expect_error(calc.stError(
    eusilc, var = "povmd60", group = c("rb090", "db040")), NA)
  expect_warning(
    calc.stError(
      eusilc, var = "povmd60", group = c("rb090", "db040"), period.mean = 4),
    "period.mean must be odd - mean over periods will not be calculated")

  diff.warning <- capture_warnings(calc.stError(
    eusilc, var = "povmd60", group = c("rb090", "db040"),
    period.diff = "2015-2008"
  ))
  expect_true(diff.warning[1] == paste0(
    "Removing 2015-2008 from period.diff -",
    " period(s) not present in column year\n"))
  expect_true(diff.warning[2] == "No differences will be calculated\n")

  expect_warning(
    calc.stError(
      eusilc, var = "povmd60", group = c("rb090", "db040"),
      period.diff = c("2015-2008", "2016-2011")),
    paste0("Removing 2015-2008 from period.diff - period(s) not",
           " present in column year"),
    fixed = TRUE)

  expect_warning(
    calc.stError(
      eusilc, var = "povmd60", group = c("rb090", "db040"),
      period.diff = c("2015-2010", "2016-2011"), period.mean = 3),
    "Cannot calculate differences between periods 2015 and 2010 over 3 periods."
  )

  expect_error(calc.stError(
    eusilc, var = "povmd60", group = c("rb090", "db040"),
    period.diff = c("2015-2011", "2016-2012"), period.mean = 3), NA)
})

test_that("test para - bias, size.limit, cv.limit, p", {

  expect_error(
    calc.stError(eusilc, var = "povmd60", group = c("rb090", "db040"),
                 bias = "FALSE"),
    "bias can only be TRUE of FALSE")
  eusilc.bias <- calc.stError(
    eusilc, weights = "db090", b.weights = paste0("w", 1:2), period = "year",
    var = "povmd60", group = c("rb090", "db040"), bias = TRUE)
  expect_true("mean_povmd60" %in% colnames(eusilc.bias$Estimates))

  expect_error(
    calc.stError(
      eusilc, var = "povmd60", group = c("rb090", "db040"), size.limit = "10"),
    "size.limit must contain one numeric value")
  expect_error(
    calc.stError(
      eusilc, var = "povmd60", group = c("rb090", "db040"), size.limit = 1:2),
    "size.limit must have length 1")
  expect_error(calc.stError(
    eusilc, var = "povmd60", group = c("rb090", "db040"), size.limit = 50), NA)

  expect_error(
    calc.stError(
      eusilc, var = "povmd60", group = c("rb090", "db040"), cv.limit = 1:10),
    "cv.limit must have length 1")
  expect_error(
    calc.stError(
      eusilc, var = "povmd60", group = c("rb090", "db040"), cv.limit = "1"),
    "cv.limit must contain one numeric value")
  expect_error(calc.stError(
    eusilc, var = "povmd60", group = c("rb090", "db040"), cv.limit = 20), NA)


  expect_error(
    calc.stError(
      eusilc, var = "povmd60", group = c("rb090", "db040"), p = ".5"),
    "p must be a numeric vector")
  expect_error(
    calc.stError(
      eusilc, var = "povmd60", group = c("rb090", "db040"), p = c(.1, .7, 1.2)),
    "Values in p must be between 0 and 1")
  expect_error(calc.stError(
    eusilc, var = "povmd60", group = c("rb090", "db040"), p = c(.1, .7, .9)),
    NA)
})

test_that("test return", {
  eusilc.est <- calc.stError(
    eusilc, var = "povmd60", group = c("rb090", "db040"))
  eusilc.comp <- rbindlist(
    list(
      eusilc[, .(V1 = weightedRatio(povmd60, db090), N_true = sum(db090)),
             by = year],
      eusilc[, .(V1 = weightedRatio(povmd60, db090), N_true = sum(db090)),
             by = list(year, rb090)],
      eusilc[, .(V1 = weightedRatio(povmd60, db090), N_true = sum(db090)),
             by = list(year, db040)]),
    use.names = TRUE, fill = TRUE
  )
  eusilc.comp <- merge(eusilc.comp, eusilc.est$Estimates[, .(year, rb090, db040,
                                                             N, val_povmd60)])
  expect_true(nrow(eusilc.comp[V1 != val_povmd60]) == 0)
  expect_true(nrow(eusilc.comp[N_true != N]) == 0)
})
