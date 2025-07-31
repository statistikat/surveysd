#################################
# test recalib()
#

context("recalib()")
library(surveysd)
library(laeken)
library(data.table)

eusilc <- surveysd:::demo.eusilc(n = 4)
eusilc <- draw.bootstrap(eusilc, REP = 2, hid = "db030", weights = "db090",
                         period = "year", strata = "db040")

# test input parameter
test_that("test para - data", {
  expect_error(
    recalib(as.matrix(eusilc), conP.var = "rb090", conH.var = "db040"),
    "dat must be a data.frame or data.table")
  expect_error(recalib(
    eusilc, conP.var = "rb090", conH.var = "db040"), NA)
  expect_error(recalib(
    eusilc, conP.var = c("rb090", "age"), conH.var = c("db040", "hsize")), NA)
})

test_that("test para - REP", {
  expect_error(
    recalib(eusilc, b.rep = "a", conP.var = "rb090", conH.var = "db040"),
    "column(s) 'a' specified in 'b.rep' not found in dat", fixed = TRUE)
  expect_error(
    recalib(eusilc, b.rep = 1:2, conP.var = "rb090", conH.var = "db040"),
    "'b.rep' must be of type character")
})

test_that("test para - hid, weights and period", {
  expect_error(
    recalib(eusilc, hid = "db030s", conP.var = "rb090", conH.var = "db040"),
    "column(s) 'db030s' specified in 'hid' not found in dat", fixed = TRUE)
  expect_error(
    recalib(eusilc, weights = "db090s", conP.var = "rb090", conH.var = "db040"),
    "column(s) 'db090s' specified in 'weights' not found in dat", fixed = TRUE)
  expect_error(
    recalib(eusilc, period = "years", conP.var = "rb090", conH.var = "db040"),
    "column(s) 'years' specified in 'years' not found in dat", fixed = TRUE)
})

test_that("test para - conP.var conH.var", {

  # check single values
  expect_error(
    recalib(eusilc, conP.var = "rb090s", conH.var = "db040"),
    "Not all elements in conP.var are column names in dat")
  expect_error(
    recalib(eusilc, conP.var = "rb090", conH.var = "db040s"),
    "Not all elements in conH.var are column names in dat")

  # check multiple values through list
  expect_error(
    recalib(eusilc, conP.var = list("rb090s", "age"), conH.var = "db040"),
    "Not all elements in conP.var are column names in dat")
  expect_error(
    recalib(eusilc, conP.var = list("rb090"),
            conH.var = list("db040s", "hsize")),
    "Not all elements in conH.var are column names in dat")

  expect_error(recalib(
    eusilc, conP.var = NULL, conH.var = "db040"), NA)
  expect_error(recalib(
    eusilc, conP.var = "rb090", conH.var = NULL), NA)
  expect_error(recalib(
      eusilc, conP.var = NULL, conH.var = NULL), NA)
  expect_error(recalib(
    eusilc, conP.var = list("rb090", "age"), conH.var = "db040"), NA)
  expect_error(recalib(
    eusilc, conP.var = list("rb090", c("age", "rb090")), conH.var = "db040"),
    NA)
  expect_error(recalib(
    eusilc, conP.var = list("rb090"), conH.var = list("db040", "hsize")), NA)
})


test_that("test para - conP conH", {

  conP1 <- xtabs(db090 ~ age + year, data = eusilc)
  conP2 <- xtabs(db090 ~ rb090 + year, data = eusilc)
  # conP2 <- xtabs(db090 ~ rb090 + db040 + year, data = eusilc)

  conH1 <- xtabs(db090 ~ hsize + year,
                 data = eusilc[!duplicated(paste(db030, year))])
  conH2 <- xtabs(db090 ~ db040 + year,
                 data = eusilc[!duplicated(paste(db030, year))])

  expect_error(
    recalib(eusilc, conP.var = "rb090", conH.var = NULL, conP = list(conP2)),
    paste("contingency table for rb090 was supplied through parameter",
          "conP AND conP.var"))
  expect_error(
    recalib(eusilc, conP.var = NULL, conH.var = "db040", conH = list(conH2)),
    paste("contingency table for db040 was supplied through parameter",
          "conH AND conH.var"))

  expect_error(recalib(
    eusilc, conP.var = NULL, conH.var = "db040", conH = list(conH1)), NA)
  expect_error(recalib(
    eusilc, conP.var = "rb090", conH.var = NULL, conP = list(conP1)), NA)
  expect_error(recalib(
    eusilc, conP.var = NULL, conH.var = NULL, conP = list(conP1),
    conH = list(conH1)), NA)
})

test_that("test return", {
  dat.calib <- recalib(
    eusilc, hid = "db030", weights = "db090", b.rep = paste0("w", 1:2), period =
      "year", conP.var = c("rb090", "age"), conH.var = c("db040", "hsize"))
  rb090.compare <- dat.calib[, lapply(.SD, sum), by = list(year, rb090),
                             .SDcols = c("db090", paste0("w", 1:2))]
  expect_true(all(!is.na(rb090.compare[, .SD, .SDcols = paste0("w", 1:2)])))
  rb090.compare <- rb090.compare[, lapply(
    .SD,
    function(z) {
      abs(db090 - z) / db090
    }
    ), .SDcols = paste0("w", 1:2)]
  expect_true(all(rb090.compare < formals(recalib)$epsP))

  db040.compare <- dat.calib[, lapply(
    .SD,
    function(z) {
      sum(z[!duplicated(db030)])
    }
    ), by = list(year, db040), .SDcols = c("db090", paste0("w", 1:2))]
  db040.compare <- db040.compare[, lapply(
    .SD,
    function(z) {
      abs(db090 - z) / db090
    }
  ), .SDcols = paste0("w", 1:2)]
  expect_true(all(db040.compare < formals(recalib)$epsH))
})
