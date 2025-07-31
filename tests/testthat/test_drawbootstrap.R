#################################
# test draw.bootstrap()
#

context("draw.bootstrap()")
library(surveysd)
library(laeken)
library(data.table)

eusilc <- surveysd:::demo.eusilc(n = 3)
eusilc[, N.households := sum(db090[!duplicated(db030)]), by = .(year, db040)]
eusilc[!duplicated(db030), N.households.error := sum(db090),
       by = .(year, db040)]
eusilc[, N.households.all := sum(db090[!duplicated(db030)]), by = .(year)]

# test input parameter
test_that("test para - data", {
  expect_error(
    draw.bootstrap(
      as.matrix(eusilc), REP = 2, hid = "db030", weights = "db090", period =
        "year", strata = "db040"),
    "dat must be a data.frame or data.table")
  expect_error(draw.bootstrap(
    eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
    strata = "db040"), NA)
})

test_that("test para - REP", {
  expect_error(
    draw.bootstrap(
      eusilc, REP = "a", hid = "db030", weights = "db090", period = "year",
      strata = "db040"),
    "'REP' must be of type numeric")
  expect_error(
    draw.bootstrap(
      eusilc, REP = 1:2, hid = "db030", weights = "db090", period = "year",
      strata = "db040"),
    "'REP' must have length 1")
})

test_that("test para - hid, weights and period", {
  expect_error(
    draw.bootstrap(
      eusilc, REP = 2, hid = "db030s", weights = "db090", period = "year",
      strata = "db040"),
    "column(s) 'db030s' specified in 'hid' not found in dat", fixed = TRUE)
  expect_error(
    draw.bootstrap(
      eusilc, REP = 2, hid = "db030", weights = "db090s", period = "year",
      strata = "db040"),
    "column(s) 'db090s' specified in 'weights' not found in dat", fixed = TRUE)
  expect_error(
    draw.bootstrap(
      eusilc, REP = 2, hid = "db030", weights = "db090", period = "years",
      strata = "db040"),
    "column(s) 'years' specified in 'period' not found in dat", fixed = TRUE)

  eusilc[, year.char := as.character(year)]
  expect_error(
    draw.bootstrap(
      eusilc, REP = 2, hid = "db030", weights = "db090", period = "year.char",
      strata = "db040"),
    "column(s) 'year.char' in parameter 'period' must correspond to numeric column(s) in dat", fixed = TRUE)
})

test_that("test para - strata, cluster and totals", {
  expect_error(
    draw.bootstrap(
      eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
      strata = "db0400"),
    "Not all elements in strata are column names in dat")
  expect_error(
    draw.bootstrap(
      eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
      strata = "db040", cluster = "a"),
    "Not all names in cluster are column names in dat")
  expect_error(
    draw.bootstrap(
      eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
      strata = c("db040", "rb090"), cluster = c("db040", "db030"),
      totals = c("N.households")),
    "totals must be specified for each stage")
  expect_error(
    draw.bootstrap(
      eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
      strata = c("db040"), cluster = c("db040", "db030")),
    paste0("strata and cluster need to have the same number of stages!\n ",
           "Please use either '1' or 'I' if there was no clustering or ",
           "stratification in one of the stages."),
    fixed = TRUE)
  expect_error(
    draw.bootstrap(
      eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
      strata = c("db040"), cluster = c("db030"), totals = c("N.something")),
    "Not all names in totals are column names in dat")
  expect_error(
    draw.bootstrap(
      eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
      strata = c("db040"), cluster = c("db030"),
      totals = c("N.households.error")),
    "Missing values found in column(s): N.households.error", fixed = TRUE)

  expect_error(draw.bootstrap(
    eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
    strata = c("db040"), cluster = c("1")), NA)
  expect_error(draw.bootstrap(
    eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
    strata = c("db040"), cluster = c("I")), NA)
  expect_error(draw.bootstrap(
    eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
    strata = NULL), NA)
  expect_error(draw.bootstrap(
    eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
    strata = "I"), NA)
  expect_error(draw.bootstrap(
    eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
    strata = "1"), NA)
})


# these are some longer tests
if (Sys.getenv("SURVEYSD_ADDITIONAL_TEST") == "TRUE") {
  test_that("test para - bootnames, split and pid", {

    expect_error(
      draw.bootstrap(
        eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
        strata = "db040", split = "FALSE"),
      "split needs to be logical")
    expect_error(
      draw.bootstrap(
        eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
        strata = "db040", split = TRUE),
      "when split is TRUE pid must be specified")

    eusilc[, rb030error := rb030[1], by = list(year, db030)]
    expect_error(
      draw.bootstrap(
        eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
        strata = "db040", split = TRUE, pid = "rb030error"),
      "pid is not unique in each household for each period")

    # create split households
    eusilc[, rb030split := rb030]
    year <- eusilc[, unique(year)]
    year <- year[-1]
    leaf_out <- c()
    for (y in year) {
      split.person <- eusilc[
        year == (y - 1) & !duplicated(db030) & !db030 %in% leaf_out,
        sample(rb030, 20)]
      overwrite.person <- eusilc[
        year == (y) & !duplicated(db030) & !db030 %in% leaf_out,
        .(rb030 = sample(rb030, 20))]
      overwrite.person[, c("rb030split", "year_curr") := .(split.person, y)]

      eusilc[overwrite.person, rb030split := i.rb030split,
             on = .(rb030, year >= year_curr)]
      leaf_out <- c(leaf_out, eusilc[rb030 %in% c(overwrite.person$rb030,
                                                  overwrite.person$rb030split),
                                     unique(db030)])
    }
    eusilc.split <- draw.bootstrap(
      eusilc, REP = 10, hid = "db030", weights = "db090", period = "year",
      strata = "db040", split = TRUE, pid = "rb030split")
    eusilc.split <- eusilc.split[, lapply(.SD, uniqueN), by = rb030split,
                                 .SDcols = paste0("w", 1:10)]
    expect_true(all(eusilc.split[, .SD, .SDcols = paste0("w", 1:10)] == 1))

    expect_error(
      draw.bootstrap(
        eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
        strata = "db040", boot.names = "1"),
      "boot.names must start with an alphabetic character")
    expect_error(draw.bootstrap(
      eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
      strata = "db040", boot.names = "weight"), NA)
  })
}


test_that("test para - single.PSU", {
  expect_warning(
    draw.bootstrap(
      eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
      strata = "db040", single.PSU = "something"),
    paste0("single.PSU was not set to either 'merge' or 'mean'!\n Bootstrap ",
           "replicates for single PSUs cases will be missing!"),
    fixed = TRUE)
})

test_that("test return", {
  dat.boot <- draw.bootstrap(
    eusilc, REP = 2, hid = "db030", weights = "db090", period = "year",
    strata = "db040")
  expect_true(ncol(dat.boot) == (2 + ncol(eusilc)))
  dat.unique <- unique(dat.boot[, mget(c("db030", paste0("w", 1:2)))])
  expect_true(nrow(dat.unique[, .N, by = db030][N > 2]) == 0)
  
  # check if any value is infinite
  expect_false(any(unlist(
    dat.boot[, lapply(
      .SD,
      function(z) {
        any(is.infinite(z))
      }),
      .SDcols = c(paste0("w", 1:2))
      ])))
  
  # check if sum of weights equals
  # numer of sampling units per strata
  tab_sums <- dat.boot[, lapply(
    .SD,
    function(w,db030) {
      s1 <- sum(w[!duplicated(db030)])
      s2 <- uniqueN(db030)
      abs(s1-s2)>1e-10
    },db030=db030),by=.(year,db040),
    .SDcols = c(paste0("w", 1:2))
  ]
  expect_false(any(unlist(tab_sums[,.SD,.SDcols=paste0("w", 1:2)])))
  
  expect_false(any(is.na(dat.boot[, .SD, .SDcols = c(paste0("w", 1:2))])))
  expect_true(all(dat.boot[, .SD, .SDcols = c(paste0("w", 1:2))] > 0))
})




