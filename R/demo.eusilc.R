#' Generate multiple years of EU-SILC data
#'
#' Create a dummy dataset to be used for demonstrating the functionalities of
#' the `surveysd` package based on [laeken::eusilc]. Please refer to the
#' documentation page of the original data for details about the variables.
#'
#' @param n Number of years to generate. Should be at least 1
#' @param prettyNames Create easy-to-read names for certain variables.
#'   Recommended for demonstration purposes. Otherwise, use the original codes
#'   documented in [laeken::eusilc].
#'
#' @details
#' If `prettyNames` is `TRUE`, the following variables will be available in an
#' easy-to-read manner.
#'
#' * `hid` Household id. Consistent with respect to the reference period
#'   (`year`)
#' * `hsize` Size of the household. derived from `hid` and `period`
#' * `region` Federal state of austria where the household is located
#' * `pid` Personal id. Consistent with respect to the reference period (`year`)
#' * `age` Age-class of the respondent
#' * `gender` A persons gender (`"male"`, `"Female"`)
#' * `ecoStat` Ecnomic status
#'   (`"part time"`, `"full time"`, `"unemployed"`, ...)
#' * `citizenship` Citizenship (`"AT"`, `"EU"`, `"other"`)
#' * `pWeight` Personal sample weight inside the reference period
#' * `year`. Simulated reference period
#' * `povertyRisk`. Logical variable determining whether a respondent is at risk
#'   of poverty
#'
#' @importFrom dplyr recode
#' @examples
#' demo.eusilc(n = 1, prettyNames = TRUE)[, c(1:8, 26, 28:30)]
#' @export
demo.eusilc <- function(n = 8, prettyNames = FALSE) {

  db030 <- rb030 <- povmd60 <- eqincome <- db090 <- eqIncome <- age <- hsize <-
    . <- povertyRisk <- ecoStat <- NULL

  data("eusilc", package = "laeken", envir = environment())
  setDT(eusilc)
  # generate yearly data for y years
  # 25% drop out from 1 year to the other
  eusilc[, year := 2010]
  eusilc.i <- copy(eusilc)
  nsamp <- round(eusilc[, uniqueN(db030)] * .25)
  hhincome <- eusilc[!duplicated(db030)][["eqIncome"]]
  nextIDs <- (1:nsamp) + eusilc[, max(db030)]
  if (n > 1)
    for (i in 1:(n - 1)) {
      eusilc.i[db030 %in% sample(unique(eusilc.i$db030), nsamp),
               c("db030", "eqIncome") := .(nextIDs[.GRP], sample(hhincome, 1L)),
               by = db030]
      eusilc.i[, year := year + 1]
      eusilc <- rbind(eusilc, eusilc.i)
      nextIDs <- (1:nsamp) + eusilc[, max(db030)]
    }

  eusilc[, rb030 := as.integer(paste0(db030, "0", 1:.N)),
         by = list(year, db030)]
  eusilc[, povmd60 := as.numeric(
    eqIncome < .6 * laeken::weightedMedian(eqIncome, w = db090)
  ), by = year]
  eusilc[, age := cut(age, c(-Inf, 16, 25, 45, 65, Inf))]
  eusilc[, hsize := cut(hsize, c(0:5, Inf))]

  if (prettyNames) {
    data.table::setnames(eusilc, "db030", "hid")
    data.table::setnames(eusilc, "db040", "region")
    data.table::setnames(eusilc, "rb030", "pid")
    data.table::setnames(eusilc, "rb090", "gender")
    data.table::setnames(eusilc, "pb220a", "citizenship")
    data.table::setnames(eusilc, "pl030", "ecoStat")
    eusilc[, ecoStat := dplyr::recode(
      ecoStat, `1` = "full time", `2` = "part time", `3` = "unemployed",
      `4` = "education", `5` = "retired", `6` = "disabled", `7` = "domestic")]
    data.table::setnames(eusilc, "rb050", "pWeight")
    data.table::setnames(eusilc, "povmd60", "povertyRisk")
    eusilc[, povertyRisk := as.logical(povertyRisk)]
  }

  return(eusilc)
}
