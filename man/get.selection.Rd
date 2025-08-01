% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.selection.R
\name{get.selection}
\alias{get.selection}
\title{Get sample selection (~deltas) from drawn bootstrap replicates}
\usage{
get.selection(
  dat,
  b.rep = attr(dat, "b.rep"),
  strata = attr(dat, "strata"),
  cluster = attr(dat, "cluster"),
  hid = attr(dat, "hid"),
  pid = attr(dat, "pid")
)
}
\arguments{
\item{dat}{either data.frame or data.table containing the survey data with
rotating panel design. Should contain only survey data from a single time period.}

\item{b.rep}{character specifying the names of the columns in \code{dat}
containing bootstrap replicates.}

\item{strata}{character vector specifying the name(s) of the column in \code{dat}
by which the population was stratified.}

\item{cluster}{character vector specifying cluster in the data.}

\item{hid}{character specifying the name of the column in \code{dat} containing
the household id. If \code{NULL} (the default), the household structure is not
regarded. \code{hid} and \code{pid} cannot both be \code{NULL}.}

\item{pid}{pid column in \code{dat} specifying the personal identifier. This
identifier needs to be unique for each person throught the whole data set.
\code{hid} and \code{pid} cannot both be \code{NULL}.}
}
\value{
Returns a list of data.tables.
The length of the list equals the number of sampling stages specified.
Each list entry contains a \code{data.table} with variables for sampling stage and/or
\code{hid}/\code{pid} as well as \code{length(attr(dat,"b.rep"))} columns each indicating if
record/cluster was drawn in the respective sampling stage for the i-th boostrap replicate.
}
\description{
Reconstruct sample selection, e.g. record was drawn or not drawn (delta = 0/1)
in each sampling stage from bootstrap replicates.
\code{get.selection()} needs the \code{cluster}, \code{strata} and \code{hid}/\code{pid} information (if not \code{NULL})
to correctly reconstruct if a record was drawn in each sampling stage for each bootstrap replicate.
Is only needed if bootstrap replicates are drawn for a survey
with existing bootstrap replicates from a previous period,
see parameter \code{already.selected} in function \code{\link[=draw.bootstrap]{draw.bootstrap()}}.
}
\examples{

library(surveysd)
library(data.table)
setDTthreads(1)
set.seed(1234)
eusilc <- demo.eusilc(n = 3, prettyNames = TRUE)

## draw replicates with stratification
dat_boot <- draw.bootstrap(eusilc[year<2012], REP = 3, weights = "pWeight",
                           strata = "region", hid = "hid",
                           period = "year")

## get selection matrix for year 2011 
dat_selection <- get.selection(dat_boot[year==2011])
print(dat_selection)

## draw bootstrap replicates for year 2012
## respecting already selected units for year 2011 ~ dat_selection
## in order to mimic rotating panel design
dat_boot_2012 <- draw.bootstrap(eusilc[year==2012], REP = 3, weights = "pWeight",
                                strata = "region", hid = "hid",
                                period = "year", 
                                already.selected = dat_selection)



}
