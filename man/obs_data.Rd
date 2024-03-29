% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{obs_data}
\alias{obs_data}
\title{Example observed data with continuous DV}
\format{
A data.table with 600 rows and 7 variables:
\describe{
  \item{ID}{Subject identifier}
  \item{TIME}{Time}
  \item{DV}{Concentration of drug}
  \item{AMT}{Amount of dosage initially administered at DV = 0, TIME = 0}
  \item{DOSE}{Dosage amount}
  \item{MDV}{Dummy indicating missing dependent variable value}
  \item{NTIME}{Nominal Time}
  \item{GENDER}{Character variable indicating subject's gender ("M", "F")}
  \item{STUDY}{Character variable indicating study type ("Study A", "Study B")}
}
}
\source{
\code{\link[vpc]{simple_data}}
}
\usage{
obs_data
}
\description{
An observed dataset from a hypothetical PK model, altered to include NTIME, GROUP, GENDER.
}
\keyword{datasets}
