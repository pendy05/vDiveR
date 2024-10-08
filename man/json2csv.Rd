% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/json2csv.R
\name{json2csv}
\alias{json2csv}
\title{JSON2CSV}
\usage{
json2csv(
  json_data,
  host_name = "unknown host",
  protein_name = "unknown protein"
)
}
\arguments{
\item{json_data}{DiMA JSON output dataframe}

\item{host_name}{name of the host species}

\item{protein_name}{name of the protein}
}
\value{
A dataframe which acts as input for the other functions in vDiveR package
}
\description{
This function converts DiMA (v5.0.9) JSON output file to a dataframe with 17
predefined columns which further acts as the input for other functions provided in this vDiveR package.
}
\examples{
inputdf<-json2csv(JSON_sample)
}
