% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ldsc_interface.R
\name{collect_ldsc}
\alias{collect_ldsc}
\title{Extract information from ldsc log file
Extracts relevant information from the .log output of LDSC. reads in all *.log
files found in dir}
\usage{
collect_ldsc(dir)
}
\arguments{
\item{dir}{filepath for directory}
}
\value{
a tibble
}
\description{
Extract information from ldsc log file
Extracts relevant information from the .log output of LDSC. reads in all *.log
files found in dir
}
\examples{
\dontrun{
collect_ldsc("/path/to/ldsc_logs")
}

}
