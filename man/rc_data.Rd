% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rc.r
\name{rc_data}
\alias{rc_data}
\title{Generate x-data for response curves}
\usage{
rc_data(x, type = c("distance", "all"), res = 200)
}
\arguments{
\item{x}{An \code{\link{mbm}} object; it is not necessary for it to be complete (i.e., fitted)}

\item{type}{The kind of response curve(s) to make}

\item{res}{How many points per response curve}
}
\value{
A list of data frames with x-values to be fit to response curves
}
\description{
Currently this is very simple, does not handle the 'all' case and it assumes covariates are 0-centered
}
