% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/citrus-functions.R
\name{PlotHierarchy}
\alias{PlotHierarchy}
\title{Plot CITRUS hierarchy}
\usage{
PlotHierarchy(clustering, regressionRes, lambda = "min", seed = 1)
}
\arguments{
\item{clustering}{A list object returned by RunCITRUS.}

\item{regressionRes}{An object returned by CITRUSRegression.}

\item{lambda}{A character to indicate lambda choice (min or 1se)}

\item{seed}{An optional seed for reproducibility.}
}
\value{
A ggplot object of the hierarchy.
}
\description{
Plot the hierarchy used in CITRUS colored by stratefying clusters.
}
