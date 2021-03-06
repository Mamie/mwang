% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/citrus-functions.R
\name{PlotRes}
\alias{PlotRes}
\title{Plot Results of CITRUS regression}
\usage{
PlotRes(clustering, cluster.attributes, outcomes, name, desc = NULL,
  ylab = "ES")
}
\arguments{
\item{clustering}{A list object returned by RunCITRUS.}

\item{cluster.attributes}{A object returned by PlotHierarchy.}

\item{name}{Name of selected channels.}

\item{desc}{Description of selected channels.}

\item{outcome}{A vector of binary outcome.}
}
\value{
A ggplot object of effect size.
}
\description{
Plot results of CITRUS regression including effect size of the
difference in abundance between groups, boxplot of abundance
between groups and distribution of differential clusters.
}
