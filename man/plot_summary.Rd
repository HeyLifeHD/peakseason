% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_summary.R
\name{plot_summary}
\alias{plot_summary}
\title{Groups samples based on signal. Useful for QC analysis}
\usage{
plot_summary(bw_summary = NULL, log_transform = TRUE, pos_const = 1,
  condition = NULL)
}
\arguments{
\item{bw_summary}{output from \code{\link{extract_summary}}}

\item{condition}{a vector of condition for each bigWig file. Samples belonging to same condition are treated as replicates. Default NULL.}

\item{fract}{Fraction of data used for PCA or clustering. Default 0.25}

\item{pcs}{Principal components to draw.}
}
\description{
Groups samples based on signal. Useful for QC analysis
}
