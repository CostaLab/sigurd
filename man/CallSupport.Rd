% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CallSupport.R
\name{CallSupport}
\alias{CallSupport}
\title{CallSupport}
\usage{
CallSupport(
  SE,
  VOI_group1,
  VOI_group2,
  group1_name = "group1",
  group2_name = "group2",
  min_mutated_reads = 3,
  min_reads = 30,
  group_factor = NULL,
  verbose = TRUE,
  return_nonsupport = FALSE
)
}
\arguments{
\item{SE}{SummarizedExperiment object.}

\item{VOI_group1}{The variants supporting the first group.}

\item{VOI_group2}{The variants supporting the second group.}

\item{group1_name}{The name used for the first group.}

\item{group2_name}{The name used for the second group.}

\item{min_mutated_reads}{The minimum number of mutated reads in a cell supporting a group.}

\item{min_reads}{Minimum number of reads per cell for a classification.}

\item{group_factor}{How much higher has the mean allele frequency to be in group 1 when compared to group 2 and vice versa? Can be a vector of length 2.}

\item{verbose}{Should the function be verbose? Default = TRUE}

\item{return_nonsupport}{Should nonsupporting cells be return too? Default FALSE.}

\item{group_of_interest}{The column data that divides the cells.}
}
\description{
Check if a cell is supported by a set of variants.
}
