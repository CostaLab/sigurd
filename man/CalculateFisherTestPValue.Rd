% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalculateFisherTestPValue.R
\name{CalculateFisherTestPValue}
\alias{CalculateFisherTestPValue}
\title{CalculateFisherTestPValue}
\usage{
CalculateFisherTestPValue(
  variant_values,
  other_mutation,
  all_variants_list,
  min_intersecting_cells = 5
)
}
\arguments{
\item{variant_values}{The fraction values you are analysing. A vector.}

\item{other_mutation}{All other variants you have. A vector of variant names.}

\item{all_variants_list}{List of fraction values for all the variants you want to compare your variant with.}

\item{min_intersecting_cells}{Minimum number of intersecting cells. Correlations with less than this will not be performed.}
}
\description{
We perform the Fisher Test of SNVs and calculate the P values.
}
