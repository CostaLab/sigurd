% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VariantWiseCorrelation.R
\name{VariantWiseCorrelation}
\alias{VariantWiseCorrelation}
\title{VariantWiseCorrelation}
\usage{
VariantWiseCorrelation(
  variants_list,
  n_cores = 1,
  p_value_adjustment = "fdr",
  verbose = TRUE
)
}
\arguments{
\item{variants_list}{List of fraction values.}

\item{n_cores}{Number of cores you want to use. Numeric.}

\item{p_value_adjustment}{Method for P value adjustment. See p.adjust for details.}

\item{verbose}{Should the function be verbose? Default = TRUE}
}
\description{
We correlate the variants with each other using the Pearson correlation.
This function calls CalculateCorrelationPValue to perform the actual correlation.
}
