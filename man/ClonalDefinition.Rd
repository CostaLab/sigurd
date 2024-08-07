% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClonalDefinition.R
\name{ClonalDefinition}
\alias{ClonalDefinition}
\title{Defining clones.}
\usage{
ClonalDefinition(
  se,
  variants_ls,
  grouping = NULL,
  identities = NULL,
  verbose = TRUE
)
}
\arguments{
\item{se}{SummarizedExperiment object.}

\item{variants_ls}{List of variants for clonal definition}

\item{grouping}{The meta data    column used to split the cells into groups. Default = NULL}

\item{identities}{Vector of groups, like samples.}

\item{verbose}{Should the function be verbose? Default = TRUE}
}
\description{
This function defines clones using a combination of variants supplied. Each set of variants can be and should be paired with a identity, like sample or cell type. The clones are then only determined in respect to this group. If this information is not provided, the clones are determined for all cells present.
}
