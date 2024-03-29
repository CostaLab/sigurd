% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CombineSEobjects.R
\name{CombineSEobjects}
\alias{CombineSEobjects}
\title{CombineSEobjects}
\usage{
CombineSEobjects(se_1, se_2, suffixes = c("_1", "_2"))
}
\arguments{
\item{se_1}{SummarizedExperiment object 1.}

\item{se_2}{SummarizedExperiment object 2.}

\item{suffixes}{The suffixes you want to add to the meta data.frame.}
}
\description{
We combine two SummarizedExperiment objects. This is originally intended to add the mitochondrial genotyping information to the somatic genotyping information.
The function can combine any two objects and it is not requried to be from two different genotyping results.
}
