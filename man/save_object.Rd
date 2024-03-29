% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_object.R
\name{save_object}
\alias{save_object}
\title{save_object}
\usage{
save_object(object, file_name, file_format = "zstd")
}
\arguments{
\item{object}{The R object to be save.}

\item{file_name}{The path were the file shall be save.}

\item{file_format}{The format of the save file. Has to be one of: zstd, lz4, gzip, bzip2, xz, nocomp.}
}
\description{
Saving function to save the RDS files quicker.
Source:https://github.com/CostaLab/CimpleG
}
