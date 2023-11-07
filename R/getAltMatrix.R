#'getAltMatrix
#'@description
#'We get the alt values from the MAEGATK results.
#'Source: https://github.com/petervangalen/MAESTER-2021
#'@importFrom SummarizedExperiment assays rowRanges
#'@param SE_object SummarizedExperiment object.
#'@param letter The base you want to use. Character.
#'@param chromosome_prefix The chromosome prefix used.
#'@export
getAltMatrix <- function(SE_object, letter, chromosome_prefix = "chrM"){
  ref_allele <- as.character(SummarizedExperiment::rowRanges(SE_object)$refAllele)
  mat <- (SummarizedExperiment::assays(SE_object)[[paste0(letter, "_counts_fw")]] + SummarizedExperiment::assays(SE_object)[[paste0(letter, "_counts_rev")]])
  rownames(mat) <- paste0(chromosome_prefix, "_", as.character(1:dim(mat)[1]), "_", toupper(ref_allele), ">", letter)
  mat <- mat[toupper(ref_allele) != letter,]
  return(mat)
}
