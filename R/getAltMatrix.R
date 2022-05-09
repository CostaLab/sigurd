#'We get the alt values from the MAEGATK results.
#'@import SummarizedExperiment
#'@param SE_object SummarizedExperiment object.
#'@param letter The base you want to use. Character.
#'@param ref_allele The reference alleles.
#'@param chromosome_prefix The chromosome prefix used.
#'@export
getAltMatrix <- function(SE_object, letter, chromosome_prefix = "chrM"){
  ref_allele <- as.character(rowRanges(SE_object)$refAllele)
  mat <- (assays(SE_object)[[paste0(letter, "_counts_fw")]] + assays(SE_object)[[paste0(letter, "_counts_rev")]])
  rownames(mat) <- paste0(chromosome_prefix, "_", as.character(1:dim(mat)[1]), "_", toupper(ref_allele), ">", letter)
  mat <- mat[toupper(ref_allele) != letter,]
  return(mat)
}
