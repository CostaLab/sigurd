#'Get the counts for a specific base over all positions.
#'@importFrom SummarizedExperiment assays rowRanges
#'@param SE SummarizedExperiment object.
#'@param letter The base for which we want the counts.
#'@param chromosome_prefix The chromosome name used as a prefix.
#'@export
getReadMatrix <- function(SE, letter, chromosome_prefix = "chrM"){
  ref_allele <- as.character(SummarizedExperiment::rowRanges(SE)$refAllele)
  mat <- (SummarizedExperiment::assays(SE)[[paste0(letter, "_counts_fw")]] + SummarizedExperiment::assays(SE)[[paste0(letter, "_counts_rev")]])
  rownames(mat) <- paste0(chromosome_prefix, "_", 1:nrow(mat), "_", toupper(ref_allele), "_", letter)
  return(mat)
}
