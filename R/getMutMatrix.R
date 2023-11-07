#'getMutMatrix
#'@description
#'This function gets the allele frequency for a specific allele. It is used in computeAFMutMatrix.
#'Source: https://github.com/petervangalen/MAESTER-2021
#'@importFrom SummarizedExperiment assays
#'@importFrom methods as
#'@param SE SummarizedExperiment object.
#'@param cov The coverage matrix from MAEGATK/MGATK.
#'@param letter The base we are interested in.
#'@param ref_allele Vector of reference alleles.
#'@param chromosome_prefix The chromosome prefix used.
#'@export
getMutMatrix <- function(SE, cov, letter, ref_allele, chromosome_prefix){
  names_rows <- paste0(chromosome_prefix, "_", 1:nrow(cov), "_", toupper(ref_allele), "_", letter)
  names_rows <- names_rows[toupper(ref_allele) != letter]
  mat_fow <- SummarizedExperiment::assays(SE)[[paste0(letter, "_counts_fw")]]
  mat_rev <- SummarizedExperiment::assays(SE)[[paste0(letter, "_counts_rev")]]
  mat <- mat_fow + mat_rev
  mat <- mat[toupper(ref_allele) != letter,]
  cov_use <- cov[toupper(ref_allele) != letter,]
  mat <- mat / cov_use
  gc()
  mat[mat > 1] <- 1
  rownames(mat) <- names_rows
  mat <- methods::as(mat, "CsparseMatrix")
  return(mat)
}
