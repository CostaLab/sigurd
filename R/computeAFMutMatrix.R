#'Calculate the allele frequency per variant.
#'This function originally written by 
#'@import SummarizedExperiment
#'@param SE SummarizedExperiment object.
#'@export
computeAFMutMatrix <- function(SE, chromosome_prefix = "chrM"){
  cov <- assays(SE)[["coverage"]] #+ 0.000001
  ref_allele <- as.character(rowRanges(SE)$refAllele)

  getMutMatrix <- function(letter){
    names_rows <- paste0(chromosome_prefix, "_", 1:nrow(cov), "_", toupper(ref_allele), "_", letter)
    names_rows <- names_rows[toupper(ref_allele) != letter]
    mat_fow <- assays(SE)[[paste0(letter, "_counts_fw")]]
    mat_rev <- assays(SE)[[paste0(letter, "_counts_rev")]]
    mat <- mat_fow + mat_rev
    mat <- mat[toupper(ref_allele) != letter,]
    cov_use <- cov[toupper(ref_allele) != letter,]
    mat <- mat / cov_use
    gc()
    mat[is.na(mat)] <- 0
    rownames(mat) <- names_rows
    mat <- as(mat, "dgCMatrix")
    return(mat)
  }
  A_matrix <- getMutMatrix("A")
  gc()
  C_matrix <- getMutMatrix("C")
  gc()
  G_matrix <- getMutMatrix("G")
  gc()
  T_matrix <- getMutMatrix("T")
  gc()
  result <- rbind(A_matrix, C_matrix, G_matrix, T_matrix)
  return(result)
}
