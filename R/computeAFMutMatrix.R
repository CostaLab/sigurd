#'Calculate the allele frequency per variant.
#'This function originally written by 
#'@import SummarizedExperiment
#'@param SE SummarizedExperiment object.
#'@export
computeAFMutMatrix <- function(SE, chromosome_prefix = "chrM"){
  cov <- assays(SE)[["coverage"]] + 0.000001
  ref_allele <- as.character(rowRanges(SE)$refAllele)

  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    mat[is.na(mat)] <- 0
    rownames(mat) <- paste0(chromosome_prefix, "_", 1:nrow(mat), "_", toupper(ref_allele), "_", letter)
    return(mat[toupper(ref_allele) != letter,])
  }

  A_matrix <- as.matrix(getMutMatrix("A"))
  C_matrix <- as.matrix(getMutMatrix("C"))
  G_matrix <- as.matrix(getMutMatrix("G"))
  T_matrix <- as.matrix(getMutMatrix("T"))
  result <- rbind(A_matrix, C_matrix, G_matrix, T_matrix)
  return(result)
}
