#'computeAFMutMatrix
#'@description
#'Calculate the allele frequency per variant.
#'Source: https://github.com/petervangalen/MAESTER-2021
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
    # We can get AF values greater than 1, which is due to uninformative reads.
    # See: https://gatk.broadinstitute.org/hc/en-us/articles/360035532252-Allele-Depth-AD-is-lower-than-expected
    # and https://github.com/caleblareau/mgatk/issues/1
    # We simply set these values to 1, since that is the actual information we have in this case.
    # This issue can be solved on the MAEGATK/GATK side.
    mat[mat > 1] <- 1
    rownames(mat) <- names_rows
    mat <- as(mat, "dgCMatrix")
    return(mat)
  }

  A_matrix <- getMutMatrix("A")
  #A_matrix <- as.matrix(A_matrix)
  gc()
  C_matrix <- getMutMatrix("C")
  #C_matrix <- as.matrix(C_matrix)
  gc()
  G_matrix <- getMutMatrix("G")
  #G_matrix <- as.matrix(G_matrix)
  gc()
  T_matrix <- getMutMatrix("T")
  #T_matrix <- as.matrix(T_matrix)
  gc()
  result <- rbind(A_matrix, C_matrix, G_matrix, T_matrix)
#  result <- as.matrix(result)
  return(result)
}
