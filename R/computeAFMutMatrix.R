#'computeAFMutMatrix
#'@description
#'Calculate the allele frequency per variant.
#'Source: https://github.com/petervangalen/MAESTER-2021
#'We can get AF values greater than 1, which is due to uninformative reads.
#'See: https://gatk.broadinstitute.org/hc/en-us/articles/360035532252-Allele-Depth-AD-is-lower-than-expected
#'and https://github.com/caleblareau/mgatk/issues/1
#'We simply set these values to 1, since that is the actual information we have in this case.
#'@importFrom SummarizedExperiment assays rowRanges
#'@param SE SummarizedExperiment object.
#'@param chromosome_prefix The prefix of the chromosome.
#'@export
computeAFMutMatrix <- function(SE, chromosome_prefix = "chrM"){
  cov <- SummarizedExperiment::assays(SE)[["coverage"]] + 0.000001
  ref_allele <- as.character(SummarizedExperiment::rowRanges(SE)$refAllele)

  A_matrix <- getMutMatrix(SE = SE, cov = cov, letter = "A", ref_allele = ref_allele, chromosome_prefix = chromosome_prefix)
  gc()
  C_matrix <- getMutMatrix(SE = SE, cov = cov, letter = "C", ref_allele = ref_allele, chromosome_prefix = chromosome_prefix)
  gc()
  G_matrix <- getMutMatrix(SE = SE, cov = cov, letter = "G", ref_allele = ref_allele, chromosome_prefix = chromosome_prefix)
  gc()
  T_matrix <- getMutMatrix(SE = SE, cov = cov, letter = "T", ref_allele = ref_allele, chromosome_prefix = chromosome_prefix)
  gc()
  result <- rbind(A_matrix, C_matrix, G_matrix, T_matrix)
  return(result)
}
