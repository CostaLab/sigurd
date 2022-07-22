#'CalculateCoverage
#'@description
#'We calculate the coverage information per variant from the MAEGATK results.
#'@import MatrixGenerics SummarizedExperiment
#'@param SE SummarizedExperiment object.
#'@param chromosome_prefix List of matrices for the alternative reads.
#'@export
CalculateCoverage <- function(SE, chromosome_prefix = "chrM"){
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  coverage <- assays(SE)[["coverage"]]
  rownames(coverage) <- paste0(chromosome_prefix, "_", 1:nrow(coverage), "_", ref_allele, "_A")
  coverage_A <- coverage[ref_allele != "A",]
  #coverage_A <- as.matrix(coverage_A)

  rownames(coverage) <- paste0(chromosome_prefix, "_", 1:nrow(coverage), "_", ref_allele, "_C")
  coverage_C <- coverage[ref_allele != "C",]
  #coverage_C <- as.matrix(coverage_C)

  rownames(coverage) <- paste0(chromosome_prefix, "_", 1:nrow(coverage), "_", ref_allele, "_G")
  coverage_G <- coverage[ref_allele != "G",]
  #coverage_G <- as.matrix(coverage_G)

  rownames(coverage) <- paste0(chromosome_prefix, "_", 1:nrow(coverage), "_", ref_allele, "_T")
  coverage_T <- coverage[ref_allele != "T",]
  #coverage_T <- as.matrix(coverage_T)

  coverage <- rbind(coverage_A, coverage_C, coverage_G, coverage_T)
  return(coverage)
}
