#'CalculateCoverage
#'@description
#'We calculate the coverage information per variant from the MAEGATK results.
# #'@import MatrixGenerics
#'@importFrom SummarizedExperiment rowRanges assays
#'@param SE SummarizedExperiment object.
#'@param chromosome_prefix List of matrices for the alternative reads.
#'@export
CalculateCoverage <- function(SE, chromosome_prefix = "chrM"){
  ref_allele <- as.character(SummarizedExperiment::rowRanges(SE)$refAllele)
  coverage <- SummarizedExperiment::assays(SE)[["coverage"]]
  rownames(coverage) <- paste0(chromosome_prefix, "_", 1:nrow(coverage), "_", ref_allele, "_A")
  coverage_A <- coverage[ref_allele != "A",]

  rownames(coverage) <- paste0(chromosome_prefix, "_", 1:nrow(coverage), "_", ref_allele, "_C")
  coverage_C <- coverage[ref_allele != "C",]

  rownames(coverage) <- paste0(chromosome_prefix, "_", 1:nrow(coverage), "_", ref_allele, "_G")
  coverage_G <- coverage[ref_allele != "G",]

  rownames(coverage) <- paste0(chromosome_prefix, "_", 1:nrow(coverage), "_", ref_allele, "_T")
  coverage_T <- coverage[ref_allele != "T",]

  coverage <- rbind(coverage_A, coverage_C, coverage_G, coverage_T)
  return(coverage)
}
