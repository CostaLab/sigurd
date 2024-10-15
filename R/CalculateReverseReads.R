#'CalculateReverseReads
#'@description
#'We calculate the number of reverse reads covering a variant.
#'@importFrom SummarizedExperiment SummarizedExperiment assays rowRanges
#'@param SE SummarizedExperiment object.
#'@param chromosome_prefix List of matrices for the alternative reads.
#'@export
CalculateReverseReads <- function(SE, chromosome_prefix = "chrM"){
  ref_allele <- as.character(SummarizedExperiment::rowRanges(SE)$refAllele)
  reads_A <- SummarizedExperiment::assays(SE)[["A_counts_rev"]]
  rownames(reads_A) <- paste0(chromosome_prefix, "_", 1:nrow(reads_A), "_", ref_allele, "_A")
  reads_A <- reads_A[ref_allele != "A",]

  reads_C <- SummarizedExperiment::assays(SE)[["C_counts_rev"]]
  rownames(reads_C) <- paste0(chromosome_prefix, "_", 1:nrow(reads_C), "_", ref_allele, "_C")
  reads_C <- reads_C[ref_allele != "C",]

  reads_G <- SummarizedExperiment::assays(SE)[["G_counts_rev"]]
  rownames(reads_G) <- paste0(chromosome_prefix, "_", 1:nrow(reads_G), "_", ref_allele, "_G")
  reads_G <- reads_G[ref_allele != "G",]

  reads_T <- SummarizedExperiment::assays(SE)[["T_counts_rev"]]
  rownames(reads_T) <- paste0(chromosome_prefix, "_", 1:nrow(reads_T), "_", ref_allele, "_T")
  reads_T <- reads_T[ref_allele != "T",]
  
  reads_rev <- rbind(reads_A, reads_C, reads_G, reads_T)
  return(reads_rev)
}
