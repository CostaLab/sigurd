#'CalculateRefReads
#'@description
#'We calculate the number of reference reads covering a variant using forward and reverse reads.
#'@importFrom SummarizedExperiment SummarizedExperiment assays
#'@importFrom Matrix Matrix
#'@param SE SummarizedExperiment object.
#'@param chromosome_prefix List of matrices for the reference reads.
#'@export
CalculateRefReads <- function(SE, chromosome_prefix = "chrM"){
  ref_allele <- as.character(SummarizedExperiment::rowRanges(SE)$refAllele)
  reads_A <- SummarizedExperiment::assays(SE)[["A_counts_fw"]] + SummarizedExperiment::assays(SE)[["A_counts_rev"]]
  reads_A <- reads_A[ref_allele == "A",]
  reads_A <- rbind(reads_A, reads_A, reads_A)
  rownames(reads_A) <- paste0(chromosome_prefix, "_", rep(which(ref_allele == "A"), 3), "_A_", rep(c("C", "G", "T"), each = sum(ref_allele == "A")))

  reads_C <- SummarizedExperiment::assays(SE)[["C_counts_fw"]] + SummarizedExperiment::assays(SE)[["C_counts_rev"]]
  reads_C <- reads_C[ref_allele == "C",]
  reads_C <- rbind(reads_C, reads_C, reads_C)
  rownames(reads_C) <- paste0(chromosome_prefix, "_", rep(which(ref_allele == "C"), 3), "_C_", rep(c("A", "G", "T"), each = sum(ref_allele == "C")))

  reads_G <- SummarizedExperiment::assays(SE)[["G_counts_fw"]] + SummarizedExperiment::assays(SE)[["G_counts_rev"]]
  reads_G <- reads_G[ref_allele == "G",]
  reads_G <- rbind(reads_G, reads_G, reads_G)
  rownames(reads_G) <- paste0(chromosome_prefix, "_", rep(which(ref_allele == "G"), 3), "_G_", rep(c("A", "C", "T"), each = sum(ref_allele == "G")))

  reads_T <- SummarizedExperiment::assays(SE)[["T_counts_fw"]] + SummarizedExperiment::assays(SE)[["T_counts_rev"]]
  reads_T <- reads_T[ref_allele == "T",]
  reads_T <- rbind(reads_T, reads_T, reads_T)
  rownames(reads_T) <- paste0(chromosome_prefix, "_", rep(which(ref_allele == "T"), 3), "_T_", rep(c("A", "C", "G"), each = sum(ref_allele == "T")))

  reads_ref <- rbind(reads_A, reads_C, reads_G, reads_T)

  N_positions <- which(ref_allele == "N")
  if(length(N_positions)){
    reads_N <- Matrix::Matrix(data = 0, nrow = length(N_positions) * 4, ncol = ncol(SE), sparse = TRUE, dimnames = list(paste0(chromosome_prefix, "_", rep(N_positions, each = 4), "_N_", c("A", "C", "G", "T")), colnames(SE)))
    reads_ref <- rbind(reads_ref, reads_N)
  }
  return(reads_ref)
}
