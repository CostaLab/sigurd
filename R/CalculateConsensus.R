#'CalculateConsensus
#'@description
#'We calculate the consensus information from the MAEGATK results.
#'We set cells that have only alternative reads to 2 (Alternative).
#'We set cells that have only reference reads to 1 (Reference).
#'We set cells that have a mixture of alternative and reference reads to 3 (Both).
#'We set cells that have no reads to 0 (NoCall).
#'
#'Please note. Cells can have reads for the reference of a specific variant and no reads for the alternative.
#'The cell can still have a reads for the other alternative alleles. Then the cell is still considered as 0 (NoCall) for this variant.
#'For example:
#'A cell has at position 3: 0 A reads, 53 T reads, 63 C reads, 148 T reads.
#'For the variant chrM_3_T_A, the cell would have 53 reference reads, but also reads for other variants at this position.
#'To make sure that there is no confusion, the cell is set to NoCall.
# #'@import MatrixGenerics
#'@importFrom SummarizedExperiment rowRanges 
#'@param SE SummarizedExperiment object.
#'@param chromosome_prefix The chromosome name used as a prefix.
#'@param verbose Should the function be verbose? Default = FALSE
#'@export
CalculateConsensus <- function(SE, chromosome_prefix = "chrM", verbose = FALSE){
  # 0 NoCall      = coverage is 0.
  # 1 Reference   = only reference reads.
  # 2 Alternative = only alternative reads of one variant.
  # 3 Both        = reads for reference and one or more variants.

  if(verbose) print("We get the read information per position.")
  letter <- c("A", "C", "G", "T")
  ref_allele <- as.character(SummarizedExperiment::rowRanges(SE)$refAllele)
  reads <- lapply(letter, getReadMatrix, SE = SE, chromosome_prefix = chromosome_prefix)
  # Since we have always the same 4 bases, we get all possible combinations by assigning numeric values.
  # A = 8, C = 4, G = 2, T = 1.
  # If there are several types of reads present at a position, we can simply add the values.
  # So, a position with A and T would have a value of 9.
  reads[[1]][reads[[1]] > 0] <- 8
  reads[[2]][reads[[2]] > 0] <- 4
  reads[[3]][reads[[3]] > 0] <- 2
  reads[[4]][reads[[4]] > 0] <- 1
  if(verbose) print("We add the values together.")
  # The row names are the names from the first matrix and not accurate any more.
  # The only relevant parts are the position and the reference base.
  variants_matrix <- reads[[1]] + reads[[2]] + reads[[3]] + reads[[4]]
  rm(reads)
  gc()

  if(verbose) print("We get the position according to their reference base.")
  # Now, we have a list for each set of position with the same base reference.
  variants_matrix_ls <- list(A = variants_matrix[grep("_A_", rownames(variants_matrix), value = TRUE), , drop = FALSE],
                             C = variants_matrix[grep("_C_", rownames(variants_matrix), value = TRUE), , drop = FALSE],
                             G = variants_matrix[grep("_G_", rownames(variants_matrix), value = TRUE), , drop = FALSE],
                             T = variants_matrix[grep("_T_", rownames(variants_matrix), value = TRUE), , drop = FALSE],
                             N = variants_matrix[grep("_N_", rownames(variants_matrix), value = TRUE), , drop = FALSE])
  # We check if the N reference is even used. If the variants_matrix_ls[["N"]] is empty (zero rows), we do not perform the consensus determination.
  n_binding <- FALSE
  if(nrow(variants_matrix_ls[["N"]]) > 0){
    n_binding <- TRUE
    variants_matrix_ls[["N"]] <- matrix(variants_matrix_ls[["N"]], nrow = 1, ncol = ncol(variants_matrix))
    colnames(variants_matrix_ls[["N"]]) <- colnames(variants_matrix)
    rownames(variants_matrix_ls[["N"]]) <- paste0(chromosome_prefix, "_3107_N_A")
  }
  rm(variants_matrix)
  gc()

  if(verbose) print("Now, we check the consensus value for all positions with the same reference base.")
  # Then we can rbind these matrices again and return one large consensus matrix in the end.
  if(verbose) print("A")
  consensus_a <- lapply(c("C", "G", "T"), get_consensus, ref_base = "A", input_matrix = as.matrix(variants_matrix_ls[[1]]), chromosome_prefix = chromosome_prefix)
  consensus_a <- do.call("rbind", consensus_a)
  if(verbose) print("C")
  consensus_c <- lapply(c("A", "G", "T"), get_consensus, ref_base = "C", input_matrix = as.matrix(variants_matrix_ls[[2]]), chromosome_prefix = chromosome_prefix)
  consensus_c <- do.call("rbind", consensus_c)
  if(verbose) print("G")
  consensus_g <- lapply(c("A", "C", "T"), get_consensus, ref_base = "G", input_matrix = as.matrix(variants_matrix_ls[[3]]), chromosome_prefix = chromosome_prefix)
  consensus_g <- do.call("rbind", consensus_g)
  if(verbose) print("T")
  consensus_t <- lapply(c("A", "C", "G"), get_consensus, ref_base = "T", input_matrix = as.matrix(variants_matrix_ls[[4]]), chromosome_prefix = chromosome_prefix)
  consensus_t <- do.call("rbind", consensus_t)
  if(n_binding){
    if(verbose) print("N")
    consensus_n <- lapply(c("A", "C", "G", "T"), get_consensus, ref_base = "N", input_matrix = variants_matrix_ls[[5]], chromosome_prefix = chromosome_prefix)
    consensus_n <- do.call("rbind", consensus_n)
  } else{
    if(verbose) print("N reference not present.")
  }
  if(verbose) print("Binding the matrices.")
  consensus <- rbind(consensus_a, consensus_c, consensus_g, consensus_t)
  if(n_binding) consensus <- rbind(consensus, consensus_n)
  return(consensus)
}
