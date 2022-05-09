#'We calculate the consensus information from the MAEGATK results.
#'@import dplyr MatrixGenerics SummarizedExperiment
#'@param SE SummarizedExperiment object.
#'@param chromosome_prefix The chromosome name used as a prefix.
#'@export
CalculateConsensus <- function(SE, chromosome_prefix = "chrM"){
  # 0 NoCall      = coverage is 0.
  # 1 Reference   = only reference reads.
  # 2 Alternative = only alternative reads of one variant.
  # 3 Both        = reads for reference and one or more variants.
  
  # We get the read information per position.
  letter <- c("A", "C", "G", "T")
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  reads <- lapply(letter, getReadMatrix, SE = SE, chromosome_prefix = chromosome_prefix)
  # Since we have always the same 4 bases, we get all possible combinations by assigning numeric values.
  # A = 8, C = 4, G = 2, T = 1.
  # If there are several types of reads present at a position, we can simply add the values.
  # So, a position with A and T would have a value of 9.
  reads[[1]][reads[[1]] > 0] <- 8
  reads[[2]][reads[[1]] > 0] <- 4
  reads[[3]][reads[[1]] > 0] <- 2
  reads[[4]][reads[[1]] > 0] <- 1
  # We add the values together.
  # The row names are the names from the first matrix and not accurate any more.
  # The only relevant parts are the position and the reference base.
  variants_matrix <- reads[[1]] + reads[[2]] + reads[[3]] + reads[[4]]
  rm(reads)
  gc()
  
  # We get the position according to their reference base.
  # Now, we have a list for each set of position with the same base reference.
  variants_matrix_ls <- list(A = variants_matrix[grep("_A_", rownames(variants_matrix)),],
                             C = variants_matrix[grep("_C_", rownames(variants_matrix)),],
                             G = variants_matrix[grep("_G_", rownames(variants_matrix)),],
                             T = variants_matrix[grep("_T_", rownames(variants_matrix)),])
  rm(variants_matrix)
  gc()
  
  # Now, we check the consensus value for all positions with the same reference base.
  # Then we can rbind these matrices again and return one large consensus matrix in the end.
  consensus_a <- lapply(c("C", "G", "T"), get_consensus, ref_base = "A", input_matrix = as.matrix(variants_matrix_ls[[1]]))
  consensus_a <- do.call("rbind", consensus_a)
  consensus_c <- lapply(c("A", "G", "T"), get_consensus, ref_base = "C", input_matrix = as.matrix(variants_matrix_ls[[2]]))
  consensus_c <- do.call("rbind", consensus_c)
  consensus_g <- lapply(c("A", "C", "T"), get_consensus, ref_base = "G", input_matrix = as.matrix(variants_matrix_ls[[3]]))
  consensus_g <- do.call("rbind", consensus_g)
  consensus_t <- lapply(c("A", "C", "G"), get_consensus, ref_base = "T", input_matrix = as.matrix(variants_matrix_ls[[4]]))
  consensus_t <- do.call("rbind", consensus_t)
  consensus <- rbind(consensus_a, consensus_c, consensus_g, consensus_t)
  return(consensus)
}
