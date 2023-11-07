#'get_consensus
#'@description
#'We get the consensus information for a specific matrix.
#'I want to remove some packages if they are not needed. See below which package apperantly wasn't needed.
#'@importFrom methods as
#'@param alt_base The alternative base.
#'@param ref_base The reference base.
#'@param input_matrix Input matrix with the present reads numerically encoded.
#'@param chromosome_prefix The chromosome name used as a prefix.
#'@export
get_consensus <- function(alt_base, ref_base, input_matrix, chromosome_prefix = "chrM"){
  base_numeric <- c("A" = 8, "C" = 4, "G" = 2, "T" = 1)
  letter_numeric <- base_numeric[alt_base]
  ref_numeric <- base_numeric[ref_base]
  #mixture_numeric <- letter_numeric + base_numeric
  #other_values <- 1:sum(base_numeric)
  #other_values <- other_values[!other_values %in% c(base_numeric[c(alt_base, ref_base)], mixture_numeric)]
  # It is possible, that a variant only has reads for a different alternative allele.
  # So, there are only C reads for the position 3107_N.
  # That means, that 3107_N>C has a value of 2 (Alt), while 3107_N>A would have a value of 3 (Both).
  # Both is not accurate in this context. Therefore, we set these cases to 0 (NoCall). 
  other_homo_values <- base_numeric[!base_numeric %in% base_numeric[c(alt_base, ref_base)]]

  output_matrix <- matrix(0, nrow = nrow(input_matrix), ncol = ncol(input_matrix))
  rownames(output_matrix) <- rownames(input_matrix)
  colnames(output_matrix) <- colnames(input_matrix)
  output_matrix[input_matrix == ref_numeric] <- 1
  output_matrix[input_matrix == letter_numeric] <- 2

  bit_and_test <- input_matrix
  bit_and_test <- matrix(bitwAnd(letter_numeric, bit_and_test), nrow = nrow(input_matrix), ncol = ncol(input_matrix))
  rownames(bit_and_test) <- rownames(input_matrix)
  colnames(bit_and_test) <- colnames(bit_and_test)
  bit_or_test <- input_matrix
  bit_or_test <- matrix(bitwOr(letter_numeric, bit_or_test), nrow = nrow(input_matrix), ncol = ncol(input_matrix))
  rownames(bit_or_test) <- rownames(input_matrix)
  colnames(bit_or_test) <- colnames(bit_and_test)
  output_matrix[(bit_and_test > 0) & (bit_or_test > letter_numeric)] <- 3

  #output_matrix[input_matrix == mixture_numeric] <- 3
  #output_matrix[input_matrix %in% other_values] <- 0
  output_matrix[input_matrix %in% other_homo_values] <- 0
  
  rownames(output_matrix) <- paste0(chromosome_prefix, "_", gsub("[^[:digit:]., ]", "", rownames(output_matrix)), "_", ref_base, "_", alt_base)
  #output_matrix <- methods::as(output_matrix, "dgCMatrix")
  output_matrix <- methods::as(output_matrix, "CsparseMatrix")
  return(output_matrix)
}
