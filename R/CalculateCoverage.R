#'We calculate the consensus information from the MAEGATK results.
#'@import MatrixGenerics SummarizedExperiment
#'@param reference_reads Reference reads matrix.
#'@param alternative_reads List of matrices for the alternative reads.
#'@export
CalculateCoverage <- function(reference_reads, alternative_reads){
  # We remove the potential N at position 3107 of the human genome.
  alternative_reads <- alternative_reads[!grepl("_N", rownames(alternative_reads)),]
  # We get the first part of the ref row name. This includes the position and the ref allele.
  rows_ref_reads <- gsub(">.*", "", rownames(reference_reads))
  # We get the first part of the alt row name. This includes the position and the ref allele.
  rows_alt_reads <- gsub(">.*", "", rownames(alternative_reads))
  # We can then subset the ref reads matrix to only include the positions from the alt matrix.
  keep <- rows_ref_reads %in% rows_alt_reads
  reference_reads <- reference_reads[keep,]
  rows_ref_reads <- rows_ref_reads[keep]
  # We match the new ref reads matrices to be the same order as the alt read matrices.
  reference_reads <- reference_reads[match(rows_alt_reads, rows_ref_reads),]
  # We divide the alt matrix by alt + ref matrix.
  coverage <- as.matrix(alternative_reads + reference_reads)
  rownames(coverage) <- gsub(">", "_", rownames(coverage))
  return(coverage)
}
