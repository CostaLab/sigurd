#'Calculating the Minor Allele Frequency.
#'@description
#'We calculate the MAF from a reference reads matrix and an alternative reads matrix.
#'This function is intended to be used with the mitochondrial genome and not with other somatic mutations.
# #'@import MatrixGenerics SummarizedExperiment
#'@param reference_reads Reference reads matrix.
#'@param alternative_reads List of matrices for the alternative reads.
#'@param pseudo_count = What is the pseudo count you want to add to the reference_reads matrix. Default = 0
#'@export
CalculateAlleleFrequency <- function(reference_reads, alternative_reads, pseudo_count = 0){
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
  allelefrequency <- as.matrix(alternative_reads / (alternative_reads + reference_reads + pseudo_count))
  allelefrequency[is.na(allelefrequency)] <- 0
  rownames(allelefrequency) <- gsub(">", "_", rownames(allelefrequency))
  return(allelefrequency)
}
