#'Supplementing scRNAseq values with Amplicon values using big.matrix
#'
#'@description
#'We replace the values from an scRNAseq experiment with values we have from an amplicon experiment.
#'
#'@import archive Matrix SummarizedExperiment VariantAnnotation bigmemory
#'@param scRNAseq The SummarizedExperiment object containing the scRNAseq data. 
#'@param amplicon The SummarizedExperiment object containing the amplicon data.
#'@export
AmpliconSupplementing_big <- function(scRNAseq, amplicon){
  # We supplement the scRNAseq data with the amplicon data.
  print("We get the new meta data.")
  new_meta_data <- merge(colData(scRNAseq), colData(amplicon), by = "Cell", all.x = TRUE, all.y = TRUE,
                         suffixes = c("scRNAseq", "Amplicon"))

  print("We get all cells and variants.")
  all_cells <- unique(c(colnames(scRNAseq), colnames(amplicon)))
  all_variants <- unique(c(rownames(scRNAseq), rownames(amplicon)))

  print("We generate our output matrices.")
  options(bigmemory.allow.dimnames=TRUE)
  consensus <- big.matrix(init = 0, ncol = length(all_cells), nrow = length(all_variants))
  rownames(consensus) <- all_variants
  colnames(consensus) <- all_cells
  fraction <- consensus
  reads <- consensus

  print("We fill the output matrices.")
  consensus[rownames(scRNAseq), colnames(scRNAseq)] <- assays(scRNAseq)$consensus[,]
  fraction[rownames(scRNAseq), colnames(scRNAseq)] <- assays(scRNAseq)$fraction[,]
  reads[rownames(scRNAseq), colnames(scRNAseq)] <- assays(scRNAseq)$coverage[,]

  print("We add the the amplicon information.")
  consensus[rownames(amplicon), colnames(amplicon)] <- assays(amplicon)$consensus[,]
  fraction[rownames(amplicon), colnames(amplicon)] <- assays(amplicon)$fraction[,]
  reads[rownames(amplicon), colnames(amplicon)] <- assays(amplicon)$coverage[,]

  #print("We add the the amplicon information.")
  #assays(scRNAseq)[["consensus"]][rownames(amplicon), colnames(amplicon)] <- as.matrix(assays(amplicon)$consensus)
  #assays(scRNAseq)[["fraction"]][rownames(amplicon), colnames(amplicon)] <- as.matrix(assays(amplicon)$fraction)
  #assays(scRNAseq)[["coverage"]][rownames(amplicon), colnames(amplicon)] <- as.matrix(assays(amplicon)$coverage)

  se <- SummarizedExperiment(assays = list(consensus = consensus, fraction = fraction, coverage = reads),
                             colData = new_meta_data)
  return(se)
}
