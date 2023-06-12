#'Supplementing scRNAseq values with Amplicon values
#'
#'@description
#'We replace the values from an scRNAseq experiment with values we have from an amplicon experiment.
#'
#'@import archive Matrix SummarizedExperiment VariantAnnotation
#'@param scRNAseq The SummarizedExperiment object containing the scRNAseq data. 
#'@param amplicon The SummarizedExperiment object containing the amplicon data.
#'@export
AmpliconSupplementing <- function(scRNAseq, amplicon){
  # We supplement the scRNAseq data with the amplicon data.
  print("We get the new meta data.")
  new_meta_data <- merge(colData(scRNAseq), colData(amplicon), by = "Cell", all.x = TRUE, all.y = TRUE,
                         suffixes = c("scRNAseq", "Amplicon"))
  new_row_data <- merge(rowData(scRNAseq), rowData(amplicon), by = "VariantName", all.x = TRUE, all.y = TRUE,
                        suffixes = c("scRNAseq", "Amplicon"))

  print("We get all cells and variants.")
  all_cells <- unique(c(colnames(scRNAseq), colnames(amplicon)))
  all_variants <- unique(c(rownames(scRNAseq), rownames(amplicon)))

  print("We generate our output matrices.")
  consensus <- matrix(0, ncol = length(all_cells), nrow = length(all_variants))
  rownames(consensus) <- all_variants
  colnames(consensus) <- all_cells
  #consensus <- sparseMatrix(i = 1, j = 1, dims = c(length(all_variants), length(all_cells)), repr = "C")
  fraction <- consensus
  reads <- consensus
  alts <- consensus
  refs <- consensus

  print("We fill the output matrices.")
  consensus[rownames(scRNAseq), colnames(scRNAseq)] <- as.matrix(assays(scRNAseq)$consensus)
  fraction[rownames(scRNAseq), colnames(scRNAseq)] <- as.matrix(assays(scRNAseq)$fraction)
  reads[rownames(scRNAseq), colnames(scRNAseq)] <- as.matrix(assays(scRNAseq)$coverage)
  alts[rownames(scRNAseq), colnames(scRNAseq)] <- as.matrix(assays(scRNAseq)$alts)
  refs[rownames(scRNAseq), colnames(scRNAseq)] <- as.matrix(assays(scRNAseq)$refs)

  print("We add the the amplicon information.")
  consensus[rownames(amplicon), colnames(amplicon)] <- as.matrix(assays(amplicon)$consensus)
  fraction[rownames(amplicon), colnames(amplicon)] <- as.matrix(assays(amplicon)$fraction)
  reads[rownames(amplicon), colnames(amplicon)] <- as.matrix(assays(amplicon)$coverage)
  alts[rownames(amplicon), colnames(amplicon)] <- as.matrix(assays(amplicon)$alts)
  refs[rownames(amplicon), colnames(amplicon)] <- as.matrix(assays(amplicon)$refs)

  #print("We add the the amplicon information.")
  #assays(scRNAseq)[["consensus"]][rownames(amplicon), colnames(amplicon)] <- as.matrix(assays(amplicon)$consensus)
  #assays(scRNAseq)[["fraction"]][rownames(amplicon), colnames(amplicon)] <- as.matrix(assays(amplicon)$fraction)
  #assays(scRNAseq)[["coverage"]][rownames(amplicon), colnames(amplicon)] <- as.matrix(assays(amplicon)$coverage)

  se <- SummarizedExperiment(assays = list(consensus = as(consensus, "dgCMatrix"), fraction = as(fraction, "dgCMatrix"), coverage = as(reads, "dgCMatrix"), alts = as(alts, "dgCMatrix"), refs = as(refs, "dgCMatrix")),
                             colData = new_meta_data, rowData = new_row_data)
  return(se)
}
