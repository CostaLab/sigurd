#'Supplementing scRNAseq values with Amplicon values
#'
#'@description
#'We replace the values from an scRNAseq experiment with values we have from an amplicon experiment.
#'@import Matrix SummarizedExperiment
#'@param scRNAseq The SummarizedExperiment object containing the scRNAseq data. 
#'@param amplicon The SummarizedExperiment object containing the amplicon data.
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
AmpliconSupplementing <- function(scRNAseq, amplicon, verbose = TRUE){
  # We supplement the scRNAseq data with the amplicon data.
  if(verbose) print("We get the new meta data.")
  new_meta_data <- S4Vectors::merge(SummarizedExperiment::colData(scRNAseq), SummarizedExperiment::colData(amplicon), by = "Cell", all.x = TRUE, all.y = TRUE,
                                    suffixes = c("scRNAseq", "Amplicon"))
  rownames(new_meta_data) <- new_meta_data$Cell
  # We add an AverageCoverage column to the new meta data.
  new_meta_data$AverageCoverage                           <- new_meta_data$AverageCoveragescRNAseq
  amplicon_value                                          <- new_meta_data$AverageCoverageAmplicon
  names(amplicon_value)                                   <- colnames(amplicon)
  amplicon_value                                          <- na.omit(amplicon_value)
  new_meta_data[names(amplicon_value), "AverageCoverage"] <- amplicon_value

  new_row_data <- merge(SummarizedExperiment::rowData(scRNAseq), SummarizedExperiment::rowData(amplicon), by = "VariantName", all.x = TRUE, all.y = TRUE,
                        suffixes = c("scRNAseq", "Amplicon"))
  rownames(new_row_data) <- new_row_data$VariantName
  if("VariantQualityscRNAseq" %in% colnames(new_row_data)){
    # We add a VariantQuality column to the row data, showing the scRNAseq quality with the supplemented amplicon quality.
    # We do the same for the concordance and the depth.
    new_row_data$VariantQuality                           <- new_row_data$VariantQualityscRNAseq
    amplicon_value                                        <- SummarizedExperiment::rowData(amplicon)$VariantQuality
    names(amplicon_value)                                 <- rownames(amplicon)
    amplicon_value                                        <- na.omit(amplicon_value)
    new_row_data[names(amplicon_value), "VariantQuality"] <- amplicon_value
  }

  if("ConcordancescRNAseq" %in% colnames(new_row_data)){
    new_row_data$Concordance                           <- new_row_data$ConcordancescRNAseq
    amplicon_value                                     <- SummarizedExperiment::rowData(amplicon)$Concordance
    names(amplicon_value)                              <- rownames(amplicon)
    amplicon_value                                     <- na.omit(amplicon_value)
    new_row_data[names(amplicon_value), "Concordance"] <- amplicon_value
  }

  new_row_data$Depth                           <- new_row_data$DepthscRNAseq
  amplicon_value                               <- SummarizedExperiment::rowData(amplicon)$Depth
  names(amplicon_value)                        <- rownames(amplicon)
  amplicon_value                               <- na.omit(amplicon_value)
  new_row_data[names(amplicon_value), "Depth"] <- amplicon_value

  if(verbose) print("We get all cells and variants.")
  all_cells     <- unique(c(colnames(scRNAseq), colnames(amplicon)))
  all_variants  <- unique(c(rownames(scRNAseq), rownames(amplicon)))
  new_meta_data <- new_meta_data[all_cells,]
  new_row_data  <- new_row_data[all_variants,]

  if(verbose) print("We generate our output matrices.")
  consensus <- matrix(0, ncol = length(all_cells), nrow = length(all_variants))
  rownames(consensus) <- all_variants
  colnames(consensus) <- all_cells
  #consensus <- Matrix::sparseMatrix(i = 1, j = 1, dims = c(length(all_variants), length(all_cells)), repr = "C")
  fraction <- consensus
  reads <- consensus
  alts <- consensus
  refs <- consensus

  if(verbose) print("We fill the output matrices.")
  consensus[rownames(scRNAseq), colnames(scRNAseq)] <- as.matrix(SummarizedExperiment::assays(scRNAseq)$consensus)
  fraction[rownames(scRNAseq), colnames(scRNAseq)]  <- as.matrix(SummarizedExperiment::assays(scRNAseq)$fraction)
  reads[rownames(scRNAseq), colnames(scRNAseq)]     <- as.matrix(SummarizedExperiment::assays(scRNAseq)$coverage)
  alts[rownames(scRNAseq), colnames(scRNAseq)]      <- as.matrix(SummarizedExperiment::assays(scRNAseq)$alts)
  refs[rownames(scRNAseq), colnames(scRNAseq)]      <- as.matrix(SummarizedExperiment::assays(scRNAseq)$refs)

  if(verbose) print("We add the the amplicon information.")
  consensus[rownames(amplicon), colnames(amplicon)] <- as.matrix(SummarizedExperiment::assays(amplicon)$consensus)
  fraction[rownames(amplicon), colnames(amplicon)]  <- as.matrix(SummarizedExperiment::assays(amplicon)$fraction)
  reads[rownames(amplicon), colnames(amplicon)]     <- as.matrix(SummarizedExperiment::assays(amplicon)$coverage)
  alts[rownames(amplicon), colnames(amplicon)]      <- as.matrix(SummarizedExperiment::assays(amplicon)$alts)
  refs[rownames(amplicon), colnames(amplicon)]      <- as.matrix(SummarizedExperiment::assays(amplicon)$refs)

  #if(verbose) print("We add the the amplicon information.")
  #SummarizedExperiment::assays(scRNAseq)[["consensus"]][rownames(amplicon), colnames(amplicon)] <- as.matrix(SummarizedExperiment::assays(amplicon)$consensus)
  #SummarizedExperiment::assays(scRNAseq)[["fraction"]][rownames(amplicon), colnames(amplicon)]  <- as.matrix(SummarizedExperiment::assays(amplicon)$fraction)
  #SummarizedExperiment::assays(scRNAseq)[["coverage"]][rownames(amplicon), colnames(amplicon)]  <- as.matrix(SummarizedExperiment::assays(amplicon)$coverage)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(consensus = as(consensus, "dgCMatrix"), fraction = as(fraction, "dgCMatrix"), coverage = as(reads, "dgCMatrix"), alts = as(alts, "dgCMatrix"), refs = as(refs, "dgCMatrix")),
                                                   colData = new_meta_data, rowData = new_row_data)
  return(se)
}
