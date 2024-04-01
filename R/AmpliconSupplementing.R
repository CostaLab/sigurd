#'Supplementing scRNAseq values with Amplicon values
#'
#'@description
#'We replace the values from an scRNAseq experiment with values we have from an amplicon experiment.
#'@importFrom S4Vectors merge
#'@importFrom SummarizedExperiment colData rowData assays SummarizedExperiment
#'@importFrom stats na.omit
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
  names(amplicon_value)                                   <- rownames(new_meta_data)
  amplicon_value                                          <- stats::na.omit(amplicon_value)
  new_meta_data[names(amplicon_value), "AverageCoverage"] <- amplicon_value
  # We add a Patient and Sample column to the new_meta_data.
  patient_values <- new_meta_data[, paste0("Patient", c("scRNAseq", "Amplicon"))]
  patient_values <- ifelse(is.na(patient_values[,1]), patient_values[,2], patient_values[,1])
  sample_values <- new_meta_data[, paste0("Sample", c("scRNAseq", "Amplicon"))]
  sample_values <- ifelse(is.na(sample_values[,1]), sample_values[,2], sample_values[,1])
  columns_to_keep <- colnames(new_meta_data)
  columns_to_keep <- columns_to_keep[!columns_to_keep %in% c("Cell", "PatientscRNAseq", "SamplecRNAseq", "PatientAmplicon", "SampleAmplicon")]
  new_meta_data <- S4Vectors::DataFrame(Cell = new_meta_data$Cell, Patient = patient_values, Sample = sample_values, new_meta_data[,columns_to_keep], row.names = new_meta_data$Cell)
  
  new_row_data <- merge(SummarizedExperiment::rowData(scRNAseq), SummarizedExperiment::rowData(amplicon), by = "VariantName", all.x = TRUE, all.y = TRUE,
                        suffixes = c("scRNAseq", "Amplicon"))
  rownames(new_row_data) <- new_row_data$VariantName
  if("VariantQualityscRNAseq" %in% colnames(new_row_data)){
    # We add a VariantQuality column to the row data, showing the scRNAseq quality with the supplemented amplicon quality.
    # We do the same for the concordance and the depth.
    new_row_data$VariantQuality                           <- new_row_data$VariantQualityscRNAseq
    amplicon_value                                        <- SummarizedExperiment::rowData(amplicon)$VariantQuality
    names(amplicon_value)                                 <- rownames(amplicon)
    amplicon_value                                        <- stats::na.omit(amplicon_value)
    new_row_data[names(amplicon_value), "VariantQuality"] <- amplicon_value
  }

  if("ConcordancescRNAseq" %in% colnames(new_row_data)){
    new_row_data$Concordance                           <- new_row_data$ConcordancescRNAseq
    amplicon_value                                     <- SummarizedExperiment::rowData(amplicon)$Concordance
    names(amplicon_value)                              <- rownames(amplicon)
    amplicon_value                                     <- stats::na.omit(amplicon_value)
    new_row_data[names(amplicon_value), "Concordance"] <- amplicon_value
  }

  new_row_data$Depth                           <- new_row_data$DepthscRNAseq
  amplicon_value                               <- SummarizedExperiment::rowData(amplicon)$Depth
  names(amplicon_value)                        <- rownames(amplicon)
  amplicon_value                               <- stats::na.omit(amplicon_value)
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

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(consensus = methods::as(consensus, "dgCMatrix"), fraction = methods::as(fraction, "dgCMatrix"), coverage = methods::as(reads, "dgCMatrix"), alts = methods::as(alts, "dgCMatrix"), refs = methods::as(refs, "dgCMatrix")),
                                                   colData = new_meta_data, rowData = new_row_data)
  return(se)
}
