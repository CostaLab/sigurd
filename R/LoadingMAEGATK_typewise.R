#'LoadingMAEGATK_typewise
#'@description
#'We load the MAEGATK output and transform it to be compatible with the VarTrix output.
#'The input file is a specifically formated csv file with all the necessary information to run the analysis.
#'Note that the source column in the input file needs to be one of the following: vartrix, mgaetk, mgatk.
#'This is hard coded and case insensitive.
#'@import Matrix SummarizedExperiment
#'@param samples_path Path to the input folder.
#'@param samples_file Path to the csv file with the samples to be loaded.
#'@param type_use The type of input. Has to be one of: scRNAseq_MT, Amplicon_MT. Only used if samples_path is not NULL.
#'@param patient The patient you want to load.
#'@param chromosome_prefix The prefix you want use. Default: "chrM"
#'@export
LoadingMAEGATK_typewise <- function(samples_file, samples_path = NULL, patient, type_use = "scRNAseq_MT", chromosome_prefix = "chrM",
                                    min_cells = 2, barcodes_path = NULL){
  if(all(!is.null(samples_path), !is.null(barcodes_path))){
    print(paste0("Loading the data for patient ", patient, "."))
    samples <- list.files(samples_path)
    samples <- grep(patient, samples, value = TRUE)
    samples_file <- data.frame(patient = patient, sample = samples, input_folder = samples_path, cells = barcodes_path)
  } else{
    print(paste0("Loading the data for patient ", patient, "."))
    print("We read in the samples file.")
    samples_file <- read.csv(samples_file)

    print("We subset to the patient of interest.")
    samples_file <- samples_file[grep("maegatk|mgatk", samples_file$source, ignore.case = TRUE),]
    samples_file <- samples_file[grep(patient, samples_file$patient),]
    samples_file <- samples_file[grep(type_use, samples_file$type),]

    print("We get the different samples.")
    samples <- samples_file$sample
  }


  print("We read in the cell barcodes output by CellRanger as a list.")
  barcodes <- lapply(samples_file$cells, read.table)
  names(barcodes) <- samples


  print("We load the MAEGATK output files.")
  se_ls <- list()
  for(i in 1:nrow(samples_file)){
    print(paste0("Loading sample ", i, " of ", nrow(samples_file)))
    input_folder_use <- samples_file$input_folder[i]
    sample_use <- samples_file$sample[i]

    # We get the final output file for either mgatk or maegatk.
    final_output_file <- list.files(paste0(input_folder_use, sample_use, "/final/"), full.names = TRUE)
    final_output_file <- grep(paste0("maegtk.rds|maegatk.rds|mgatk.rds|", sample_use, ".rds"), final_output_file, value = TRUE)
    se_ls[[sample_use]] <- load_object(final_output_file)
    colnames(se_ls[[sample_use]]) <- paste0(sample_use, "_", colnames(se_ls[[sample_use]]))
    barcodes_use <- paste0(sample_use, "_", barcodes[[sample_use]][,1])
    barcodes_use <- barcodes_use[barcodes_use %in% colnames(se_ls[[sample_use]])]
    se_ls[[sample_use]] <- se_ls[[sample_use]][,barcodes_use]
  }


  print("We merge the samples.")
  se_merged <- do.call("cbind", se_ls)
  rm(se_ls)
  gc()


  print("We get the allele frequency.")
  fraction <- computeAFMutMatrix(se_merged, chromosome_prefix = chromosome_prefix)


  print("We get the coverage information.")
  coverage <- CalculateCoverage(SE = se_merged, chromosome_prefix = chromosome_prefix)
  if(!all(rownames(fraction) == rownames(coverage))){
    coverage <- coverage[match(rownames(fraction), rownames(coverage)),]
  }


  print("We get the number of alternative reads per variant.")
  reads_alt <- CalculateAltReads(SE = se_merged, chromosome_prefix = chromosome_prefix)


  print("We get the quality information.")
  variant_quality <- CalculateQuality(SE = se_merged, variants = rownames(reads_alt), chromosome_prefix = chromosome_prefix)


  print("We get the number of reference reads.")
  reads_ref <- coverage - reads_alt

  print("Calculating the strand concordance.")
  concordance <- CalculateStrandCorrelation(SE = se_merged, chromosome_prefix = chromosome_prefix)


  print("We calculate the consensus information.")
  consensus <- CalculateConsensus(SE = se_merged, chromosome_prefix = chromosome_prefix)
  # We order the consensus matrix like the coverage matrix.
  if(!all(rownames(fraction) == rownames(consensus))){
    consensus <- consensus[match(rownames(fraction), rownames(consensus)),]
  }


  print("We perform some filtering to reduce the memory needed.")
  print(paste0("We remove variants, which are not covered in at least ", min_cells, " cells ."))
  keep_variants <- rowSums(consensus >= 1)
  keep_variants <- keep_variants >= min_cells
  # If we only have one cell or one variant, we loose the matrix.
  cell_ids <- colnames(consensus)
  variant_names <- names(keep_variants[keep_variants])
  # consensus <- consensus[keep_variants,]
  # coverage <- coverage[keep_variants,]
  # fraction <- fraction[keep_variants,]
  # concordance <- concordance[keep_variants]
  # reads_alt <- reads_alt[keep_variants,]
  # reads_ref <- reads_ref[keep_variants,]
  consensus <- consensus[keep_variants,]
  consensus <- suppressWarnings(matrix(consensus, nrow = length(variant_names), ncol = length(cell_ids)))
  colnames(consensus) <- cell_ids
  rownames(consensus) <- variant_names
  consensus <- as(consensus, "dgCMatrix")
  coverage <- coverage[keep_variants,]
  coverage <- suppressWarnings(matrix(coverage, nrow = length(variant_names), ncol = length(cell_ids)))
  colnames(coverage) <- cell_ids
  rownames(coverage) <- variant_names
  coverage <- as(coverage, "dgCMatrix")
  fraction <- fraction[keep_variants,]
  fraction <- suppressWarnings(matrix(fraction, nrow = length(variant_names), ncol = length(cell_ids)))
  colnames(fraction) <- cell_ids
  rownames(fraction) <- variant_names
  fraction <- as(fraction, "dgCMatrix")
  concordance <- concordance[keep_variants]
  variant_quality <- variant_quality[keep_variants]
  reads_alt <- reads_alt[keep_variants,]
  reads_alt <- suppressWarnings(matrix(reads_alt, nrow = length(variant_names), ncol = length(cell_ids)))
  colnames(reads_alt) <- cell_ids
  rownames(reads_alt) <- variant_names
  reads_alt <- as(reads_alt, "dgCMatrix")
  reads_ref <- reads_ref[keep_variants,]
  reads_ref <- suppressWarnings(matrix(reads_ref, nrow = length(variant_names), ncol = length(cell_ids)))
  colnames(reads_ref) <- cell_ids
  rownames(reads_ref) <- variant_names
  reads_ref <- as(reads_ref, "dgCMatrix")


  print("We remove cells that are always NoCall.")
  consensus_test <- consensus > 0
  keep_cells <- colSums(consensus_test) > 0
  # If we only have one cell or one variant, we loose the matrix.
  cell_ids <- colnames(consensus)
  variant_names <- names(keep_variants[keep_variants])
  # consensus <- consensus[,keep_cells]
  # coverage <- coverage[,keep_cells]
  # fraction <- fraction[,keep_cells]
  # reads_alt <- reads_alt[,keep_cells]
  # reads_ref <- reads_ref[,keep_cells]
  consensus <- consensus[,keep_cells]
  consensus <- suppressWarnings(matrix(consensus, nrow = length(variant_names), ncol = length(cell_ids)))
  colnames(consensus) <- cell_ids
  rownames(consensus) <- variant_names
  consensus <- as(consensus, "dgCMatrix")
  coverage <- coverage[,keep_cells]
  coverage <- suppressWarnings(matrix(coverage, nrow = length(variant_names), ncol = length(cell_ids)))
  colnames(coverage) <- cell_ids
  rownames(coverage) <- variant_names
  coverage <- as(coverage, "dgCMatrix")
  fraction <- fraction[,keep_cells]
  fraction <- suppressWarnings(matrix(fraction, nrow = length(variant_names), ncol = length(cell_ids)))
  colnames(fraction) <- cell_ids
  rownames(fraction) <- variant_names
  fraction <- as(fraction, "dgCMatrix")
  reads_alt <- reads_alt[,keep_cells]
  reads_alt <- suppressWarnings(matrix(reads_alt, nrow = length(variant_names), ncol = length(cell_ids)))
  colnames(reads_alt) <- cell_ids
  rownames(reads_alt) <- variant_names
  reads_alt <- as(reads_alt, "dgCMatrix")
  reads_ref <- reads_ref[,keep_cells]
  reads_ref <- suppressWarnings(matrix(reads_ref, nrow = length(variant_names), ncol = length(cell_ids)))
  colnames(reads_ref) <- cell_ids
  rownames(reads_ref) <- variant_names
  reads_ref <- as(reads_ref, "dgCMatrix")


  # We check if the matrices are empty (0 cells, 0 variants). Then we simply return NULL.
  dim_test <- dim(coverage)
  if(any(dim_test == 0)){
    print(paste0("The filtering left ", dim_test[1], " variants and ", dim_test[2], "cells."))
    print("Returning NULL.")
    return(NULL)
  } else{
    print("We add the information to the merged matrices.")
    coverage_depth_per_cell <- rownames(coverage)
    coverage_depth_per_cell <- gsub("_._.$", "", coverage_depth_per_cell)
    coverage_depth_per_cell <- !duplicated(coverage_depth_per_cell)
    coverage_depth_per_cell <- coverage[coverage_depth_per_cell,]
    cell_ids <- colnames(coverage_depth_per_cell)
    variant_names <- rownames(coverage_depth_per_cell)
    coverage_depth_per_cell <- suppressWarnings(matrix(coverage, nrow = length(variant_names), ncol = length(cell_ids)))
    colnames(coverage_depth_per_cell) <- cell_ids
    rownames(coverage_depth_per_cell) <- variant_names
    coverage_depth_per_variant <- rowMeans(coverage)
    coverage_depth_per_cell <- colMeans(coverage_depth_per_cell)
    meta_data_col <- data.frame(Cell = colnames(consensus), AverageCoverage = coverage_depth_per_cell)
    rownames(meta_data_col) <- meta_data_col$Cell
    meta_data_row <- data.frame(VariantName = rownames(consensus), Concordance = concordance, VariantQuality = variant_quality, Depth = coverage_depth_per_variant)
    rownames(meta_data_row) <- meta_data_row$VariantName
    se_output <- SummarizedExperiment(assays = list(consensus = consensus, fraction = fraction, coverage = coverage, alts = reads_alt, refs = reads_ref),
                                      colData = meta_data_col, rowData = meta_data_row)
    return(se_output)
  }
}
