#'LoadingMGATK_typewise
#'@description
#'We load the MGATK output and transform it to be compatible with the VarTrix output.
#'The input file is a specifically formatted csv file with all the necessary information to run the analysis.
#'Note that the source column in the input file needs to be mgaetk or mgatk for this function. This is case insensitive.
#'If you want to only load a single sample without the use of an input file, you have to set the following variables.
#' \enumerate{
#'   \item samples_path
#'   \item barcodes_path
#'   \item patient
#'   \item samples_file = NULL
#' }
#'@importFrom utils read.csv read.table
#'@importFrom SummarizedExperiment SummarizedExperiment
#'@importFrom Matrix rowSums colSums
#'@param samples_path Path to the input folder.
#'@param samples_file Path to the csv file with the samples to be loaded.
#'@param type_use The type of input. Only rows that have the specified type will be loaded.
#'@param patient The patient you want to load.
#'@param patient_column The column that contains the patient information. Use merge, if all samples should be merged.
#'@param chromosome_prefix The prefix you want use. Default: "chrM"
#'@param min_cells The minimum number of cells with coverage for a variant. Variants with coverage in less than this amount of cells are removed. Default = 2
#'@param cells_include A vector of cell barcodes. Only these cells will be retained. 
#'@param cells_exclude A vector of cell barcodes. These cells will be removed from the output.
#'@param barcodes_path Path to the barcodes file tsv. Default = NULL
#'@param cellbarcode_length The length of the cell barcode. This should be the length of the actual barcode plus two for the suffix (-1). Default = 18
#'@param ignore_quality Should the quality assay be ignore? This can be useful if the quality was not reported in mgatk.
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
LoadingMGATK_typewise <- function(samples_file, patient, samples_path = NULL, patient_column = "patient", type_use = "scRNAseq_MT", chromosome_prefix = "chrM",
                                  min_cells = 2, cells_include = NULL, cells_exclude = NULL, barcodes_path = NULL, cellbarcode_length = 18, ignore_quality = TRUE, verbose = TRUE){
  if(all(!is.null(samples_path), !is.null(barcodes_path))){
    if(verbose) print(paste0("Loading the data for patient ", patient, "."))
    samples <- patient
    samples_file <- data.frame(patient = patient, sample = samples, input_path = samples_path, cells = barcodes_path)
  } else{
    if(verbose) print(paste0("Loading the data for patient ", patient, "."))
    if(verbose) print("We read in the design matrix.")
    samples_file <- utils::read.csv(samples_file)
    if(!patient_column %in% colnames(samples_file) & patient_column != "merge"){
      stop(paste0("Error: the column ", patient_column, " is not in your design matrix."))
    }

    if(verbose) print("We subset to the relevant files.")
    samples_file <- samples_file[grep("mgatk", samples_file$source, ignore.case = TRUE),]
    if(patient_column != "merge") samples_file <- samples_file[samples_file[,patient_column] == patient,]
    samples_file <- samples_file[samples_file$type == type_use,]

    if(verbose) print("We get the different samples.")
    samples <- samples_file$sample
  }


  if(verbose) print("We read in the cell barcodes output by CellRanger as a list.")
  barcodes        <- lapply(samples_file$cells, utils::read.table)
  names(barcodes) <- samples


  if(verbose) print("We load the MAEGATK output files.")
  se_ls <- list()
  for(i in 1:nrow(samples_file)){
    if(verbose) print(paste0("Loading sample ", i, " of ", nrow(samples_file)))
    input_file_use <- samples_file$input_path[i]
    sample_use     <- samples_file$sample[i]

    # We check if the file exists.
    if(!file.exists(input_file_use)){
      stop(paste0("Error: the file ", input_file_use, " does not exist."))
    }

    # We get the final output file for either mgatk or maegatk.
    se_ls[[sample_use]]           <- sigurd::load_object(input_file_use)
    colnames(se_ls[[sample_use]]) <- paste0(sample_use, "_", colnames(se_ls[[sample_use]]))
    barcodes_use                  <- paste0(sample_use, "_", barcodes[[sample_use]][,1])
    barcodes_use                  <- barcodes_use[barcodes_use %in% colnames(se_ls[[sample_use]])]
    se_ls[[sample_use]]           <- se_ls[[sample_use]][, barcodes_use]
  }


  if(verbose) print("We merge the samples.")
  se_merged <- do.call("cbind", se_ls)
  rm(se_ls)
  gc()


  if(!is.null(cells_include)){
    if(verbose) "We remove all cells not in the allow list."
    # Check if any cells would be left.
    check_leftover <- any(colnames(se_merged) %in% cells_include)
    if(!check_leftover) stop("No cells are in your supplied list.")
    se_merged <- se_merged[, colnames(se_merged) %in% cells_include]
  }


  if(!is.null(cells_exclude)){
    if(verbose) "We remove all cells in the exclusion list."
    # Check if any cells would be left.
    check_leftover <- all(colnames(se_merged) %in% cells_exclude)
    if(check_leftover) stop("All cells are in your exclusion list.")
    se_merged <- se_merged[, !colnames(se_merged) %in% cells_exclude]
  }


  if(verbose) print("We get the allele frequency.")
  fraction <- computeAFMutMatrix(SE = se_merged, chromosome_prefix = chromosome_prefix)


  if(verbose) print("We get the coverage information.")
  coverage <- CalculateCoverage(SE = se_merged, chromosome_prefix = chromosome_prefix)
  if(!all(rownames(fraction) == rownames(coverage))){
    coverage <- coverage[match(rownames(fraction), rownames(coverage)),]
  }


  if(verbose) print("We get the number of alternative reads per variant.")
  reads_alt <- CalculateAltReads(SE = se_merged, chromosome_prefix = chromosome_prefix)


  if(verbose) print("We get the number of forward alternative reads per variant.")
  reads_forward <- CalculateForwardReads(SE = se_merged, chromosome_prefix = chromosome_prefix)


  if(verbose) print("We get the number of reverse alternative reads per variant.")
  reads_reverse <- CalculateReverseReads(SE = se_merged, chromosome_prefix = chromosome_prefix)


  if(!ignore_quality){
    if(verbose) print("We get the quality information.")
    # Checking if the quality assays are present.
    quality_assay_check <- sum(grepl("_qual_", names(SummarizedExperiment::assays(se_merged)))) # This has to be 8 (1 Forward and 1 Reverse for 4 bases.)
    if(quality_assay_check != 8){
      quality_assays_present <- paste0(grep("_qual", names(SummarizedExperiment::assays(se_merged)), value = TRUE), collapse = ", ")
      stop(paste0("Your quality assays are", quality_assays_present, ". Check if you forgot to emit the quality in mgatk."))
    }
    variant_quality <- CalculateQuality(SE = se_merged, variants = rownames(reads_alt), chromosome_prefix = chromosome_prefix)
  }

  if(verbose) print("We get the number of reference reads.")
  reads_ref <- CalculateRefReads(SE = se_merged, chromosome_prefix = chromosome_prefix)
  reads_ref <- reads_ref[rownames(coverage), , drop = FALSE]


  if(verbose) print("Calculating the strand concordance.")
  concordance <- CalculateStrandCorrelation(SE = se_merged, chromosome_prefix = chromosome_prefix)


  if(verbose) print("We calculate the consensus information.")
  consensus <- CalculateConsensus(SE = se_merged, chromosome_prefix = chromosome_prefix)
  # We order the consensus matrix like the coverage matrix.
  if(!all(rownames(fraction) == rownames(consensus))){
    consensus <- consensus[match(rownames(fraction), rownames(consensus)),]
  }


  if(verbose) print("We perform some filtering to reduce the memory needed.")
  if(verbose) print(paste0("We remove variants, which are not covered in at least ", min_cells, " cells ."))
  keep_variants <- Matrix::rowSums(consensus >= 1)
  keep_variants <- keep_variants >= min_cells
  consensus     <- consensus[keep_variants, , drop = FALSE]
  coverage      <- coverage[keep_variants, ,  drop = FALSE]
  fraction      <- fraction[keep_variants, ,  drop = FALSE]
  concordance   <- concordance[keep_variants]
  if(!ignore_quality) variant_quality <- variant_quality[keep_variants]
  reads_alt     <- reads_alt[keep_variants, , drop = FALSE]
  reads_ref     <- reads_ref[keep_variants, , drop = FALSE]
  reads_forward <- reads_forward[keep_variants, , drop = FALSE]
  reads_reverse <- reads_reverse[keep_variants, , drop = FALSE]


  if(verbose) print("We remove cells that are always NoCall.")
  consensus_test <- consensus > 0
  keep_cells     <- Matrix::colSums(consensus_test) > 0
  consensus      <- consensus[, keep_cells, drop = FALSE]
  coverage       <- coverage[,  keep_cells, drop = FALSE]
  fraction       <- fraction[,  keep_cells, drop = FALSE]
  reads_alt      <- reads_alt[, keep_cells, drop = FALSE]
  reads_ref      <- reads_ref[, keep_cells, drop = FALSE]
  reads_forward  <- reads_forward[, keep_cells, drop = FALSE]
  reads_reverse  <- reads_reverse[, keep_cells, drop = FALSE]


  # We check if the matrices are empty (0 cells, 0 variants). Then we simply return NULL.
  dim_test <- dim(coverage)
  if(any(dim_test == 0)){
    if(verbose) print(paste0("The filtering left ", dim_test[1], " variants and ", dim_test[2], "cells."))
    if(verbose) print("Returning NULL.")
    return(NULL)
  } else{
    if(verbose) print("We add the information to the merged matrices.")
    coverage_depth_per_cell    <- rownames(coverage)
    coverage_depth_per_cell    <- gsub("_._.$", "", coverage_depth_per_cell)
    coverage_depth_per_cell    <- !duplicated(coverage_depth_per_cell)
    coverage_depth_per_cell    <- coverage[coverage_depth_per_cell,]
    coverage_depth_per_variant <- rowMeans(coverage)
    coverage_depth_per_cell    <- colMeans(coverage_depth_per_cell)
    meta_data_col              <- data.frame(Cell = colnames(consensus), Patient = patient, Sample = substr(x = colnames(consensus), start = 1, stop = nchar(colnames(consensus))-(cellbarcode_length+1)), AverageCoverage = coverage_depth_per_cell)
    rownames(meta_data_col)    <- meta_data_col$Cell
    if(ignore_quality){
      meta_data_row <- data.frame(VariantName = rownames(consensus), Concordance = concordance, Depth = coverage_depth_per_variant)
    } else{
      meta_data_row <- data.frame(VariantName = rownames(consensus), Concordance = concordance, VariantQuality = variant_quality, Depth = coverage_depth_per_variant)
    }
    rownames(meta_data_row) <- meta_data_row$VariantName

    se_output <- SummarizedExperiment::SummarizedExperiment(assays = list(consensus = as(consensus, "CsparseMatrix"), fraction = as(fraction, "CsparseMatrix"), coverage = as(coverage, "CsparseMatrix"), alts = as(reads_alt, "CsparseMatrix"), refs = as(reads_ref, "CsparseMatrix"),
                                                                          forward = as(reads_forward, "CsparseMatrix"), reverse = as(reads_reverse, "CsparseMatrix")),
                                                            colData = meta_data_col, rowData = meta_data_row)
    return(se_output)
  }
}
