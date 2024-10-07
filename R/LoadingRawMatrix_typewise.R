#'LoadingRawMatrix_typewise
#'@description
#'A function to load a dense matrix of genotyping information and reformat it. Each matrix should have a column with the cell barcodes, one with the reference reads and one with the alternative reads.
#'The type in the design matrix is used as the variant name. 
#'@importFrom utils read.table
#'@importFrom SummarizedExperiment SummarizedExperiment
#'@importFrom Matrix colMeans rowMeans
#'@param samples_path Path to the input folder. Must include a barcodes file.
#'@param samples_file Path to the csv file with the samples to be loaded.
#'@param patient The patient you want to load.
#'@param patient_column The column that contains the patient information.
#'@param variant_use The variant of the respective matrices.
#'@param cells_include A vector of cell barcodes. Only these cells will be retained. 
#'@param cells_exclude A vector of cell barcodes. These cells will be removed from the output.
#'@param cellbarcode_length The length of the cell barcode. This should be the length of the actual barcode plus two for the suffix (-1). Default = 18
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
LoadingRawMatrix_typewise <- function(samples_file, patient, patient_column = "patient", variant_use = NULL, cells_include = NULL, cells_exclude = NULL, matrix_column_separator = ",", matrix_header = TRUE, cellbarcode_index = 1, ref_index = 2, alt_index = 3, cellbarcode_length = 18, verbose = TRUE){
  if(is.null(variant_use)) stop ("You have to supply a variant.")
  if(cellbarcode_index == 0) stop("Your cellbarcode_index cannot be 0.")
  if(ref_index == 0) stop("Your ref_index cannot be 0.")
  if(alt_index == 0) stop("Your alt_index cannot be 0.")

  if(verbose) print("We load the design matrix.")
  samples_file <- read.table(samples_file, sep = ",", header = TRUE)
  if(!patient_column %in% colnames(samples_file)) stop("Your patient_column is not in the design matrix.")
  if(!patient %in% samples_file$patient) stop("Your patient is not in the design matrix.")
  samples_file <- subset(samples_file, patient == patient)
  samples_file <- samples_file[grep("matrix", samples_file$source, ignore.case = TRUE),]
  if(any(!variant_use %in% samples_file$type)){
    stop("Your variant in not in the design matrix.")
  }
  samples_file <- samples_file[samples_file$type == variant_use, , drop = FALSE]
  samples <- samples_file$sample
  variant_use <- samples_file$type


  if(verbose) print("We load the data.")
  loaded_matrix_ls <- c()
  cellbarcodes_total <- c()
  for(i in 1:nrow(samples_file)){
    loaded_matrix <- read.table(samples_file$input_path[i], sep = matrix_column_separator, header = matrix_header)
    if(cellbarcode_index > ncol(loaded_matrix)) stop(paste0("Your cellbarcode_index is outside the matrix dimensions for the ", i , " sample."))
    if(ref_index > ncol(loaded_matrix)) stop(paste0("Your ref_index is outside the matrix dimensions for the ", i , " sample."))
    if(alt_index > ncol(loaded_matrix)) stop(paste0("Your alt_index is outside the matrix dimensions for the ", i , " sample."))
    loaded_matrix[, cellbarcode_index] <- paste0(samples[i], "_", loaded_matrix[, cellbarcode_index])
    if(!is.null(cells_include)){
      if(verbose) paste0("We remove all cells not in the white list for sample ", i, ".")
      # Check if any cells would be left.
      check_leftover <- any(loaded_matrix[, cellbarcode_index] %in% cells_include)
      if(!check_leftover) stop("No cells are in your supplied list for the ", i, " sample.")
      loaded_matrix <- loaded_matrix[colnames(coverage_matrix_total) %in% cells_include, , drop = FALSE]
    }
    if(!is.null(cells_exclude)){
      if(verbose) paste0("We remove all cells in the exclusion list for sample ", i, ".")
      # Check if any cells would be left.
      check_leftover <- all(loaded_matrix[, cellbarcode_index] %in% cells_exclude)
      if(check_leftover) stop(paste0("All cells are in your exclusion list for sample ", i, "."))
      loaded_matrix <- loaded_matrix[!loaded_matrix[, cellbarcode_index] %in% cells_exclude, , drop = FALSE]
    }
    cellbarcodes_total <- unique(c(cellbarcodes_total, loaded_matrix[, cellbarcode_index]))
    loaded_matrix_ls[[paste0(samples[i], "_", samples_file$type[i])]] <- loaded_matrix
  }


  alts <- refs <- matrix(0, ncol = length(cellbarcodes_total), nrow = length(variant_use), dimnames = list(variant_use, cellbarcodes_total))
  for(i in 1:length(loaded_matrix_ls)){
    alts[i, loaded_matrix_ls[[i]][, cellbarcode_index]] <- loaded_matrix_ls[[i]][, alt_index]
    refs[i, loaded_matrix_ls[[i]][, cellbarcode_index]] <- loaded_matrix_ls[[i]][, ref_index]
  }
  coverage <- alts + refs


  ref_construction <- refs > 0
  alt_construction <- alts > 0
  alt_construction[alt_construction] <- 2
  consensus <- ref_construction + alt_construction
  fraction <- alts / coverage


  coverage <- as(coverage, "CsparseMatrix")
  alts <- as(alts, "CsparseMatrix")
  refs <- as(refs, "CsparseMatrix")
  consensus <- as(consensus, "CsparseMatrix")
  fraction <- as(fraction, "CsparseMatrix")


  if(verbose) print("We generate a SummarizedExperiment object containing the fraction and the consensus matrices.")
  # We want an assay for the Consensus information and for the fraction.
  # As meta data we add a data frame showing the cell id, the associated patient and the sample.
  average_depth_per_cell <- Matrix::colMeans(coverage)
  coverage_depth_per_variant <- Matrix::rowMeans(coverage)
  meta_data <- data.frame(Cell = colnames(consensus), Patient = patient, Sample = substr(x = colnames(consensus), start = 1, stop = nchar(colnames(consensus))-(cellbarcode_length+1)), AverageCoverage = average_depth_per_cell)
  rownames(meta_data) <- meta_data$Cell
  meta_row <- data.frame(VariantName = rownames(consensus), Depth = coverage_depth_per_variant)
  rownames(meta_row) <- meta_row$VariantName
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(consensus = consensus, fraction = fraction, coverage = coverage, alts = alts, refs = refs), colData = meta_data, rowData = meta_row)
  return(se)
}

