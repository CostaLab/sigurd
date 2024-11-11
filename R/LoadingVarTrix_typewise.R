#'LoadingVarTrix_typewise
#'@description
#'When we load all the different types of results (scRNAseq/amplicon and MT/amplicon),
#'we might need extreme amounts of memory. To solve this issue, I will load each type separately.
#'In a following function (AmpliconSupplementing), we can add the amplicon information to the
#'scRNAseq information.
#'The input file is a specifically formatted csv file with all the necessary information to run the analysis.
#'Note that the source column in the input file needs to be vartrix for this function.
#'This is hard coded and case insensitive.
#'@importFrom utils read.csv read.table
#'@importFrom VariantAnnotation readVcf info
#'@importFrom SummarizedExperiment SummarizedExperiment
#'@importFrom Matrix readMM colSums rowSums colMeans rowMeans
#'@param samples_path Path to the input folder. Must include a barcodes file.
#'@param samples_file Path to the csv file with the samples to be loaded.
#'@param vcf_path Path to the VCF file with the variants.
#'@param snp_path Path to the SNP file used for VarTrix (SNV.loci.txt). 
#'@param patient The patient you want to load.
#'@param patient_column The column that contains the patient information. Use merge, if all samples should be merged.
#'@param type_use The type of input. Only rows that have the specified type will be loaded.
#'@param min_reads The minimum number of reads we want. Otherwise we treat this as a NoCall. Default = NULL.
#'@param min_cells The minimum number of cells for a variant. Otherwise, we will remove a variant. Default = 2.
#'@param cells_include A vector of cell barcodes. Only these cells will be retained. 
#'@param cells_exclude A vector of cell barcodes. These cells will be removed from the output.
#'@param barcodes_path The path to the cell barcodes tsv. Default = NULL
#'@param cellbarcode_length The length of the cell barcode. This should be the length of the actual barcode plus two for the suffix (-1). Default = 18
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
LoadingVarTrix_typewise <- function(samples_file, samples_path = NULL, barcodes_path = NULL, snp_path = NULL, vcf_path, patient, patient_column = "patient", type_use = "scRNAseq_Somatic", min_reads = NULL, min_cells = 2, cells_include = NULL, cells_exclude = NULL, cellbarcode_length = 18, verbose = TRUE){
  if(all(!is.null(samples_path), !is.null(barcodes_path), !is.null(snp_path))){
    if(verbose) print(paste0("Loading the data for sample ", patient, "."))
    samples_file <- data.frame(patient = patient, sample = patient, input_path = samples_path, cells = barcodes_path)
    samples <- samples_file$sample
  } else{
    if(verbose) print(paste0("Loading the data for patient ", patient, "."))
    if(verbose) print("We read in the design matrix.")
    samples_file <- utils::read.csv(samples_file, stringsAsFactors = FALSE)
    if(!patient_column %in% colnames(samples_file) & patient_column != "merge"){
      stop(paste0("Error: the column ", patient_column, " is not in your design matrix."))
    }

    if(verbose) print("We subset to the relevant files.")
    samples_file <- samples_file[grep("vartrix", samples_file$source, ignore.case = TRUE),]
    if(patient_column != "merge") samples_file <- samples_file[samples_file[,patient_column] == patient,]
    samples_file <- samples_file[samples_file$type == type_use,]

    if(verbose) print("We get the different samples.")
    samples <- samples_file$sample
  }


  if(verbose) print("We load the SNV files.")
  if(!is.null(snp_path)){
    path_snps <- snp_path
  } else{
    path_snps <- paste0(samples_file$input_path, "/SNV.loci.txt")
  }


  if(verbose) print("We read the variants.")
  snps_list <- lapply(path_snps, utils::read.table, header = FALSE)
  names(snps_list) <- samples


  if(verbose) print("We read in the cell barcodes output by CellRanger as a list.")
  barcodes <- lapply(samples_file$cells, utils::read.table)
  names(barcodes) <- samples


  if(verbose) print("We read in the vcf file.")
  vcf         <- VariantAnnotation::readVcf(vcf_path)
  vcf_info    <- VariantAnnotation::info(vcf)


  if(verbose) print("We generate more accessible names.")
  if(all(c("GENE", "AA", "CDS") %in% colnames(vcf_info))){
    new_names <- paste0(vcf_info$GENE, "_", vcf_info$AA, "_", vcf_info$CDS)
  } else{
    new_names <- rownames(vcf_info)
    new_names <- gsub(":|\\/|\\?", "_", new_names)
  }

  if(verbose) print("We read in the different sparse genotype matrices as a list.")
  if(verbose) print("We have a slot per type of input data.")
  coverage_matrices <- list()
  ref_matrices <- list()
  consensus_matrices <- list()
  for(i in 1:nrow(samples_file)){
    if(verbose) print(paste0("Loading sample ", i, " of ", nrow(samples_file)))
    input_folder_use <- samples_file$input_path[i]
    sample_use <- samples_file$sample[i]

    # The cell barcodes and variants.
    cellbarcodes_use <- barcodes[[sample_use]]
    variants_use     <- snps_list[[sample_use]]

    # We load the different matrices.
    coverage_matrices[[sample_use]]  <- Matrix::readMM(paste0(input_folder_use, sample_use, "/out_matrix_coverage.mtx"))
    ref_matrices[[sample_use]]       <- Matrix::readMM(paste0(input_folder_use, sample_use, "/ref_matrix_coverage.mtx"))
    consensus_matrices[[sample_use]] <- Matrix::readMM(paste0(input_folder_use, sample_use, "/out_matrix_consensus.mtx"))


    # We rename the rows and columns.
    colnames(coverage_matrices[[sample_use]])  <- paste0(sample_use, "_", cellbarcodes_use$V1)
    colnames(ref_matrices[[sample_use]])       <- paste0(sample_use, "_", cellbarcodes_use$V1)
    colnames(consensus_matrices[[sample_use]]) <- paste0(sample_use, "_", cellbarcodes_use$V1)
    rownames(coverage_matrices[[sample_use]])  <- variants_use$V1
    rownames(ref_matrices[[sample_use]])       <- variants_use$V1
    rownames(consensus_matrices[[sample_use]]) <- variants_use$V1
  }


  if(verbose) print("We merge the matrices.")
  coverage_matrix_total  <- do.call("cbind", coverage_matrices)
  ref_matrix_total       <- do.call("cbind", ref_matrices)
  consensus_matrix_total <- do.call("cbind", consensus_matrices)


  if(verbose) print("We remove the matrix lists.")
  rm(coverage_matrices, ref_matrices, consensus_matrices)
  gc()


  if(!is.null(cells_include)){
    if(verbose) "We remove all cells not in the allow list."
    # Check if any cells would be left.
    check_leftover <- any(colnames(coverage_matrix_total) %in% cells_include)
    if(!check_leftover) stop("No cells are in your supplied list.")
    coverage_matrix_total <- coverage_matrix_total[, colnames(coverage_matrix_total) %in% cells_include, drop = FALSE]
    ref_matrix_total <- ref_matrix_total[, colnames(ref_matrix_total) %in% cells_include, drop = FALSE]
    consensus_matrix_total <- consensus_matrix_total[, colnames(consensus_matrix_total) %in% cells_include, drop = FALSE]
  }


  if(!is.null(cells_exclude)){
    if(verbose) "We remove all cells in the exclusion list."
    # Check if any cells would be left.
    check_leftover <- all(colnames(coverage_matrix_total) %in% cells_exclude)
    if(check_leftover) stop("All cells are in your exclusion list.")
    coverage_matrix_total <- coverage_matrix_total[, !colnames(coverage_matrix_total) %in% cells_exclude]
    ref_matrix_total <- ref_matrix_total[, !colnames(ref_matrix_total) %in% cells_exclude]
    consensus_matrix_total <- consensus_matrix_total[, !colnames(consensus_matrix_total) %in% cells_exclude]
  }


  if(!is.null(min_reads)){
    if(verbose) print(paste0("We set read values below the threshold of ", min_reads, " to 0."))
    if(verbose) print("We then generate the consensus matrix again.")
    ref_matrix_total@x[ref_matrix_total@x < min_reads] <- 0
    coverage_matrix_total@x[coverage_matrix_total@x < min_reads] <- 0

    reference_construction <- ref_matrix_total
    reference_construction@x[reference_construction@x > 0] <- 1

    coverage_construction <- coverage_matrix_total
    coverage_construction@x[coverage_construction@x > 0] <- 2

    consensus_matrix_total <- reference_construction + coverage_construction
    rm(reference_construction, coverage_construction)
  }


  # We check if number of rows of the matrices are the same as the length of the new names.
  if(all(nrow(coverage_matrix_total) != length(new_names))){
    input_rows <- nrow(coverage_matrix_total)
    new_names_length <- length(new_names)
    stop(paste0("Error: you have ", input_rows, " variants in you matrix and ", new_names_length, " actual variant names."))
  }
  
  rownames(coverage_matrix_total)  <- new_names
  rownames(ref_matrix_total)       <- new_names
  rownames(consensus_matrix_total) <- new_names


  if(verbose) print(paste0("We remove variants, that are not detected in at least ", min_cells, " cells."))
  keep_variants <- Matrix::rowSums(consensus_matrix_total >= 1)
  keep_variants <- keep_variants >= min_cells
  # If we only have one cell or one variant, we loose the matrix.
  consensus_matrix_total <- consensus_matrix_total[keep_variants, , drop = FALSE]
  coverage_matrix_total <- coverage_matrix_total[keep_variants, , drop = FALSE]
  ref_matrix_total <- ref_matrix_total[keep_variants, , drop = FALSE]


  if(verbose) print("We remove cells that are always NoCall.")
  consensus_test <- consensus_matrix_total > 0
  keep_cells <- Matrix::colSums(consensus_test) > 0
  # If we only have one cell or one variant, we loose the matrix.
  consensus_matrix_total <- consensus_matrix_total[, keep_cells, drop = FALSE]
  coverage_matrix_total <- coverage_matrix_total[, keep_cells, drop = FALSE]
  ref_matrix_total <- ref_matrix_total[, keep_cells, drop = FALSE]


  if(verbose) print(paste0(type_use, " Variants: ", nrow(consensus_matrix_total)))
  if(verbose) print(paste0(type_use, " Cells: ", ncol(consensus_matrix_total)))

  rm(consensus_test, keep_variants, keep_cells)
  gc(verbose = FALSE)


  if(verbose) print("We transform the sparse matrices to matrices, so we can calculate the fraction.")
  reads_total                           <- coverage_matrix_total + ref_matrix_total
  fraction_total                        <- coverage_matrix_total / reads_total
  fraction_total[is.na(fraction_total)] <- 0
  gc(verbose = FALSE)


  # We check if the matrices are empty (0 cells, 0 variants). Then we simply return NULL.
  dim_test <- dim(reads_total)
  if(any(dim_test == 0)){
    print(paste0("The filtering left ", dim_test[1], " variants and ", dim_test[2], "cells."))
    print("Returning NULL.")
    return(NULL)
  } else {
    if(verbose) print("We generate a SummarizedExperiment object containing the fraction and the consensus matrices.")
    # We want an assay for the Consensus information and for the fraction.
    # As meta data we add a data frame showing the cell id, the associated patient and the sample.
    coverage_depth_per_cell <- Matrix::colMeans(reads_total)
    coverage_depth_per_variant <- Matrix::rowMeans(reads_total)
    meta_data <- data.frame(Cell = colnames(consensus_matrix_total), Patient = patient, Sample = substr(x = colnames(consensus_matrix_total), start = 1, stop = nchar(colnames(consensus_matrix_total))-(cellbarcode_length+1)), Type = type_use, AverageCoverage = coverage_depth_per_cell)
    rownames(meta_data) <- meta_data$Cell
    meta_row <- data.frame(VariantName = rownames(consensus_matrix_total), Depth = coverage_depth_per_variant)
    rownames(meta_row) <- meta_row$VariantName
    se_merged <- SummarizedExperiment::SummarizedExperiment(assays = list(consensus = as(consensus_matrix_total, "CsparseMatrix"), fraction = as(fraction_total, "CsparseMatrix"), coverage = as(reads_total, "CsparseMatrix"), alts = as(coverage_matrix_total, "CsparseMatrix"), refs = as(ref_matrix_total, "CsparseMatrix")),
                                                            colData = meta_data, rowData = meta_row)
    return(se_merged)
  }
}

