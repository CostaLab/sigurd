#'LoadingVCF_typewise
#'@description
#' We load a cellwise pileup result from a VCF file.
#' If you want to only load a single sample without the use of an input file, you have to set the following variables.
#' \enumerate{
#'   \item samples_path
#'   \item patient
#'   \item samples_file = NULL
#' }
#' Note that the source column in the input file needs to be vcf for this function. This is case insensitive.
#'
#' It has happened that reads with an N allele were aligned. This can cause problems since these variants are typically not in variants lists.
#' We can remove all of these variants by setting remove_N_alternative to TRUE (the default).
#' Set this option to FALSE, if you really want to retain these variants.
#'@importFrom GenomeInfoDb seqnames
#'@importFrom BiocGenerics start
#'@importFrom utils read.table read.csv
#'@importFrom VariantAnnotation readVcf geno info ref alt ScanVcfParam
#'@importFrom SummarizedExperiment SummarizedExperiment
#'@importFrom Matrix rowSums colSums rowMeans colMeans
#'@param samples_path Path to the input folder.
#'@param samples_file Path to the csv file with the samples to be loaded.
#'@param vcf_path Path to the VCF file with the variants.
#'@param patient The patient you want to load.
#'@param patient_column The column that contains the patient information. Use merge, if all samples should be merged.
#'@param type_use The type of input. Only rows that have the specified type will be loaded.
#'@param min_reads The minimum number of reads we want. Otherwise we treat this as a NoCall. Default = NULL.
#'@param min_cells The minimum number of cells for a variant. Otherwise, we will remove a variant. Default = 2.
#'@param cells_include A vector of cell barcodes. Only these cells will be retained. 
#'@param cells_exclude A vector of cell barcodes. These cells will be removed from the output.
#'@param remove_N_alternative Remove all variants that have N as an alternative, see Description. Default = TRUE
#'@param cellbarcode_length The length of the cell barcode. This should be the length of the actual barcode plus two for the suffix (-1). Default = 18
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
LoadingVCF_typewise <- function(samples_file, samples_path = NULL, vcf_path, patient, patient_column = "patient", type_use = "scRNAseq_Somatic", min_reads = NULL, min_cells = 2, cells_include = NULL, cells_exclude = NULL, remove_N_alternative = TRUE, cellbarcode_length = 18, verbose = TRUE){
  if(!is.null(samples_path)){
    if(verbose) print(paste0("Loading the data for sample ", patient, "."))
    samples_file <- data.frame(patient = patient, sample = patient, input_path = samples_path)
    samples <- samples_file$sample
  } else{
    if(verbose) print(paste0("Loading the data for patient ", patient, "."))
    if(verbose) print("We read in the central input file.")
    samples_file <- utils::read.csv(samples_file, stringsAsFactors = FALSE)
    if(!patient_column %in% colnames(samples_file) & patient_column != "merge"){
      stop(paste0("Error: the column ", patient_column, " is not in your central input file."))
    }

    if(verbose) print("We subset to the relevant files.")
    samples_file <- samples_file[grep("vcf", samples_file$source, ignore.case = TRUE),]
    if(patient_column != "merge") samples_file <- samples_file[samples_file[,patient_column] == patient,]
    samples_file <- samples_file[samples_file$type == type_use,]

    if(verbose) print("We get the different samples.")
    samples <- samples_file$sample
  }


  if(verbose) print("We read in the vcf file.")
  vcf      <- VariantAnnotation::readVcf(vcf_path)
  vcf_info <- VariantAnnotation::info(vcf)


  if(verbose) print("We load the samples VCF file.")
  reads_matrix_total <- c() # The total number of reads
  coverage_matrix_total <- c() # The alternative reads
  ref_matrix_total <- c() # The reference reads
  consensus_matrix_total <- c()
  for(i in 1:length(samples)){
    if(verbose) print(paste0("Loading sample ", i, " of ", nrow(samples_file)))
    input_folder_use <- samples_file$input_path[i]
    sample_use <- samples_file$sample[i]

    # We load the VCF file.
    depth_to_add                      <- VariantAnnotation::readVcf(file = input_folder_use, param = VariantAnnotation::ScanVcfParam(geno = "DP"))
    depth_to_add                      <- VariantAnnotation::geno(depth_to_add)$DP
    depth_to_add[is.na(depth_to_add)] <- 0
    rownames(depth_to_add)            <- make.names(rownames(depth_to_add))
    colnames(depth_to_add)            <- paste0(sample_use, "_", colnames(depth_to_add))
    depth_to_add                      <- methods::as(depth_to_add, "CsparseMatrix")
    reads_matrix_total                <- cbind(reads_matrix_total, depth_to_add)

    alts_to_add                       <- VariantAnnotation::readVcf(file = input_folder_use, param = VariantAnnotation::ScanVcfParam(geno = "AD"))
    alts_to_add                       <- VariantAnnotation::geno(alts_to_add)$AD
    alts_to_add[is.na(alts_to_add)]   <- 0
    rownames(alts_to_add)             <- make.names(rownames(alts_to_add))
    colnames(alts_to_add)             <- paste0(sample_use, "_", colnames(alts_to_add))
    alts_to_add                       <- methods::as(alts_to_add, "CsparseMatrix")
    coverage_matrix_total             <- cbind(coverage_matrix_total, alts_to_add)

    consensus_to_add                  <- VariantAnnotation::readVcf(file = input_folder_use, param = VariantAnnotation::ScanVcfParam(geno = "GT"))
    consensus_to_add                  <- VariantAnnotation::geno(consensus_to_add)$GT
    consensus_to_add                  <- matrix(sapply(consensus_to_add, char_to_numeric), nrow = nrow(consensus_to_add), dimnames = list(make.names(rownames(consensus_to_add)), paste0(sample_use, "_", colnames(consensus_to_add))))
    consensus_to_add                  <- methods::as(consensus_to_add, "CsparseMatrix")
    consensus_matrix_total            <- cbind(consensus_matrix_total, consensus_to_add)
  }
  ref_matrix_total <- reads_matrix_total - coverage_matrix_total
  rm(consensus_to_add, alts_to_add, depth_to_add)


  if(!is.null(cells_include)){
    if(verbose) "We remove all cells not in the allow list."
    # Check if any cells would be left.
    check_leftover <- any(colnames(coverage_matrix_total) %in% cells_include)
    if(!check_leftover) stop("No cells are in your supplied list.")
    coverage_matrix_total <- coverage_matrix_total[, colnames(coverage_matrix_total) %in% cells_include, drop = FALSE]
    ref_matrix_total <- ref_matrix_total[, colnames(ref_matrix_total) %in% cells_include, drop = FALSE]
    reads_matrix_total <- reads_matrix_total[, colnames(reads_matrix_total) %in% cells_include, drop = FALSE]
  }


  if(!is.null(cells_exclude)){
    if(verbose) "We remove all cells in the exclusion list."
    # Check if any cells would be left.
    check_leftover <- all(colnames(coverage_matrix_total) %in% cells_exclude)
    if(!check_leftover) stop("All cells are in your supplied list.")
    coverage_matrix_total <- coverage_matrix_total[, !colnames(coverage_matrix_total) %in% cells_exclude]
    ref_matrix_total <- ref_matrix_total[, !colnames(ref_matrix_total) %in% cells_exclude]
    reads_matrix_total <- reads_matrix_total[, !colnames(reads_matrix_total) %in% cells_exclude]
  }


  # We can get the N allele as an alternative allele. This happened in a visium data set.
  # We remove all variants with the N allele as alternative.
  if(remove_N_alternative){
    ref_matrix_total_n     <- substr(rownames(ref_matrix_total), start = nchar(rownames(ref_matrix_total)), stop = nchar(rownames(ref_matrix_total)))
    ref_matrix_total_n     <- ref_matrix_total_n != "N"
    ref_matrix_total       <- ref_matrix_total[ref_matrix_total_n, , drop = FALSE]
    reads_matrix_total     <- reads_matrix_total[ref_matrix_total_n, , drop = FALSE]
    coverage_matrix_total  <- coverage_matrix_total[ref_matrix_total_n, , drop = FALSE]
    consensus_matrix_total <- consensus_matrix_total[ref_matrix_total_n, , drop = FALSE]
    rm(ref_matrix_total_n)
  } else{
    print("We keep all variants with an N as alternative allele. Please ensure that these variants are in your variant VCF file.")
  }


  if(verbose) print("We generate more accessible names.")
  if(all(c("GENE", "AA", "CDS") %in% colnames(vcf_info))){
    new_names <- paste0(vcf_info$GENE, "_", vcf_info$AA, "_", vcf_info$CDS)
    names(new_names) <- make.names(paste0(as.character(rep(GenomeInfoDb::seqnames(vcf)@values, GenomeInfoDb::seqnames(vcf)@lengths)), ".", BiocGenerics::start(vcf), "_", as.character(VariantAnnotation::ref(vcf)), ".", as.character(unlist(VariantAnnotation::alt(vcf)))))
    new_names <- new_names[rownames(ref_matrix_total)]
  } else{
    new_names <- rownames(vcf_info)
    new_names <- gsub(":|\\/|\\?", "_", new_names)
    names(new_names) <-	make.names(paste0(as.character(rep(GenomeInfoDb::seqnames(vcf)@values, GenomeInfoDb::seqnames(vcf)@lengths)), ".", BiocGenerics::start(vcf), "_", as.character(VariantAnnotation::ref(vcf)), ".", as.character(unlist(VariantAnnotation::alt(vcf)))))
    # new_names <- new_names[names(new_names) %in% rownames(ref_matrix_total)]
    new_names <- new_names[rownames(ref_matrix_total)]
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
  cell_ids <- colnames(consensus_matrix_total)
  variant_names <- names(keep_variants[keep_variants])
  consensus_matrix_total <- consensus_matrix_total[keep_variants, , drop = FALSE]
  coverage_matrix_total <- coverage_matrix_total[keep_variants, , drop = FALSE]
  ref_matrix_total <- ref_matrix_total[keep_variants, , drop = FALSE]


  if(verbose) print("We remove cells that are always NoCall.")
  consensus_test <- consensus_matrix_total > 0
  keep_cells <- Matrix::colSums(consensus_test) > 0
  # If we only have one cell or one variant, we loose the matrix.
  cell_ids <- names(keep_cells[keep_cells])
  variant_names <- rownames(consensus_matrix_total)
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
    if(verbose) print(paste0("The filtering left ", dim_test[1], " variants and ", dim_test[2], "cells."))
    if(verbose) print("Returning NULL.")
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

