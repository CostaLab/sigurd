#'LoadingVCF_typewise
#'@description
#'We load a cellwise pileup result from a VCF file.
#'
#'@import Matrix SummarizedExperiment VariantAnnotation
#'@param samples_path Path to the input folder. Must include a barcodes file.
#'@param samples_file Path to the csv file with the samples to be loaded.
#'@param vcf_path Path to the VCF file with the variants.
#'@param patient The patient you want to load.
#'@param type_use The type of input. Has to be one of: scRNAseq_Somatic, Amplicon_Somatic, scRNAseq_MT, Amplicon_MT.
#'@param min_reads The minimum number of reads we want. Otherwise we treat this as a NoCall. Default = NULL.
#'@param min_cells The minimum number of cells for a variant. Otherwise, we will remove a variant. Default = 2.
#'@export
LoadingVCF_typewise <- function(samples_file, samples_path = NULL, barcodes_path = NULL, vcf_path, patient, sample = NULL, type_use = "scRNAseq_Somatic", min_reads = NULL, min_cells = 2){
  if(all(!is.null(samples_path), !is.null(barcodes_path), !is.null(sample))){
    print(paste0("Loading the data for sample ", sample, "."))
    samples_file <- data.frame(patient = patient, sample = sample, input_folder = samples_path, cells = barcodes_path)
    samples <- samples_file$sample
  } else{
    print(paste0("Loading the data for patient ", patient, "."))
    print("We read in the samples file.")
    samples_file <- read.csv(samples_file, stringsAsFactors = FALSE)


    print("We subset to the patient of interest.")
    samples_file <- samples_file[grep("vcf", samples_file$source, ignore.case = TRUE),]
    samples_file <- samples_file[samples_file$patient == patient,]
    samples_file <- samples_file[samples_file$type == type_use,]


    print("We get the different samples.")
    samples <- samples_file$sample
  }


  # A function to convert the heterozyguous/homozyguous information from the VCF to the consensus information from VarTrix.
  # We define this function here, since we will not use it anywhere else.
  # Maybe move it to a separate function.
  char_to_numeric <- function(char_value) {
    if(char_value == "1/1") return(2)
    if(char_value %in% c("1/0", "0/1")) return(2)
    if(char_value == "0/0") return(1)
    return(0)
  }


  print("We read in the cell barcodes output by CellRanger as a list.")
  barcodes <- lapply(samples_file$cells, read.table)
  names(barcodes) <- samples


  print("We read in the vcf file.")
  vcf                <- readVcf(vcf_path)
  vcf_info           <- info(vcf)


  print("We load the VCF file.")
  reads_matrix_total <- c() # The total number of reads
  coverage_matrix_total <- c() # The alternative reads
  ref_matrix_total <- c() # The reference reads
  consensus_matrix_total <- c()
  for(i in 1:length(samples)){
    print(paste0("Loading sample ", i, " of ", nrow(samples_file)))
    input_folder_use <- samples_file$input_folder[i]
    sample_use <- samples_file$sample[i]

    # The cell barcodes and variants.
    cellbarcodes_use <- barcodes[[sample_use]]

    # We load the VCF file.
    vcf_data <- paste0(input_folder_use, sample_use, "/cellSNP.cells.sorted.vcf.gz")
    depth_to_add                      <- VariantAnnotation::readGeno(vcf_data, "DP")
    depth_to_add[is.na(depth_to_add)] <- 0
    rownames(depth_to_add)            <- make.names(rownames(depth_to_add))
    colnames(depth_to_add)            <- paste0(sample_use, "_", colnames(depth_to_add))
    depth_to_add                      <- as(depth_to_add, "sparseMatrix")
    reads_matrix_total                <- cbind(reads_matrix_total, depth_to_add)

    alts_to_add                       <- readGeno(vcf_data, "AD")
    alts_to_add[is.na(alts_to_add)]   <- 0
    rownames(alts_to_add)             <- make.names(rownames(alts_to_add))
    colnames(alts_to_add)             <- paste0(sample_use, "_", colnames(alts_to_add))
    alts_to_add                       <- as(alts_to_add, "sparseMatrix")
    coverage_matrix_total             <- cbind(coverage_matrix_total, alts_to_add)

    consensus_to_add                  <- readGeno(vcf_data, "GT")
    consensus_to_add                  <- matrix(sapply(consensus_to_add, char_to_numeric), nrow = nrow(consensus_to_add), dimnames = list(make.names(rownames(consensus_to_add)), paste0(sample_use, "_", colnames(consensus_to_add))))
    consensus_to_add                  <- as(consensus_to_add, "sparseMatrix")
    consensus_matrix_total            <- cbind(consensus_matrix_total, consensus_to_add)
  }
  ref_matrix_total <- reads_matrix_total - coverage_matrix_total
  rm(consensus_to_add, alts_to_add, depth_to_add)


  print("We generate more accessible names.")
  if(all(c("GENE", "AA", "CDS") %in% colnames(vcf_info))){
    new_names <- paste0(vcf_info$GENE, "_", vcf_info$AA, "_", vcf_info$CDS)
    names(new_names) <- make.names(paste0(as.character(rep(seqnames(vcf)@values, seqnames(vcf)@lengths)), ".", start(vcf), "_", as.character(ref(vcf)), ".", as.character(unlist(alt(vcf)))))
    new_names <- new_names[rownames(ref_matrix_total)]
  } else{
    new_names <- rownames(vcf_info)
    new_names <- gsub(":|\\/|\\?", "_", new_names)
    names(new_names) <-	make.names(paste0(as.character(rep(seqnames(vcf)@values, seqnames(vcf)@lengths)), ".", start(vcf), "_", as.character(ref(vcf)), ".", as.character(unlist(alt(vcf)))))
    new_names <- new_names[rownames(ref_matrix_total)]
  }


  if(!is.null(min_reads)){
    print(paste0("We set read values below the threshold of ", min_reads, " to 0."))
    print("We then generate the consensus matrix again.")
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


  print(paste0("We remove variants, that are not detected in at least ", min_cells, " cells."))
  keep_variants <- rowSums(consensus_matrix_total >= 1)
  keep_variants <- keep_variants >= min_cells
  # If we only have one cell or one variant, we loose the matrix.
  cell_ids <- colnames(consensus_matrix_total)
  variant_names <- names(keep_variants[keep_variants])
  consensus_matrix_total <- consensus_matrix_total[keep_variants, ]
  consensus_matrix_total <- matrix(consensus_matrix_total, nrow = length(variant_names), ncol = length(cell_ids))
  colnames(consensus_matrix_total) <- cell_ids
  rownames(consensus_matrix_total) <- variant_names
  consensus_matrix_total <- as(consensus_matrix_total, "dgCMatrix")
  coverage_matrix_total <- coverage_matrix_total[keep_variants, ]
  coverage_matrix_total <- matrix(coverage_matrix_total, nrow = length(variant_names), ncol = length(cell_ids))
  colnames(coverage_matrix_total) <- cell_ids
  rownames(coverage_matrix_total) <- variant_names
  coverage_matrix_total <- as(coverage_matrix_total, "dgCMatrix")
  ref_matrix_total <- ref_matrix_total[keep_variants, ]
  ref_matrix_total <- matrix(ref_matrix_total, nrow = length(variant_names), ncol = length(cell_ids))
  colnames(ref_matrix_total) <- cell_ids
  rownames(ref_matrix_total) <- variant_names
  ref_matrix_total <- as(ref_matrix_total, "dgCMatrix")


  print("We remove cells that are always NoCall.")
  consensus_test <- consensus_matrix_total > 0
  keep_cells <- colSums(consensus_test) > 0
  # If we only have one cell or one variant, we loose the matrix.
  cell_ids <- names(keep_cells[keep_cells])
  variant_names <- rownames(consensus_matrix_total)
  consensus_matrix_total <- consensus_matrix_total[, keep_cells]
  consensus_matrix_total <- matrix(consensus_matrix_total, nrow = length(variant_names), ncol = length(cell_ids))
  colnames(consensus_matrix_total) <- cell_ids
  rownames(consensus_matrix_total) <- variant_names
  consensus_matrix_total <- as(consensus_matrix_total, "dgCMatrix")
  coverage_matrix_total <- coverage_matrix_total[, keep_cells]
  coverage_matrix_total <- matrix(coverage_matrix_total, nrow = length(variant_names), ncol = length(cell_ids))
  colnames(coverage_matrix_total) <- cell_ids
  rownames(coverage_matrix_total) <- variant_names
  coverage_matrix_total <- as(coverage_matrix_total, "dgCMatrix")
  ref_matrix_total <- ref_matrix_total[, keep_cells]
  ref_matrix_total <- matrix(ref_matrix_total, nrow = length(variant_names), ncol = length(cell_ids))
  colnames(ref_matrix_total) <- cell_ids
  rownames(ref_matrix_total) <- variant_names
  ref_matrix_total <- as(ref_matrix_total, "dgCMatrix")


  print(paste0(type_use, " Variants: ", nrow(consensus_matrix_total)))
  print(paste0(type_use, " Cells: ", ncol(consensus_matrix_total)))

  rm(consensus_test, keep_variants, keep_cells)
  gc()


  print("We transform the sparse matrices to matrices, so we can calculate the fraction.")
  coverage_matrix_total                             <- as.matrix(coverage_matrix_total)
  ref_matrix_total                                  <- as.matrix(ref_matrix_total)
  consensus_matrix_total                            <- as.matrix(consensus_matrix_total)
  reads_total                                       <- coverage_matrix_total + ref_matrix_total
  fraction_total                                    <- coverage_matrix_total / reads_total
  fraction_total[is.na(fraction_total)] <- 0
  gc()


  # We check if the matrices are empty (0 cells, 0 variants). Then we simply return NULL.
  dim_test <- dim(reads_total)
  if(any(dim_test == 0)){
    print(paste0("The filtering left ", dim_test[1], " variants and ", dim_test[2], "cells."))
    print("Returning NULL.")
    return(NULL)
  } else {
    print("We generate a SummarizedExperiment object containing the fraction and the consensus matrices.")
    # We want an assay for the Consensus information and for the fraction.
    # As meta data we add a data frame showing the cell id, the associated patient and the sample.
    coverage_depth_per_cell <- colMeans(reads_total)
    coverage_depth_per_variant <- rowMeans(reads_total)
    meta_data <- data.frame(Cell = colnames(consensus_matrix_total), Type = type_use, AverageCoverage = coverage_depth_per_cell)
    rownames(meta_data) <- meta_data$Cell
    meta_row <- data.frame(VariantName = rownames(consensus_matrix_total), Depth = coverage_depth_per_variant)
    rownames(meta_row) <- meta_row$VariantName
    #se_merged <- SummarizedExperiment(assays = list(consensus = as(consensus_matrix_total, "dgCMatrix"), fraction = as(fraction_total, "dgCMatrix"), coverage = as(reads_total, "dgCMatrix")),
    #                                  colData = meta_data)
    se_merged <- SummarizedExperiment(assays = list(consensus = as(consensus_matrix_total, "CsparseMatrix"), fraction = as(fraction_total, "CsparseMatrix"), coverage = as(reads_total, "CsparseMatrix"),
                                                    alts = as(coverage_matrix_total, "CsparseMatrix"), refs = as(ref_matrix_total, "CsparseMatrix")),
                                      colData = meta_data, rowData = meta_row)
    return(se_merged)
  }
}
