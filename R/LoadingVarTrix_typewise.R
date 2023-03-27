#'LoadingVarTrix_typewise
#'@description
#'When we load all the different types of results (scRNAseq/amplicon and MT/amplicon),
#'we might need extreme amounts of memory. To solve this issue, I will load each type separately.
#'In a following function (AmpliconSupplementing), we can add the amplicon information to the
#'scRNAseq information.
#'The input file is a specifically formated csv file with all the necessary information to run the analysis.
#'Note that the source column in the input file needs to be one of the following: vartrix, mgaetk, mgatk.
#'This is hard coded and case insensitive.
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
LoadingVarTrix_typewise <- function(samples_file, samples_path = NULL, barcodes_path = NULL, snp_path = NULL, vcf_path, patient, sample = NULL, type_use = "scRNAseq_Somatic", min_reads = NULL, min_cells = 2){
  if(all(!is.null(samples_path), !is.null(barcodes_path), !is.null(sample), !is.null(snp_path))){
    #samples <- list.files(samples_path)
    #samples <- grep(patient, samples, value = TRUE)
    
    #barcodes_files <- list.files(path = samples_path, pattern = "barcodes")
    #barcodes_files <- unlist(lapply(paste0(samples_path, samples, "/"), list.files, pattern = "barcodes", full.names = TRUE))
    

    #samples_file <- data.frame(patient = patient, sample = samples, input_folder = samples_path, cells = barcodes_files)
    samples_file <- data.frame(patient = patient, sample = sample, input_folder = samples_path, cells = barcodes_path)
    samples <- samples_file$sample
  } else{
    print("We read in the samples file.")
    samples_file <- read.csv(samples_file, stringsAsFactors = FALSE)
    
    
    print("We subset to the patient of interest.")
    samples_file <- samples_file[grep("vartrix", samples_file$source, ignore.case = TRUE),]
    samples_file <- samples_file[samples_file$patient == patient,]
    samples_file <- samples_file[samples_file$type == type_use,]
    
    
    print("We get the different samples.")
    samples <- samples_file$sample
  }


  print("We load the SNV files.")
  if(!is.null(snp_path)){
    path_snps <- snp_path
  } else{
    path_snps <- paste0(samples_file$input_folder, "/SNV.loci.txt")
  }


  print("We read the variants.")
  snps_list <- lapply(path_snps, read.table, header = FALSE)
  names(snps_list) <- samples


  print("We read in the cell barcodes output by CellRanger as a list.")
  barcodes <- lapply(samples_file$cells, read.table)
  names(barcodes) <- samples


  print("We read in the vcf file.")
  vcf         <- readVcf(vcf_path)
  vcf_info    <- info(vcf)


  print("We generate more accessible names.")
  if(all(c("GENE", "AA", "CDS") %in% colnames(vcf_info))){
    new_names <- paste0(vcf_info$GENE, "_", vcf_info$AA, "_", vcf_info$CDS)
  } else{
    new_names <- rownames(vcf_info)
    new_names <- gsub(":|\\/|\\?", "_", new_names)
  }

  print("We read in the different sparse genotype matrices as a list.")
  print("We have a slot per type of input data.")
  coverage_matrices <- list()
  ref_matrices <- list()
  consensus_matrices <- list()
  for(i in 1:nrow(samples_file)){
    print(paste0("Loading sample ", i, " of ", nrow(samples_file)))
    input_folder_use <- samples_file$input_folder[i]
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


  print("We generate a large data.frame of all the snv matrices.")
  coverage_matrix_total  <- do.call("cbind", coverage_matrices)
  ref_matrix_total       <- do.call("cbind", ref_matrices)
  consensus_matrix_total <- do.call("cbind", consensus_matrices)


  print("We remove the matrix lists.")
  rm(coverage_matrices, ref_matrices, consensus_matrices)
  gc()


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
  keep_variants <- rowSums(consensus_matrix_total >= 2)
  keep_variants <- keep_variants >= min_cells
  # If we only have 1 variant left, we would loose the matrix.
  if(sum(keep_variants) > 1){
    consensus_matrix_total <- consensus_matrix_total[keep_variants,]
    coverage_matrix_total <- coverage_matrix_total[keep_variants,]
    ref_matrix_total <- ref_matrix_total[keep_variants,]
  } else if(sum(keep_variants) == 1){
    cell_ids <- colnames(consensus_matrix_total)
    variant_names <- names(keep_variants[keep_variants])
    consensus_matrix_total <- consensus_matrix_total[keep_variants,]
    consensus_matrix_total <- matrix(consensus_matrix_total, nrow = 1, ncol = length(consensus_matrix_total))
    rownames(consensus_matrix_total) <- variant_names
    colnames(consensus_matrix_total) <- cell_ids
    consensus_matrix_total <- as(consensus_matrix_total, "dgCMatrix")
    coverage_matrix_total <- coverage_matrix_total[keep_variants,]
    coverage_matrix_total <- matrix(coverage_matrix_total, nrow = 1, ncol = length(coverage_matrix_total))
    rownames(coverage_matrix_total) <- variant_names
    colnames(coverage_matrix_total) <- cell_ids
    coverage_matrix_total <- as(coverage_matrix_total, "dgCMatrix")
    ref_matrix_total <- ref_matrix_total[keep_variants,]
    ref_matrix_total <- matrix(ref_matrix_total, nrow = 1, ncol = length(ref_matrix_total))
    rownames(ref_matrix_total) <- variant_names
    colnames(ref_matrix_total) <- cell_ids
    ref_matrix_total <- as(ref_matrix_total, "dgCMatrix")
    rm(cell_ids, variant_names)
  }


  print("We remove cells that are always NoCall.")
  consensus_test <- consensus_matrix_total > 0
  keep_cells <- colSums(consensus_test) > 0
  # If we only have one cell, we loose the matrix.
  if(sum(keep_cells) > 1){
    consensus_matrix_total <- consensus_matrix_total[, keep_cells]
    coverage_matrix_total <- coverage_matrix_total[, keep_cells][]
    ref_matrix_total <- ref_matrix_total[, keep_cells]
  } else if(sum(keep_cells) == 1){
    cell_ids <- names(keep_cells[keep_cells])
    variant_names <- rownames(consensus_matrix_total)
    consensus_matrix_total <- consensus_matrix_total[,keep_cells]
    consensus_matrix_total <- matrix(consensus_matrix_total, nrow = length(variant_names), ncol = 1)
    rownames(consensus_matrix_total) <- variant_names
    colnames(consensus_matrix_total) <- cell_ids
    consensus_matrix_total <- as(consensus_matrix_total, "dgCMatrix")
    coverage_matrix_total <- coverage_matrix_total[,keep_cells]
    coverage_matrix_total <- matrix(coverage_matrix_total, nrow = length(variant_names), ncol = 1)
    rownames(coverage_matrix_total) <- variant_names
    colnames(coverage_matrix_total) <- cell_ids
    coverage_matrix_total <- as(coverage_matrix_total, "dgCMatrix")
    ref_matrix_total <- ref_matrix_total[, keep_cells]
    ref_matrix_total <- matrix(ref_matrix_total, nrow = length(variant_names), ncol = 1)
    rownames(ref_matrix_total) <- variant_names
    colnames(ref_matrix_total) <- cell_ids
    ref_matrix_total <- as(ref_matrix_total, "dgCMatrix")
    rm(cell_ids, variant_names)
  }

  print(paste0(type_use, " Variants: ", nrow(consensus_matrix_total)))
  print(paste0(type_use, " Cells: ", ncol(consensus_matrix_total)))

  rm(consensus_test, keep_variants, keep_cells)
  gc()


  print("We transform the sparse matrices to matrices, so we can calculate the fraction.")
  # For test purposes
  #coverage_matrix_total_ori <- coverage_matrix_total
  #ref_matrix_total_ori <- ref_matrix_total
  #consensus_matrix_total_ori <- consensus_matrix_total
  #coverage_matrix_total <- coverage_matrix_total_ori
  #ref_matrix_total <- ref_matrix_total_ori
  #consensus_matrix_total <- consensus_matrix_total_ori

  #coverage_matrix_total <- coverage_matrix_total[1:5000,1:5000]
  #ref_matrix_total <- ref_matrix_total[1:5000,1:5000]
  #consensus_matrix_total <- coverage_matrix_total[1:5000,1:5000]
  #colnames(coverage_matrix_total) <- make.names(colnames(coverage_matrix_total))
  #rownames(coverage_matrix_total) <- make.names(rownames(coverage_matrix_total))
  #colnames(ref_matrix_total) <- make.names(colnames(ref_matrix_total))
  #rownames(ref_matrix_total) <- make.names(rownames(ref_matrix_total))
  #colnames(consensus_matrix_total) <- make.names(colnames(consensus_matrix_total))
  #rownames(consensus_matrix_total) <- make.names(rownames(consensus_matrix_total))
  coverage_matrix_total                             <- as.matrix(coverage_matrix_total)
  ref_matrix_total                                  <- as.matrix(ref_matrix_total)
  consensus_matrix_total                            <- as.matrix(consensus_matrix_total)
  reads_total                                       <- coverage_matrix_total + ref_matrix_total
  fraction_total                                    <- coverage_matrix_total / reads_total
  fraction_total[is.na(fraction_total)] <- 0
  #fraction_total <- sdiv(X = coverage_matrix_total, Y = reads_total,
  #                       names = dimnames(coverage_matrix_total))
  #rm(coverage_matrix_total, ref_matrix_total)
  gc()


  print("We generate a SummarizedExperiment object containing the fraction and the consensus matrices.")
  # We want an assay for the Consensus information and for the fraction.
  # As meta data we add a data frame showing the cell id, the associated patient and the sample.
  coverage_depth_per_cell <- colMeans(reads_total)
  meta_data <- data.frame(Cell = colnames(consensus_matrix_total), Type = type_use, AverageCoverage = coverage_depth_per_cell)
  #se_merged <- SummarizedExperiment(assays = list(consensus = as(consensus_matrix_total, "dgCMatrix"), fraction = as(fraction_total, "dgCMatrix"), coverage = as(reads_total, "dgCMatrix")),
  #                                  colData = meta_data)
  se_merged <- SummarizedExperiment(assays = list(consensus = as(consensus_matrix_total, "CsparseMatrix"), fraction = as(fraction_total, "CsparseMatrix"), coverage = as(reads_total, "CsparseMatrix"),
                                                  alts = as(coverage_matrix_total, "CsparseMatrix"), refs = as(ref_matrix_total, "CsparseMatrix")),
                                    colData = meta_data)
  return(se_merged)
}

