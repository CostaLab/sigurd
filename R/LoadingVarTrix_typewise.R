# When we load all the different types of results (scRNAseq/amplicon and MT/amplicon),
# we might need extreme amounts of memory. To solve this issue, I will load each type separately.
# In a following function (AmpliconSupplementing), we can add the amplicon information to the
# scRNAseq information.
#'@import Matrix SummarizedExperiment VariantAnnotation
#'@importfrom Matrix readMM
#'@param samples_path vcf_path patient type_use 
LoadingVarTrix_typewise <- function(samples_path, vcf_path, patient, type_use = "scRNAseq_Somatic"){
  print("We read in the samples file.")
  samples_file <- read.csv(samples_path)


  print("We subset to the patient of interest.")
  samples_file <- samples_file[grep("vartrix", samples_file$source, ignore.case = TRUE),]
  samples_file <- samples_file[grep(patient, samples_file$patient),]
  samples_file <- samples_file[grep(type_use, samples_file$type),]


  print("We get the different samples.")
  samples <- samples_file$sample


  print("We load the SNV files.")
  path_snps <- paste0(samples_file$input_folder, "/SNV.loci.txt")


  print("We read the variants.")
  snps_list <- lapply(path_snps, read.table, header = FALSE)
  names(snps_list) <- samples


  print("We read in the cell barcodes output by CellRanger as a list.")
  barcodes <- lapply(samples_file$cells, read.table)
  names(barcodes) <- samples


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
    coverage_matrices[[sample_use]]  <- readMM(paste0(input_folder_use, sample_use, "/out_matrix_coverage.mtx"))
    ref_matrices[[sample_use]]       <- readMM(paste0(input_folder_use, sample_use, "/ref_matrix_coverage.mtx"))
    consensus_matrices[[sample_use]] <- readMM(paste0(input_folder_use, sample_use, "/out_matrix_consensus.mtx"))

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


  print("We read in the vcf file so we can subset it together with the total matrices.")
  vcf         <- readVcf(vcf_path)
  vcf_info    <- info(vcf)


  print("We change the names of the mutations, to make it more accessible.")
  if(all(c("GENE", "AA", "CDS") %in% colnames(vcf_info))){
    new_names <- paste0(vcf_info$GENE, "_", vcf_info$AA, "_", vcf_info$CDS)
  } else{
    new_names <- rownames(vcf_info)
    new_names <- gsub(":|\\/|\\?", "_", new_names)
  }
  print(paste0(type_use, " Variants: ", nrow(consensus_matrix_total)))
  print(paste0(type_use, " Cells: ", ncol(consensus_matrix_total)))

  rownames(coverage_matrix_total)  <- new_names
  rownames(ref_matrix_total)       <- new_names
  rownames(consensus_matrix_total) <- new_names

  keep_variants <- rowSums(coverage_matrix_total >= 2)
  keep_variants <- keep_variants >= 2
  consensus_matrix_total <- consensus_matrix_total[keep_variants,]
  coverage_matrix_total <- coverage_matrix_total[keep_variants,]
  ref_matrix_total <- ref_matrix_total[keep_variants,]

  consensus_test <- consensus_matrix_total > 0
  keep_cells <- colSums(consensus_test) > 0
  consensus_matrix_total <- consensus_matrix_total[,keep_cells]
  coverage_matrix_total <- coverage_matrix_total[,keep_cells]
  ref_matrix_total <- ref_matrix_total[,keep_cells]

  print(paste0(type_use, " Variants: ", nrow(consensus_matrix_total)))
  print(paste0(type_use, " Cells: ", ncol(consensus_matrix_total)))

  rm(consensus_test, keep_variants, keep_cells)
  gc()


  print("We transform the sparse matrices to matrices, so we can calculate the fraction.")
  #coverage_matrix_total_ori <- coverage_matrix_total
  #ref_matrix_total_ori <- ref_matrix_total
  #consensus_matrix_total_ori <- consensus_matrix_total
  #coverage_matrix_total <- coverage_matrix_total_ori
  #ref_matrix_total <- ref_matrix_total_ori
  #consensus_matrix_total <- consensus_matrix_total_ori

  #coverage_matrix_total <- coverage_matrix_total[1:5000,1:5000]
  #ref_matrix_total <- ref_matrix_total[1:5000,1:5000]
  #consensus_matrix_total <- coverage_matrix_total[1:5000,1:5000]
#  colnames(coverage_matrix_total) <- make.names(colnames(coverage_matrix_total))
#  rownames(coverage_matrix_total) <- make.names(rownames(coverage_matrix_total))
#  colnames(ref_matrix_total) <- make.names(colnames(ref_matrix_total))
#  rownames(ref_matrix_total) <- make.names(rownames(ref_matrix_total))
#  colnames(consensus_matrix_total) <- make.names(colnames(consensus_matrix_total))
#  rownames(consensus_matrix_total) <- make.names(rownames(consensus_matrix_total))
  coverage_matrix_total                             <- as.matrix(coverage_matrix_total)
  ref_matrix_total                                  <- as.matrix(ref_matrix_total)
  consensus_matrix_total                            <- as.matrix(consensus_matrix_total)
  reads_total                                       <- coverage_matrix_total + ref_matrix_total
  fraction_total                                    <- coverage_matrix_total / reads_total
  fraction_total[is.na(fraction_total)] <- 0
#  fraction_total <- sdiv(X = coverage_matrix_total, Y = reads_total,
#                         names = dimnames(coverage_matrix_total))
  rm(coverage_matrix_total, ref_matrix_total)
  gc()


  print("We generate a SummarizedExperiment object containing the fraction and the consensus matrices.")
  # We want an assay for the Consensus information and for the fraction.
  # As meta data we add a data frame showing the cell id, the associated patient and the sample.
  coverage_depth_per_cell <- colMeans(reads_total)
  meta_data <- data.frame(Cell = colnames(consensus_matrix_total), Type = type_use, AverageCoverage = coverage_depth_per_cell)
  se_merged <- SummarizedExperiment(assays = list(consensus = consensus_matrix_total, fraction = fraction_total, coverage = reads_total),
                                    colData = meta_data)
  return(se_merged)
}
