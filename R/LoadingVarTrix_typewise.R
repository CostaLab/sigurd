#'When we load all the different types of results (scRNAseq/amplicon and MT/amplicon),
#'we might need extreme amounts of memory. To solve this issue, I will load each type separately.
#'In a following function (AmpliconSupplementing), we can add the amplicon information to the
#'scRNAseq information.
#'@import Matrix SummarizedExperiment VariantAnnotation
#'@param samples_path Path to the csv file with the samples to be loaded.
#'@param vcf_path Path to the VCF file with the variants.
#'@param patient The patient you want to load.
#'@param type_use The type of input. Has to be one of: scRNAseq_Somatic, Amplicon_Somatic, scRNAseq_MT, Amplicon_MT.
#'@export
LoadingVarTrix_typewise <- function(samples_path, vcf_path, patient, type_use = "scRNAseq_Somatic", chromosome_prefix = "chrM"){
  library(Matrix)
  library(SummarizedExperiment)
  library(VariantAnnotation)
  print("TEST")
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


  print(type_use)
  if(type_use %in% c("scRNAseq_MT", "Amplicon_MT")){
    # Since we have multiple mitochondria per cell, we need to take all variants per position into account.
    print("We load the coverage (alt) and ref matrices.")
    coverage_matrices <- list()
    ref_matrices <- list()
    
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

      
      # We rename the rows and columns.
      colnames(coverage_matrices[[sample_use]])  <- paste0(sample_use, "_", cellbarcodes_use$V1)
      colnames(ref_matrices[[sample_use]])       <- paste0(sample_use, "_", cellbarcodes_use$V1)
      rownames(coverage_matrices[[sample_use]])  <- variants_use$V1
      rownames(ref_matrices[[sample_use]])       <- variants_use$V1
    }


    print("We generate a large data.frame of all the snv matrices.")
    coverage_matrix_total  <- do.call("cbind", coverage_matrices)
    coverage_matrix_total  <- as.matrix(coverage_matrix_total)
    ref_matrix_total       <- do.call("cbind", ref_matrices)
    ref_matrix_total       <- as.matrix(ref_matrix_total)


    print("We remove the matrix lists.")
    rm(coverage_matrices, ref_matrices)
    gc()

    print("We change the names.")
    rownames(coverage_matrix_total)  <- new_names
    rownames(ref_matrix_total)       <- new_names


    print("We get the positions.")
    positions <- start(vcf)
    positions <- unique(positions)
    ref_alleles <- data.frame(as.numeric(start(vcf)), gsub(":", "_", rownames(vcf_info)))
    ref_alleles <- ref_alleles[order(ref_alleles[,1]),][,2]
    ref_alleles <- gsub("/.*", "", ref_alleles)
    ref_alleles <- unique(ref_alleles)


    print("We split the alternative matrix by alternative allele.")
    # Now, we have a list of matrices, where each cell is the number of alternative reads for this position.
    coverage_matrix_total_list <- list(A = coverage_matrix_total[grep("_A$", rownames(coverage_matrix_total), value = TRUE),],
                                       C = coverage_matrix_total[grep("_C$", rownames(coverage_matrix_total), value = TRUE),],
                                       G = coverage_matrix_total[grep("_G$", rownames(coverage_matrix_total), value = TRUE),],
                                       T = coverage_matrix_total[grep("_T$", rownames(coverage_matrix_total), value = TRUE),])

    # We generate a matrix with all positions 
    alt_read_matrix <- matrix(data = 0, nrow = length(positions), ncol = ncol(coverage_matrix_total))
    rownames(alt_read_matrix) <- unique(gsub("_.$", "", rownames(coverage_matrix_total)))
    colnames(alt_read_matrix) <- colnames(coverage_matrix_total)
    #alt_read_matrix <- sparseMatrix(i = length(positions), j = ncol(coverage_matrix_total), x = 0,
    #                                dimnames = list(unique(gsub("_.$", "", rownames(coverage_matrix_total))), colnames(coverage_matrix_total)))


    alt_read_matrix[gsub("_A$", "", rownames(coverage_matrix_total_list[["A"]])),] <- alt_read_matrix[gsub("_A$", "", rownames(coverage_matrix_total_list[["A"]])),] + coverage_matrix_total_list[["A"]]
    alt_read_matrix[gsub("_C$", "", rownames(coverage_matrix_total_list[["C"]])),] <- alt_read_matrix[gsub("_C$", "", rownames(coverage_matrix_total_list[["C"]])),] + coverage_matrix_total_list[["C"]]
    alt_read_matrix[gsub("_G$", "", rownames(coverage_matrix_total_list[["G"]])),] <- alt_read_matrix[gsub("_G$", "", rownames(coverage_matrix_total_list[["G"]])),] + coverage_matrix_total_list[["G"]]
    alt_read_matrix[gsub("_T$", "", rownames(coverage_matrix_total_list[["T"]])),] <- alt_read_matrix[gsub("_T$", "", rownames(coverage_matrix_total_list[["T"]])),] + coverage_matrix_total_list[["T"]]


    # We get the new reference matrix.
    new_ref_matrix_total <- ref_matrix_total

    # We get the read matrix for the >A variants.
    #time1_start <- Sys.time()
    read_matrix_subset <- alt_read_matrix
    rownames(read_matrix_subset) <- paste0(rownames(read_matrix_subset), "_A")
    rownames_subset <- rownames(read_matrix_subset)[rownames(read_matrix_subset) %in% rownames(new_ref_matrix_total)]
    read_matrix_subset <- read_matrix_subset[rownames_subset,]
    new_ref_matrix_total[rownames_subset, ] <- new_ref_matrix_total[rownames_subset, ] + read_matrix_subset
    #time1_stop <- Sys.time()
    #time1_stop - time1_start


    #time2_start <- Sys.time()
    read_matrix_subset <- alt_read_matrix
    rownames(read_matrix_subset) <- paste0(rownames(read_matrix_subset), "_C")
    rownames_subset <- rownames(read_matrix_subset)[rownames(read_matrix_subset) %in% rownames(new_ref_matrix_total)]
    read_matrix_subset <- read_matrix_subset[rownames_subset,]
    new_ref_matrix_total[rownames_subset, ] <- new_ref_matrix_total[rownames_subset, ] + read_matrix_subset
    #time2_stop <- Sys.time()
    #time2_stop - time2_start

    
    #time3_start <- Sys.time()
    read_matrix_subset <- alt_read_matrix
    rownames(read_matrix_subset) <- paste0(rownames(read_matrix_subset), "_G")
    rownames_subset <- rownames(read_matrix_subset)[rownames(read_matrix_subset) %in% rownames(new_ref_matrix_total)]
    read_matrix_subset <- read_matrix_subset[rownames_subset,]
    new_ref_matrix_total[rownames_subset, ] <- new_ref_matrix_total[rownames_subset, ] + read_matrix_subset
    #time3_stop <- Sys.time()
    #time3_stop - time3_start

  
    #time4_start <- Sys.time()
    read_matrix_subset <- alt_read_matrix
    rownames(read_matrix_subset) <- paste0(rownames(read_matrix_subset), "_T")
    rownames_subset <- rownames(read_matrix_subset)[rownames(read_matrix_subset) %in% rownames(new_ref_matrix_total)]
    read_matrix_subset <- read_matrix_subset[rownames_subset,]
    new_ref_matrix_total[rownames_subset, ] <- new_ref_matrix_total[rownames_subset, ] + read_matrix_subset
    #time4_stop <- Sys.time()
    #time4_stop - time4_start

    rm(read_matrix_subset)
    rm(alt_read_matrix)
    gc()


    print("We calculate the fraction.")
    reads_total <- coverage_matrix_total + new_ref_matrix_total
    fraction_total <- coverage_matrix_total / reads_total
    fraction_total[is.na(fraction_total)] <- 0
    rm(new_ref_matrix_total)
    gc()


    #print("We get a matrix with the reference reads per position.")
    #ref_matrix <- matrix(data = 0, nrow = length(positions), ncol = ncol(coverage_matrix_total))
    #rownames(ref_matrix) <- unique(gsub("_.$", "", rownames(coverage_matrix_total)))
    #colnames(ref_matrix) <- colnames(coverage_matrix_total)
    #read_matrix <- ref_matrix

    print("We split the reference matrix by reference allele.")
    ref_matrix_total_list <- list(A = ref_matrix_total[grep("_A_T", rownames(ref_matrix_total), value = TRUE),],
                                  C = ref_matrix_total[grep("_C_A", rownames(ref_matrix_total), value = TRUE),],
                                  G = ref_matrix_total[grep("_G_A", rownames(ref_matrix_total), value = TRUE),],
                                  T = ref_matrix_total[grep("_T_A", rownames(ref_matrix_total), value = TRUE),])#,
#                                  N = matrix(ref_matrix_total[grep("_N_A", rownames(ref_matrix_total), value = TRUE),], ncol = ncol(ref_matrix_total)))
    rownames(ref_matrix_total_list[["A"]]) <- gsub("_T$", "", rownames(ref_matrix_total_list[["A"]]))
    rownames(ref_matrix_total_list[["C"]]) <- gsub("_A$", "", rownames(ref_matrix_total_list[["C"]]))
    rownames(ref_matrix_total_list[["G"]]) <- gsub("_A$", "", rownames(ref_matrix_total_list[["G"]]))
    rownames(ref_matrix_total_list[["T"]]) <- gsub("_A$", "", rownames(ref_matrix_total_list[["T"]]))
    #colnames(ref_matrix_total_list[["N"]]) <- colnames(ref_matrix_total)
    #rownames(ref_matrix_total_list[["N"]]) <- paste0(chromosome_prefix, "_3107_N")


    print("We get the reads per type (A/C/G/T) per position.")
    ###
    ### Check the results.
    ###
    # VarTrix does not return the same reference values for the different variants.
    # Example:
    #                 Minus_ET_83497_AAACCCAAGCTCCATA-1 Minus_ET_83497_AAACCCAAGGAACTAT-1 Minus_ET_83497_AAACCCATCAATCTTC-1
    # chrM_16043_A_C                                  0                                 1                                 1
    # chrM_16043_A_G                                  0                                 1                                 1
    # chrM_16043_A_T                                  0                                 0                                 1
    # For now, we will ignore this issue. Hopefully this was solved in the newest version of VarTrix.
    reads_matrix <- matrix(0, nrow = length(ref_alleles), ncol = ncol(ref_matrix_total))
    rownames(reads_matrix) <- ref_alleles
    colnames(reads_matrix) <- colnames(ref_matrix_total)
    reads_matrix_ls <- list(A = reads_matrix, C = reads_matrix, G = reads_matrix, T = reads_matrix)#, N = reads_matrix)
    reads_matrix_ls[["A"]][rownames(ref_matrix_total_list[["A"]]),] <- ref_matrix_total_list[["A"]]
    reads_matrix_ls[["C"]][rownames(ref_matrix_total_list[["C"]]),] <- ref_matrix_total_list[["C"]]
    reads_matrix_ls[["G"]][rownames(ref_matrix_total_list[["G"]]),] <- ref_matrix_total_list[["G"]]
    reads_matrix_ls[["T"]][rownames(ref_matrix_total_list[["T"]]),] <- ref_matrix_total_list[["T"]]
    #reads_matrix_ls[["N"]][rownames(ref_matrix_total_list[["N"]]),] <- ref_matrix_total_list[["N"]]
    
    reads_matrix_ls[["A"]][gsub("_A$", "", rownames(coverage_matrix_total_list[["A"]])),] <- coverage_matrix_total_list[["A"]]
    reads_matrix_ls[["C"]][gsub("_C$", "", rownames(coverage_matrix_total_list[["C"]])),] <- coverage_matrix_total_list[["C"]]
    reads_matrix_ls[["G"]][gsub("_G$", "", rownames(coverage_matrix_total_list[["G"]])),] <- coverage_matrix_total_list[["G"]]
    reads_matrix_ls[["T"]][gsub("_T$", "", rownames(coverage_matrix_total_list[["T"]])),] <- coverage_matrix_total_list[["T"]]
    
    reads_matrix_ls[["A"]][reads_matrix_ls[["A"]] > 0] <- 8
    reads_matrix_ls[["C"]][reads_matrix_ls[["C"]] > 0] <- 4
    reads_matrix_ls[["G"]][reads_matrix_ls[["G"]] > 0] <- 2
    reads_matrix_ls[["T"]][reads_matrix_ls[["T"]] > 0] <- 1
    reads_matrix <- reads_matrix_ls[["A"]] + reads_matrix_ls[["C"]] + reads_matrix_ls[["G"]] + reads_matrix_ls[["T"]]
    rm(reads_matrix_ls)
    gc()
    reads_matrix_ls <- list(A = reads_matrix[grep("_A$", rownames(reads_matrix)),],
                            C = reads_matrix[grep("_C$", rownames(reads_matrix)),],
                            G = reads_matrix[grep("_G$", rownames(reads_matrix)),],
                            T = reads_matrix[grep("_T$", rownames(reads_matrix)),],
                            N = reads_matrix[grep("_N$", rownames(reads_matrix)),])
    reads_matrix_ls[["N"]] <- matrix(reads_matrix_ls[["N"]], nrow = 1, ncol = length(reads_matrix_ls[["N"]]))
    colnames(reads_matrix_ls[["N"]]) <- colnames(ref_matrix_total)
    rownames(reads_matrix_ls[["N"]]) <- paste0(chromosome_prefix, "_3107_N")
    rm(reads_matrix)
    gc()


    consensus_a <- lapply(c("C", "G", "T"), get_consensus, ref_base = "A", input_matrix = reads_matrix_ls[["A"]])
    consensus_a <- do.call("rbind", consensus_a)
    consensus_c <- lapply(c("A", "G", "T"), get_consensus, ref_base = "C", input_matrix = reads_matrix_ls[["C"]])
    consensus_c <- do.call("rbind", consensus_c)
    consensus_g <- lapply(c("A", "C", "T"), get_consensus, ref_base = "G", input_matrix = reads_matrix_ls[["G"]])
    consensus_g <- do.call("rbind", consensus_g)
    consensus_t <- lapply(c("A", "C", "G"), get_consensus, ref_base = "T", input_matrix = reads_matrix_ls[["T"]])
    consensus_t <- do.call("rbind", consensus_t)
    consensus_n <- lapply(c("A", "C", "G", "T"), get_consensus, ref_base = "N", input_matrix = reads_matrix_ls[["N"]])
    consensus_n <- do.call("rbind", consensus_n)
    consensus <- rbind(consensus_a, consensus_c, consensus_g, consensus_t, consensus_n)
    rm(consensus_a, consensus_c, consensus_g, consensus_t, consensus_n)
    gc()
    consensus <- consensus[match(rownames(ref_matrix_total), rownames(consensus)),]


    print("We generate a SummarizedExperiment object containing the fraction and the consensus matrices.")
    # We want an assay for the Consensus information and for the fraction.
    # As meta data we add a data frame showing the cell id, the associated patient and the sample.
    coverage_depth_per_cell <- colMeans(reads_total)
    meta_data <- data.frame(Cell = colnames(consensus), Type = type_use, AverageCoverage = coverage_depth_per_cell)
    se_merged <- SummarizedExperiment(assays = list(consensus = consensus, fraction = fraction_total, coverage = reads_total),
                                      colData = meta_data)
    return(se_merged)

  } else if(type_use %in% c("scRNAseq_somatic", "Amplicon_somatic.")){
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


    rownames(coverage_matrix_total)  <- new_names
    rownames(ref_matrix_total)       <- new_names
    rownames(consensus_matrix_total) <- new_names


    # We remove variants, that are not detected in at least 2 cells.
    keep_variants <- rowSums(consensus_matrix_total >= 2)
    keep_variants <- keep_variants >= 2
    consensus_matrix_total <- consensus_matrix_total[keep_variants,]
    coverage_matrix_total <- coverage_matrix_total[keep_variants,]
    ref_matrix_total <- ref_matrix_total[keep_variants,]


    # We remove cells that are always NoCall.
    consensus_test <- consensus_matrix_total > 0
    keep_cells <- colSums(consensus_test) > 0
    consensus_matrix_total <- consensus_matrix_total[, keep_cells]
    coverage_matrix_total <- coverage_matrix_total[, keep_cells]
    ref_matrix_total <- ref_matrix_total[, keep_cells]

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
}
