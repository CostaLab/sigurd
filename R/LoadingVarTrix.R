LoadingVarTrix <- function(samples_path, vcf_path, vcf_path_MT, patient){
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(VariantAnnotation))
  suppressPackageStartupMessages(library(archive))
  suppressPackageStartupMessages(library(SummarizedExperiment))

  print("We read in the samples file.")
  samples_file <- read.csv(samples_path)


  print("We subset to the patient of interest.")
  samples_file <- samples_file[grep("vartrix", samples_file$source, ignore.case = TRUE),]
  samples_file <- samples_file[grep(patient, samples_file$patient),]


  print("We get the different samples.")
  samples <- samples_file$sample
  samples_unique <- unique(samples)
  types <- samples_file$type
  types_unique <- unique(types)


  print("We load the SNV files.")
  path_snps <- paste0(samples_file$input_folder, "/SNV.loci.txt")


  print("We read the variants.")
  snps_list <- lapply(path_snps, read.table, header = FALSE)
  names(snps_list) <- paste0(types, "_", samples)


  print("We read in the cell barcodes output by CellRanger as a list.")
  barcodes <- lapply(samples_file$cells, read.table)
  names(barcodes) <- paste0(types, "_", samples)


  print("We read in the different sparse genotype matrices as a list.")
  print("We have a slot per type of input data.")
  coverage_matrices <- list()
  ref_matrices <- list()
  consensus_matrices <- list()
  for(i in 1:nrow(samples_file)){
    print(paste0("Loading sample ", i, " of ", nrow(samples_file)))
    type_use <- samples_file$type[i]
    input_folder_use <- samples_file$input_folder[i]
    sample_use <- samples_file$sample[i]

    # The cell barcodes and variants.
    cellbarcodes_use <- barcodes[[paste0(type_use, "_", sample_use)]]
    variants_use     <- snps_list[[paste0(type_use, "_", sample_use)]]

    # We load the different matrices.
    coverage_matrices[[type_use]][[sample_use]]  <- readMM(paste0(input_folder_use, sample_use, "/out_matrix_coverage.mtx"))
    ref_matrices[[type_use]][[sample_use]]       <- readMM(paste0(input_folder_use, sample_use, "/ref_matrix_coverage.mtx"))
    consensus_matrices[[type_use]][[sample_use]] <- readMM(paste0(input_folder_use, sample_use, "/out_matrix_consensus.mtx"))

    # We rename the rows and columns.
    colnames(coverage_matrices[[type_use]][[sample_use]])  <- paste0(sample_use, "_", cellbarcodes_use$V1)
    colnames(ref_matrices[[type_use]][[sample_use]])       <- paste0(sample_use, "_", cellbarcodes_use$V1)
    colnames(consensus_matrices[[type_use]][[sample_use]]) <- paste0(sample_use, "_", cellbarcodes_use$V1)
    rownames(coverage_matrices[[type_use]][[sample_use]])  <- variants_use$V1
    rownames(ref_matrices[[type_use]][[sample_use]])       <- variants_use$V1
    rownames(consensus_matrices[[type_use]][[sample_use]]) <- variants_use$V1
  }


  print("We generate a large data.frame of all the snv matrices.")
  print("Since the different types are in the same list, we get all the used types and loop over them.")
  coverage_matrix_total  <- lapply(coverage_matrices, do.call, what = "cbind")
  ref_matrix_total       <- lapply(ref_matrices, do.call, what = "cbind")
  consensus_matrix_total <- lapply(consensus_matrices, do.call, what = "cbind")


  print("We remove the matrix lists.")
  rm(coverage_matrices, ref_matrices, consensus_matrices)
  gc()


  print("We read in the vcf file so we can subset it together with the total matrices.")
  vcf         <- readVcf(vcf_path)
  vcf_mt      <- readVcf(vcf_path_MT)
  vcf_info    <- info(vcf)
  vcf_mt_info <- info(vcf_mt)


  print("We change the names of the mutations, to make it more accessible.")
  if(all(c("GENE", "AA", "CDS") %in% colnames(vcf_info))){
    new_names <- paste0(vcf_info$GENE, "_", vcf_info$AA, "_", vcf_info$CDS)
  } else{
    new_names <- rownames(vcf_info)
    new_names <- gsub(":|\\/|\\?", "_", new_names)
  }
  if(all(c("GENE", "AA", "CDS") %in% colnames(vcf_mt_info))){
    new_names_mt <- paste0(vcf_mt_info$GENE, "_", vcf_mt_info$AA, "_", vcf_mt_info$CDS)
  } else{
    new_names_mt <- rownames(vcf_mt_info)
    new_names_mt <- gsub(":|\\/|\\?", "_", new_names_mt)
  }
  for(i in 1:length(types_unique)){
    type_use <- types_unique[i]
    if(type_use %in% c("scRNAseq_Somatic", "Amplicon_Somatic")){
      print(paste0(type_use, " Variants: ", nrow(consensus_matrix_total[[type_use]])))
      print(paste0(type_use, " Cells: ", ncol(consensus_matrix_total[[type_use]])))
      rownames(coverage_matrix_total[[type_use]])  <- new_names
      rownames(ref_matrix_total[[type_use]])       <- new_names
      rownames(consensus_matrix_total[[type_use]]) <- new_names
      consensus_test <- consensus_matrix_total[[type_use]] > 0
      keep_variants <- rowSums(consensus_test) > 0
      keep_cells <- colSums(consensus_test) > 0
      consensus_matrix_total[[type_use]] <- consensus_matrix_total[[type_use]][keep_variants,keep_cells]
      coverage_matrix_total[[type_use]] <- coverage_matrix_total[[type_use]][keep_variants,keep_cells]
      ref_matrix_total[[type_use]] <- ref_matrix_total[[type_use]][keep_variants,keep_cells]
      print(paste0(type_use, " Variants: ", nrow(consensus_matrix_total[[type_use]])))
      print(paste0(type_use, " Cells: ", ncol(consensus_matrix_total[[type_use]])))
    } else if(type_use %in% c("scRNAseq_MT", "Amplicon_MT")){
      print(paste0(type_use, " Variants: ", nrow(consensus_matrix_total[[type_use]])))
      print(paste0(type_use, " Cells: ", ncol(consensus_matrix_total[[type_use]])))
      rownames(coverage_matrix_total[[type_use]])  <- new_names_mt
      rownames(ref_matrix_total[[type_use]])       <- new_names_mt
      rownames(consensus_matrix_total[[type_use]]) <- new_names_mt
      consensus_test <- consensus_matrix_total[[type_use]] > 0
      keep_variants <- rowSums(consensus_test) > 0
      keep_cells <- colSums(consensus_test) > 0
      consensus_matrix_total[[type_use]] <- consensus_matrix_total[[type_use]][keep_variants,keep_cells]
      coverage_matrix_total[[type_use]] <- coverage_matrix_total[[type_use]][keep_variants,keep_cells]
      ref_matrix_total[[type_use]] <- ref_matrix_total[[type_use]][keep_variants,keep_cells]
      print(paste0(type_use, " Variants: ", nrow(consensus_matrix_total[[type_use]])))
      print(paste0(type_use, " Cells: ", ncol(consensus_matrix_total[[type_use]])))
    }
  }
  rm(consensus_test, keep_variants, keep_cells)
  gc()


  print("We transform the sparse matrices to matrices, so we can calculate the fraction.")
  fraction_total <- list()
  reads_total <- list()
  #coverage_matrix_total_ori <- coverage_matrix_total
  #ref_matrix_total_ori <- ref_matrix_total
  #consensus_matrix_total_ori <- consensus_matrix_total
  #coverage_matrix_total <- coverage_matrix_total_ori
  #ref_matrix_total <- ref_matrix_total_ori
  #consensus_matrix_total <- consensus_matrix_total_ori
  for(i in 1:length(types_unique)){
    type_use <- types_unique[i]
    print(type_use)
    #coverage_matrix_total[[type_use]] <- coverage_matrix_total[[type_use]][1:5000,1:5000]
    #ref_matrix_total[[type_use]] <- ref_matrix_total[[type_use]][1:5000,1:5000]
    #consensus_matrix_total[[type_use]] <- coverage_matrix_total[[type_use]][1:5000,1:5000]
#    colnames(coverage_matrix_total[[type_use]]) <- make.names(colnames(coverage_matrix_total[[type_use]]))
#    rownames(coverage_matrix_total[[type_use]]) <- make.names(rownames(coverage_matrix_total[[type_use]]))
#    colnames(ref_matrix_total[[type_use]]) <- make.names(colnames(ref_matrix_total[[type_use]]))
#    rownames(ref_matrix_total[[type_use]]) <- make.names(rownames(ref_matrix_total[[type_use]]))
#    colnames(consensus_matrix_total[[type_use]]) <- make.names(colnames(consensus_matrix_total[[type_use]]))
#    rownames(consensus_matrix_total[[type_use]]) <- make.names(rownames(consensus_matrix_total[[type_use]]))
    coverage_matrix_total[[type_use]]                             <- as.matrix(coverage_matrix_total[[type_use]])
    ref_matrix_total[[type_use]]                                  <- as.matrix(ref_matrix_total[[type_use]])
    consensus_matrix_total[[type_use]]                            <- as.matrix(consensus_matrix_total[[type_use]])
    reads_total[[type_use]]                                       <- coverage_matrix_total[[type_use]] + ref_matrix_total[[type_use]]
    fraction_total[[type_use]]                                    <- coverage_matrix_total[[type_use]] / reads_total[[type_use]]
    fraction_total[[type_use]][is.na(fraction_total[[type_use]])] <- 0
#    fraction_total[[type_use]] <- sdiv(X = coverage_matrix_total[[type_use]], Y = reads_total[[type_use]],
#                                       names = dimnames(coverage_matrix_total[[type_use]]))
  }
  rm(coverage_matrix_total, ref_matrix_total)
  gc()


  print("We generate 2 maximal matrices with all the somatic/MT variants and all the cells.")
  all_cells <- lapply(consensus_matrix_total, colnames)
  all_cells <- unique(unlist(all_cells))
  if(any(names(consensus_matrix_total) %in% c("scRNAseq_Somatic", "Amplicon_Somatic"))){
    somatic_variants <- consensus_matrix_total[grepl("scRNAseq_Somatic|Amplicon_Somatic", names(consensus_matrix_total))]
    somatic_variants <- lapply(somatic_variants, rownames)
    somatic_variants <- unique(unlist(somatic_variants))
  }
  if(any(names(consensus_matrix_total) %in% c("scRNAseq_MT", "Amplicon_MT"))){
    MT_variants <- consensus_matrix_total[grepl("scRNAseq_MT|Amplicon_MT", names(consensus_matrix_total))]
    MT_variants <- lapply(MT_variants, rownames)
    MT_variants <- unique(unlist(MT_variants))
  }
  consensus_matrix_total_merged_somatic <- matrix(0, ncol = length(all_cells), nrow = length(somatic_variants))
  rownames(consensus_matrix_total_merged_somatic) <- somatic_variants
  colnames(consensus_matrix_total_merged_somatic) <- all_cells
  fraction_total_merged_somatic <- consensus_matrix_total_merged_somatic
  reads_total_merged_somatic <- consensus_matrix_total_merged_somatic
  consensus_matrix_total_merged_MT <- matrix(0, ncol = length(all_cells), nrow = length(MT_variants))
  rownames(consensus_matrix_total_merged_MT) <- MT_variants
  colnames(consensus_matrix_total_merged_MT) <- all_cells
  fraction_total_merged_MT <- consensus_matrix_total_merged_MT
  reads_total_merged_MT <- consensus_matrix_total_merged_MT


  print("We add the information to the merged matrices.")
  meta_data_somatic <- data.frame(Cell = colnames(consensus_matrix_total_merged_somatic))
  meta_data_MT      <- data.frame(Cell = colnames(consensus_matrix_total_merged_MT))
  if(any(names(consensus_matrix_total) == "scRNAseq_Somatic")){
    cells_somatic_scRNAseq <- colnames(consensus_matrix_total[["scRNAseq_Somatic"]])
    consensus_matrix_total_merged_somatic[rownames(consensus_matrix_total[["scRNAseq_Somatic"]]), colnames(consensus_matrix_total[["scRNAseq_Somatic"]])] <- consensus_matrix_total[["scRNAseq_Somatic"]]
    fraction_total_merged_somatic[rownames(fraction_total[["scRNAseq_Somatic"]]), colnames(fraction_total[["scRNAseq_Somatic"]])] <- fraction_total[["scRNAseq_Somatic"]]
    reads_total_merged_somatic[rownames(reads_total[["scRNAseq_Somatic"]]), colnames(reads_total[["scRNAseq_Somatic"]])] <- reads_total[["scRNAseq_Somatic"]]
    meta_data_somatic <- data.frame(meta_data_somatic, scRNAseq = colnames(consensus_matrix_total_merged_somatic) %in% cells_somatic_scRNAseq)
  }
  if(any(names(consensus_matrix_total) == "Amplicon_Somatic")){
    cells_somatic_amplicon <- colnames(consensus_matrix_total[["Amplicon_Somatic"]])
    consensus_matrix_total_merged_somatic[rownames(consensus_matrix_total[["Amplicon_Somatic"]]), colnames(consensus_matrix_total[["Amplicon_Somatic"]])] <- consensus_matrix_total[["Amplicon_Somatic"]]
    fraction_total_merged_somatic[rownames(fraction_total[["Amplicon_Somatic"]]), colnames(fraction_total[["Amplicon_Somatic"]])] <- fraction_total[["Amplicon_Somatic"]]
    reads_total_merged_somatic[rownames(reads_total[["Amplicon_Somatic"]]), colnames(reads_total[["Amplicon_Somatic"]])] <- reads_total[["Amplicon_Somatic"]]
    meta_data_somatic <- data.frame(meta_data_somatic, Amplicon = colnames(consensus_matrix_total_merged_somatic) %in% cells_somatic_amplicon)
  }
  if(any(names(consensus_matrix_total) == "scRNAseq_MT")){
    cells_MT_scRNAseq <- colnames(consensus_matrix_total[["scRNAseq_MT"]])
    consensus_matrix_total_merged_MT[rownames(consensus_matrix_total[["scRNAseq_MT"]]), colnames(consensus_matrix_total[["scRNAseq_MT"]])] <- consensus_matrix_total[["scRNAseq_MT"]]
    fraction_total_merged_MT[rownames(fraction_total[["scRNAseq_MT"]]), colnames(fraction_total[["scRNAseq_MT"]])] <- fraction_total[["scRNAseq_MT"]]
    reads_total_merged_MT[rownames(reads_total[["scRNAseq_MT"]]), colnames(reads_total[["scRNAseq_MT"]])] <- reads_total[["scRNAseq_MT"]]
    meta_data_MT <- data.frame(meta_data_MT, scRNAseq = colnames(consensus_matrix_total_merged_MT) %in% cells_MT_scRNAseq)
  }
  if(any(names(consensus_matrix_total) == "Amplicon_MT")){
    cells_MT_amplicon <- colnames(consensus_matrix_total[["Amplicon_MT"]])
    consensus_matrix_total_merged_MT[rownames(consensus_matrix_total[["Amplicon_MT"]]), colnames(consensus_matrix_total[["Amplicon_MT"]])] <- consensus_matrix_total[["Amplicon_MT"]]
    fraction_total_merged_MT[rownames(fraction_total[["Amplicon_MT"]]), colnames(fraction_total[["Amplicon_MT"]])] <- fraction_total[["Amplicon_MT"]]
    reads_total_merged_MT[rownames(reads_total[["Amplicon_MT"]]), colnames(reads_total[["Amplicon_MT"]])] <- reads_total[["Amplicon_MT"]]
    meta_data_MT <- data.frame(meta_data_MT, Amplicon = colnames(consensus_matrix_total_merged_MT) %in% cells_MT_amplicon)
  }


  print("We generate a SummarizedExperiment object containing the fraction and the consensus matrices.")
  # We want an assay for the Consensus information and for the fraction.
  # As meta data we add a data frame showing the cell id, the associated patient and the sample.
  consensus <- rbind(consensus_matrix_total_merged_somatic, consensus_matrix_total_merged_MT)
  fraction <- rbind(fraction_total_merged_somatic, fraction_total_merged_MT)
  reads <- rbind(reads_total_merged_somatic, reads_total_merged_MT)
  coverage_depth_per_cell <- colMeans(reads)
  meta_data <- cbind(merge(x = meta_data_somatic, y = meta_data_MT, by = "Cell"), coverage_depth_per_cell)
  colnames(meta_data) <-  c("Cell", paste0("Somatic_", colnames(meta_data_somatic)[-1]), paste0("MT_", colnames(meta_data_MT)[-1]), "depth")
  se_merged <- SummarizedExperiment(assays = list(consensus = consensus, fraction = fraction, coverage = reads), colData = meta_data)
  return(se_merged)
}
