#'We load the MAEGATK output and transform it to be compatible with the VarTrix output.
#'@import Matrix SummarizedExperiment
#'@param samples_path Path to the input folder.
#'@param samples_file Path to the csv file with the samples to be loaded.
#'@param type_use The type of input. Has to be one of: scRNAseq_MT, Amplicon_MT. Only used if samples_path is not NULL.
#'@param patient The patient you want to load.
#'@param chromosome_prefix The prefix you want use. Default: "chrM"
#'@export
LoadingMAEGATK_typewise <- function(samples_file, samples_path = NULL, patient, type_use = "scRNAseq_MT", chromosome_prefix = "chrM",
                                    min_cells = 2){
  if(!is.null(samples_path)){
    samples <- list.files(samples_path)
    samples <- grep(patient, samples, value = TRUE)
    samples_file <- data.frame(patient = patient, sample = samples, input_folder = samples_path)
  } else{
    print("We read in the samples file.")
    samples_file <- read.csv(samples_file)

    print("We subset to the patient of interest.")
    samples_file <- samples_file[grep("maegatk|mgatk", samples_file$source, ignore.case = TRUE),]
    samples_file <- samples_file[grep(patient, samples_file$patient),]
    samples_file <- samples_file[grep(type_use, samples_file$type),]

    print("We get the different samples.")
    samples <- samples_file$sample
  }


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
  }


  print("We merge the samples.")
  se_merged <- do.call("cbind", se_ls)
  rm(se_ls)
  gc()


  print("We get the allele frequency.")
  fraction <- computeAFMutMatrix(se_merged, chromosome_prefix = chromosome_prefix)
  #fraction <- data.matrix(fraction)
  #fraction <- fraction[!rownames(fraction) %in% paste0(chromosome_prefix, "_", c("3107_N_A", "3107_N_C", "3107_N_G", "3107_N_T")),]


  print("We get the coverage information.")
  coverage <- CalculateCoverage(SE = se_merged, chromosome_prefix = chromosome_prefix)
  coverage <- coverage[match(rownames(fraction), rownames(coverage)),]
  

  print("We calculate the consensus information.")
  consensus <- CalculateConsensus(SE = se_merged, chromosome_prefix = chromosome_prefix)
  # We order the consensus matrix like the coverage matrix.
  consensus <- consensus[match(rownames(fraction), rownames(consensus)),]


  print("We perform some filtering to reduce the memory needed.")
  print(paste0("We remove variants, which are not detected in at least ", min_cells, "."))
  keep_variants <- rowSums(consensus >= 2)
  keep_variants <- keep_variants >= min_cells
  consensus <- consensus[keep_variants,]
  coverage <- coverage[keep_variants,]
  fraction <- fraction[keep_variants,]

  print("We remove cells that are always NoCall.")
  consensus_test <- consensus > 0
  keep_cells <- colSums(consensus_test) > 0
  consensus <- consensus[,keep_cells]
  coverage <- coverage[,keep_cells]
  fraction <- fraction[,keep_cells]


  print("We add the information to the merged matrices.")
  coverage_depth_per_cell <- rownames(coverage)
  coverage_depth_per_cell <- gsub("_._.$", "", coverage_depth_per_cell)
  coverage_depth_per_cell <- !duplicated(coverage_depth_per_cell)
  coverage_depth_per_cell <- coverage[coverage_depth_per_cell,]
  coverage_depth_per_cell <- colMeans(coverage_depth_per_cell)
  meta_data <- data.frame(Cell = colnames(consensus), AverageCoverage = coverage_depth_per_cell)
  se_output <- SummarizedExperiment(assays = list(consensus = consensus, fraction = fraction, coverage = coverage),
                                    colData = meta_data)
  return(se_output)
}
