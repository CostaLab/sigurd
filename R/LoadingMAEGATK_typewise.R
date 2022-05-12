#'We load the MAEGATK output and transform it to be compatible with the VarTrix output.
#'@import Matrix SummarizedExperiment VariantAnnotation
#'@param samples_path Path to the csv file with the samples to be loaded.
#'@param vcf_path Path to the VCF file with the somatic variants.
#'@param vcf_path_MT Path to the VCF file with the MT variants.
#'@param type_use The type of input. Has to be one of: scRNAseq_MT, Amplicon_MT.
#'@param patient The patient you want to load.
#'@export
LoadingMAEGATK_typewise <- function(samples_path, vcf_path, patient, type_use = "scRNAseq_MT", chromosome_prefix = "chrM"){
  print("We read in the samples file.")
  samples_file <- read.csv(samples_path)


  print("We subset to the patient of interest.")
  samples_file <- samples_file[grep("maegatk|mgatk", samples_file$source, ignore.case = TRUE),]
  samples_file <- samples_file[grep(patient, samples_file$patient),]
  samples_file <- samples_file[grep(type_use, samples_file$type),]


  print("We get the different samples.")
  samples <- samples_file$sample


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


  print("We calculate the coverage information.")
  # We first get the ref reads per sample.
  ref_reads <- rbind(getRefMatrix(SE_object = se_merged, letter = "A", chromosome_prefix = chromosome_prefix), 
                     getRefMatrix(SE_object = se_merged, letter = "C", chromosome_prefix = chromosome_prefix),
                     getRefMatrix(SE_object = se_merged, letter = "G", chromosome_prefix = chromosome_prefix), 
                     getRefMatrix(SE_object = se_merged, letter = "T", chromosome_prefix = chromosome_prefix))
  # We then get the alt reads per sample and a specific base.
  alt_reads <- list(getAltMatrix(SE_object = se_merged, letter = "A", chromosome_prefix = chromosome_prefix),
                    getAltMatrix(SE_object = se_merged, letter = "C", chromosome_prefix = chromosome_prefix),
                    getAltMatrix(SE_object = se_merged, letter = "G", chromosome_prefix = chromosome_prefix),
                    getAltMatrix(SE_object = se_merged, letter = "T", chromosome_prefix = chromosome_prefix))
  coverage <- lapply(alt_reads, CalculateCoverage, ref_reads = ref_reads)
  coverage <- do.call("rbind", coverage)
  rm(ref_reads, alt_reads)


  print("We calculate the consensus information.")
  consensus <- CalculateConsensus(SE = se_merged, chromosome_prefix = chromosome_prefix)
  # We order the consensus matrix like the coverage matrix.
  consensus <- consensus[match(rownames(coverage), rownames(consensus)),]


  print("We get the allele frequency.")
  fraction <- computeAFMutMatrix(se_merged, chromosome_prefix = "chrM")
  fraction <- data.matrix(fraction)
  fraction <- fraction[!rownames(fraction) %in% paste0(chromosome_prefix, "_", c("3107_N_A", "3107_N_C", "3107_N_G", "3107_N_T")),]
  
  
  print("We perform some filtering to reduce the memory needed.")
  keep_variants <- rowSums(consensus >= 2)
  keep_variants <- keep_variants >= 2
  consensus <- consensus[keep_variants,]
  coverage <- coverage[keep_variants,]
  fraction <- fraction[keep_variants,]

  # We remove cells that are always NoCall.
  consensus_test <- consensus > 0
  keep_cells <- colSums(consensus_test) > 0
  consensus <- consensus[,keep_cells]
  coverage <- coverage[,keep_cells]
  fraction <- fraction[,keep_cells]


  print("We add the information to the merged matrices.")
  coverage_depth_per_cell <- colMeans(coverage)
  meta_data <- data.frame(Cell = colnames(consensus), AverageCoverage = coverage_depth_per_cell)
  se_output <- SummarizedExperiment(assays = list(consensus = consensus, fraction = fraction, coverage = coverage),
                                    colData = meta_data)
  return(se_output)
}
