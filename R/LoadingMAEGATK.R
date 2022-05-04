samples_path <- "~/labcluster_data/scRNA/SingleCell_Variant_Correlation_Test/JAK2Amplicon/MPN_Test_Data_local.csv"
patient <- "ET_83497"
vcf_path <- "~/labcluster_data/scRNA/SingleCell_Variant_Correlation_Test/VariantsOfInterest/ALFA_subset_MAF2_prefix.vcf"
vcf_path_MT <- "~/labcluster_data/mg000001/MPN/Mutations/Mitochondrial/results/chrM_Input_VCF_NoMAF_Filtering.vcf"

source("~/labcluster_data/mg000001/sigurd/R/LoadingVarTrix.R")
source("~/labcluster_data/mg000001/sigurd/R/save_loader_helper.R")
source("~/labcluster_data/mg000001/sigurd/R/%fin%.R")
source("~/labcluster_data/mg000001/sigurd/R/getAltMatrix.R")
source("~/labcluster_data/mg000001/sigurd/R/getRefMatrix.R")
source("~/labcluster_data/mg000001/sigurd/R/CalculateConsensus.R")
source("~/labcluster_data/mg000001/sigurd/R/Merging_SE_list.R")
source("~/labcluster_data/mg000001/sigurd/R/computeAFMutMatrix.R")

LoadingMAEGATK <- function(samples_path, vcf_path, vcf_path_MT, patient){
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(VariantAnnotation))
  suppressPackageStartupMessages(library(archive))
  suppressPackageStartupMessages(library(SummarizedExperiment))

  print("We read in the samples file.")
  samples_file <- read.csv(samples_path)


  print("We subset to the patient of interest.")
  samples_file <- samples_file[grep("maegatk|mgatk", samples_file$source, ignore.case = TRUE),]
  samples_file <- samples_file[grep(patient, samples_file$patient),]


  print("We get the different samples.")
  samples <- samples_file$sample
  samples_unique <- unique(samples)
  types <- samples_file$type
  types_unique <- unique(types)


  print("We read in the cell barcodes output by CellRanger as a list.")
  barcodes <- lapply(samples_file$cells, read.table)
  names(barcodes) <- paste0(types, "_", samples)


  print("We load the MAEGATK output files.")
  se_ls <- list()
  for(i in 1:nrow(samples_file)){
    print(paste0("Loading sample ", i, " of ", nrow(samples_file)))
    type_use <- samples_file$type[i]
    input_folder_use <- samples_file$input_folder[i]
    sample_use <- samples_file$sample[i]

    # We get the final output file for either mgatk or maegatk.
    final_output_file <- list.files(paste0(input_folder_use, sample_use, "/final/"), full.names = TRUE)
    final_output_file <- grep(paste0("maegtk.rds|maegatk.rds|mgatk.rds|", sample_use, ".rds"), final_output_file, value = TRUE)
    se_ls[[type_use]][[sample_use]] <- load_object(final_output_file)
    colnames(se_ls[[type_use]][[sample_use]]) <- paste0(sample_use, "_", colnames(se_ls[[type_use]][[sample_use]]))
  }


  print("We merge the samples.")
  se_merged_ls <- lapply(se_ls, Merging_SE_list)
  rm(se_ls)
  gc()


  input_base <- c("A", "C", "G", "T")
  print(paste0("Consensus: ", patient))
  consensus <- lapply(se_merged_ls, function(x) lapply(input_base, CalculateConsensus, input_se = x))
  consensus <- lapply(consensus, do.call, what = "rbind")


  print("We get the allele frequency.")
  fraction <- lapply(se_merged_ls, computeAFMutMatrix)
  fraction <- lapply(fraction, data.matrix)
  fraction <- lapply(fraction, function(x){
    rownames(x) <- gsub(">", ".", rownames(x))
    return(x)
    })
  # The fractions include the position 3107_N. Since this is always "mutated", we remove it.
  for(i in 1:length(fraction)) fraction[[i]] <- fraction[[i]][!rownames(fraction[[i]]) %fin% c("3107_N.A", "3107_N.C", "3107_N.G", "3107_N.T"),]

  print("We generate 2 maximal matrices with all the somatic/MT variants and all the cells.")
  all_cells <- lapply(consensus, colnames)
  all_cells <- unique(unlist(all_cells))
  MT_variants <- rownames(consensus[[1]])

  consensus_merged <- matrix(0, ncol = length(all_cells), nrow = length(MT_variants))
  rownames(consensus_merged) <- MT_variants
  colnames(consensus_merged) <- all_cells
  fraction_merged <- consensus_merged


  print("We add the information to the merged matrices.")
  meta_data_MT <- data.frame(Cell = colnames(consensus_merged))
  if(any(names(consensus) == "scRNAseq_MT")){
    cells_MT_scRNAseq <- colnames(consensus[["scRNAseq_MT"]])
    consensus_merged[rownames(consensus[["scRNAseq_MT"]]), colnames(consensus[["scRNAseq_MT"]])] <- consensus[["scRNAseq_MT"]]
    fraction_merged[rownames(fraction[["scRNAseq_MT"]]), colnames(fraction[["scRNAseq_MT"]])] <- fraction[["scRNAseq_MT"]]
    meta_data_MT <- data.frame(meta_data_MT, scRNAseq = colnames(consensus_merged) %in% cells_MT_scRNAseq)
  }
  if(any(names(consensus) == "Amplicon_MT")){
    cells_MT_amplicon <- colnames(consensus[["Amplicon_MT"]])
    consensus_merged[rownames(consensus[["Amplicon_MT"]]), colnames(consensus[["Amplicon_MT"]])] <- consensus[["Amplicon_MT"]]
    fraction_merged[rownames(fraction[["Amplicon_MT"]]), colnames(fraction[["Amplicon_MT"]])] <- fraction[["Amplicon_MT"]]
    meta_data_MT <- data.frame(meta_data_MT, Amplicon = colnames(consensus_merged) %in% cells_MT_amplicon)
  }


  print("We generate a SummarizedExperiment object containing the fraction and the consensus matrices.")
  # We want an assay for the Consensus information and for the fraction.
  # As meta data we add a data frame showing the cell id, the associated patient and the sample.
  se_MT <- SummarizedExperiment(assays = list(consensus_MT = consensus_merged, fraction_MT = fraction_merged),
                                colData = meta_data_MT)
  return(se_MT)
  #print("We save the unfiltered consensus matrix.")
  #save_object(consensus_matrix_total_merged_somatic, paste0(output, "/", patient, "/", patient, "_unfiltered_consensus_somatic_lz4.Rds"), file_format = "lz4")
  #save_object(fraction_total_merged_somatic,         paste0(output, "/", patient, "/", patient, "_unfiltered_fraction_somatic_lz4.Rds"),  file_format = "lz4")
  #save_object(consensus_matrix_total_merged_MT,      paste0(output, "/", patient, "/", patient, "_unfiltered_consensus_MT_lz4.Rds"),      file_format = "lz4")
  #save_object(fraction_total_merged_MT,              paste0(output, "/", patient, "/", patient, "_unfiltered_fraction_MT_lz4.Rds"),       file_format = "lz4")
  #save_object(se_somatic,                            paste0(output, "/", patient, "/", patient, "_unfiltered_SE_somatic_lz4.Rds"),        file_format = "lz4")
  #save_object(se_MT,                                 paste0(output, "/", patient, "/", patient, "_unfiltered_SE_MT_lz4.Rds"),             file_format = "lz4")
}
