#'We combine two SummarizedExperiment objects with big matrices.
#'@import SummarizedExperiment BiocGenerics bigmemory
#'@param se_somatic SummarizedExperiment object for the somatic variants.
#'@param se_MT SummarizedExperiment object for the MT variants.
#'@param suffixes The suffixes you want to add to the meta data.frame.
#'@export
CombineSEobjects_big <- function(se_somatic, se_MT, suffixes = c("_somatic", "_MT")){
  features <- combine_NAMES(names(se_somatic), names(se_MT))
  cells <- combine_NAMES(colnames(se_somatic), colnames(se_MT))

  meta_data_somatic <- colData(se_somatic)
  meta_data_MT      <- colData(se_MT)
  meta_data <- merge(meta_data_somatic, meta_data_MT, by = "Cell", all = TRUE, suffixes = suffixes)
  meta_data <- meta_data[match(cells, meta_data$Cell),]

  assays_somatic <- assays(se_somatic)
  assays_somatic <- lapply(assays_somatic, as.matrix)
  assays_MT <- assays(se_MT)
  assays_MT <- lapply(assays_MT, as.matrix)

  assays_combined <- S4Vectors::mendoapply(BiocGenerics::combine, assays_somatic, assays_MT)
  assays_combined[["consensus"]][is.na(assays_combined[["consensus"]])] <- 0
  assays_combined[["fraction"]][is.na(assays_combined[["fraction"]])] <- 0
  assays_combined[["coverage"]][is.na(assays_combined[["coverage"]])] <- 0
  assays_combined <- lapply(assays_combined, as.big.matrix)
  se_combined <- SummarizedExperiment(assays = assays_combined,
                                      colData = meta_data)
  return(se_combined)
}

