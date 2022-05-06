#'We combine two SummarizedExperiment objects.
#'@import SummarizedExperiment
#'@param se_somatic se_MT suffixes min_intersecting_cells
CombineSEobjects <- function(se_somatic, se_MT, suffixes = c("_somatic", "_MT")){
  features <- combine_NAMES(names(se_somatic), names(se_MT))
  cells <- combine_NAMES(colnames(se_somatic), colnames(se_MT))

  meta_data_somatic <- colData(se_somatic)
  meta_data_MT      <- colData(se_MT)
  meta_data <- merge(meta_data_somatic, meta_data_MT, by = "Cell", all = TRUE, suffixes = suffixes)
  meta_data <- meta_data[match(cells, meta_data$Cell),]
  
  assays_combined <- mendoapply(combine, assays(se_somatic), assays(se_MT))
  se_combined <- SummarizedExperiment(assays = assays_combined,
                                      colData = meta_data)
  return(se_combined)
}

