#'CombineSEobjects
#'@description
#'We combine two SummarizedExperiment objects.
#'@import SummarizedExperiment BiocGenerics
#'@param se_somatic SummarizedExperiment object for the somatic variants.
#'@param se_MT SummarizedExperiment object for the MT variants.
#'@param suffixes The suffixes you want to add to the meta data.frame.
#'@export
CombineSEobjects <- function(se_somatic, se_MT, suffixes = c("_somatic", "_MT")){
  # We check if the assays are equally named.
  assay_names_somatic <- names(assays(se_somatic))
  assay_names_MT <- names(assays(se_MT))
  if(!all(assay_names_somatic == assay_names_MT)){
    stop("Your assays are not equally named or ordered.")
  }
  features <- combine_NAMES(names(se_somatic), names(se_MT))
  cells <- combine_NAMES(colnames(se_somatic), colnames(se_MT))

  meta_data_somatic <- colData(se_somatic)
  meta_data_MT      <- colData(se_MT)
  meta_data <- merge(meta_data_somatic, meta_data_MT, by = "Cell", all = TRUE, suffixes = suffixes)
  meta_data <- meta_data[match(cells, meta_data$Cell),]

  #assays_somatic <- assays(se_somatic)
  #assays_somatic <- lapply(assays_somatic, as.matrix)
  #assays_MT <- assays(se_MT)
  #assays_MT <- lapply(assays_MT, as.matrix)
  #assays_combined <- S4Vectors::mendoapply(BiocGenerics::combine, assays_somatic, assays_MT)
  #assays_combined[[1]] <- as(assays_combined[[1]], "dgCMatrix")
  #assays_combined[[2]] <- as(assays_combined[[2]], "dgCMatrix")
  #assays_combined[[3]] <- as(assays_combined[[3]], "dgCMatrix")
  #assays_combined[["consensus"]]@x[is.na(assays_combined[["consensus"]]@x)] <- 0
  #assays_combined[["fraction"]]@x[is.na(assays_combined[["fraction"]]@x)] <- 0
  #assays_combined[["coverage"]]@x[is.na(assays_combined[["coverage"]]@x)] <- 0

  assays_combined <- lapply(assay_names_somatic, function(x){
    result <- combine_SparseMatrix(assays(se_somatic)[[x]], assays(se_MT)[[x]])
  })
  names(assays_combined) <- assay_names_somatic

  se_combined <- SummarizedExperiment(assays = assays_combined,
                                      colData = meta_data)
  return(se_combined)
}

