#'GetVariantInfo
#'@description
#'We add the genotyping information for a set of variants to a Seurat object.
#'The function returns a matrix with the values from the specified assay.
#'@import SummarizedExperiment
#'@param SE SummarizedExperiment object.
#'@param Seurat The Seurat object.
#'@param information The assay with the desired information. Default: consensus
#'@param variants A vector of variants.
#'@export
SetVariantInfo <- function(SE, Seurat, information = "consensus", variants = NULL){
  # We check if a vector of variants has been supplied.
  if(is.null(variants)){
    stop("Error: You forgot to supply a vector of variants.")
  }
  # We check if the variants are actually in a vector.
  if(!is.vector(variants)){
    stop("Error: Your variants are not in a vector.")
  }
  # We check if all variants are in the SE object.
  variants_check <- variants %in% rownames(SE)
  if(!all(variants_check)){
    variants_check <- sum(variants_check)
    stop(paste0("Only ", variants_check, " of ", length(variants), " are in the SE object."))
  }
  # We check if the requested assay is actually present.
  assay_check <- information %in% names(assays(SE))
  if(!assay_check){
    stop("The assay you wants is not present in the object.")
  }
  res <- t(assays(SE)[[information]][variants, , drop = FALSE])
  # We check if all the cells are actually in the Seurat object.
  # If not, we only add the information for the ones present. 
  # We execute an error function if there are zero cells present.
  cell_check_SE <- rownames(res) %in% colnames(Seurat)
  if(sum(cell_check_SE) == 0){
    stop("Error: 0 cells from the SE object are in the Seurat object.")
  }
  res <- res[rownames(res) %in% colnames(Seurat), , drop = FALSE]
  Seurat <- AddMetaData(Seurat, res, col.name = colnames(res))
  return(Seurat)
}
