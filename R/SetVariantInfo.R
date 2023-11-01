#'GetVariantInfo
#'@description
#'We add the genotyping information for a set of variants to a Seurat object.
#'The function returns a matrix with the values from the specified assay.
#'@import SummarizedExperiment Seurat
#'@param SE SummarizedExperiment object.
#'@param seurat_object The Seurat object.
#'@param information The assay with the desired information. Default: consensus
#'@param variants A vector of variants.
#'@export
SetVariantInfo <- function(SE, seurat_object, information = "consensus", variants = NULL){
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
  assay_check <- information %in% names(SummarizedExperiment::assays(SE))
  if(!assay_check){
    stop("The assay you wants is not present in the object.")
  }
  res <- t(SummarizedExperiment::assays(SE)[[information]][variants, , drop = FALSE])
  # We check if all the cells are actually in the Seurat object.
  # If not, we only add the information for the ones present. 
  # We execute an error function if there are zero cells present.
  cell_check_SE <- rownames(res) %in% colnames(seurat_object)
  if(sum(cell_check_SE) == 0){
    stop("Error: 0 cells from the SE object are in the Seurat object.")
  }
  res <- res[rownames(res) %in% colnames(seurat_object), , drop = FALSE]
  seurat_object <- Seurat::AddMetaData(seurat_object, res, col.name = paste0(colnames(res), "_", information))
  return(seurat_object)
}
