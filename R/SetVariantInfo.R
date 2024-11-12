#'SetVariantInfo
#'@description
#'We add the genotyping information for a set of variants to a Seurat object.
#'The function returns a Seurat object with the values from the specified assay added as meta data.
#'@importFrom SummarizedExperiment assays
#'@importFrom Matrix t
#'@importFrom Seurat AddMetaData
#'@param SE SummarizedExperiment object.
#'@param seurat_object The Seurat object.
#'@param information The assay with the desired information. Default: consensus
#'@param variants A vector of variants.
#'@param consensus_character Should the consensus information be save as a character value?
#'@param consensus_group Should the consensus groups Alt and Both be merged? Has no effect if information is not consensus.
#'@export
SetVariantInfo <- function(SE, seurat_object, information = "consensus", variants = NULL, consensus_character = TRUE, consensus_group = TRUE){
  # We check if a vector of variants has been supplied.
  if(is.null(variants)){
    print("You did not supply a vector of variants. All variants will be used.")
    variants <- rownames(SE)
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
  if(information == "consensus"){
    if(consensus_group){
     res@x[res@x == 3] <- 2
    }
    if(consensus_character){
      res <- as.matrix(res)
      res[res == 0] <- "NoCall"
      res[res == 1] <- "Ref"
      res[res == 2] <- "Alt"
      res[res == 3] <- "Both"
    }
  }
  seurat_object <- Seurat::AddMetaData(seurat_object, res, col.name = paste0(make.names(colnames(res)), "_", information))
  return(seurat_object)
}
