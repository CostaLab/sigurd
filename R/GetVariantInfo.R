#'GetVariantInfo
#'@description
#'We get the genotyping information for a set of variants.
#'The function returns a matrix with the values from the specified assay.
#'@importFrom SummarizedExperiment assays
#'@param SE SummarizedExperiment object.
#'@param information The assay with the desired information. Default: consensus
#'@param variants A vector of variants.
#'@param cells A vector of cell IDs. On default all cells are returned. Default: NULL
#'@export
GetVariantInfo <- function(SE, information = "consensus", variants = NULL, cells = NULL){
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
  res <- SummarizedExperiment::assays(SE)[[information]][variants, , drop = FALSE]
  # We subset the result to only include the cells of interest.
  # We check if the cells vector is not NULL.
  if(!is.null(cells)){
    # We check if all cell IDs are in the SE object.
    cell_check <- cells %in% colnames(SE)
    if(!all(cell_check)){
      stop(paste0("Error: Only ", sum(cell_check), " cells are present in the object."))
    } else{
      res <- res[, cells, drop = FALSE]
    }
  }
  return(res)
}
