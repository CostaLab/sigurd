#'CalculateQuality
#'@description
#'We calculate the quality per variant.
# #'@import MatrixGenerics
#'@importFrom SummarizedExperiment assays
#'@importFrom Matrix rowSums
#'@importFrom utils tail
#'@param SE SummarizedExperiment object.
#'@param variants The variants you want to get the quality for.
#'@param chromosome_prefix List of matrices for the alternative reads.
#'@export
CalculateQuality <- function(SE, variants, chromosome_prefix = "chrM"){
  variants <- gsub(paste0(chromosome_prefix, "_"), "", variants)
  qualities <- lapply(c("A", "T", "C", "G"), function(x){
    fwrev <- cbind(SummarizedExperiment::assays(SE)[[paste0(x, "_counts_fw")]], SummarizedExperiment::assays(SE)[[paste0(x, "_counts_rev")]])
    qualities_fwrev <- cbind(SummarizedExperiment::assays(SE)[[paste0(x, "_qual_fw")]], SummarizedExperiment::assays(SE)[[paste0(x, "_qual_rev")]])
    variants_use <- strsplit(variants, "")
    variants_use <- sapply(variants_use, utils::tail, n = 1)
    variants_use <- variants_use == x
    variants_use_names <- variants[variants_use]
    variants_use <- as.numeric(gsub("_.*", "", variants_use_names))
    variants_use_names <- paste0(chromosome_prefix, "_", variants_use_names)
    fwrev <- fwrev[variants_use,]
    rownames(fwrev) <- variants_use_names
    qualities_fwrev <- qualities_fwrev[variants_use,]
    rownames(qualities_fwrev) <- variants_use_names
    fwrev <- fwrev > 0
    qualities_fwrev <- qualities_fwrev * fwrev
    qualities <- apply(qualities_fwrev, 1, sum) / Matrix::rowSums(fwrev > 0)
    qualities[qualities == 0] <- NA
    return(qualities)
  })
  qualities <- unlist(qualities)
  qualities <- qualities[paste0(chromosome_prefix, "_", variants)]
  return(qualities)  
}
