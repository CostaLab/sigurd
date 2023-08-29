#'CalculateQuality
#'@description
#'We calculate the quality per variant.
#'@import MatrixGenerics SummarizedExperiment
#'@param SE SummarizedExperiment object.
#'@param chromosome_prefix List of matrices for the alternative reads.
#'@export
CalculateQuality <- function(SE, variants = rownames(reads_alt), chromosome_prefix = "chrM"){
  variants <- gsub(paste0(chromosome_prefix, "_"), "", variants)
  qualities <- lapply(c("A", "T", "C", "G"), function(x){
    fwrev <- cbind(assays(SE)[[paste0(x, "_counts_fw")]], assays(SE)[[paste0(x, "_counts_rev")]])
    qualities_fwrev <- cbind(assays(SE)[[paste0(x, "_qual_fw")]], assays(SE)[[paste0(x, "_qual_rev")]])
    variants_use <- strsplit(variants, "")
    variants_use <- sapply(variants_use, tail, n = 1)
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
    qualities <- apply(qualities_fwrev, 1, sum) / rowSums(fwrev > 0)
    qualities[qualities == 0] <- NA
    return(qualities)
  })
  qualities <- unlist(qualities)
  qualities <- qualities[paste0(chromosome_prefix, "_", variants)]
  return(qualities)  
}
