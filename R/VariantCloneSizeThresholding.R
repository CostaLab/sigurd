#'VariantCloneSizeThresholding
#'@description
#'We get variants of interest using a clone size thresholding.
#'Source: https://github.com/petervangalen/MAESTER-2021
#'@importFrom SummarizedExperiment assays
#'@importFrom Matrix rowMeans
#'@param SE SummarizedExperiment object.
#'@param min_coverage Minimum coverage a variant needs to have.
#'@param fraction_negative_cells The fraction of negative cells needed. 
#'@param min_clone_size minimum number of cells.
#'@param vaf_threshold Variant Allele Threshold. Cells above this threshold are considered mutated.
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
VariantCloneSizeThresholding <- function(SE, min_coverage = 2, fraction_negative_cells = 0.9, min_clone_size = 10, vaf_threshold = 0.5, verbose = TRUE){
  if(verbose) print("Get the mean allele frequency and coverage.")
  mean_fraction <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["fraction"]], na.rm = TRUE)
  mean_coverage <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["coverage"]], na.rm = TRUE)

  if(verbose) print("Collect all information in a tibble")
  vars <- do.call(cbind, c(list(mean_fraction), list(mean_coverage)))
  colnames(vars) <- c("mean_fraction", "mean_coverage")
  
  if(verbose) print("We add the number of cells that exceed the VAF thresholds.")
  n_cells <- apply(SummarizedExperiment::assays(SE)[["fraction"]], 1, function(x) sum(x > vaf_threshold, na.rm = TRUE))
  VAF_threshold <- apply(SummarizedExperiment::assays(SE)[["fraction"]], 1, function(x) sum(x > vaf_threshold, na.rm = TRUE))
  vars <- cbind(vars, n_cells, VAF_threshold)
  
  if(verbose) print("Thresholding using the clone size approach.")
  voi_ch <- subset(vars, vars$mean_coverage > min_coverage & 
                   vars$n_cells > ceiling(fraction_negative_cells * ncol(SE)) &
                   vars$VAF_threshold > min_clone_size)
  voi_ch <- rownames(voi_ch)
  return(voi_ch)
}
