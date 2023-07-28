#'VariantQuantileThresholding
#'@description
#'We get variants of interest using the quantile thresholding.
#'This function is adapted from the Peter van Galen.
#'Source: https://github.com/petervangalen/MAESTER-2021
#'@import dplyr SummarizedExperiment
#'@param SE SummarizedExperiment object.
#'@param min_coverage Minimum coverage needed.
#'@param quantiles The lower and upper quantile you want to use. 
#'@param thresholds The VAF thresholds you want to use for the quantiles.
#'@param min_quality The minimum quality you want for the Variants of Interest. Can be ignored by setting it to NULL.
#'@export
VariantQuantileThresholding <- function(SE, min_coverage = 2, quantiles = c(0.1, 0.9), thresholds = c(0.1, 0.9), min_quality = 30){
  print("Get the mean allele frequency and coverage.")
  mean_af <- rowMeans(assays(SE)[["fraction"]], na.rm = TRUE)
  mean_cov <- rowMeans(assays(SE)[["coverage"]], na.rm = TRUE)


  print("Get the quantiles of the VAFs of each variant.")
  quantiles <- lapply(quantiles, function(x) apply(assays(SE)[["fraction"]], 1, quantile, x, na.rm = TRUE))


  print("Collect all information in a tibble")
  vars <- do.call(cbind, c(list(mean_af), list(mean_cov), list(rowData(SE)$VariantQuality), quantiles))


  print("Thresholding using the quantile approach.")
  voi_ch <- subset(vars, vars[,2] > min_coverage & vars[,4] < thresholds[1] & vars[,5] > thresholds[2])
  if(!is.null(min_quality)){
    voi_ch <- voi_ch[!is.na(voi_ch[,3]),]
    voi_ch <- subset(voi_ch, voi_ch[,3] > min_quality)
  }
  voi_ch <- rownames(voi_ch)
  return(voi_ch)
}
