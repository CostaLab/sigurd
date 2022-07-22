#'VariantQuantileThresholding
#'@description
#'We get variants of interest using the quantile thresholding.
#'Source: https://github.com/petervangalen/MAESTER-2021
#'@import dplyr SummarizedExperiment
#'@param SE SummarizedExperiment object.
#'@param min_coverage Minimum coverage needed.
#'@param quantiles The lower and upper quantile you want to use. 
#'@param thresholds The VAF thresholds you want to use for the quantiles.
#'@export
VariantQuantileThresholding <- function(SE, min_coverage = 2, quantiles = c(0.1, 0.9), thresholds = c(0.1, 0.9)){
  # This function is adapted from the Peter van Galen.
  print("Get the mean allele frequency and coverage.")
  mean_af <- rowMeans(assays(SE)[["fraction"]], na.rm = TRUE)
  mean_cov <- rowMeans(assays(SE)[["coverage"]], na.rm = TRUE)


  print("Get the quantiles of the VAFs of each variant.")
  #quantiles <- c("q01" = 0.01, "q10" = 0.1, "q50" = 0.5, "q90" = 0.9, "q99" = 0.99)
  start_time <- Sys.time()
  quantiles <- lapply(quantiles, function(x) apply(assays(SE)[["fraction"]], 1, quantile, x, na.rm = TRUE))
  Sys.time() - start_time


  print("Collect all information in a tibble")
  vars_tib <- as_tibble(do.call(cbind, c(list(mean_af), list(mean_cov), quantiles)), rownames = "var")
  #colnames(vars_tib)[2] <- "mean_af"
  #colnames(vars_tib)[3] <- "mean_cov"


  print("Thresholding using the quantile approach.")
  #voi_ch <- subset(vars_tib, mean_cov > min_coverage & V3 < thresholds[1] & V4 > thresholds[2])$var
  voi_ch <- subset(vars_tib, V2 > min_coverage & V3 < thresholds[1] & V4 > thresholds[2])$var
  return(voi_ch)
}
