#'VariantSelection_Quantile
#'@description
#'We get variants of interest using the quantile thresholding.
#'This function is adapted from the Peter van Galen.
#'Source: https://github.com/petervangalen/MAESTER-2021
#'@importFrom SummarizedExperiment assays colData rowData
#'@importFrom stats quantile
#'@importFrom Matrix rowMeans rowSums
#'@param SE SummarizedExperiment object.
#'@param min_coverage Minimum coverage needed.
#'@param quantiles The lower and upper quantile you want to use. 
#'@param thresholds The VAF thresholds you want to use for the quantiles.
#'@param min_quality The minimum quality you want for the Variants of Interest. Can be ignored by setting it to NULL.
#'@param mean_allele_frequency The minimum mean allele frequency. Default = 0
#'@param remove_nocall Should NoCall cells (consensus = 0) be disregarded during the analysis?
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
VariantSelection_Quantile <- function(SE, min_coverage = 2, quantiles = c(0.1, 0.9), thresholds = c(0.1, 0.9), min_quality = NULL, mean_allele_frequency = 0, remove_nocall = FALSE, verbose = TRUE){
  if(remove_nocall){
    if(verbose) print("We set NoCall cells as NA.")
    nocall_check <- SummarizedExperiment::assays(SE)[["consensus"]] == 0
    SummarizedExperiment::assays(SE)[["fraction"]][nocall_check] <- NA
    SummarizedExperiment::assays(SE)[["coverage"]][nocall_check] <- NA
  }

  if(verbose) print("Get the mean allele frequency and coverage.")
  mean_af  <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["fraction"]], na.rm = TRUE)
  mean_cov <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["coverage"]], na.rm = TRUE)

  if(verbose) print("Get the quantiles of the VAFs of each variant.")
  quantiles <- lapply(quantiles, function(x) apply(SummarizedExperiment::assays(SE)[["fraction"]], 1, quantile, x, na.rm = TRUE))

  # Apply min_quality filtering
  if(!is.null(min_quality)){
    vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, VariantQuality = SummarizedExperiment::rowData(SE)$VariantQuality, Quantile1 = quantiles[[1]], Quantile2 = quantiles[[2]])
    vars <- subset(vars, !is.na(vars$VariantQuality) & vars$VariantQuality > min_quality)
  } else{
    vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, Quantile1 = quantiles[[1]], Quantile2 = quantiles[[2]])
  }

  if(verbose) print("Thresholding using the quantile approach.")
  vois <- subset(vars, vars$Mean_AF > mean_allele_frequency & vars$Mean_Cov > min_coverage & vars$Quantile1 <= thresholds[1] & vars$Quantile2 >= thresholds[2])

  return(rownames(vois))
}
