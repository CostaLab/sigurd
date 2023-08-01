#'VariantQuantileThresholding
#'@description
#'We get variants of interest using the quantile thresholding.
#'If you use top_cells and top_VAF, you have to only supply one quantil value (quantiles = 0.9, thresholds = 0).
#'This function is adapted from the Peter van Galen.
#'Source: https://github.com/petervangalen/MAESTER-2021
#'@import dplyr SummarizedExperiment
#'@param SE SummarizedExperiment object.
#'@param min_coverage Minimum coverage needed.
#'@param quantiles The lower and upper quantile you want to use. 
#'@param thresholds The VAF thresholds you want to use for the quantiles.
#'@param min_quality The minimum quality you want for the Variants of Interest. Can be ignored by setting it to NULL.
#'@param top_cells The number of cells with at least top_VAF percent for a variant.
#'@param top_VAF The VAF for the top cells.
#'@export
VariantQuantileThresholding <- function(SE, min_coverage = 2, quantiles = c(0.1, 0.9), thresholds = c(0.1, 0.9), top_cells = NULL, top_VAF = NULL, min_quality = 30){
  print("Get the mean allele frequency and coverage.")
  mean_af <- rowMeans(assays(SE)[["fraction"]], na.rm = TRUE)
  mean_cov <- rowMeans(assays(SE)[["coverage"]], na.rm = TRUE)

  if(any(is.null(top_cells), is.null(top_VAF))){
    print("Get the quantiles of the VAFs of each variant.")
    quantiles <- lapply(quantiles, function(x) apply(assays(SE)[["fraction"]], 1, quantile, x, na.rm = TRUE))

    print("Collect all information in a tibble")
    vars <- do.call(cbind, c(list(mean_af), list(mean_cov), list(rowData(SE)$VariantQuality), quantiles))
    vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, Quality = rowData(SE)$VariantQuality, Quantile1 = quantiles[[1]], Quantile2 = quantiles[[2]])

    print("Thresholding using the quantile approach.")
    if(length(quantiles) != 2) stop("Your quantiles are not of length 2.")
    if(length(thresholds) != 2) stop("Your thresholds are not of length 2.")
    voi_ch <- subset(vars, vars[,2] > min_coverage & vars[,4] < thresholds[1] & vars[,5] > thresholds[2])
    if(!is.null(min_quality)){
      voi_ch <- voi_ch[!is.na(voi_ch[,3]),]
      voi_ch <- subset(voi_ch, voi_ch[,3] > min_quality)
      }
    } else{
      print("Get the quantile of the VAF of each variant.")
      if(length(quantiles) > 1) stop("You are providing more than 1 quantile. You should only provide 1.")
      if(length(thresholds) > 1) stop("You are providing more than 1 threshold. You should only provide 1.")
      quantiles <- lapply(quantiles, function(x) apply(assays(SE)[["fraction"]], 1, quantile, x, na.rm = TRUE))
      quantiles <- quantiles[[1]]
      top_cells_values <- assays(SE)[["fraction"]]
      top_cells_values <- top_cells_values > top_VAF
      top_cells_values <- rowSums(top_cells_values)
      vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, Quality = rowData(SE)$VariantQuality, Quantile = quantiles, TopCells = top_cells_values)
      voi_ch <- subset(vars, Mean_Cov > min_coverage & Quantile <= thresholds[1] & TopCells >= top_cells)
      voi_ch <- voi_ch[!is.na(voi_ch$Quality),]
      voi_ch <- subset(voi_ch, Quality > min_quality)
    }
  voi_ch <- rownames(voi_ch)
  return(voi_ch)
}
