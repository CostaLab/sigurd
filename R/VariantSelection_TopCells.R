#'VariantSelection_TopCells
#'@description
#'We get variants of interest by selecting a small number of cells with a high VAF for a variant.
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
#'@param remove_nocall Should NoCall cells (consensus = 0) be disregarded during the analysis?
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
VariantSelection_TopCells <- function(SE, min_coverage = 2, quantiles = 0.9, thresholds = 0, top_cells = 10, top_VAF = 0.5, min_quality = NULL, remove_nocall = TRUE, verbose = TRUE){
  if(remove_nocall){
    if(verbose) print("We set NoCall cells as NA.")
    nocall_check <- SummarizedExperiment::assays(SE)[["consensus"]] == 0
    SummarizedExperiment::assays(SE)[["fraction"]][nocall_check] <- NA
    SummarizedExperiment::assays(SE)[["coverage"]][nocall_check] <- NA
  }

  if(verbose) print("Get the mean allele frequency and coverage.")
  mean_af  <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["fraction"]], na.rm = TRUE)
  mean_cov <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["coverage"]], na.rm = TRUE)

  if(verbose) print("Get the quantile of the VAF of each variant.")
  quantiles <- apply(SummarizedExperiment::assays(SE)[["fraction"]], 1, quantile, quantiles, na.rm = TRUE)

  top_cells_values <- SummarizedExperiment::assays(SE)[["fraction"]] > top_VAF
  top_cells_values <- Matrix::rowSums(top_cells_values, na.rm = TRUE)

  # Apply min_quality filtering
  if(!is.null(min_quality)){
    vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, VariantQuality = SummarizedExperiment::rowData(SE)$VariantQuality, Quantile = quantiles, TopCells = top_cells_values)
    vars <- subset(vars, !is.na(vars$VariantQuality) & vars$VariantQuality > min_quality)
  } else{
    vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, Quantile = quantiles, TopCells = top_cells_values)
  }

  vois <- subset(vars, vars$Mean_Cov > min_coverage & vars$Quantile <= thresholds & vars$TopCells >= top_cells)

  return(rownames(vois))
}
