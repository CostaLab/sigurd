#'VariantSelection_Group
#'@description
#'We get variants of interest by selecting variants with a high VAF difference between two groups.
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
#'@param group_of_interest The column data that divides the cells.
#'@param group1 The first group of interest. If set, the quantiles are only calculated for this group.
#'@param group2 The second group of interest.
#'@param group_factor How much higher has the mean allele frequency to be in group 1 when compared to group 2?
#'@param remove_nocall Should NoCall cells (consensus = 0) be disregarded during the analysis?
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
VariantSelection_Group <- function(SE, min_coverage = 2, quantiles = c(0.1, 0.9), thresholds = c(0.1, 0.9), min_quality = NULL, mean_allele_frequency = 0,
                                   group_of_interest = NULL, group1 = NULL, group2 = NULL, group_factor = 5, remove_nocall = TRUE, verbose = TRUE){
  if(remove_nocall){
    if(verbose) print("We set NoCall cells as NA.")
    nocall_check <- SummarizedExperiment::assays(SE)[["consensus"]] == 0
    SummarizedExperiment::assays(SE)[["fraction"]][nocall_check] <- NA
    SummarizedExperiment::assays(SE)[["coverage"]][nocall_check] <- NA
  }
  
  if(verbose) print("Get the mean allele frequency and coverage.")
  mean_af  <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["fraction"]], na.rm = TRUE)
  mean_cov <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["coverage"]], na.rm = TRUE)

  if(!group_of_interest %in% colnames(SummarizedExperiment::colData(SE))) stop("Error: Your group_of_interest is not in the colData.")
  if(!group1 %in% SummarizedExperiment::colData(SE)[, group_of_interest]) stop("Error: Your group1 is not in the group_of_interest.")
  if(!group2 %in% SummarizedExperiment::colData(SE)[, group_of_interest]) stop("Error: Your group2 is not in the group_of_interest.")
  
  # Processing group1 and group2
  cells_group1 <- SummarizedExperiment::colData(SE)[, group_of_interest, drop = FALSE]
  cells_group1 <- cells_group1[cells_group1[, group_of_interest] == group1, , drop = FALSE]
  cells_group2 <- SummarizedExperiment::colData(SE)[, group_of_interest, drop = FALSE]
  cells_group2 <- cells_group2[cells_group2[, group_of_interest] == group2, , drop = FALSE]

  mean_af_group1 <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["fraction"]][, rownames(cells_group1)], na.rm = TRUE)
  mean_af_group2 <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["fraction"]][, rownames(cells_group2)], na.rm = TRUE)

  mean_af_group_check <- mean_af_group1 > (group_factor * mean_af_group2)

  if(verbose) print("Get the quantiles of the VAFs of each variant.")
  quantiles_group1 <- lapply(quantiles, function(x) apply(SummarizedExperiment::assays(SE)[["fraction"]][, rownames(cells_group1)], 1, stats::quantile, x, na.rm = TRUE))

  # Apply min_quality filtering
  if(!is.null(min_quality)){
    vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, VariantQuality = SummarizedExperiment::rowData(SE)$VariantQuality,
                       Quantile1 = quantiles_group1[[1]], Quantile2 = quantiles_group1[[2]])
    vars <- subset(vars, !is.na(vars$VariantQuality) & vars$VariantQuality > min_quality)
  } else{
    vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, Quantile1 = quantiles_group1[[1]], Quantile2 = quantiles_group1[[2]])
  }

  mean_af_group_check <- mean_af_group_check[rownames(vars)]
  vars <- vars[mean_af_group_check, , drop = FALSE]

  if(verbose) print("Thresholding using the quantile approach.")
  vois <- subset(vars, vars$Mean_AF > mean_allele_frequency & vars$Mean_Cov > min_coverage & vars$Quantile1 < thresholds[1] & vars$Quantile2 > thresholds[2])

  return(rownames(vois))
}
