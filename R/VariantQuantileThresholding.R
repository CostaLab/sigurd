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
#'@param mean_allele_frequency The minimum mean allele frequency. Default = 0
#'@param group_of_interest The column data that divides the cells.
#'@param group1 The first group of interest.
#'@param group2 The second group of interest.
#'@param group_factor How much higher has the mean allele frequency to be in group 1 when compared to group 2?
#'@export
VariantQuantileThresholding <- function(SE, min_coverage = 2, quantiles = c(0.1, 0.9), thresholds = c(0.1, 0.9), top_cells = NULL, top_VAF = NULL, min_quality = 30, mean_allele_frequency = 0,
                                        group_of_interest = NULL, group1 = NULL, group2 = NULL, group_factor = NULL){
  print("Get the mean allele frequency and coverage.")
  mean_af <- rowMeans(assays(SE)[["fraction"]], na.rm = TRUE)
  mean_cov <- rowMeans(assays(SE)[["coverage"]], na.rm = TRUE)
  if(all(!is.null(group_of_interest), !is.null(group1), !is.null(group2))){
    if(!group_of_interest %in% colnames(colData(SE))) stop("Error: Your group_of_interest is not in the colData.")
    if(!group1 %in% colData(SE)[,group_of_interest]) stop("Error: Your group1 is not in the group_of_interest.")
    if(!group2 %in% colData(SE)[,group_of_interest]) stop("Error: Your group2 is not in the group_of_interest.")
    cells_group1 <- colData(SE)[,group_of_interest, drop = FALSE]
    cells_group1 <- cells_group1[cells_group1[, group_of_interest] == group1, , drop = FALSE]
    cells_group2 <- colData(SE)[,group_of_interest, drop = FALSE]
    cells_group2 <- cells_group2[cells_group2[, group_of_interest] == group2, , drop = FALSE]
    mean_af_group1 <- rowMeans(assays(SE)[["fraction"]][,rownames(cells_group1)], na.rm = TRUE)
    mean_af_group2 <- rowMeans(assays(SE)[["fraction"]][,rownames(cells_group2)], na.rm = TRUE)
    mean_af_group_check <- mean_af_group1 > (group_factor * mean_af_group2)

    print("Get the quantiles of the VAFs of each variant.")
    quantiles <- lapply(quantiles, function(x) apply(assays(SE)[["fraction"]], 1, quantile, x, na.rm = TRUE))
    # vars <- do.call(cbind, c(list(mean_af), list(mean_cov), list(rowData(SE)$VariantQuality), quantiles))
    vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, Quality = rowData(SE)$VariantQuality, Quantile1 = quantiles[[1]], Quantile2 = quantiles[[2]])
    vars <- vars[mean_af_group_check,]

    print("Thresholding using the quantile approach.")
    if(length(quantiles) != 2) stop("Your quantiles are not of length 2.")
    if(length(thresholds) != 2) stop("Your thresholds are not of length 2.")
    #voi_ch <- subset(vars, vars[,1] > mean_allele_frequency & vars[,2] > min_coverage & vars[,4] < thresholds[1] & vars[,5] > thresholds[2])
    voi_ch <- subset(vars, Mean_AF > mean_allele_frequency & Mean_Cov > min_coverage & Quantile1 < thresholds[1] & Quantile2 > thresholds[2])
    if(!is.null(min_quality)){
      voi_ch <- voi_ch[!is.na(voi_ch$Quality),]
      voi_ch <- subset(voi_ch, Quality > min_quality)
    }
  } else if(any(is.null(top_cells), is.null(top_VAF))){
    print("Get the quantiles of the VAFs of each variant.")
    quantiles <- lapply(quantiles, function(x) apply(assays(SE)[["fraction"]], 1, quantile, x, na.rm = TRUE))

    # vars <- do.call(cbind, c(list(mean_af), list(mean_cov), list(rowData(SE)$VariantQuality), quantiles))
    vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, Quality = rowData(SE)$VariantQuality, Quantile1 = quantiles[[1]], Quantile2 = quantiles[[2]])

    print("Thresholding using the quantile approach.")
    if(length(quantiles) != 2) stop("Your quantiles are not of length 2.")
    if(length(thresholds) != 2) stop("Your thresholds are not of length 2.")
    #voi_ch <- subset(vars, vars[,1] > mean_allele_frequency & vars[,2] > min_coverage & vars[,4] < thresholds[1] & vars[,5] > thresholds[2])
    voi_ch <- subset(vars, Mean_AF > mean_allele_frequency & Mean_Cov > min_coverage & Quantile1 < thresholds[1] & Quantile2 > thresholds[2])
    if(!is.null(min_quality)){
      voi_ch <- voi_ch[!is.na(voi_ch$Quality),]
      voi_ch <- subset(voi_ch, Quality > min_quality)
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
    if(!is.null(min_quality)){
      voi_ch <- voi_ch[!is.na(voi_ch$Quality),]
      voi_ch <- subset(voi_ch, Quality > min_quality)
    }
  }
  voi_ch <- rownames(voi_ch)
  return(voi_ch)
}
