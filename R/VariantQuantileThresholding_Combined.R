#'VariantQuantileThresholding_Combined
#'@description
#'We get variants of interest using the quantile thresholding. This function combines the functions VariantSelection_Quantile, VariantSelection_Group and VariantSelection_TopCells.
#'If you use top_cells and top_VAF, you have to only supply one quantil value (quantiles = 0.9, thresholds = 0).
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
#'@param top_cells The number of cells with at least top_VAF percent for a variant.
#'@param top_VAF The VAF for the top cells.
#'@param mean_allele_frequency The minimum mean allele frequency. Default = 0
#'@param group_of_interest The column data that divides the cells.
#'@param group1 The first group of interest. If set, the quantiles are only calculated for this group.
#'@param group2 The second group of interest.
#'@param group_factor How much higher has the mean allele frequency to be in group 1 when compared to group 2?
#'@param remove_nocall Should NoCall cells (consensus = 0) be disregarded during the analysis?
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
VariantQuantileThresholding_Combined <- function(SE, min_coverage = 2, quantiles = c(0.1, 0.9), thresholds = c(0.1, 0.9), top_cells = NULL, top_VAF = NULL, min_quality = NULL, mean_allele_frequency = 0,
                                        group_of_interest = NULL, group1 = NULL, group2 = NULL, group_factor = NULL, remove_nocall = TRUE, verbose = TRUE){
  if(remove_nocall){
    if(verbose) print("We set NoCall cells as NA.")
    nocall_check <- SummarizedExperiment::assays(SE)[["consensus"]] == 0
    SummarizedExperiment::assays(SE)[["fraction"]][nocall_check] <- NA
    SummarizedExperiment::assays(SE)[["coverage"]][nocall_check] <- NA
  }

  if(verbose) print("Get the mean allele frequency and coverage.")
  mean_af  <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["fraction"]], na.rm = TRUE)
  mean_cov <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["coverage"]], na.rm = TRUE)
  if(all(is.null(min_quality), is.numeric(min_quality))) stop("Error: Your minimum quality is not either NULL or a numeric.")
  if(all(is.null(mean_allele_frequency), is.numeric(mean_allele_frequency))) stop("Error: Your mean allele frequency is not either NULL or a numeric.")

  if(all(!is.null(group_of_interest), !is.null(group1), !is.null(group2))){
    if(!group_of_interest %in% colnames(SummarizedExperiment::colData(SE))) stop("Error: Your group_of_interest is not in the colData.")
    if(!group1 %in% SummarizedExperiment::colData(SE)[, group_of_interest]) stop("Error: Your group1 is not in the group_of_interest.")
    if(!group2 %in% SummarizedExperiment::colData(SE)[, group_of_interest]) stop("Error: Your group2 is not in the group_of_interest.")
    cells_group1 <- SummarizedExperiment::colData(SE)[, group_of_interest, drop = FALSE]
    cells_group1 <- cells_group1[cells_group1[, group_of_interest] == group1, , drop = FALSE]
    cells_group2 <- SummarizedExperiment::colData(SE)[, group_of_interest, drop = FALSE]
    cells_group2 <- cells_group2[cells_group2[, group_of_interest] == group2, , drop = FALSE]
    mean_af_group1 <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["fraction"]][, rownames(cells_group1)], na.rm = TRUE)
    mean_af_group2 <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["fraction"]][, rownames(cells_group2)], na.rm = TRUE)
    mean_af_group_check <- mean_af_group1 > (group_factor * mean_af_group2)

    if(verbose) print("Get the quantiles of the VAFs of each variant.")
    quantiles_group1 <- lapply(quantiles, function(x) apply(SummarizedExperiment::assays(SE)[["fraction"]][, rownames(cells_group1)], 1, stats::quantile, x, na.rm = TRUE))
    if(!is.null(min_quality)){
      vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, VariantQuality = SummarizedExperiment::rowData(SE)$VariantQuality, Quantile1 = quantiles_group1[[1]], Quantile2 = quantiles_group1[[2]])
      vars <- vars[!is.na(vars$VariantQuality), ]
      vars <- subset(vars, vars$VariantQuality > min_quality)
    } else{
      vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, Quantile1 = quantiles_group1[[1]], Quantile2 = quantiles_group1[[2]])
    }
    mean_af_group_check <- mean_af_group_check[rownames(vars)]
    vars <- vars[mean_af_group_check, , drop = FALSE]

    if(verbose) print("Thresholding using the quantile approach.")
    if(length(quantiles)  != 2) stop("Your quantiles are not of length 2.")
    if(length(thresholds) != 2) stop("Your thresholds are not of length 2.")
    vois <- subset(vars, vars$Mean_AF > mean_allele_frequency & vars$Mean_Cov > min_coverage & vars$Quantile1 < thresholds[1] & vars$Quantile2 > thresholds[2])

  } else if(any(is.null(top_cells), is.null(top_VAF))){
    if(verbose) print("Get the quantiles of the VAFs of each variant.")
    quantiles <- lapply(quantiles, function(x) apply(SummarizedExperiment::assays(SE)[["fraction"]], 1, quantile, x, na.rm = TRUE))

    if(!is.null(min_quality)){
      vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, VariantQuality = SummarizedExperiment::rowData(SE)$VariantQuality, Quantile1 = quantiles[[1]], Quantile2 = quantiles[[2]])
      vars <- vars[!is.na(vars$VariantQuality), ]
      vars <- subset(vars, vars$VariantQuality > min_quality)
    } else{
      vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, Quantile1 = quantiles[[1]], Quantile2 = quantiles[[2]])
    }

    if(verbose) print("Thresholding using the quantile approach.")
    if(length(quantiles) != 2) stop("Your quantiles are not of length 2.")
    if(length(thresholds) != 2) stop("Your thresholds are not of length 2.")
    vois <- subset(vars, vars$Mean_AF > mean_allele_frequency & vars$Mean_Cov > min_coverage & vars$Quantile1 < thresholds[1] & vars$Quantile2 > thresholds[2])
  } else{
    if(verbose) print("Get the quantile of the VAF of each variant.")
    if(length(quantiles) > 1) stop("You are providing more than 1 quantile. You should only provide 1.")
    if(length(thresholds) > 1) stop("You are providing more than 1 threshold. You should only provide 1.")
    quantiles <- lapply(quantiles, function(x) apply(SummarizedExperiment::assays(SE)[["fraction"]], 1, quantile, x, na.rm = TRUE))
    quantiles <- quantiles[[1]]
    top_cells_values <- SummarizedExperiment::assays(SE)[["fraction"]]
    top_cells_values <- top_cells_values > top_VAF
    top_cells_values <- Matrix::rowSums(top_cells_values, na.rm = TRUE)
    if(!is.null(min_quality)){
      vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, VariantQuality = SummarizedExperiment::rowData(SE)$VariantQuality, Quantile = quantiles, TopCells = top_cells_values)
      vars <- vars[!is.na(vars$VariantQuality), ]
      vars <- subset(vars, vars$VariantQuality > min_quality)
    } else{
      vars <- data.frame(Mean_AF = mean_af, Mean_Cov = mean_cov, Quantile = quantiles, TopCells = top_cells_values)
    }
    vois <- subset(vars, vars$Mean_Cov > min_coverage & vars$Quantile <= thresholds[1] & vars$TopCells >= top_cells)
  }
  vois <- rownames(vois)
  return(vois)
}
