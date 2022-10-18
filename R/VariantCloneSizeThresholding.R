#'VariantCloneSizeThresholding
#'@description
#'We get variants of interest using a clone size thresholding.
#'Source: https://github.com/petervangalen/MAESTER-2021
#'@import dplyr Matrix SummarizedExperiment tidyverse
#'@param se SummarizedExperiment object.
#'@param min_coverage Minimum coverage a variant needs to have.
#'@param fraction_negative_cells The fraction of negative cells needed. 
#'@param min_clone_size minimum number of cells.
#'@param vaf_threshold Variant Allele Threshold. Cells above this threshold are considered mutated.
#'@export
VariantCloneSizeThresholding <- function(se, min_coverage = 2, fraction_negative_cells = 0.9, min_clone_size = 10, vaf_threshold = 0.5){
  # This function is adapted from the Peter van Galen.
  print("Get the mean allele frequency and coverage.")
  mean_af <- rowMeans(assays(se)[["fraction"]], na.rm = TRUE)
  mean_cov <- rowMeans(assays(se)[["coverage"]], na.rm = TRUE)

  print("Collect all information in a tibble")
  #vars_tib <- as_tibble(do.call(cbind, c(list(mean_af), list(mean_cov))), rownames = "var")
  vars <- do.call(cbind, c(list(mean_af), list(mean_cov)))
  #colnames(vars_tib)[2] <- "mean_af"
  #colnames(vars_tib)[3] <- "mean_cov"
  colnames(vars) <- c("mean_af", "mean_cov")
  
  print("We add the number of cells that exceed the VAF thresholds.")
  #vars_tib <- vars_tib %>% 
  #  mutate(n0 = apply(assays(se)[["fraction"]], 1, function(x) sum(x > vaf_threshold, na.rm = TRUE))) %>%
  #  mutate(VAF_threshold = apply(assays(se)[["fraction"]], 1, function(x) sum(x > vaf_threshold, na.rm = TRUE)))
  n0 <- apply(assays(se)[["fraction"]], 1, function(x) sum(x > vaf_threshold, na.rm = TRUE))
  VAF_threshold <- apply(assays(se)[["fraction"]], 1, function(x) sum(x > vaf_threshold, na.rm = TRUE))
  vars <- cbind(vars, n0, VAF_threshold)
  
  print("Thresholding using the clone size approach.")
  #voi_ch <- subset(vars_tib, mean_cov > min_coverage & 
  #                           n0 > ceiling(fraction_negative_cells * ncol(se)) &
  #                           VAF_threshold > min_clone_size)$var
  voi_ch <- subset(vars, mean_cov > min_coverage & 
                         n0 > ceiling(fraction_negative_cells * ncol(se)) &
                         VAF_threshold > min_clone_size)
  voi_ch <- rownames(voi_ch)
  return(voi_ch)
}
