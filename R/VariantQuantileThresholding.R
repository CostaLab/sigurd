#'We get variants of interest using the quantile thresholding.
#'@import SummarizedExperiment tidyverse
#'@param se SummarizedExperiment object.
#'@param min_coverage Minimum coverage needed.
#'@export
VariantQuantileThresholding <- function(se, min_coverage = 2){
  # This function is adapted from the Peter van Galen.
  print("Get the mean allele frequency and coverage.")
  mean_af <- rowMeans(assays(se)[["fraction"]], na.rm = TRUE)
  mean_cov <- rowMeans(assays(se)[["coverage"]], na.rm = TRUE)

  print("Get the quantiles of the VAFs of each variant.")
  quantiles <- c("q01" = 0.01, "q10" = 0.1, "q50" = 0.5, "q90" = 0.9, "q99" = 0.99)
  start_time <- Sys.time()
  quantiles <- lapply(quantiles, function(x) apply(assays(se)[["fraction"]], 1, quantile, x, na.rm = TRUE))
  Sys.time() - start_time


  print("Collect all information in a tibble")
  vars_tib <- as_tibble(do.call(cbind, c(list(mean_af), list(mean_cov), quantiles)), rownames = "var")
  colnames(vars_tib)[2] <- "mean_af"
  colnames(vars_tib)[3] <- "mean_cov"

  
  print("Thresholding using the quantile approach.")
  voi_ch <- vars_tib %>% filter(.[,"mean_cov"] > min_coverage,
                                .[,"q10"] < 0.1,
                                .[,"q90"] > 0.9) %>% .$var
  return(voi_ch)
}
