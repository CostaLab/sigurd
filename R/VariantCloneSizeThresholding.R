#'@import Matrix SummarizedExperiment tidyverse
#'@param se min_coverage fraction_negative_cells min_clone_size vaf_threshold
VariantCloneSizeThresholding <- function(se, min_coverage = 2, fraction_negative_cells = 0.9, min_clone_size = 10, vaf_threshold = 0.5){
  # This function is adapted from the Peter van Galen.
  print("Get the mean allele frequency and coverage.")
  mean_af <- rowMeans(assays(se)[["fraction"]], na.rm = TRUE)
  mean_cov <- rowMeans(assays(se)[["coverage"]], na.rm = TRUE)

  print("Collect all information in a tibble")
  vars_tib <- as_tibble(do.call(cbind, c(list(mean_af), list(mean_cov))), rownames = "var")
  colnames(vars_tib)[2] <- "mean_af"
  colnames(vars_tib)[3] <- "mean_cov"

  print("We add the number of cells that exceed the VAF thresholds.")
  vars_tib <- vars_tib %>% 
    mutate(n0  = apply(assays(se)[["fraction"]], 1, function(x) sum(x > vaf_threshold, na.rm = TRUE))) %>%
    mutate(VAF_threshold  = apply(assays(se)[["fraction"]], 1, function(x) sum(x > vaf_threshold, na.rm = TRUE)))

  print("Thresholding using the clone size approach.")
  voi_ch <- vars_tib %>% filter(.[,"mean_cov"] > min_coverage,
                                .[,"n0"] > ceiling(fraction_negative_cells * ncol(se)),
                                .[,paste0("VAF_threshold")] > min_clone_size) %>% .$var
  return(voi_ch)
}
