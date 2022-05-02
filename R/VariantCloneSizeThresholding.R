VariantQuantileThresholding <- function(se, min_coverage = 2, fraction_negative_cells = 0.9, min_clone_size = 10, vaf_threshold = 0.5){
  # This function is adapted from the Peter van Galen.
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(SummarizedExperiment))

  print("Get the mean allele frequency and coverage.")
  mean_af <- rowMeans(assays(se)[["fraction"]])
  mean_cov <- rowMeans(assays(se)[["coverage"]])

  print("Collect all information in a tibble")
  vars_tib <- as_tibble(do.call(cbind, c(list(mean_af), list(mean_cov))), rownames = "var")
  colnames(vars_tib)[2] <- "mean_af"
  colnames(vars_tib)[3] <- "mean_cov"

  print("We add the number of cells that exceed the VAF thresholds.")
  vars_tib <- vars_tib %>%
      mutate(n0  = apply(assays(se)[["fraction"]], 1, function(x) sum(x == 0)))    %>%
      mutate(n1  = apply(assays(se)[["fraction"]], 1, function(x) sum(x >  0.01))) %>%
      mutate(n5  = apply(assays(se)[["fraction"]], 1, function(x) sum(x >  0.05))) %>%
      mutate(n10 = apply(assays(se)[["fraction"]], 1, function(x) sum(x >  0.1)))  %>%
      mutate(n50 = apply(assays(se)[["fraction"]], 1, function(x) sum(x >  0.5)))  %>%
      mutate(n75 = apply(assays(se)[["fraction"]], 1, function(x) sum(x >  0.75))) %>%
      mutate(n90 = apply(assays(se)[["fraction"]], 1, function(x) sum(x >  0.9)))


  print("Thresholding using the quantile approach.")
  voi_ch <- vars_tib %>% filter(.[,"mean_cov"] > min_coverage,
                                .[,"n0"] > ceiling(fraction_negative_cells * ncol(se)),
                                .[,paste0("n", vaf_threshold*100)] > min_clone_size) %>% .$var
  return(voi_ch)
}
