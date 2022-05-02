VariantQuantileThresholding <- function(se, min_coverage = 2){
  # This function is adapted from the Peter van Galen.
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(SummarizedExperiment))

  print("Get the mean allele frequency and coverage.")
  mean_af <- rowMeans(assays(se)[["fraction"]])
  mean_cov <- rowMeans(assays(se)[["coverage"]])

  print("Get the quantiles of the VAFs of each variant.")
  quantiles <- c("q01" = 0.01, "q10" = 0.1, "q50" = 0.5, "q90" = 0.9, "q99" = 0.99)
  start_time <- Sys.time()
  quantiles <- lapply(quantiles, function(x) apply(assays(se)[["fraction"]], 1, quantile, x))
  Sys.time() - start_time


  print("Collect all information in a tibble")
  vars_tib <- as_tibble(do.call(cbind, c(list(mean_af), list(mean_cov), quantiles)), rownames = "var")
  colnames(vars_tib)[2] <- "mean_af"
  colnames(vars_tib)[3] <- "mean_cov"

  print("We add the number of cells that exceed the VAF thresholds.")
  vars_tib <- vars_tib %>%
      mutate(n0  = apply(assays(se)[["fraction"]], 1, function(x) sum(x == 0))) %>%
      mutate(n1  = apply(assays(se)[["fraction"]], 1, function(x) sum(x >  0.01))) %>%
      mutate(n5  = apply(assays(se)[["fraction"]], 1, function(x) sum(x >  0.05))) %>%
      mutate(n10 = apply(assays(se)[["fraction"]], 1, function(x) sum(x >  0.1))) %>%
      mutate(n50 = apply(assays(se)[["fraction"]], 1, function(x) sum(x >  0.5)))

  print("Thresholding using the quantile approach.")
  voi_ch <- vars_tib %>% filter(.[,"mean_cov"] > min_coverage,
                                .[,"q10"] < 0.1,
                                .[,"q90"] > 0.9) %>% .$var
  return(voi_ch)
}
