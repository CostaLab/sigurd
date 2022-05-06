#'We get the variant information per cell.
#'@import dplyr SummarizedExperiment tibble tidyverse
#'@param se, voi_ch voi_ch
GetCellInfoPerVariant <- function(se, voi_ch){
  print("Generate matrices with coverage, allele frequency and reference / variant reads")
  cov_voi_mat <- assays(se)[["coverage"]][voi_ch,]
  af_voi_mat  <- assays(se)[["fraction"]][voi_ch,]

  print("Add coverage and allele frequency info from variants of interest to cells_tib.")
  cells_tib <- tibble(cell = colnames(se),
                      Mean_Cov = se$depth)
  for(voi in voi_ch){
    cells_tib <- cells_tib %>%
      left_join(as_tibble(assays(se)[["coverage"]][voi,], rownames = "cell"), by = "cell") %>%
      left_join(as_tibble(assays(se)[["fraction"]][voi,], rownames = "cell"), by = "cell")
    colnames(cells_tib) <- gsub("value.x", paste0("cov_", voi), colnames(cells_tib))
    colnames(cells_tib) <- gsub("value.y", paste0("af_", voi), colnames(cells_tib))
  }
  return(cells_tib)
}
