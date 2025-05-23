#'We get the variant information per cell.
# #'@import dplyr SummarizedExperiment tibble tidyverse
#'@importFrom SummarizedExperiment assays
#'@importFrom dplyr left_join %>%
#'@importFrom tibble tibble as_tibble
#'@param se SummarizedExperiment object.
#'@param voi_ch Variants of interest.
#'@param verbose Should the function be verbose? Default = FALSE
#'@export
GetCellInfoPerVariant <- function(se, voi_ch, verbose = FALSE){
  if(verbose) print("Generate matrices with coverage, allele frequency and reference / variant reads")
  cov_voi_mat <- SummarizedExperiment::assays(se)[["coverage"]][voi_ch,]
  af_voi_mat  <- SummarizedExperiment::assays(se)[["fraction"]][voi_ch,]

  if(verbose) print("Add coverage and allele frequency info from variants of interest to cells_tib.")
  cells_tib <- tibble::tibble(cell = colnames(se), Mean_Cov = se$depth)
  for(voi in voi_ch){
    cells_tib <- cells_tib %>%
      dplyr::left_join(tibble::as_tibble(SummarizedExperiment::assays(se)[["coverage"]][voi,], rownames = "cell"), by = "cell") %>%
      dplyr::left_join(tibble::as_tibble(SummarizedExperiment::assays(se)[["fraction"]][voi,], rownames = "cell"), by = "cell")
    colnames(cells_tib) <- gsub("value.x", paste0("cov_", voi), colnames(cells_tib))
    colnames(cells_tib) <- gsub("value.y", paste0("af_", voi), colnames(cells_tib))
  }
  return(cells_tib)
}
