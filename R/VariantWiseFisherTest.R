#'VariantWiseFisherTest
#'@description
#'We perform the Fisher test to determine which variants are associated.
#'This function calls CalculateFisherTestPValue to perform the actual testing.
#'@importFrom parallel mclapply
#'@importFrom stats p.adjust
#'@param variants_list List of fraction values.
#'@param n_cores Number of cores you want to use. Numeric.
#'@param p_value_adjustment Method for P value adjustment. See p.adjust for details.
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
VariantWiseFisherTest <- function(variants_list, n_cores = 1, p_value_adjustment = "fdr", verbose = TRUE){
  # We correlate the somatic variants with each other and the MT variants.
  # Since we have tens of thousands of MT variants, we do not correlate them with each other.
  variants <- names(variants_list)
  MT_check <- grepl("^chrM_|^chrMT|^MT|^chrMt|^chrM|^Mt", variants)
  variants <- variants[!MT_check]

  results_total <- c()
  number_of_variants <- length(variants)
  for(i in 1:number_of_variants){
    variant_use <- variants[i]
    if(verbose) print(paste0("Testing Variant: ", variant_use, ", ", i, " out of ", number_of_variants))
    variants_values_use <- variants_list[[variant_use]]
    variants_list_use <- variants_list[names(variants_list) != variant_use]
    all_variants <- names(variants_list_use)
    results <- parallel::mclapply(X = all_variants, CalculateFisherTestPValue, variant_values = variants_values_use, all_variants_list = variants_list_use, mc.cores = n_cores)
    results <- do.call("rbind", results)
    results <- data.frame(Variant1 = variant_use, Variant2 = all_variants, P = results[,1], OddsRatio = results[,2],
                          Cells_Alt_1_2 = results[,3], Cells_Alt_1_Ref_2 = results[,4], Cells_Alt_2_Ref_1 = results[,5], Cells_Ref_1_2 = results[,6])
    results_total <- rbind(results_total, results)
  }

  if(verbose) print("We remove the NA P values.")
  results_total <- results_total[!is.na(results_total$P),]

  if(verbose) print(paste0("Adjusting P values using ", p_value_adjustment, "."))
  results_total$P_adj <- stats::p.adjust(results_total$P, method = p_value_adjustment)
  rownames(results_total) <- NULL
  return(results_total)
}
