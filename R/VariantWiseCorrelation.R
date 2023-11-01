#'VariantWiseCorrelation
#'@description
#'We correlate the variants with each other using the Pearson correlation.
#'This function calls CalculateCorrelationPValue to perform the actual correlation.
#'Packages I want to remove.
#'SummarizedExperiment
#'@import Matrix parallel
#'@param variants_list List of fraction values.
#'@param n_cores Number of cores you want to use. Numeric.
#'@param p_value_adjustment Method for P value adjustment. See p.adjust for details.
#'@export
VariantWiseCorrelation <- function(variants_list, n_cores = 1, p_value_adjustment = "fdr"){
  # We correlate the somatic variants with each other and the MT variants.
  # Since we have tens of thousands of MT variants, we do not correlate them with each other.
  variants <- names(variants_list)
  MT_check <- grepl("^chrM_|^chrMT_|^MT_|^chrMt_|^chrM_|^Mt_", variants)
  variants <- variants[!MT_check]

  results_total <- c()
  number_of_variants <- length(variants)
  for(i in 1:length(variants)){
    variant_use <- variants[i]
    print(paste0("Correlating Variant: ", variant_use, ", ", i, " out of ", number_of_variants))
    variants_values_use <- variants_list[[variant_use]]
    variants_list_use <- variants_list[names(variants_list) != variant_use]
    all_variants <- names(variants_list_use)
    results <- parallel::mclapply(X = all_variants, CalculateCorrelationPValue, variant_values = variants_values_use, all_variants_list = variants_list_use, mc.cores = n_cores)
    results <- do.call("rbind", results)
    results <- data.frame(Variant1 = variant_use, Variant2 = all_variants, P = results[,1], Corr = results[,2],
                          Cells_1_Alt = results[,3], Cells_1_Ref = results[,4], Cells_2_Alt = results[,5], Cells_2_Ref = results[,6])
    results_total <- rbind(results_total, results)
  }

  print("We remove the NA P values.")
  results_total <- results_total[!is.na(results_total$P),]

  print("We remove the negative corrlated SNPs.")
  results_total <- subset(results_total, Corr > 0)

  print(paste0("Adjusting P values using ", p_value_adjustment, "."))
  results_total$P_adj <- p.adjust(results_total$P, method = p_value_adjustment)
  return(results_total)
}
