#'VariantFisherTestHeatmap
#'@description
#'We generate a heatmap showing the Fisher test of somatic variants with the MT variants.
#'Packages I want to remove.
#'@importFrom ComplexHeatmap columnAnnotation rowAnnotation Heatmap
#'@importFrom circlize colorRamp2
#'@importFrom grid gpar
#'@importFrom stats na.omit
#'@param fisher_results Data.frame with the correlation results.
#'@param patient The patient for this heatmap.
#'@param min_alt_cells Minimum number of mutated cells needed, otherwise an association will not be plotted.
#'@param min_oddsratio Minimum correlation needed. 
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
VariantFisherTestHeatmap <- function(fisher_results, patient, min_alt_cells = 5, min_oddsratio = 1, verbose = TRUE){
  fisher_results$P_adj_logged <- -log10(fisher_results$P_adj)
  fisher_results <- subset(fisher_results, fisher_results$P_adj_logged > -log10(0.05))
  fisher_results <- subset(fisher_results, fisher_results$Cells_Alt_1_2 >= min_alt_cells)
  fisher_results <- subset(fisher_results, fisher_results$OddsRatio > min_oddsratio)


  if(verbose) print("We get the unique variants.")
  somatic_uniques <- unique(fisher_results$Variant1)
  mt_uniques      <- unique(fisher_results$Variant2)


  if(verbose) print("Getting the maximum P value.")
  pvalue_max <- as.numeric(stats::na.omit(fisher_results$P_adj_logged))
  if(length(pvalue_max) > 1){
    pvalue_max <- pvalue_max[pvalue_max != Inf]
    if(length(pvalue_max) >= 1){
      pvalue_max <- max(pvalue_max[!is.na(pvalue_max)])
      pvalue_max <- ifelse(pvalue_max == 0, 100, pvalue_max)
    } else{
      pvalue_max <- 100
    }
  } else{
    pvalue_max <- max(pvalue_max, 100)
  }
  fisher_results$P_adj_logged[fisher_results$P_adj_logged == Inf] <- pvalue_max
  col_fun <- circlize::colorRamp2(c(0,pvalue_max), c("white", "red"))
  
  
  if(verbose) print("We set insignificant P values to NA.")
  fisher_results$P_adj_logged <- ifelse(fisher_results$P_adj_logged > -log10(0.05), fisher_results$P_adj_logged, NA)
  
  
  if(verbose) print("We generate a matrix with the adjusted P values.")
  p_values <- matrix(NA, nrow = length(somatic_uniques), ncol = length(mt_uniques))
  rownames(p_values) <- somatic_uniques
  colnames(p_values) <- mt_uniques
  for(i in 1:length(somatic_uniques)){
    fisher_results_subset <- subset(fisher_results, fisher_results$Variant1 == somatic_uniques[i])
    p_values_use <- fisher_results_subset$P_adj_logged
    names(p_values_use) <- fisher_results_subset$Variant2
    p_values[somatic_uniques[i], names(p_values_use)] <- p_values_use
  }
  
  
  if(verbose) print("Setting the column and row annotations for the heat map.")
  annotation_top <- ComplexHeatmap::columnAnnotation(Mutations = mt_uniques,
                                                     show_legend = FALSE, show_annotation_name = FALSE)
  annotation_left <- ComplexHeatmap::rowAnnotation(Mutations = somatic_uniques,
                                                   show_legend = FALSE, show_annotation_name = FALSE)
  
  if(verbose) print("Since we might have no results left after the subsetting, we check if the P value matrix has values.")
  if(all(dim(p_values) > 0)){
    if(verbose) print("Generating the actual heat map.")
    p <- ComplexHeatmap::Heatmap(p_values, name = "-log10(P)",
                                 column_title = paste0("Patient ", patient, "\nLogged adj. P values between the variants"),
                                 row_title = "", show_row_names = TRUE, show_column_names = TRUE,
                                 col = col_fun, left_annotation = annotation_left, top_annotation = annotation_top,
                                 column_names_rot = 45, row_names_side = "left",
                                 column_names_gp = grid::gpar(hjust = 1),
                                 cluster_columns = FALSE, cluster_rows = FALSE, use_raster = FALSE, show_row_dend = FALSE, show_column_dend = FALSE)
    return(p)
  } else{
    return(NULL)
  }
}
