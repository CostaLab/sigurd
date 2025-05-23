#'VariantCorrelationHeatmap
#'@description
#'We generate a heatmap showing the correlation of somatic variants with the MT variants.
#'Packages I want to remove. I cannot see where they are used.
#'ggplot2 parallel rcompanion tidyr
#'@importFrom circlize colorRamp2
#'@importFrom ComplexHeatmap columnAnnotation rowAnnotation Heatmap draw
#'@importFrom grid gpar
#'@importFrom stats na.omit
#'@importFrom grDevices png dev.off
#'@param correlation_results Data.frame with the correlation results.
#'@param output_path Path to the output folder.
#'@param patient The patient for this heatmap.
#'@param min_alt_cells Minimum number of mutated cells needed, otherwise a correlation will not be plotted.
#'@param min_correlation Minimum correlation needed. 
#'@param width_use Width of the heatmap in px.
#'@param height_use Height of the heatmap in px.
#'@param padding_use Space around the heatmap in mm. If this is to low, the variant names might be cut off.
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
VariantCorrelationHeatmap <- function(correlation_results, output_path = NULL, patient, min_alt_cells = 5, min_correlation = 0.5,
                                      width_use = 2000, height_use = 2000, padding_use = c(165,165,2,2), verbose = TRUE){
  correlation_results$P_adj_logged <- -log10(correlation_results$P_adj)
  correlation_results <- subset(correlation_results, correlation_results$P_adj_logged > -log10(0.05))
  correlation_results <- subset(correlation_results, correlation_results$Cells_1_Alt >= min_alt_cells & correlation_results$Cells_2_Alt >= min_alt_cells)
  correlation_results <- subset(correlation_results, correlation_results$Corr > min_correlation)
  
  
  if(verbose) print("We get the unique variants.")
  somatic_uniques <- unique(correlation_results$Variant1)
  mt_uniques      <- unique(correlation_results$Variant2)
  
  
  if(verbose) print("Getting the maximum P value.")
  pvalue_max <- as.numeric(stats::na.omit(correlation_results$P_adj_logged))
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
  correlation_results$P_adj_logged[correlation_results$P_adj_logged == Inf] <- pvalue_max
  col_fun <- circlize::colorRamp2(c(0,pvalue_max), c("white", "red"))
  
  
  if(verbose) print("We set insignificant P values to NA.")
  correlation_results$P_adj_logged <- ifelse(correlation_results$P_adj_logged > -log10(0.05), correlation_results$P_adj_logged, NA)
  
  
  if(verbose) print("We generate a matrix with the adjusted P values.")
  p_values <- matrix(NA, nrow = length(somatic_uniques), ncol = length(mt_uniques))
  rownames(p_values) <- somatic_uniques
  colnames(p_values) <- mt_uniques
  for(i in 1:length(somatic_uniques)){
    correlation_results_subset <- subset(correlation_results, correlation_results$Variant1 == somatic_uniques[i])
    p_values_use <- correlation_results_subset$P_adj_logged
    names(p_values_use) <- correlation_results_subset$Variant2
    p_values[somatic_uniques[i],names(p_values_use)] <- p_values_use
  }
  
  
  if(verbose) print("Setting the column and row annotations for the heat map.")
  annotation_top <- ComplexHeatmap::columnAnnotation(Mutations = mt_uniques,
                                                     show_legend = FALSE, show_annotation_name = FALSE)
  annotation_left <- ComplexHeatmap::rowAnnotation(Mutations = somatic_uniques,
                                                   show_legend = FALSE, show_annotation_name = FALSE)
  
  if(verbose) print("Since we can have no results left after the subsetting, we check if the P value matrix has values.")
  if(all(dim(p_values) > 0)){
    if(verbose) print("Generating the actual heat map.")
    p1 <- ComplexHeatmap::Heatmap(p_values, name = "-log10(P)",
                                  column_title = paste0("Patient ", patient, "\nLogged adj. P values between the mutations"),
                                  row_title = "", show_row_names = TRUE, show_column_names = TRUE,
                                  col = col_fun, left_annotation = annotation_left, top_annotation = annotation_top,
                                  column_title_gp = grid::gpar(fontsize = 40), row_title_gp = grid::gpar(fontsize = 40),
                                  column_names_gp = grid::gpar(fontsize = 40), row_names_gp = grid::gpar(fontsize = 40),
                                  column_names_rot = 45,
                                  row_names_side = "left",
                                  heatmap_legend_param = list(labels_gp = grid::gpar(fontsize = 40), title_gp = grid::gpar(fontsize = 40, fontface = "bold")),
                                  cluster_columns = FALSE, cluster_rows = FALSE, use_raster = FALSE, show_row_dend = FALSE, show_column_dend = FALSE)
    
    
    if(!is.null(output_path)){
      if(verbose) print("Saving the png.")
      grDevices::png(paste0(output_path, "Correlation_Pvalue_", patient, ".png"), width = width_use, height = height_use, units = "px", type = "cairo", antialias = "none")
      ComplexHeatmap::draw(p1, padding = unit(padding_use, "mm"))
      grDevices::dev.off()
    } else{
      return(p1)
    }
  }
}
