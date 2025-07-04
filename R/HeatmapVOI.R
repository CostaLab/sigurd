#'HeatmapVoi
#'@description
#'We plot a heatmap of a set of Variants Of Interest using the Variant Allele Frequency values of a SummarizedExperiment object.
#'@importFrom ComplexHeatmap columnAnnotation Heatmap
#'@importFrom SummarizedExperiment assays colData 
#'@importFrom grid gpar unit
#'@importFrom circlize colorRamp2
#'@importFrom scales hue_pal
#'@importFrom Matrix colSums
#'@param SE SummarizedExperiment object.
#'@param voi Variants Of Interest.
#'@param annotation_trait Cell Annotation at the bottom of the heat map.
#'@param column_title The title of the heat map. Default = NULL
#'@param minimum_coverage The minimum coverage per cell to be plotted.
#'@param sort_cells Should the cells be sorted by ordering the cells according the the largest clones? FALSE uses default complete clustering with euclidean distance.
#'@param cluster_variants Should the variants be clustered? Default FALSE.
#'@param cluster_variants_distance The distance for clustering the variants. Default euclidean.
#'@param cluster_variants_method The distance for the variant clustering. Default complete
#'@param remove_empty_cells Should cells that have a fraction of 0 for all variants be removed? Default = FALSE
#'@param minimum_allele_freq Minimum allele frequency to include a cell.
#'@export
HeatmapVoi <- function(SE, voi, annotation_trait = NULL, column_title = NULL, minimum_coverage = 0, sort_cells = FALSE, remove_empty_cells = FALSE, minimum_allele_freq = 0, cluster_variants = FALSE, cluster_variants_distance = "euclidean", cluster_variants_method = "complete"){
  if(minimum_coverage > 0){
    coverage_test <- SummarizedExperiment::assays(SE)[["coverage"]][voi,]
    coverage_test <- coverage_test > minimum_coverage
    coverage_test <- Matrix::colSums(coverage_test)
    coverage_test <- coverage_test == length(voi)
    SE <- SE[,coverage_test]
  }
  fraction <- SummarizedExperiment::assays(SE)[["fraction"]][voi,]
  fraction[is.na(fraction)] <- 0
  if(length(voi) == 1){
    fraction <- t(as.matrix(fraction))
    rownames(fraction) <- voi
  } else if(length(voi) > 1){
    fraction <- as.matrix(fraction)
  }
  # We remove cells that are negative for all variants.
  if(remove_empty_cells){
    cell_check <- Matrix::colSums(fraction > minimum_allele_freq) > 0
    fraction <- fraction[, cell_check, drop = FALSE]
    SE <- SE[, cell_check]
  }

  # We get a different column title.
  if(is.null(column_title)){
    column_title <- "Cells"
  }

  if(sort_cells){
    cell_order <- fraction
    for(variant in rev(rownames(fraction))){
      cell_order <- cell_order[, order(cell_order[variant,], decreasing = TRUE), drop = FALSE]
    }
    cell_order <- colnames(cell_order)
    fraction <- fraction[, cell_order, drop = FALSE]
    SE <- SE[, cell_order]
  }

  if(!is.null(annotation_trait)){
    colours_use <- scales::hue_pal()(length(unique(SummarizedExperiment::colData(SE)[,annotation_trait])))
    names(colours_use) <- unique(SummarizedExperiment::colData(SE)[,annotation_trait])
    ha <- ComplexHeatmap::columnAnnotation(annotation_trait = SummarizedExperiment::colData(SE)[,annotation_trait],
                                           col = list(annotation_trait = colours_use))
  } else if(is.null(annotation_trait)){
    ha <- NULL
  }

  heatmap_voi <- ComplexHeatmap::Heatmap(fraction,
                                         column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                                         row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                                         row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                                         col = circlize::colorRamp2(seq(0, ceiling(max(fraction, na.rm = TRUE)), length.out = 9),
                                                                    c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
                                         show_row_names = TRUE, show_column_names = FALSE, cluster_columns = ifelse(sort_cells, FALSE, TRUE), cluster_rows = cluster_variants, name = "VAF",
					 clustering_distance_rows = cluster_variants_distance, clustering_method_rows = cluster_variants_method,
                                         heatmap_legend_param = list(border = "#000000"),
                                         bottom_annotation = ha, border = TRUE, use_raster = FALSE,
                                         column_title = column_title,
                                         row_title = "Variants")
  return(heatmap_voi)
}
