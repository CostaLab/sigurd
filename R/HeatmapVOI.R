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
#'@param remove_empty_cells Should cells that have a fraction of 0 for all variants be removed? Default = FALSE
#'@export
HeatmapVoi <- function(SE, voi, annotation_trait = NULL, column_title = NULL, remove_empty_cells = FALSE){

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
    cell_check <- Matrix::colSums(fraction > 0) > 0
    fraction <- fraction[,cell_check, drop = FALSE]
  }
  if(!is.null(annotation_trait)){
    colours_use <- scales::hue_pal()(length(unique(SummarizedExperiment::colData(SE)[,annotation_trait])))
    names(colours_use) <- unique(SummarizedExperiment::colData(SE)[,annotation_trait])
    ha <- ComplexHeatmap::columnAnnotation(annotation_trait = SummarizedExperiment::colData(SE)[,annotation_trait],
                                           col = list(annotation_trait = colours_use))
  } else if(is.null(annotation_trait)){
    ha <- NULL
  }

  # We get a different column title.
  if(is.null(column_title)){
    column_title <- "Cells"
  }

  heatmap_voi <- ComplexHeatmap::Heatmap(fraction,
                                         column_title_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                                         row_title_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                                         row_names_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                                         col = circlize::colorRamp2(seq(0, round(max(fraction, na.rm = TRUE)), length.out = 9),
                                                                    c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
                                         show_row_names = T, show_column_names = F, cluster_columns = T, clustering_method_columns = "complete", cluster_rows = F, name = "VAF",
                                         heatmap_legend_param = list(border = "#000000", grid_height = grid::unit(10, "mm")),
                                         bottom_annotation = ha, border = T, use_raster = T,
                                         column_title = column_title,
                                         row_title = "Variants")
  return(heatmap_voi)
}
