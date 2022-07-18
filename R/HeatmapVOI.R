#'We plot a heatmap of a set of Variants Of Interest using the Variant Allele Frequency values of a SummarizedExperiment object.
#'@import ComplexHeatmap SummarizedExperiment circlize ggsci Seurat scales
#'@param SE SummarizedExperiment object.
#'@param voi Variants Of Interest.
#'@param annotation_trait Cell Annotation at the bottom of the heat map. 
#'@export
HeatmapVoi <- function(SE, voi, annotation_trait = NULL){
  #colours_list <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
  #                  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
  #                  "#ffff99", "deeppink", "green", "blue", "gold", "indianred3",
  #                  "firebrick4", "#FF6F00FF", "#C71000FF", "#008EA0FF", "#8A4198FF",
  #                  "#5A9599FF", "#FF6348FF", "#84D7E1FF", "#FF95A8FF", "#3D3B25FF",
  #                  "#ADE2D0FF", "#1A5354FF", "#3F4041FF", "chocolate4",  "cornsilk3",
  #                  "cyan", "darkgreen", "gray0", "dodgerblue4", "forestgreen",
  #                  "darkviolet", "indianred3", "midnightblue", "mintcream",
  #                  "mediumaquamarine", "slateblue4"
  #)

  fraction <- assays(SE)[["fraction"]][voi,]
  fraction[is.na(fraction)] <- 0
  if(!is.null(annotation_trait)){
    colours_use <- hue_pal(length(unique(colData(SE)[,annotation_trait])))
    #colours_use <- colours_list[1:length(unique(colData(SE)[,annotation_trait]))]
    #colours_use <- DiscretePalette(length(unique(colData(SE)[,annotation_trait])))
    names(colours_use) <- unique(colData(SE)[,annotation_trait])
    ha <- ComplexHeatmap::columnAnnotation(annotation_trait = colData(SE)[,annotation_trait],
                                           col = list(annotation_trait = colours_use))
  } else if(is.null(annotation_trait)){
    ha <- NULL
  }
  fraction <- as.matrix(fraction)
  heatmap_voi <- Heatmap(fraction,
                         column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                         row_title_gp = gpar(fontsize = 20, fontface = "bold"),
                         row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                         col = colorRamp2(seq(0, round(max(fraction, na.rm = TRUE)), length.out = 9),
                                          c("#FCFCFC","#FFEDB0","#FFDF5F","#FEC510","#FA8E24","#F14C2B","#DA2828","#BE2222","#A31D1D")),
                         show_row_names = T, show_column_names = F, cluster_columns = T, cluster_rows = F, name = "VAF",
                         heatmap_legend_param = list(border = "#000000", grid_height = unit(10, "mm")),
                         bottom_annotation = ha, border = T, use_raster = T,
                         column_title = "Cells",
                         row_title = "Variants")
  return(heatmap_voi)
}
