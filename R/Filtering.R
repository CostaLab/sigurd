#'Filtering the loaded genotyping data.
#'
#'@description
#'We filter a SummarizedExperiment object to exclude variants and cells.
#'
#'@details
#'We do this for one sample at a time.
#'We want to remove:
#' \enumerate{
#'   \item all cells that are blacklisted,
#'   \item all cells that are not in a Seurat object,
#'   \item all cells that do not have at least one variant with >1 (Reference),
#'   \item all variants that are for alternative transcripts,
#'   \item all variants that are always NoCall,
#'   \item set variants with a VAF below a threshold to reference.
#' }
#'@import fastmatch Matrix Seurat SummarizedExperiment
#'@param se SummarizedExperiment object.
#'@param blacklisted_barcodes_path Barcodes you want to remove. Path to a file with one column without header.
#'@param fraction_threshold Variants with an VAF below this threshold are set to 0. Numeric. Default = NULL.
#'@param alts_threshold Variants with a number of alt reads less than this threshold are set to 0. Numeric. Default = NULL.
#'@param min_cells_per_variant In how many cells should a variant be present to be included? Numeric. Default = 2.
#'@param path_seurat Path to a Seurat object. Cells not present in the object will be removed.
#'@param min_variants_per_cell How many variants should be covered in a cell have to be included? Default = 1.
#'@examples
#' \dontrun{
#'   # Removing all variants that are not detected in at least 2 cells.
#'   # Before we remove the variants, we set fraction value of variants below 0.05 to 0.
#'   se_geno <- Filtering(se_geno, min_cells_per_variant = 2, fraction_threshold = 0.05)
#' }
#'@export
Filtering <- function(se, blacklisted_barcodes_path = NULL, fraction_threshold = NULL, alts_threshold = NULL, path_seurat = NULL, min_cells_per_variant = 2, min_variants_per_cell = 1){
  if(!is.null(blacklisted_barcodes_path)){
    print("We remove the unwanted cell barcodes.")
    blacklisted_barcodes <- read.table(blacklisted_barcodes_path, header = FALSE)
    blacklisted_barcodes <- blacklisted_barcodes[,1]
    barcodes_keep <- colnames(se)
    barcodes_keep <- barcodes_keep[!barcodes_keep %in% blacklisted_barcodes]
    se <- se[, barcodes_keep]
  }

  if(!is.null(path_seurat)){
    print("We load the Seurat object.")
    scrna <- load_object(path_seurat)

    print("We remove all cells that are not in the Seurat object.")
    seu_cells <- colnames(se)
    seu_cells <- seu_cells[seu_cells %in% colnames(scrna)]
    se <- se[,seu_cells]
  }

  #print("We remove all variants for alternative transcripts.")
  #keep_variants <- !grepl("_ENST", rownames(se))
  #se <- se[keep_variants,]


  if(!is.null(fraction_threshold)){
    if(any(fraction_threshold >= 1, fraction_threshold <= 0)){
      stop("Your fraction threshold is not 0 < x < 1.")
    }
    print("We assume that cells with a fraction smaller than our threshold are actually reference.")
    print("We set consensus values to 1 (Ref) and fraction values to 0.")
    print(paste0("We do not set fractions between ", fraction_threshold, " and 1 to 1."))
    print("This way, we retain the heterozygous information.")
    # Filtering using sparse matrices.
    consensus_matrix <- assays(se)$consensus
    fraction_matrix <- assays(se)$fraction
    position_matrix <- summary(fraction_matrix)
    position_matrix <- subset(position_matrix, x > 0 & x < fraction_threshold)
    # If no elements fall between 0 and the fraction_threshold, we do not have to change the matrices.
    if(nrow(position_matrix) > 0){
      ij <- as.matrix(position_matrix[, 1:2])
      consensus_matrix[ij] <- 1
      fraction_matrix[ij] <- 0
      assays(se)$consensus <- consensus_matrix
      assays(se)$fraction <- fraction_matrix
    }
  }


  if(!is.null(alts_threshold)){
    if(any(!is.numeric(alts_threshold), is.infinite(alts_threshold))){
      stop("Your alts_threshold should be a numeric value.")
    }
    print("We assume that cells with a number of alternative reads smaller than our threshold are actually reference.")
    print("We set consensus values to 1 (Ref), fraction values to 0 and alts to 0.")
    # Filtering using sparse matrices.
    consensus_matrix <- assays(se)$consensus
    fraction_matrix <- assays(se)$fraction
    coverage_matrix <- assays(se)$coverage
    alts_matrix <- assays(se)$alts
    position_matrix <- summary(alts_matrix)
    position_matrix <- subset(position_matrix, x < alts_threshold)
    # If no elements fall between 0 and the alts_threshold, we do not have to change the matrices.
    if(nrow(position_matrix) > 0){
      ij <- as.matrix(position_matrix[, 1:2])
      consensus_matrix[ij] <- 1
      fraction_matrix[ij] <- 0
      coverage_matrix[ij] <- coverage_matrix[ij] - alts_matrix[ij]
      alts_matrix[ij] <- 0
      assays(se)$consensus <- consensus_matrix
      assays(se)$fraction <- fraction_matrix
      assays(se)$coverage <- coverage_matrix
      assays(se)$alts <- alts_matrix
    }
  }


  print("We remove all the variants that are always NoCall.")
  consensus_test <- assays(se)$consensus > 0
  keep_variants <- rowSums(consensus_test) > 0
  se <- se[keep_variants,]

  print(paste0("We remove variants, that are not at least detected in ", min_cells_per_variant, " cells."))
  keep_variants <- rowSums(assays(se)$consensus >= 1)
  keep_variants <- keep_variants >= min_cells_per_variant
  se <- se[keep_variants,]


  print(paste0("We remove all cells that are not >= 1 (Ref) for at least ", min_variants_per_cell, " variant."))
  consensus_test <- assays(se)$consensus >= 1
  keep_cells <- colSums(consensus_test) > min_variants_per_cell
  se <- se[,keep_cells]
  return(se)
}
