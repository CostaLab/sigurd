#'UMIQC
#'@description
#'We load the molecule_info.h5 file from CellRanger and plot the molecules by rank.
#'A molecule is defined as a combination of UMI and CB.
#'The molecules are ranked by their number of replicates.
#'@importFrom DropletUtils read10xMolInfo
#'@param central_input_file Path to the csv file with the samples to be loaded.
#'@param molecule_info_path Path to a single molecule_info file.
#'@param subset_by_cbs Should the list of cell barcodes from the molecule_info.h5 file be subset to a list of approved cell barcodes? Default: TRUE
#'@param cellbarcodes_path Path to a cell barcodes file. Not used if the central input file is used.
#'@param umi_length The length of the UMI. Default = 12
#'@param min_reads The minimum number of reads per molecule.
#'@param sep The separator used in your samples_file. Default: ","
#'@param samples_column The column for the samples in your samples_file. Default: "sample"
#'@param molecule_info_column The column for the paths to the molecule_info files in your samples_file. Default: "molecule_info"
#'@param cellbarcode_column The column for the paths to the cell barcodes files in your central input file. Default: "cells"
#'@param add_10x_suffix Should a suffix be added to the cell barcodes? Default to the 10X standard: "-1"
#'@export
UMIQC <- function(central_input_file = NULL, molecule_info_path = NULL, subset_by_cbs = TRUE, cellbarcodes_path = NULL, umi_length = 12, min_reads = 100, sep = ",", samples_column = "sample", molecule_info_column = "molecule_info", cellbarcode_column = "cells", add_10x_suffix = "-1"){
  if(!is.null(molecule_info_path)){
    if(!file.exists(molecule_info_path)){
      stop("Your molecule_info_path does not point to a file.")
    } else{
      molecule_info <- DropletUtils::read10xMolInfo(molecule_info_path)
      molecule_info <- molecule_info$data
      molecule_info$umi <- UMI_seq(number = molecule_info$umi, umi_length = umi_length)
      if(!is.null(add_10x_suffix)){
        molecule_info[,"cell"] <- paste0(molecule_info[,"cell"], add_10x_suffix)
      }
      molecule_info <- data.frame(cell = molecule_info$cell, umi = molecule_info$umi, molecule = paste0(molecule_info$cell, "_", molecule_info$umi), reads = molecule_info$reads)
      molecule_info <- molecule_info[base::order(molecule_info$reads, decreasing = TRUE),]
      molecule_info$ranks <- 1:nrow(molecule_info)
      molecule_info <- subset(molecule_info, reads >= min_reads)
      if(!is.null(cellbarcodes_path) & !file.exists(cellbarcodes_path)){
        stop("Your cellbarcodes_path does not point to a file.")
      } else{
        cellbarcodes <- read.table(cellbarcodes_path, header = FALSE)
        molecule_info <- subset(molecule_info, cell %in% cellbarcodes[,1])
      }
      cells_retained <- length(unique(molecule_info[,1]))
      umis_per_cell <- lapply(molecule_info[,1], function(x){
        res <- subset(molecule_info, molecule_info$cell == x)
        umis <- nrow(res)
        names(umis) <- x
        return(umis)
      })
      umis_per_cell <- unlist(umis_per_cell)
      sample_summary <- data.frame(Cells = cells_retained, Average_UMIs_Per_Cells = mean(umis_per_cell), Average_Reads_Per_UMI = mean(molecule_info[,"reads"]))
      molecule_info <- list(sample_summary = sample_summary, molecule_info = molecule_info)
    }
  }
  if(!is.null(central_input_file)){
    if(file.exists(central_input_file)){
      samples_file <- read.table(central_input_file, header = TRUE, sep = sep)
      if(!molecule_info_column %in% colnames(samples_file)){
        stop("Your molecule_info_column is not in the columns of your samples file.")
      } else{
        molecule_info <- lapply(1:nrow(samples_file), function(x){
          mol_info <- DropletUtils::read10xMolInfo(samples_file[x, molecule_info_column])$data
          mol_info$umi <- UMI_seq(number = mol_info$umi, umi_length = umi_length)
          if(!is.null(add_10x_suffix)){
            mol_info[,"cell"] <- paste0(mol_info[,"cell"], add_10x_suffix)
          }
          mol_info <- data.frame(Sample = samples_file[x, samples_column], cell = mol_info$cell, umi = mol_info$umi, molecule = paste0(mol_info$cell, "_", mol_info$umi), reads = mol_info$reads)
          mol_info <- mol_info[base::order(mol_info$reads, decreasing = TRUE),]
          mol_info$ranks <- 1:nrow(mol_info)
          mol_info <- subset(mol_info, mol_info$reads >= min_reads)
          if(subset_by_cbs){
            if(!cellbarcode_column %in% colnames(samples_file)) stop("Your cellbarcode_column is not in the central input file.")
            if(!file.exists(samples_file[x, cellbarcode_column])){
              stop("Your cellbarcodes_path does not point to a file.")
            } else{
              cellbarcodes <- read.table(samples_file[x, cellbarcode_column], header = FALSE)
              mol_info <- subset(mol_info, cell %in% cellbarcodes[,1])
            }
          }
          cells_retained <- length(unique(mol_info[,"cell"]))
          umis_per_cell <- lapply(mol_info[,"cell"], function(x){
            res <- subset(mol_info, mol_info$cell == x)
            umis <- nrow(res)
            names(umis) <- x
            return(umis)
          })
          umis_per_cell <- unlist(umis_per_cell)
          sample_summary <- data.frame(Sample = samples_file[x, samples_column], Cells = cells_retained, Average_UMIs_Per_Cells = mean(umis_per_cell), Average_Reads_Per_UMI = mean(mol_info[,"reads"]))
          return(list(sample_summary = sample_summary, molecule_info = mol_info))
        })
        molecule_info <- list(sample_summary = do.call("rbind", lapply(molecule_info, "[[", 1)), 
                              molecule_info = do.call("rbind", lapply(molecule_info, "[[", 2)))
      }
    }
  }
  return(molecule_info)
}

