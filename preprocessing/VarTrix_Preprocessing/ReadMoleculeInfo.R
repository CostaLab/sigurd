# We read in the molecule_info.h5 file per sample.
# We only use filtered cells.
# We then only retain CB UB combinations with more than min_reads_per_umi occurences.
# We save the combinations in a file and use them to subset the modified bam files.

# Libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(sigurd))

# Variables.
option_list <- list(
  optparse::make_option("--molecule_input_path",  type = "character", default = "", help = "The path to the molecule_input.h5 file.", metavar = "character"),
  optparse::make_option("--output",               type = "character", default = "", help = "The output path.", metavar = "character"),
  optparse::make_option("--sample",               type = "character", default = "", help = "The sample analysed.", metavar = "character"),
  optparse::make_option("--cell_barcodes",        type = "character", default = "", help = "The cell barcodes.", metavar = "character"),
  optparse::make_option("--min_reads_per_umi",    type = "character", default = 10, help = "The minimum number of reads a UMI needs.", metavar = "numeric"),
  optparse::make_option("--umi_length",           type = "character", default = 12, help = "The length of a UMI.", metavar = "numeric"),
  optparse::make_option("--samples_column",       type = "character", default = "", help = "The column of the individual samples in the central input file.", metavar = "character"),
  optparse::make_option("--molecule_info_column", type = "character", default = "", help = "The column containing the path to the molecule_info.h5 file.", metavar = "character")
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
molecule_input_path <- opt$molecule_input_path
output <- opt$output
sample_use <- opt$sample
cell_barcodes <- opt$cell_barcodes
min_reads_per_umi <- as.numeric(opt$min_reads_per_umi)
umi_length <- as.numeric(opt$umi_length)
samples_column <- opt$samples_column
molecule_info_column <- opt$molecule_info_column

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


# We get the UMI QC values and save the overview and the molecule specific information.
umi_qc_values <- UMIQC(
  molecule_info_path = molecule_input_path,
  cellbarcodes_path = cell_barcodes,
  umi_length = umi_length,
  min_reads = min_reads_per_umi,
  samples_column = samples_column,
  molecule_info_column = molecule_info_column)
write.table(umi_qc_values[[1]], file.path(output, paste0("UMIQC_Summary_", sample_use, ".csv")), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(umi_qc_values[[2]], file.path(output, paste0("UMIQC_Molecule_Info_", sample_use, ".csv")), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

# We plot the ordered molecules and save the retained molecules for the subsetting of the BAM files.
number_of_molecules_above_threshold <- nrow(umi_qc_values[[2]])
write.table(umi_qc_values[[2]][, "molecule", drop = FALSE], file.path(output, paste0("Molecules_", sample_use, "_Filtered.tsv")), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

cells_above_threshold <- data.frame(CBs = unique(gsub("_.*", "", umi_qc_values[[2]])))
write.table(cells_above_threshold, file.path(output, paste0("CBs_", sample_use, "_Filtered.tsv")), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

number_of_cells_above_threshold <- nrow(cells_above_threshold)
p <- ggplot2::ggplot(umi_qc_values[[2]], ggplot2::aes(x = ranks, y = reads)) +
  ggplot2::geom_point(color = "blue") + ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
  ggplot2::geom_hline(yintercept = min_reads_per_umi, col = "red", linewidth = 2) +
  ggplot2::ylab("Number of Reads") + ggplot2::xlab("Rank of Molecule") +
  ggplot2::ggtitle(paste0(sample_use, "\nRetained Cells: ", number_of_cells_above_threshold, ", Retained Molecules: ", number_of_molecules_above_threshold))
ggplot2::ggsave(file.path(output, paste0(sample_use, "_MinReads", min_reads_per_umi, ".png")), p, width = 4, height = 4, units = "in", dpi = 60)
