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
  optparse::make_option("--samples_column",       type = "character", default = "sample", help = "The column of the individual samples in the central input file.", metavar = "character"),
  optparse::make_option("--molecule_info_column", type = "character", default = "molecule_info", help = "The column containing the path to the molecule_info.h5 file.", metavar = "character")
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

# We get the UMI QC values and save the overview and the molecule specific information.
umi_qc_values <- sigurd::UMIQC(
  molecule_info_path = molecule_input_path,
  cellbarcodes_path = cell_barcodes,
  umi_length = umi_length,
  min_reads = min_reads_per_umi,
  samples_column = samples_column,
  molecule_info_column = molecule_info_column)
write.table(umi_qc_values[[1]], file.path(output, paste0("UMIQC_Summary_", sample_use, ".csv")), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(umi_qc_values[[2]], file.path(output, paste0("UMIQC_Molecule_Info_", sample_use, ".csv")), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(umi_qc_values[[2]][, "molecule", drop = FALSE], file.path(output, paste0("Molecules_", sample_use, "_Filtered.tsv")), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
cells_above_threshold <- data.frame(CBs = unique(umi_qc_values[[2]][,"cell"]))
write.table(cells_above_threshold, file.path(output, paste0("CBs_", sample_use, "_Filtered.tsv")), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# We plot the ordered molecules and save the retained molecules for the subsetting of the BAM files.
p <- ggplot2::ggplot(umi_qc_values[[2]], ggplot2::aes(x = ranks, y = reads)) +
  ggplot2::geom_point(color = "blue") + ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
  ggplot2::geom_hline(yintercept = min_reads_per_umi, col = "red", linewidth = 2) +
  ggplot2::ylab("Number of Reads") + ggplot2::xlab("Rank of Molecule") +
  ggplot2::ggtitle(paste0(sample_use, "\nRetained Cells: ", nrow(cells_above_threshold), ", Retained Molecules: ", nrow(umi_qc_values[[2]])))
ggplot2::ggsave(file.path(output, paste0(sample_use, "_MinReads", min_reads_per_umi, ".png")), p, width = 4, height = 4, units = "in", dpi = 60)
