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
min_reads_per_umi <- as.numeric(opt$min_reads_per_umi)
umi_length <- as.numeric(opt$umi_length)
samples_column <- opt$samples_column
molecule_info_column <- opt$molecule_info_column

# We get the UMI QC values and save the overview and the molecule specific information.
umi_qc_values <- sigurd::UMIQC(molecule_info_path = molecule_input_path, umi_length = umi_length, min_reads = min_reads_per_umi, samples_column = samples_column, molecule_info_column = molecule_info_column)
write.table(umi_qc_values[[1]], file.path(output, paste0("UMIQC_Summary_", sample_use, ".csv")), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(umi_qc_values[[2]], file.path(output, paste0("UMIQC_Molecule_Info_", sample_use, ".csv")), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

# We plot the ordered molecules and save the retained molecules for the subsetting of the BAM files.
number_of_molecules_above_threshold <- sum(umi_qc_values[[2]]$reads > min_reads)
molecules_above_threshold <- subset(umi_qc_values[[2]], reads > min_reads)
molecules_above_threshold <- data.frame(molecules_above_threshold$molecule)
write.table(molecules_above_threshold, file.path(output, paste0("CBs_", sample_use, "_Filtered.tsv")), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

number_of_cells_above_threshold <- length(unique(gsub("_.*", "", molecules_above_threshold[,1])))
p <- ggplot2::ggplot(umi_qc_values[[2]], ggplot2::aes(x = ranks, y = reads)) +
  ggplot2::geom_point(color = "blue") + ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
  ggplot2::geom_hline(yintercept = min_reads, col = "red", linewidth = 2) +
  ggplot2::ylab("Number of Reads") + ggplot2::xlab("Rank of Molecule") +
  ggplot2::ggtitle(paste0(sample_use, "\nRetained Cells: ", number_of_cells_above_threshold, ", Retained Molecules: ", number_of_molecules_above_threshold))
ggplot2::ggsave(file.path(output, paste0(sample_use, "_MinReads", min_reads, ".png")), p, width = 4, height = 4, units = "in", dpi = 60)
