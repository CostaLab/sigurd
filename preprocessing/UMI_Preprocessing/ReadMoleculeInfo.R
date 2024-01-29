# We read in the molecule_info.h5 file per sample.
# We only use filtered cells.
# We then only retain CB UB combinations with more than min_reads_per_umi occurences.
# We save the combinations in a file and use them to subset the modified bam files.

print("Libraries.")
suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(optparse))


print("Variables.")
option_list = list(
  make_option("--molecule_info",     type = "character", default = "", help = "The molecule info input file.",            metavar = "character"),
  make_option("--barcodes_path",     type = "character", default = "", help = "The barcodes.",                            metavar = "character"),
  make_option("--output",            type = "character", default = "", help = "The output path.",                         metavar = "character"),
  make_option("--sample",            type = "character", default = "", help = "The sample used.",                         metavar = "character"),
  make_option("--min_reads_per_umi", type = "character", default = 10, help = "The minimum number of reads a UMI needs.", metavar = "numeric"),
  make_option("--umi_length",        type = "character", default = 12, help = "The length of a UMI.",                     metavar = "numeric")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

molecule_info_path <- opt$molecule_info
barcodes_path <- opt$barcodes_path
output <- opt$output
sample <- opt$sample
min_reads_per_umi <- as.numeric(opt$min_reads_per_umi)
umi_length <- as.numeric(opt$umi_length)

print("Custom function.")
# This function is from 10X.
# https://kb.10xgenomics.com/hc/en-us/articles/360004105372-How-do-you-decompress-the-2-bit-barcode-sequences-in-molecule-info-h5-file-
UMI_seq <- function(number, umi_length) {
  # Binary string needs to be zero-padded to ensure consistent length
  string_length <- umi_length * 2
  binary_number <- R.utils::intToBin(number)
  if(nchar(binary_number) < string_length){
    zeros_for_padding <- string_length - nchar(binary_number)
    zeros_for_padding <- paste0(base::rep(0, zeros_for_padding), collapse = "")
    binary_number <- paste0(zeros_for_padding, binary_number)
  }
  # Split binary_str into chunks of 2.
  nuc_list <- substring(binary_number, seq(1, nchar(binary_number), by = 2), seq(2, nchar(binary_number), by = 2))
  UMI_seq <- ""
  for (i in nuc_list) {
    if (i == "00") {
      UMI_seq <- paste0(UMI_seq, "A")
    } else if (i == "01") {
      UMI_seq <- paste0(UMI_seq, "C")
    } else if (i == "10") {
      UMI_seq <- paste0(UMI_seq, "G")
    } else {
      UMI_seq <- paste0(UMI_seq, "T")
    }
  }
  return(UMI_seq)
}


print("We read the molecule_info.h5 file.")
molecule_info <- read10xMolInfo(molecule_info_path)
molecule_info <- molecule_info$data
molecule_info$cell <- paste0(molecule_info$cell, "-1")


print("We load the barcodes.")
barcodes <- read.table(barcodes_path, header = FALSE, sep = "\t")[,1]


print("We get the number of UGs (CB+UB combinations) that appear often enough.")
molecule_info <- subset(molecule_info, cell %in% barcodes)
molecule_info <- subset(molecule_info, reads >= min_reads_per_umi)
molecule_info$umi <- unlist(lapply(molecule_info$umi, UMI_seq, umi_length = umi_length))
UGs <- data.frame(paste0(molecule_info$cell, "_", molecule_info$umi))


print("We save the UGs.")
dir.create(paste0(output, sample), showWarnings = FALSE, recursive = TRUE)
write.table(UGs, paste0(output, sample, "/UGs_Filtered.csv"), sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
